import ROOT
import correctionlib.convert as convert
import json
import argparse
import os
import concurrent.futures
from array import array
from tqdm import tqdm


# ---------------------------------------------------------------------
# Analysis efficiency-map binning
# ---------------------------------------------------------------------
# BTV recommends a pT granularity along the lines of
#   [20, 30, 50, 70, 100, 140, 200, 300, 600, 1000] GeV.
# NIsoMuon only uses AK4 jets with pT > 30 GeV, so the efficiency map starts
# at 30 GeV.  The last [600,1000] histogram bin is used as the >=600 GeV bin
# in the correction because the correctionlib node uses flow='clamp'.  Any ROOT
# pT overflow is explicitly merged into that last bin before the efficiency is
# calculated.
TARGET_PT_EDGES = [30.0, 50.0, 70.0, 100.0, 140.0, 200.0, 300.0, 600.0, 1000.0]
TARGET_ETA_EDGES = [0.0, 2.4]


def axis_edges(axis):
    """Return all visible bin edges of a ROOT TAxis as Python floats."""
    return [float(axis.GetBinLowEdge(1))] + [
        float(axis.GetBinUpEdge(i)) for i in range(1, axis.GetNbins() + 1)
    ]


def find_matching_edge(value, edges, tolerance=1.0e-5):
    for edge in edges:
        if abs(float(edge) - float(value)) < tolerance:
            return float(edge)
    return None


def rebin_efficiency_axes(hist, new_name):
    """
    Merge the existing TH2 to the NIsoMuon efficiency-map binning.

    Target binning:
      |eta| : [0.0, 2.4]             (one inclusive analysis bin)
      pT    : [30, 50, 70, 100, 140, 200, 300, 600, inf) GeV

    correctionlib needs a finite upper edge, so the stored pT edges end at
    1000 GeV and flow='clamp' is used.  All visible ROOT bins above 1000 GeV
    and the ROOT pT overflow are explicitly folded into the last [600,1000]
    bin before the efficiency is formed.  Thus that last efficiency value is
    the one used for every pT >= 600 GeV.

    Numerator and denominator are rebinned independently and divided only
    afterwards, so every output efficiency is

        sum(N_pass) / sum(N_total),

    never an average of per-bin efficiencies.
    """
    if not hist or not hist.InheritsFrom('TH2'):
        raise RuntimeError(f'{new_name}: expected a TH2 histogram')

    source_eta_edges = axis_edges(hist.GetXaxis())
    source_pt_edges = axis_edges(hist.GetYaxis())

    # The requested pT operation must be a pure merge: every finite target
    # edge must already exist in the source histogram.
    resolved_pt_edges = []
    for requested in TARGET_PT_EDGES:
        matched = find_matching_edge(requested, source_pt_edges)
        if matched is None:
            raise RuntimeError(
                f'{new_name}: requested pT edge {requested:g} GeV is not present '
                f'in source histogram edges {source_pt_edges}'
            )
        resolved_pt_edges.append(matched)

    # Likewise, the one eta bin must align with the source analysis range.
    eta_low = find_matching_edge(TARGET_ETA_EDGES[0], source_eta_edges)
    eta_high = find_matching_edge(TARGET_ETA_EDGES[1], source_eta_edges)
    if eta_low is None or eta_high is None:
        raise RuntimeError(
            f'{new_name}: requested |eta| range [0,2.4] is not aligned with '
            f'source histogram edges {source_eta_edges}'
        )

    out = ROOT.TH2D(
        new_name, hist.GetTitle(),
        1, array('d', [eta_low, eta_high]),
        len(resolved_pt_edges) - 1, array('d', resolved_pt_edges),
    )
    out.SetDirectory(0)
    out.Sumw2()

    # Sum every visible source eta bin contained in 0 <= |eta| < 2.4.
    for ix in range(1, hist.GetNbinsX() + 1):
        x_low = float(hist.GetXaxis().GetBinLowEdge(ix))
        x_high = float(hist.GetXaxis().GetBinUpEdge(ix))

        if x_high <= eta_low + 1.0e-7:
            continue
        if x_low >= eta_high - 1.0e-7:
            continue
        if x_low < eta_low - 1.0e-5 or x_high > eta_high + 1.0e-5:
            raise RuntimeError(
                f'{new_name}: source eta bin [{x_low},{x_high}] straddles the '
                f'requested [0,2.4] range; refusing to split bins'
            )

        target_x = 1

        # Visible pT bins.
        for iy in range(1, hist.GetNbinsY() + 1):
            low = float(hist.GetYaxis().GetBinLowEdge(iy))
            high = float(hist.GetYaxis().GetBinUpEdge(iy))
            centre = float(hist.GetYaxis().GetBinCenter(iy))

            if high <= resolved_pt_edges[0] + 1.0e-7:
                continue

            # Any visible source bin starting at/above 1000 GeV contributes to
            # the final >=600 GeV efficiency bin.
            if low >= resolved_pt_edges[-1] - 1.0e-7:
                target_y = out.GetNbinsY()
            else:
                target_y = out.GetYaxis().FindFixBin(centre)
                if target_y < 1:
                    continue
                if target_y > out.GetNbinsY():
                    target_y = out.GetNbinsY()

            old_content = float(out.GetBinContent(target_x, target_y))
            old_error = float(out.GetBinError(target_x, target_y))
            add_content = float(hist.GetBinContent(ix, iy))
            add_error = float(hist.GetBinError(ix, iy))
            out.SetBinContent(target_x, target_y, old_content + add_content)
            out.SetBinError(
                target_x,
                target_y,
                (old_error * old_error + add_error * add_error) ** 0.5,
            )

        # ROOT pT overflow is also part of the final >=600 GeV bin.
        iy_overflow = hist.GetNbinsY() + 1
        add_content = float(hist.GetBinContent(ix, iy_overflow))
        add_error = float(hist.GetBinError(ix, iy_overflow))
        if add_content != 0.0 or add_error != 0.0:
            target_y = out.GetNbinsY()
            old_content = float(out.GetBinContent(target_x, target_y))
            old_error = float(out.GetBinError(target_x, target_y))
            out.SetBinContent(target_x, target_y, old_content + add_content)
            out.SetBinError(
                target_x,
                target_y,
                (old_error * old_error + add_error * add_error) ** 0.5,
            )

    return out

# ---------------------------------------------------------------------
# 1. Parse histogram keys from a ROOT file
# ---------------------------------------------------------------------
def histParser(root_file):
    """
    Parse all histogram keys in the given ROOT file.
    Expected key format (example):
      tagging#b##era#2022EE##tagger#deepJet##working_point#M##flavor#5##systematic#central##num
      or for the denominator:
      tagging#b##era#2022EE##flavor#5##systematic#central##den
    Returns a list of dictionaries with extracted fields and the full key.
    """
    hist_list = []
    for keyObj in root_file.GetListOfKeys():
        key = keyObj.GetName()
        if "##" not in key:
            continue
        segments = key.split("##")
        info = {}
        # Process each segment: if it contains '#' then split into field and value.
        for seg in segments:
            if '#' not in seg:
                # Fallback: sometimes the type might be appended at the end.
                info["hist_type"] = seg.strip()
                continue
            parts = seg.split("#", 1)
            if len(parts) < 2:
                continue
            field, value = parts
            info[field] = value
        # Try to set hist_type if not set: look at the key ending.
        if "hist_type" not in info:
            if key.endswith("num"):
                info["hist_type"] = "num"
            elif key.endswith("den"):
                info["hist_type"] = "den"
            elif key.endswith("eff"):
                info["hist_type"] = "eff"
        info["hist_key"] = key
        hist_list.append(info)
    return hist_list

# ---------------------------------------------------------------------
# 2. Create efficiency histograms by matching numerator and denominator histograms
# ---------------------------------------------------------------------
def makeTempEffHist(in_file):
    """
    From the input ROOT file (in_file), find matching numerator and denominator histograms.
    
    Note that the numerator histogram key contains extra fields (tagger and working_point)
    while the denominator is shared among taggers/working points. We match based on:
    
      - tagging
      - era
      - flavor
      - systematic
    
    For each numerator histogram, the corresponding denominator is the one whose key
    matches these fields and has hist_type "den".
    
    The resulting efficiency histogram (with key "eff") is saved to "output.root".
    """
    hist_list = histParser(in_file)
    
    # Separate numerator and denominator histograms
    num_list = [h for h in hist_list if h.get("hist_type", "") == "num"]
    den_list = [h for h in hist_list if h.get("hist_type", "") == "den"]

    eff_hists = []
    # Matching criteria for numerator vs. denominator:
    # Only match on these fields: tagging, era, flavor, systematic.
    match_fields = ["tagging", "era", "flavor", "systematic"]
    for num_info in num_list:
        # For each numerator, try to find a denominator that matches on the common fields.
        matching_den = None
        for den_info in den_list:
            print(f"Will find denominator for numerator: {[num_info[f] for f in match_fields]}")
            print("Total list of denominator: ")
            for den in den_list:
                print(f"Denominator: {[den[f] for f in match_fields]}")
            if all(num_info.get(f) == den_info.get(f) for f in match_fields):
                matching_den = den_info
                break
        if not matching_den:
            # If no matching denominator is found, skip this numerator.
            print(f"WARNING: No denominator found for numerator: {num_info['hist_key']}")
            continue

        num_hist = in_file.Get(num_info["hist_key"])
        den_hist = in_file.Get(matching_den["hist_key"])
        if not num_hist or not den_hist:
            print(f"WARNING: Missing histogram object for {num_info['hist_key']} or {matching_den['hist_key']}")
            continue

        # Rebin numerator and denominator BEFORE calculating the efficiency.
        # Merge |eta| to one inclusive [0,2.4] bin and pT to
        # [30, 50, 70, 100, 140, 200, 300, 600, inf] GeV.
        eff_name = num_info["hist_key"].replace("num", "eff")
        num_rebinned = rebin_efficiency_axes(num_hist, eff_name + "__num_rebinned")
        den_rebinned = rebin_efficiency_axes(den_hist, eff_name + "__den_rebinned")

        eff_hist = num_rebinned.Clone(eff_name)
        eff_hist.SetDirectory(0)
        eff_hist.Divide(den_rebinned)
        eff_hists.append(eff_hist)
        print(f"Created efficiency histogram: {eff_name}")
        print(f"  pT edges: {TARGET_PT_EDGES[:-1]} + [inf]")
        print(f"  |eta| edges: {TARGET_ETA_EDGES}")

    # Write all efficiency histograms to a new ROOT file "output.root"
    out_file = ROOT.TFile("output.root", "RECREATE")
    out_file.cd()
    for hist in eff_hists:
        hist.Write()
    out_file.Close()

# ---------------------------------------------------------------------
# 3. Worker function to convert an efficiency histogram to correction data
# ---------------------------------------------------------------------
def process_histogram(hist_info):
    """
    Reopen the "output.root" file and convert the efficiency histogram given by hist_info.
    Returns a tuple: (tagger, systematic, working_point, flavor, conversion_data)
    """
    hist_path = f"output.root:{hist_info['hist_key']}"
    try:
        conv = convert.from_uproot_THx(hist_path)
        conv_data = conv.dict()['data']
        # Adjust conversion dictionary as in your original Code 2.
        conv_data['inputs'] = ['abseta', 'pt']
        conv_data['flow'] = 'clamp'
    except Exception as e:
        raise RuntimeError(f"Conversion failed for {hist_info['hist_key']}: {e}")

    # For grouping, extract tagger, working_point, and flavor.
    # (For the denominator these fields are missing, but we only process efficiency histograms,
    #  which come from numerators and thus include tagger and working_point.)
    tagger = hist_info.get("tagger")
    systematic = hist_info.get("systematic")
    working_point = hist_info.get("working_point")
    try:
        flavor = int(hist_info.get("flavor"))
    except (TypeError, ValueError):
        flavor = None

    return (tagger, systematic, working_point, flavor, conv_data)

# ---------------------------------------------------------------------
# 4. Build the JSON corrections file dynamically
# ---------------------------------------------------------------------
def makingJson(era, tagging_mode):
    """
    Open the merged efficiency ROOT file ("output.root"), parse its histogram keys,
    filter for a given era and tagging mode ("b" or "c"), convert each histogram in parallel,
    and group the results into a nested JSON structure.
    """
    out_file = ROOT.TFile("output.root", "READ")
    hist_list = histParser(out_file)
    out_file.Close()

    # Filter for efficiency histograms matching the era and tagging mode.
    filtered = [h for h in hist_list 
                if h.get("hist_type", "") == "eff" and
                   h.get("era") == era and
                   h.get("tagging") == tagging_mode]

    conversion_results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_histogram, h) for h in filtered]
        for future in tqdm(concurrent.futures.as_completed(futures),
                           total=len(futures),
                           desc="Converting histograms"):
            try:
                conversion_results.append(future.result())
            except Exception as e:
                print(f"Error converting histogram: {e}")

    if not filtered:
        raise RuntimeError(
            f"No efficiency histograms found for era={era}, tagging={tagging_mode}"
        )

    if len(conversion_results) != len(filtered):
        raise RuntimeError(
            f"Histogram conversion failed for era={era}, tagging={tagging_mode}: "
            f"{len(conversion_results)}/{len(filtered)} converted successfully."
        )

    # Group the conversion results by tagger → systematic → working_point → flavor.
    grouped = {}
    for tagger, systematic, wp, flavor, data in conversion_results:
        if None in (tagger, systematic, wp, flavor):
            continue
        grouped.setdefault(tagger, {}).setdefault(systematic, {}).setdefault(wp, {})[flavor] = data

    corrections = []
    for tagger, syst_dict in grouped.items():
        corr_entry = {
            "name": tagger,
            "version": 0,
            "inputs": [
                {"name": "systematic", "type": "string"},
                {"name": "working_point", "type": "string", "description": "b-tagging working point"},
                {"name": "flavor", "type": "int", "description": "hadron flavor definition: 5=b, 4=c, 0=udsg"},
                {"name": "abseta", "type": "real"},
                {"name": "pt", "type": "real"}
            ],
            "output": {"name": "eff", "type": "real"},
            "data": {
                "nodetype": "category",
                "input": "systematic",
                "content": []
            }
        }
        for systematic, wp_dict in syst_dict.items():
            syst_entry = {
                "key": systematic,
                "value": {
                    "nodetype": "category",
                    "input": "working_point",
                    "content": []
                }
            }
            for wp, flav_dict in wp_dict.items():
                wp_entry = {
                    "key": wp,
                    "value": {
                        "nodetype": "category",
                        "input": "flavor",
                        "content": []
                    }
                }
                for flav, data in flav_dict.items():
                    wp_entry["value"]["content"].append({"key": flav, "value": data})
                wp_entry["value"]["content"].sort(key=lambda x: x["key"])
                syst_entry["value"]["content"].append(wp_entry)
            syst_entry["value"]["content"].sort(key=lambda x: x["key"])
            corr_entry["data"]["content"].append(syst_entry)
        corr_entry["data"]["content"].sort(key=lambda x: x["key"])
        corrections.append(corr_entry)

    if not corrections:
        raise RuntimeError(
            f"No corrections were built for era={era}, tagging={tagging_mode}"
        )

    main_json = {
        "schema_version": 2,
        "description": "This json file contains the b-tagging efficiency corrections",
        "corrections": corrections
    }
    return main_json

# ---------------------------------------------------------------------
# 5. Main script: process the input file and write JSON files
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Create b-tagging efficiency JSON file dynamically")
    parser.add_argument('--input_folder', dest='input_folder',
                        help="Folder containing ROOT files",
                        default=os.path.join(os.environ.get('SKNANO_OUTPUT', ''), 'MeasureJetTaggingEff'))
    parser.add_argument('--input', dest='input',
                        help="Input ROOT file name", default='TTLJ_powheg.root')
    parser.add_argument('--out_name_str', dest='out_name_str',
                        help="Output JSON file name prefix", default='')
    args = parser.parse_args()

    # List of eras to process (adjust as needed)
    totalEras = ['2016preVFP', '2016postVFP', '2017', '2018', '2022', '2022EE', '2023', '2023BPix']
    for era in totalEras:
        out_dir = os.path.join(os.environ.get('SKNANO_DATA', ''), era, 'BTV')
        os.makedirs(out_dir, exist_ok=True)
        input_path = os.path.join(args.input_folder, era, args.input)
        in_file = ROOT.TFile(input_path)
        if not in_file or in_file.IsZombie():
            raise RuntimeError(f"Cannot open file: {input_path}")

        # Create efficiency histograms by matching numerator and denominator histograms.
        makeTempEffHist(in_file)
        in_file.Close()

        # Process both tagging modes: "b" and "c"
        #for tagging_mode, json_suffix in zip(["b", "c"], ["btaggingEff.json", "ctaggingEff.json"]):
        for tagging_mode, json_suffix in zip(["b"], ["btaggingEff.json"]):
            main_json = makingJson(era, tagging_mode)
            out_path = os.path.join(out_dir, args.out_name_str + json_suffix)
            with open(out_path, "w") as fout:
                json.dump(main_json, fout, indent=4)
            print(f"Wrote JSON file: {out_path}")

if __name__ == "__main__":
    main()
