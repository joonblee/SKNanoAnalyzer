#!/bin/bash

##############
### setups ###
##############

# run setups
RUN_DT=true
RUN_MC=true
RUN_QCDonly=false   # Available only with 'RUN_MC=true'
RUN_SIG=false

UseSkim=true
UseRun3SignalSkim=false  # false: use raw Sample/ForSNU metadata for Run-3 signals
# NOTE: TriggerEfficiency requires TrigObj branches in the input.
#       If UseSkim=true, the Skim_NIsoMuon output must retain TrigObj_* branches.

# -----------------------------------------------------------------------------
# NIsoMuon analysis modes / systematic modes
#
# Each entry is run INDEPENDENTLY.  The entries are NOT combined with each other.
#
#   ""                    : nominal NIsoDimuon analysis
#   "RunSyst"             : NIsoDimuon experimental-systematic variations
#   "RunXSecSyst"         : NIsoDimuon generator PDF/scale/alpha_s variations
#   "MuonIDEfficiency"    : muon-ID tag-and-probe mode (central only)
#   "TriggerEfficiency"   : HighPtMuon-trigger tag-and-probe mode (central only)
#
# Examples:
#
#   # nominal + both systematic productions
#   flags=("" "RunSyst" "RunXSecSyst")
#
#   # only the two efficiency measurements
#   flags=("MuonIDEfficiency" "TriggerEfficiency")
#
#   # everything
#   flags=("" "RunSyst" "RunXSecSyst" "MuonIDEfficiency" "TriggerEfficiency")
#
# Current selection:
flags=("" "RunSyst" "RunXSecSyst" "MuonIDEfficiency" "TriggerEfficiency")
flags=("" "RunSyst" "RunXSecSyst")

# analysis setups
analysis="NIsoMuon"
skim="Skim_NIsoMuon"
TriggerSets=("HighPtMuon")

signalset_Run2="SampleLists/Run2signal.txt"
signalset_Run3="SampleLists/Run3signal.txt"

mcset_Run2="SampleLists/Run2mc.txt"
#mcset_Run3="SampleLists/Run3mc_all.txt"
mcset_Run3="SampleLists/Run3mc.txt"

qcdset_Run2="SampleLists/Run2qcd.txt"
qcdset_Run3="SampleLists/Run3qcd.txt"

systset_Run2="SampleLists/Run2Syst.txt"
systset_Run3="SampleLists/Run3Syst.txt"

effset_Run3="SampleLists/Run3eff.txt"

Eras=(2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix)
Eras=(2022 2022EE 2023 2023BPix)

nBatch=50
batchname=""       # e.g. "NIsoMuon"
extra_args=()      # e.g. ("--no_exec") or ("--memory" "4000")

mkdir -p log

list_to_input() {
  local listfile="$1"
  local prefix="$2"

  awk -v prefix="$prefix" '
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*$/ {next}
    {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
      print prefix $0
    }
  ' "$listfile" | paste -sd, -
}

prefix_input() {
  local samples="$1"
  local prefix="$2"
  local out=""
  local sample

  IFS=',' read -ra sample_array <<< "$samples"
  for sample in "${sample_array[@]}"; do
    if [[ -n "$out" ]]; then
      out+=","
    fi
    out+="${prefix}${sample}"
  done
  echo "$out"
}

is_efficiency_mode() {
  local flag="$1"
  [[ "$flag" == "MuonIDEfficiency" || "$flag" == "TriggerEfficiency" ]]
}

is_systematic_mode() {
  local flag="$1"
  [[ "$flag" == "RunSyst" || "$flag" == "RunXSecSyst" ]]
}

echo ""
echo "// ------------------------------------------------------ //"
echo "// ----------------- Run SKNanoAnalyzer ----------------- //"
echo "// ------------------------------------------------------ //"
echo ""
echo -n "// time: "
date
echo ""

for trig in "${TriggerSets[@]}"; do

  for era in "${Eras[@]}"; do

    if [[ "$era" == "2016preVFP" || "$era" == "2016postVFP" || "$era" == "2017" || "$era" == "2018" ]]; then
      dataset="SingleMuon"
      signalset="$signalset_Run2"
      mcset="$mcset_Run2"
      qcdset="$qcdset_Run2"
      systset="$systset_Run2"
    elif [[ "$era" == "2022" || "$era" == "2022EE" ]]; then
      dataset="Muon"
      signalset="$signalset_Run3"
      mcset="$mcset_Run3"
      qcdset="$qcdset_Run3"
      systset="$systset_Run3"
    elif [[ "$era" == "2023" || "$era" == "2023BPix" ]]; then
      dataset="Muon0,Muon1"
      signalset="$signalset_Run3"
      mcset="$mcset_Run3"
      qcdset="$qcdset_Run3"
      systset="$systset_Run3"
    fi

    for flag in "${flags[@]}"; do

      cmd_common=(SKNano.py -a "$analysis" -e "$era" -n "$nBatch" --nmax 4000)

      echo " - Era: $era"
      echo " - Trigger: $trig"
      echo " - Dataset: $dataset"
      echo " - MCset: $mcset"

      if [[ -n "$flag" ]]; then
        cmd_common+=(--userflags "$flag")
        echo " - flag: $flag"
        if is_efficiency_mode "$flag"; then
          echo " - mode: central-only efficiency measurement"
        elif [[ "$flag" == "RunSyst" ]]; then
          echo " - mode: NIsoDimuon experimental systematics"
        elif [[ "$flag" == "RunXSecSyst" ]]; then
          echo " - mode: NIsoDimuon generator theory systematics"
        fi
      else
        echo " - flag: <none>"
        echo " - mode: nominal NIsoDimuon"
      fi
      echo ""

      mcset_this="$mcset"
      
      if $RUN_QCDonly; then
        mcset_this="$qcdset"
      fi
      
      if is_systematic_mode "$flag"; then
        mcset_this="$systset"
      fi
      
      if is_efficiency_mode "$flag"; then
        if [[ "$era" == "2022" || "$era" == "2022EE" || \
              "$era" == "2023" || "$era" == "2023BPix" ]]; then
          mcset_this="$effset_Run3"
        fi
      fi

      # MuonIDEfficiency and TriggerEfficiency are central-only modes.
      # They use the ordinary nominal MC list above, not the RunXSecSyst list.

      if [[ -n "$batchname" ]]; then
        cmd_common+=(--batchname "$batchname")
      fi

      if ((${#extra_args[@]})); then
        cmd_common+=("${extra_args[@]}")
      fi

      if $RUN_DT && ! is_systematic_mode "$flag"; then
        if $UseSkim; then
          data_input=$(prefix_input "$dataset" "${skim}_")
        else
          data_input="$dataset"
        fi
      fi

      if $RUN_MC; then
        if $UseSkim; then
          mc_input=$(list_to_input "$mcset_this" "${skim}_")
        else
          mc_input=$(list_to_input "$mcset_this" "")
        fi
      fi

      if $RUN_SIG && ! is_efficiency_mode "$flag" && [[ "$flag" != "RunXSecSyst" ]]; then
        if [[ "$era" == "2022" || "$era" == "2022EE" || \
              "$era" == "2023" || "$era" == "2023BPix" ]]; then
          if $UseRun3SignalSkim; then
            sig_input=$(list_to_input "$signalset" "${skim}_")
          else
            sig_input=$(list_to_input "$signalset" "")
          fi
        elif $UseSkim; then
          sig_input=$(list_to_input "$signalset" "${skim}_")
        else
          sig_input=$(list_to_input "$signalset" "")
        fi
      fi

      if $RUN_DT && ! is_systematic_mode "$flag"; then
        "${cmd_common[@]}" \
          -i "$data_input" \
          &> "log/submit_${era}_${dataset//,/_}_${trig}${flag:+__${flag}}.log"
        echo "[SKNano.py] Run analyzer for data: $dataset trigger: $trig era: $era"
      else
        echo "[SKNano.py] Do not make DATA samples"
      fi

      if $RUN_MC; then
        "${cmd_common[@]}" \
          -i "$mc_input" \
          &> "log/submit_${era}_MC_${trig}${flag:+__${flag}}.log"
        echo "[SKNano.py] Run analyzer for MC backgrounds: $mcset_this trigger: $trig era: $era"
      else
        echo "[SKNano.py] Do not make MC samples"
      fi

      if $RUN_SIG && ! is_efficiency_mode "$flag" && [[ "$flag" != "RunXSecSyst" ]]; then
        "${cmd_common[@]}" \
          -i "$sig_input" \
          &> "log/submit_${era}_sig_${trig}${flag:+__${flag}}.log"
        echo "[SKNano.py] Run analyzer for signals: $signalset trigger: $trig era: $era"
      elif $RUN_SIG && is_efficiency_mode "$flag"; then
        echo "[SKNano.py] Skip signal samples for efficiency mode: $flag"
      else
        echo "[SKNano.py] Do not make signal samples"
      fi

    done
  done
done
