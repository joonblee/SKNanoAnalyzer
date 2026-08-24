#include "MeasureJetTaggingEff.h"
#include "JetTaggingParameter.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

constexpr float kLeadingMuonPt = 52.0f;
constexpr float kSubleadingMuonPt = 10.0f;
constexpr float kMuonEta = 2.4f;

constexpr float kJetPt = 30.0f;
constexpr float kJetEta = 2.4f;
constexpr float kMuonJetDR = 0.3f;
constexpr float kPairMassCut = 2.0f;

bool HasPretagOSDimuonCandidate(
    const RVec<Muon> &muons,
    const RVec<Jet> &jets
) {
    if (muons.size() < 2 || jets.empty()) return false;
    if (!(muons.front().Pt() > kLeadingMuonPt)) return false;

    for (const auto &jet : jets) {
        for (std::size_t iMuon = 0; iMuon < muons.size(); ++iMuon) {
            if (!(muons[iMuon].Pt() > kLeadingMuonPt)) break;
            if (!(jet.DeltaR(muons[iMuon]) < kMuonJetDR)) continue;

            for (std::size_t jMuon = iMuon + 1; jMuon < muons.size(); ++jMuon) {
                if (muons[iMuon].Charge() * muons[jMuon].Charge() >= 0) continue;
                if (!(jet.DeltaR(muons[jMuon]) < kMuonJetDR)) continue;

                const float dimuonMass = (muons[iMuon] + muons[jMuon]).M();
                if (!(dimuonMass > kPairMassCut)) continue;

                return true;
            }
        }
    }

    return false;
}

} // namespace


MeasureJetTaggingEff::MeasureJetTaggingEff()
{
}

MeasureJetTaggingEff::~MeasureJetTaggingEff()
{
}


void MeasureJetTaggingEff::initializeAnalyzer()
{
    if (IsDATA) {
        throw std::runtime_error(
            "MeasureJetTaggingEff is intended for MC only."
        );
    }

    // ------------------------------------------------------------------
    // Branch selection
    // ------------------------------------------------------------------
    //
    // IMPORTANT:
    // "Jet_*" does NOT include "nJet", and "Muon_*" does NOT include
    // "nMuon".  GetAllJets()/GetAllMuons() need these counters.
    //
    // GetAllJets() also re-applies the nominal JEC and therefore needs:
    //   - run
    //   - rho
    //
    // We do not call GetEvent() in this analyzer.  That avoids unnecessarily
    // requiring MET, GenMET, trigger, luminosity-block, and other branches.
    // ------------------------------------------------------------------
    fChain->SetBranchStatus("*", 0);

    auto RequireAndEnable = [&](const TString &branchName) {
        if (fChain->GetBranch(branchName) == nullptr) {
            throw std::runtime_error(
                std::string("MeasureJetTaggingEff: required branch is missing: ") +
                branchName.Data()
            );
        }
        fChain->SetBranchStatus(branchName, 1);
    };

    RequireAndEnable("nJet");
    fChain->SetBranchStatus("Jet_*", 1);

    RequireAndEnable("nMuon");
    fChain->SetBranchStatus("Muon_*", 1);

    // GetAllJets() performs nominal JER smearing for MC via
    //     SmearJets(Jets, GetAllGenJets())
    // so the GenJet collection and MET_pt (used as the random seed)
    // must remain active even though this analyzer does not access them
    // explicitly.
    RequireAndEnable("nGenJet");
    fChain->SetBranchStatus("GenJet_*", 1);
    RequireAndEnable("MET_pt");

    RequireAndEnable("nGenPart");
    fChain->SetBranchStatus("GenPart_*", 1);

    RequireAndEnable("genWeight");
    RequireAndEnable("Pileup_nTrueInt");
    RequireAndEnable("run");

    if (Run == 3) {
        RequireAndEnable("Rho_fixedGridRhoFastjetAll");
    } else if (Run == 2) {
        RequireAndEnable("fixedGridRhoFastjetAll");
    } else {
        throw std::runtime_error(
            "MeasureJetTaggingEff: unsupported Run value."
        );
    }

    myCorr = new MyCorrection(
        DataEra,
        DataPeriod,
        MCSample,
        false
    );

    // NIsoMuon uses DeepJet fixed-WP decisions only.
    // Loose  : b-tag veto in the LightJet control region.
    // Medium : independent b-tagged tag jet in the BJet region.
    Taggers.clear();
    WPs.clear();

    Taggers.push_back(JetTagging::JetFlavTagger::DeepJet);
    WPs.push_back(JetTagging::JetFlavTaggerWP::Loose);
    WPs.push_back(JetTagging::JetFlavTaggerWP::Medium);

    // epsilon_MC(|eta|, pT, flavour)
    vec_etabins = {
        0.0, 0.8, 1.6, 2.0, 2.4
    };

    vec_ptbins = {
        30.0, 40.0, 50.0, 70.0, 100.0,
        140.0, 200.0, 300.0, 600.0, 1000.0
    };

    PtMax = vec_ptbins.back();
    NEtaBin = static_cast<int>(vec_etabins.size()) - 1;
    NPtBin = static_cast<int>(vec_ptbins.size()) - 1;

    etabins = new float[NEtaBin + 1];
    for (int i = 0; i <= NEtaBin; ++i) {
        etabins[i] = vec_etabins.at(i);
    }

    ptbins = new float[NPtBin + 1];
    for (int i = 0; i <= NPtBin; ++i) {
        ptbins[i] = vec_ptbins.at(i);
    }
}


void MeasureJetTaggingEff::executeEvent()
{
    static unsigned long debugEvent = 0;
    const bool debug = (debugEvent < 3);
    if (debug) std::cerr << "[BTagEffDBG] event " << debugEvent
                         << " : enter executeEvent" << std::endl;
    // ------------------------------------------------------------------
    // Unweighted diagnostic cut flow.
    //
    // This is intentionally filled before the tight NIsoMuon selection.
    // If the job processes even one event, hists_*.root should therefore
    // contain CutFlow.
    //
    // bins:
    //   0 : processed event
    //   1 : >= 2 POG-Medium muons with pT > 10 GeV
    //   2 : leading selected muon pT > 52 GeV
    //   3 : >= 1 tight AK4 jet with pT > 30 GeV and |eta| < 2.4
    //   4 : OS dimuon candidate inside one selected jet, m(mumu) > 2 GeV
    // ------------------------------------------------------------------
    if (debug) std::cerr << "[BTagEffDBG] before initial CutFlow" << std::endl;
    FillHist(
        "CutFlow",
        0.5,
        1.0,
        5,
        0.0,
        5.0
    );
    if (debug) std::cerr << "[BTagEffDBG] after initial CutFlow" << std::endl;

    // GetAllJets() needs nJet, Jet_*, run, rho, GenJet_*, and MET_pt.
    // GetAllMuons() needs nMuon and Muon_*.
    if (debug) {
        std::cerr << "[BTagEffDBG] before GetAllJets"
                  << " nJet=" << nJet
                  << " nGenJet=" << nGenJet
                  << " rho=" << fixedGridRhoFastjetAll
                  << " MET_pt=" << MET_pt
                  << " run=" << RunNumber
                  << std::endl;
    }
    const RVec<Jet> allJets = GetAllJets();
    if (debug) std::cerr << "[BTagEffDBG] after GetAllJets size="
                         << allJets.size() << std::endl;

    if (debug) std::cerr << "[BTagEffDBG] before GetAllMuons nMuon="
                         << nMuon << std::endl;
    const RVec<Muon> allMuons = GetAllMuons();
    if (debug) std::cerr << "[BTagEffDBG] after GetAllMuons size="
                         << allMuons.size() << std::endl;

    // ------------------------------------------------------------------
    // NIsoMuon pretag phase space
    // ------------------------------------------------------------------
    //
    // No DeepJet requirement is imposed before the denominator is filled.
    // ------------------------------------------------------------------
    RVec<Muon> muons = SelectMuons(
        allMuons,
        Muon::MuonID::POG_MEDIUM,
        kSubleadingMuonPt,
        kMuonEta
    );

    std::sort(
        muons.begin(),
        muons.end(),
        [](const Muon &a, const Muon &b) {
            return a.Pt() > b.Pt();
        }
    );

    if (muons.size() < 2) {
        if (debug) std::cerr << "[BTagEffDBG] <2 selected muons; return" << std::endl;
        ++debugEvent;
        return;
    }

    FillHist(
        "CutFlow",
        1.5,
        1.0,
        5,
        0.0,
        5.0
    );

    if (!(muons.front().Pt() > kLeadingMuonPt)) {
        if (debug) std::cerr << "[BTagEffDBG] leading muon below threshold; return" << std::endl;
        ++debugEvent;
        return;
    }

    FillHist(
        "CutFlow",
        2.5,
        1.0,
        5,
        0.0,
        5.0
    );

    RVec<Jet> jets = SelectJets(
        allJets,
        Jet::JetID::TIGHT,
        kJetPt,
        kJetEta
    );

    std::sort(
        jets.begin(),
        jets.end(),
        [](const Jet &a, const Jet &b) {
            return a.Pt() > b.Pt();
        }
    );

    if (jets.empty()) {
        if (debug) std::cerr << "[BTagEffDBG] no selected jets; return" << std::endl;
        ++debugEvent;
        return;
    }

    FillHist(
        "CutFlow",
        3.5,
        1.0,
        5,
        0.0,
        5.0
    );

    if (!HasPretagOSDimuonCandidate(muons, jets)) {
        if (debug) std::cerr << "[BTagEffDBG] no candidate; return" << std::endl;
        ++debugEvent;
        return;
    }

    FillHist(
        "CutFlow",
        4.5,
        1.0,
        5,
        0.0,
        5.0
    );

    // ------------------------------------------------------------------
    // MC weight
    // ------------------------------------------------------------------
    //
    // MCweight() gives the sample normalization in 1/pb.
    // That is sufficient when several MC samples are hadd'ed before taking
    // numerator/denominator.  Multiplication by the era luminosity would be
    // a common factor for all samples in the same era and cancels in the
    // efficiency ratio, so GetEvent()/GetTriggerLumi("Full") is unnecessary.
    //
    // PU reweighting is retained because tagging performance depends on PU.
    // ------------------------------------------------------------------
    const float wGenNorm = MCweight();
    const float wPU = myCorr->GetPUWeight(Pileup_nTrueInt);
    const float weight = wGenNorm * wPU;

    // ------------------------------------------------------------------
    // Fill denominator and Loose/Medium numerators
    // ------------------------------------------------------------------
    for (const auto &jet : jets) {
        TString flavour = "0";

        if (std::abs(jet.hadronFlavour()) == 4) {
            flavour = "4";
        }

        if (std::abs(jet.hadronFlavour()) == 5) {
            flavour = "5";
        }

        const float absEta = std::fabs(jet.Eta());

        const float pt =
            jet.Pt() < PtMax
                ? jet.Pt()
                : PtMax - 1.0f;

        // Denominator: all eligible jets, before any DeepJet decision.
        FillHist(
            std::string("tagging#b") +
                "##era#" + DataEra.Data() +
                "##flavor#" + std::string(flavour.Data()) +
                "##systematic#central##den",
            absEta,
            pt,
            weight,
            NEtaBin,
            etabins,
            NPtBin,
            ptbins
        );

        FillHist(
            "DeepJetBTaggingScore" + flavour,
            jet.GetBTaggerResult(
                JetTagging::JetFlavTagger::DeepJet
            ),
            weight,
            100,
            0.0,
            1.0
        );

        for (const auto tagger : Taggers) {
            for (const auto wp : WPs) {

                // Use the direct overload.  No global tagging state and no
                // c-tagging payload are needed.
                const float bTagCut =
                    myCorr->GetBTaggingWP(
                        tagger,
                        wp
                    );

                if (!(jet.GetBTaggerResult(tagger) > bTagCut)) {
                    continue;
                }

                FillHist(
                    std::string("tagging#b") +
                        "##era#" + DataEra.Data() +
                        "##tagger#" +
                        JetTagging::GetTaggerCorrectionLibStr(tagger).Data() +
                        "##working_point#" +
                        JetTagging::GetTaggerCorrectionWPStr(wp).Data() +
                        "##flavor#" +
                        std::string(flavour.Data()) +
                        "##systematic#central##num",
                    absEta,
                    pt,
                    weight,
                    NEtaBin,
                    etabins,
                    NPtBin,
                    ptbins
                );
            }
        }
    }

    if (debug) std::cerr << "[BTagEffDBG] normal end executeEvent" << std::endl;
    ++debugEvent;
}
