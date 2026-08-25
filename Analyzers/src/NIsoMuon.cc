#include "NIsoMuon.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {
constexpr float kLeadingMuonPt = 52.0f;
constexpr float kSubleadingMuonPt = 10.0f;
constexpr float kMuonEta = 2.4f;
constexpr float kJetPt = 30.0f;
constexpr float kJetEta = 2.4f;
constexpr float kMuonJetDR = 0.3f;
constexpr float kTagJetDR = 0.4f;
constexpr float kPairMassCut = 1.8f;
constexpr float kHistogramMassCut = 2.0f;

const TString kMuonIDSFKey = "NUM_MediumID_DEN_TrackerMuons";
const std::string kRun3HighPtTriggerSFKey =
    "NUM_HLT_DEN_HighPtLooseRelIsoProbes";

constexpr double kMuonIDSFMinPt = 15.0;
constexpr double kMuonIDSFMaxPt = 200.0;

// BTV fixed-WP uncertainty source. The direction is carried by the existing
// ExperimentalSyst::BTagUp/Down enum while the source is kept separately so
// NIsoMuon.h does not need an analysis-specific enum expansion.
TString gActiveBTagSource = "central";

// L1 ECAL prefiring weight variation. 0=nominal, +1=up, -1=down.  Keeping the
// state here likewise avoids changing the shared ExperimentalSyst enum.
int gActiveL1PrefireVariation = 0;

double ClampSFPt(double pt, double minPt, double maxPt) {
    if (pt < minPt) {
        return minPt;
    }

    // correctionlib uses [low, high), so stay infinitesimally
    // below the upper boundary.
    if (pt >= maxPt) {
        return std::nextafter(maxPt, minPt);
    }

    return pt;
}
}

NIsoMuon::NIsoMuon()
    : analysisMode(AnalysisMode::NIsoDimuon),
      runSyst(false),
      runXSecSyst(false),
      deepJetMediumWP(0.0f),
      deepJetLooseWP(0.0f) {}

NIsoMuon::~NIsoMuon() {}

TString NIsoMuon::CategoryName(Category category) {
    return category == Category::BJet ? "BJet" : "LightJet";
}

TString NIsoMuon::SignName(DileptonSign sign) {
    return sign == DileptonSign::OS ? "OS" : "SS";
}

void NIsoMuon::AddTriggerIfAvailable(
    const TString &name,
    bool coveredByTriggerSF
) {
    if (TriggerMap.find(name) != TriggerMap.end()) {
        highPtMuonTriggers.push_back(name);
        if (coveredByTriggerSF) triggerSFCoveredPaths.push_back(name);
        return;
    }
    std::cout << "[NIsoMuon::ConfigureTriggers] HLT path unavailable in this era/input: "
              << name << std::endl;
}

void NIsoMuon::ConfigureTriggers() {
    highPtMuonTriggers.clear();
    triggerSFCoveredPaths.clear();
    highPtEfficiencyTriggers.clear();
    isolatedMuonTriggers.clear();

    auto addTo = [this](const TString &name, RVec<TString> &collection) {
        if (TriggerMap.find(name) != TriggerMap.end()) {
            collection.push_back(name);
        } else {
            std::cout << "[NIsoMuon::ConfigureTriggers] HLT path unavailable in this era/input: "
                      << name << std::endl;
        }
    };

    if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
        AddTriggerIfAvailable("HLT_Mu50", true);
        AddTriggerIfAvailable("HLT_TkMu50", true);

        if (analysisMode != AnalysisMode::NIsoDimuon) {
            addTo("HLT_Mu50", highPtEfficiencyTriggers);
            addTo("HLT_TkMu50", highPtEfficiencyTriggers);
        }
        if (analysisMode == AnalysisMode::TriggerEfficiency) {
            addTo("HLT_IsoMu24", isolatedMuonTriggers);
            addTo("HLT_IsoTkMu24", isolatedMuonTriggers);
        }
    } else if (DataEra == "2017") {
        AddTriggerIfAvailable("HLT_Mu50", true);
        AddTriggerIfAvailable("HLT_OldMu100", true);
        AddTriggerIfAvailable("HLT_TkMu100", true);

        if (analysisMode != AnalysisMode::NIsoDimuon) {
            addTo("HLT_Mu50", highPtEfficiencyTriggers);
            addTo("HLT_OldMu100", highPtEfficiencyTriggers);
            addTo("HLT_TkMu100", highPtEfficiencyTriggers);
        }
        if (analysisMode == AnalysisMode::TriggerEfficiency) {
            addTo("HLT_IsoMu27", isolatedMuonTriggers);
        }
    } else if (DataEra == "2018") {
        AddTriggerIfAvailable("HLT_Mu50", true);
        AddTriggerIfAvailable("HLT_OldMu100", true);
        AddTriggerIfAvailable("HLT_TkMu100", true);

        if (analysisMode != AnalysisMode::NIsoDimuon) {
            addTo("HLT_Mu50", highPtEfficiencyTriggers);
            addTo("HLT_OldMu100", highPtEfficiencyTriggers);
            addTo("HLT_TkMu100", highPtEfficiencyTriggers);
        }
        if (analysisMode == AnalysisMode::TriggerEfficiency) {
            addTo("HLT_IsoMu24", isolatedMuonTriggers);
        }
    } else if (DataEra == "2022" || DataEra == "2022EE" ||
               DataEra == "2023" || DataEra == "2023BPix") {
        AddTriggerIfAvailable("HLT_Mu50", true);
        AddTriggerIfAvailable("HLT_CascadeMu100", true);
        AddTriggerIfAvailable("HLT_HighPtTkMu100", true);

        // Kept only in the main NIsoDimuon event-level OR. It is deliberately
        // excluded from the dedicated measured HighPtMuon OR below.
        AddTriggerIfAvailable("HLT_Mu50_L1SingleMuShower", false);

        if (analysisMode != AnalysisMode::NIsoDimuon) {
            addTo("HLT_Mu50", highPtEfficiencyTriggers);
            addTo("HLT_CascadeMu100", highPtEfficiencyTriggers);
            addTo("HLT_HighPtTkMu100", highPtEfficiencyTriggers);
        }
        if (analysisMode == AnalysisMode::TriggerEfficiency) {
            addTo("HLT_IsoMu24", isolatedMuonTriggers);
        }
    } else {
        std::cerr << "[NIsoMuon::ConfigureTriggers] Unsupported era: "
                  << DataEra << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (highPtMuonTriggers.empty()) {
        std::cerr << "[NIsoMuon::ConfigureTriggers] No requested high-pT trigger is available"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (analysisMode != AnalysisMode::NIsoDimuon &&
        highPtEfficiencyTriggers.empty()) {
        std::cerr << "[NIsoMuon::ConfigureTriggers] No target HighPtMuon trigger is available"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (analysisMode == AnalysisMode::TriggerEfficiency &&
        isolatedMuonTriggers.empty()) {
        std::cerr << "[NIsoMuon::ConfigureTriggers] No isolated reference trigger is available"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    auto printTriggers = [](const char *label, const RVec<TString> &collection) {
        std::cout << label;
        for (const auto &name : collection) std::cout << " " << name;
        std::cout << std::endl;
    };

    printTriggers("[NIsoMuon::ConfigureTriggers] main HighPtMuon OR:",
                  highPtMuonTriggers);
    if (analysisMode != AnalysisMode::NIsoDimuon) {
        printTriggers("[NIsoMuon::ConfigureTriggers] measured HighPtMuon OR:",
                      highPtEfficiencyTriggers);
    }
    if (analysisMode == AnalysisMode::TriggerEfficiency) {
        printTriggers("[NIsoMuon::ConfigureTriggers] isolated reference OR:",
                      isolatedMuonTriggers);
    }
}

void NIsoMuon::ConfigureMuonScaleFactors() {
    run2TriggerSFKey = "";
    run3MuonZCorrections.reset();
    run3HighPtMuonCorrections.reset();

    if (IsDATA) return;

    if (Run == 2) {
        if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
            run2TriggerSFKey =
                "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose";
        } else if (DataEra == "2017" || DataEra == "2018") {
            run2TriggerSFKey =
                "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose";
        } else {
            std::cerr << "[NIsoMuon::ConfigureMuonScaleFactors] Unsupported Run-2 era: "
                      << DataEra << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::cout << "[NIsoMuon::ConfigureMuonScaleFactors] Medium-ID SF key: "
                  << kMuonIDSFKey << std::endl;
        std::cout << "[NIsoMuon::ConfigureMuonScaleFactors] High-pT trigger SF key: "
                  << run2TriggerSFKey << std::endl;
        return;
    }

    TString jsonFolder;
    if (DataEra == "2022") {
        jsonFolder = "2022_Summer22";
    } else if (DataEra == "2022EE") {
        jsonFolder = "2022_Summer22EE";
    } else if (DataEra == "2023") {
        jsonFolder = "2023_Summer23";
    } else if (DataEra == "2023BPix") {
        jsonFolder = "2023_Summer23BPix";
    } else {
        std::cerr << "[NIsoMuon::ConfigureMuonScaleFactors] Unsupported Run-3 era: "
                  << DataEra << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const char *jsonPogPath = std::getenv("JSONPOG_REPO_PATH");
    if (jsonPogPath == nullptr) {
        std::cerr << "[NIsoMuon::ConfigureMuonScaleFactors] JSONPOG_REPO_PATH is not set"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const std::string muonZFile =
        std::string(jsonPogPath) + "/POG/MUO/" +
        std::string(jsonFolder.Data()) + "/muon_Z.json.gz";
    const std::string highPtFile =
        std::string(jsonPogPath) + "/POG/MUO/" +
        std::string(jsonFolder.Data()) + "/muon_HighPt.json.gz";

    try {
        run3MuonZCorrections =
            correction::CorrectionSet::from_file(muonZFile);
        run3HighPtMuonCorrections =
            correction::CorrectionSet::from_file(highPtFile);
        (void)run3MuonZCorrections->at(std::string(kMuonIDSFKey.Data()));
        (void)run3HighPtMuonCorrections->at(kRun3HighPtTriggerSFKey);
    } catch (const std::exception &error) {
        std::cerr << "[NIsoMuon::ConfigureMuonScaleFactors] Failed to load Run-3 muon SF payloads"
                  << std::endl;
        std::cerr << "  ID file:      " << muonZFile << std::endl;
        std::cerr << "  ID key:       " << kMuonIDSFKey << std::endl;
        std::cerr << "  trigger file: " << highPtFile << std::endl;
        std::cerr << "  trigger key:  " << kRun3HighPtTriggerSFKey << std::endl;
        std::cerr << "  error: " << error.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::cout << "[NIsoMuon::ConfigureMuonScaleFactors] Medium-ID SF file: "
              << muonZFile << std::endl;
    std::cout << "[NIsoMuon::ConfigureMuonScaleFactors] Medium-ID SF key: "
              << kMuonIDSFKey << std::endl;
    std::cout << "[NIsoMuon::ConfigureMuonScaleFactors] High-pT trigger SF file: "
              << highPtFile << std::endl;
    std::cout << "[NIsoMuon::ConfigureMuonScaleFactors] High-pT trigger SF key: "
              << kRun3HighPtTriggerSFKey << std::endl;

    if (TriggerMap.find("HLT_Mu50_L1SingleMuShower") != TriggerMap.end()) {
        std::cout << "[NIsoMuon::ConfigureMuonScaleFactors] WARNING: "
                  << "HLT_Mu50_L1SingleMuShower is retained in the event OR, "
                  << "but it is not covered by the selected Muon POG trigger SF. "
                  << "Events firing only this path receive trigger SF = 1."
                  << std::endl;
    }
}

void NIsoMuon::initializeAnalyzer() {
    const bool runMuonIDMode = HasFlag("MuonIDEfficiency");
    const bool runTriggerMode = HasFlag("TriggerEfficiency");

    if (runMuonIDMode && runTriggerMode) {
        std::cerr << "[NIsoMuon::initializeAnalyzer] Select only one efficiency mode: "
                  << "MuonIDEfficiency or TriggerEfficiency" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    analysisMode = AnalysisMode::NIsoDimuon;
    if (runMuonIDMode) analysisMode = AnalysisMode::MuonIDEfficiency;
    if (runTriggerMode) analysisMode = AnalysisMode::TriggerEfficiency;

    runSyst = HasFlag("RunSyst");
    runXSecSyst = HasFlag("RunXSecSyst");

    if (IsDATA && (runSyst || runXSecSyst)) {
        std::cerr << "[NIsoMuon::initializeAnalyzer] RunSyst and RunXSecSyst are "
                  << "MC-only modes; refusing to run on data." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (analysisMode != AnalysisMode::NIsoDimuon &&
        (runSyst || runXSecSyst)) {
        std::cerr << "[NIsoMuon::initializeAnalyzer] Efficiency modes are central-only. "
                  << "Do not combine MuonIDEfficiency/TriggerEfficiency with "
                  << "RunSyst or RunXSecSyst." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const char *modeName = "NIsoDimuon";
    if (analysisMode == AnalysisMode::MuonIDEfficiency) {
        modeName = "MuonIDEfficiency";
    } else if (analysisMode == AnalysisMode::TriggerEfficiency) {
        modeName = "TriggerEfficiency";
    }

    std::cout << "[NIsoMuon::initializeAnalyzer] Analysis mode = "
              << modeName << std::endl;
    std::cout << "[NIsoMuon::initializeAnalyzer] RunSyst = "
              << runSyst << std::endl;
    std::cout << "[NIsoMuon::initializeAnalyzer] RunXSecSyst = "
              << runXSecSyst << std::endl;

    myCorr = new MyCorrection(
        DataEra,
        DataPeriod,
        IsDATA ? DataStream : MCSample,
        IsDATA
    );

    // Keep the existing NIsoDimuon setup unchanged. The efficiency modes do
    // not use b-tagging or the scale factors that they are intended to measure.
    if (analysisMode == AnalysisMode::NIsoDimuon) {
        deepJetMediumWP = myCorr->GetBTaggingWP(
            JetTagging::JetFlavTagger::DeepJet,
            JetTagging::JetFlavTaggerWP::Medium
        );
        deepJetLooseWP = myCorr->GetBTaggingWP(
            JetTagging::JetFlavTagger::DeepJet,
            JetTagging::JetFlavTaggerWP::Loose
        );
    }

    ConfigureTriggers();

    if (analysisMode == AnalysisMode::NIsoDimuon) {
        ConfigureMuonScaleFactors();
    }
}

bool NIsoMuon::PassDileptonSign(
    DileptonSign sign,
    const Muon &muon1,
    const Muon &muon2
) const {
    const int chargeProduct = muon1.Charge() * muon2.Charge();
    return sign == DileptonSign::OS ? chargeProduct < 0 : chargeProduct > 0;
}

bool NIsoMuon::FindCandidate(
    Category category,
    DileptonSign sign,
    const RVec<Muon> &muons,
    const RVec<Jet> &allJets,
    const RVec<Jet> &tagJets,
    Muon &leadingMuon,
    Muon &subleadingMuon,
    Jet &dimuonJet,
    Jet &tagJet
) const {
    if (muons.size() < 2 || allJets.empty()) return false;
    if (!(muons[0].Pt() > kLeadingMuonPt)) return false;

    for (const auto &jet : allJets) {
        for (std::size_t iMuon = 0; iMuon < muons.size(); ++iMuon) {
            if (!(muons[iMuon].Pt() > kLeadingMuonPt)) break;
            if (!(jet.DeltaR(muons[iMuon]) < kMuonJetDR)) continue;

            for (std::size_t jMuon = iMuon + 1; jMuon < muons.size(); ++jMuon) {
                if (!PassDileptonSign(sign, muons[iMuon], muons[jMuon])) continue;
                if (!(jet.DeltaR(muons[jMuon]) < kMuonJetDR)) continue;

                const float mass = (muons[iMuon] + muons[jMuon]).M();
                if (!(mass > kPairMassCut)) continue;

                if (category == Category::BJet) {
                    bool foundTagJet = false;
                    for (const auto &candidateTagJet : tagJets) {
                        if (!(candidateTagJet.DeltaR(jet) > kTagJetDR)) continue;
                        tagJet = candidateTagJet;
                        foundTagJet = true;
                        break;
                    }
                    if (!foundTagJet) continue;
                }

                leadingMuon = muons[iMuon];
                subleadingMuon = muons[jMuon];
                dimuonJet = jet;
                return true;
            }
        }
    }

    return false;
}

float NIsoMuon::CalculateBTagWeight(
    const RVec<Jet> &jets,
    JetTagging::JetFlavTaggerWP workingPoint,
    ExperimentalSyst syst
) {
    using Variation = MyCorrection::variation;
    const auto tagger = JetTagging::JetFlavTagger::DeepJet;
    const auto method = JetTagging::JetTaggingSFMethod::comb;

    // BTV fixed-WP uncertainties are split into heavy-flavour (b/c) and
    // light-flavour components.  Explicitly split the input jets so that an
    // hf_* source changes only b/c jets and an lf_* source changes only
    // light-flavour jets.  This also avoids passing an hf source to the
    // light-flavour correction set, or vice versa.
    RVec<Jet> heavyJets;
    RVec<Jet> lightJets;
    heavyJets.reserve(jets.size());
    lightJets.reserve(jets.size());

    for (const auto &jet : jets) {
        const int flavour = std::abs(jet.hadronFlavour());
        if (flavour == 4 || flavour == 5) {
            heavyJets.push_back(jet);
        } else {
            lightJets.push_back(jet);
        }
    }

    auto evaluate = [&](const RVec<Jet> &selectedJets,
                        Variation variation,
                        const TString &source) -> float {
        if (selectedJets.empty()) return 1.0f;
        return myCorr->GetBTaggingReweightMethod1a(
            selectedJets,
            tagger,
            workingPoint,
            method,
            variation,
            source
        );
    };

    const float centralHF = evaluate(heavyJets, Variation::nom, "central");
    const float centralLF = evaluate(lightJets, Variation::nom, "central");
    const float central = centralHF * centralLF;

    if (syst != ExperimentalSyst::BTagUp &&
        syst != ExperimentalSyst::BTagDown) {
        return central;
    }

    const Variation variation =
        (syst == ExperimentalSyst::BTagUp) ? Variation::up : Variation::down;

    if (gActiveBTagSource == "hf_corr" ||
        gActiveBTagSource == "hf_uncorr") {
        return evaluate(heavyJets, variation, gActiveBTagSource) * centralLF;
    }

    if (gActiveBTagSource == "lf_corr" ||
        gActiveBTagSource == "lf_uncorr") {
        return centralHF * evaluate(lightJets, variation, gActiveBTagSource);
    }

    std::cerr << "[NIsoMuon::CalculateBTagWeight] Invalid BTV uncertainty source: "
              << gActiveBTagSource << std::endl;
    std::exit(EXIT_FAILURE);
    return central;
}

float NIsoMuon::CalculateMuonIDSF(
    const RVec<Muon> &muons,
    ExperimentalSyst syst
) const {
    if (IsDATA) return 1.0f;

    using Variation = MyCorrection::variation;
    Variation variation = Variation::nom;
    if (syst == ExperimentalSyst::MuonIDSFUp) {
        variation = Variation::up;
    } else if (syst == ExperimentalSyst::MuonIDSFDown) {
        variation = Variation::down;
    }

    //if (Run == 2) {
    //    return myCorr->GetMuonIDSF(kMuonIDSFKey, muons, variation);
    //}
    if (Run == 2) {
        RVec<Muon> sfMuons = muons;
    
        for (auto &muon : sfMuons) {
            const double ptForSF = ClampSFPt(
                static_cast<double>(muon.OriginalPt()),
                kMuonIDSFMinPt,
                kMuonIDSFMaxPt
            );
    
            muon.SetOriginalPt(static_cast<float>(ptForSF));
        }
    
        return myCorr->GetMuonIDSF(
            kMuonIDSFKey,
            sfMuons,
            variation
        );
    }

    if (!run3MuonZCorrections) {
        std::cerr << "[NIsoMuon::CalculateMuonIDSF] Run-3 muon-Z payload is not loaded"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string variationName = "nominal";
    if (variation == Variation::up) variationName = "systup";
    if (variation == Variation::down) variationName = "systdown";

    try {
        const auto correction = run3MuonZCorrections->at(
            std::string(kMuonIDSFKey.Data())
        );

        float weight = 1.0f;
        for (const auto &muon : muons) {
            // Run-3 2023 Z-based Muon POG SFs use signed eta.  The 2022
            // payload accepts the same eta input and internally applies the
            // required absolute-value treatment.
            const double eta = std::max(
                -2.4 + 1.0e-6,
                std::min(static_cast<double>(muon.Eta()), 2.4 - 1.0e-6)
            );
            //const double pt = std::max(
            //    2.0 + 1.0e-3,
            //    static_cast<double>(muon.Pt())
            //);
            const double pt = ClampSFPt(
                static_cast<double>(muon.Pt()),
                kMuonIDSFMinPt,
                kMuonIDSFMaxPt
            );

            weight *= static_cast<float>(correction->evaluate({
                eta,
                pt,
                variationName
            }));
        }
        return weight;
    } catch (const std::exception &error) {
        std::cerr << "[NIsoMuon::CalculateMuonIDSF] Failed to evaluate Run-3 Medium-ID SF"
                  << std::endl;
        std::cerr << "  era = " << DataEra
                  << ", variation = " << variationName << std::endl;
        std::cerr << "  error: " << error.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return 1.0f;
}

float NIsoMuon::CalculateMuonTriggerSF(
    const Event &ev,
    const Muon &leadingMuon,
    ExperimentalSyst syst
) const {
    if (IsDATA) return 1.0f;

    // No official SF is assigned to an event that fires only a path outside
    // the correction payload, notably HLT_Mu50_L1SingleMuShower in Run 3.
    if (triggerSFCoveredPaths.empty() ||
        !ev.PassTrigger(triggerSFCoveredPaths)) {
        return 1.0f;
    }

    using Variation = MyCorrection::variation;
    Variation variation = Variation::nom;
    if (syst == ExperimentalSyst::MuonTriggerSFUp) {
        variation = Variation::up;
    } else if (syst == ExperimentalSyst::MuonTriggerSFDown) {
        variation = Variation::down;
    }

    //if (Run == 2) {
    //    return myCorr->GetMuonIDSF( run2TriggerSFKey, leadingMuon, variation );
    //}

    if (Run == 2) {
        Muon sfMuon = leadingMuon;
    
        const double ptForSF = std::max(
            52.0 + 1.0e-3,
            static_cast<double>(leadingMuon.OriginalPt())
        );
    
        sfMuon.SetOriginalPt(
            static_cast<float>(ptForSF)
        );
    
        return myCorr->GetMuonIDSF(
            run2TriggerSFKey,
            sfMuon,
            variation
        );
    }

    if (!run3HighPtMuonCorrections) {
        std::cerr << "[NIsoMuon::CalculateMuonTriggerSF] Run-3 trigger SF payload is not loaded"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string variationName = "nominal";
    if (variation == Variation::up) variationName = "systup";
    if (variation == Variation::down) variationName = "systdown";

    // The official high-pT HLT payload is tabulated from 50 GeV and up to
    // 1000 GeV. The analysis threshold is above 52 GeV; only the upper edge
    // therefore needs a practical clamp for exceptional events.
    const double eta = std::min(
        std::abs(static_cast<double>(leadingMuon.Eta())),
        2.4 - 1.0e-6
    );
    const double pt = std::max(
        50.0 + 1.0e-3,
        std::min(static_cast<double>(leadingMuon.Pt()), 1000.0 - 1.0e-3)
    );

    try {
        const auto correction =
            run3HighPtMuonCorrections->at(kRun3HighPtTriggerSFKey);
        return static_cast<float>(correction->evaluate({
            eta,
            pt,
            variationName
        }));
    } catch (const std::exception &error) {
        std::cerr << "[NIsoMuon::CalculateMuonTriggerSF] Failed to evaluate high-pT trigger SF"
                  << std::endl;
        std::cerr << "  era = " << DataEra
                  << ", eta = " << eta
                  << ", pt = " << pt
                  << ", variation = " << variationName << std::endl;
        std::cerr << "  error: " << error.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return 1.0f;
}

float NIsoMuon::CalculateWeight(
    const Event &ev,
    Category category,
    ExperimentalSyst syst,
    const RVec<Muon> &selectedMuons,
    const RVec<Jet> &eligibleJets,
    const Jet &tagJet,
    float theoryWeight
) {
    if (IsDATA) return 1.0f;

    using Variation = MyCorrection::variation;

    Variation puVariation = Variation::nom;
    if (syst == ExperimentalSyst::PUUp) puVariation = Variation::up;
    if (syst == ExperimentalSyst::PUDown) puVariation = Variation::down;

    Variation l1PrefireVariation = Variation::nom;
    if (gActiveL1PrefireVariation > 0) l1PrefireVariation = Variation::up;
    if (gActiveL1PrefireVariation < 0) l1PrefireVariation = Variation::down;

    float weight = MCweight();
    weight *= ev.GetTriggerLumi("Full");
    weight *= GetL1PrefireWeight(l1PrefireVariation);
    weight *= myCorr->GetPUWeight(ev.nTrueInt(), puVariation);
    weight *= CalculateMuonIDSF(selectedMuons, syst);
    weight *= CalculateMuonTriggerSF(
        ev,
        selectedMuons.at(0),
        syst
    );

    if (category == Category::BJet) {
        const RVec<Jet> tagJets = {tagJet};
        weight *= CalculateBTagWeight(
            tagJets,
            JetTagging::JetFlavTaggerWP::Medium,
            syst
        );
    } else {
        weight *= CalculateBTagWeight(
            eligibleJets,
            JetTagging::JetFlavTaggerWP::Loose,
            syst
        );
    }

    weight *= theoryWeight;
    return weight;
}

void NIsoMuon::RunVariation(
    const Event &ev,
    const RVec<Muon> &allMuons,
    const RVec<Jet> &allJets,
    const RVec<GenJet> &allGenJets,
    ExperimentalSyst syst,
    const TString &systSuffix,
    float theoryWeight
) {
    using Variation = MyCorrection::variation;

    RVec<Muon> shiftedMuons = allMuons;
    RVec<Jet> shiftedJets = allJets;

    if (!IsDATA) {
        if (syst == ExperimentalSyst::JetResUp) {
            shiftedJets = SmearJets(allJets, allGenJets, Variation::up);
        } else if (syst == ExperimentalSyst::JetResDown) {
            shiftedJets = SmearJets(allJets, allGenJets, Variation::down);
        } else if (syst == ExperimentalSyst::JetEnUp) {
            shiftedJets = ScaleJets(allJets, Variation::up);
        } else if (syst == ExperimentalSyst::JetEnDown) {
            shiftedJets = ScaleJets(allJets, Variation::down);
        } else if (syst == ExperimentalSyst::MuonEnUp) {
            shiftedMuons = ScaleMuons(allMuons, "up");
        } else if (syst == ExperimentalSyst::MuonEnDown) {
            shiftedMuons = ScaleMuons(allMuons, "down");
        }
    }

    RVec<Muon> muons = SelectMuons(
        shiftedMuons,
        Muon::MuonID::POG_MEDIUM,
        kSubleadingMuonPt,
        kMuonEta
    );
    std::sort(muons.begin(), muons.end(), PtComparing);
    if (muons.size() < 2) return;
    if (!(muons[0].Pt() > kLeadingMuonPt)) return;

    RVec<Jet> jets = SelectJets(
        shiftedJets,
        Jet::JetID::TIGHT,
        kJetPt,
        kJetEta
    );
    std::sort(jets.begin(), jets.end(), PtComparing);
    if (jets.empty()) return;

    RVec<Jet> mediumBJets;
    bool hasLooseBJet = false;
    for (const auto &jet : jets) {
        const float score = jet.GetBTaggerResult(
            JetTagging::JetFlavTagger::DeepJet
        );
        if (score > deepJetMediumWP) mediumBJets.push_back(jet);
        if (score > deepJetLooseWP) hasLooseBJet = true;
    }

    const std::array<Category, 2> categories = {
        Category::BJet,
        Category::LightJet
    };
    const std::array<DileptonSign, 2> signs = {
        DileptonSign::OS,
        DileptonSign::SS
    };

    for (const auto category : categories) {
        if (category == Category::BJet && mediumBJets.empty()) continue;
        if (category == Category::LightJet && hasLooseBJet) continue;

        for (const auto sign : signs) {
            if ((runSyst || runXSecSyst) &&
                (category != Category::BJet || sign != DileptonSign::OS)) {
                continue;
            }
            if (category == Category::LightJet && sign == DileptonSign::SS) {
                continue;
            }
            if (MCSample.Contains("Zp") &&
                (category == Category::LightJet || sign == DileptonSign::SS)) {
                continue;
            }

            Muon leadingMuon;
            Muon subleadingMuon;
            Jet dimuonJet;
            Jet tagJet;

            if (!FindCandidate(
                    category,
                    sign,
                    muons,
                    jets,
                    mediumBJets,
                    leadingMuon,
                    subleadingMuon,
                    dimuonJet,
                    tagJet)) {
                continue;
            }

            const float dimuonMass = (leadingMuon + subleadingMuon).M();
            if (!(dimuonMass > kHistogramMassCut)) continue;

            const bool outsideUpsilon =
                dimuonMass < 9.0f || dimuonMass > 10.4f;
            if (!(sign == DileptonSign::SS ||
                  category == Category::LightJet ||
                  outsideUpsilon)) {
                continue;
            }

            const RVec<Muon> selectedMuons = {
                leadingMuon,
                subleadingMuon
            };

            const float weight = CalculateWeight(
                ev,
                category,
                syst,
                selectedMuons,
                jets,
                tagJet,
                theoryWeight
            );

            const TString baseRegion =
                SignName(sign) +
                "_POGMedium_tight_" +
                CategoryName(category);
            const TString region =
                baseRegion + systSuffix + "_NIsoDimuon";
            const TString histogram =
                region + "/Dilepton_Mass___" + region;

            // Fill hisgrams //
            FillHist(histogram, dimuonMass, weight, 7500, 0.0, 150.0);
  					
            if (!(10.4 < dimuonMass && dimuonMass < 80.)) continue;

            for (unsigned int i = 0; i < selectedMuons.size(); ++i) {
                const TString index = TString::Itoa(i, 10);
            
                FillHist(region + "/Lepton_" + index + "_Pt___"  + region,
                         selectedMuons.at(i).Pt(),  weight, 500, 0.0, 5000.0);
                FillHist(region + "/Lepton_" + index + "_Eta___" + region,
                         selectedMuons.at(i).Eta(), weight, 600, -3.0, 3.0);
                FillHist(region + "/Lepton_" + index + "_Phi___" + region,
                         selectedMuons.at(i).Phi(), weight, 640, -3.2, 3.2);
            }
            
            // Dimuon jet = Jet_0
            FillHist(region + "/Jet_0_Pt___"  + region,
                     dimuonJet.Pt(),  weight, 100, 0.0, 500.0);
            FillHist(region + "/Jet_0_Eta___" + region,
                     dimuonJet.Eta(), weight, 60, -3.0, 3.0);
            FillHist(region + "/Jet_0_Phi___" + region,
                     dimuonJet.Phi(), weight, 60, -3.0, 3.0);
            
            // Separate b-tagged tag jet = Jet_1
            if (category == Category::BJet) {
                FillHist(region + "/Jet_1_Pt___"  + region,
                         tagJet.Pt(),  weight, 100, 0.0, 500.0);
                FillHist(region + "/Jet_1_Eta___" + region,
                         tagJet.Eta(), weight, 60, -3.0, 3.0);
                FillHist(region + "/Jet_1_Phi___" + region,
                         tagJet.Phi(), weight, 60, -3.0, 3.0);
            }
        }
    }
}

// =============================================================================
// Dedicated muon-ID and HighPtMuon-trigger efficiency modes
// NIsoMuon efficiency modes patch 20260819_v2
// =============================================================================

float NIsoMuon::CalculateEfficiencyWeight(const Event &ev) {
    if (IsDATA) return 1.0f;

    using Variation = MyCorrection::variation;
    float weight = MCweight();
    weight *= ev.GetTriggerLumi("Full");
    weight *= GetL1PrefireWeight(Variation::nom);
    weight *= myCorr->GetPUWeight(ev.nTrueInt(), Variation::nom);
    return weight;
}

bool NIsoMuon::MatchMuonToTriggerBit(
    const Muon &muon,
    const RVec<TrigObj> &triggerObjects,
    int bit,
    float minimumObjectPt
) const {
    for (const auto &triggerObject : triggerObjects) {
        if (!triggerObject.isMuon()) continue;
        if (!(triggerObject.DeltaR(muon) < 0.3)) continue;
        if (!(triggerObject.Pt() > minimumObjectPt)) continue;
        if (triggerObject.hasBit(bit)) return true;
    }
    return false;
}

bool NIsoMuon::MatchMuonToIsolatedTriggerObject(
    const Muon &muon,
    const RVec<TrigObj> &triggerObjects
) const {
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
        // In 2016 NanoAOD, bit 1 represents IsoMu and bit 3 represents IsoTkMu.
        // The trigger-object pT threshold follows the corresponding HLT path.
        return MatchMuonToTriggerBit(muon, triggerObjects, 1, 24.0f) ||
               MatchMuonToTriggerBit(muon, triggerObjects, 3, 24.0f);
    }
    if (DataEra == "2017") {
        return MatchMuonToTriggerBit(muon, triggerObjects, 3, 27.0f);
    }
    return MatchMuonToTriggerBit(muon, triggerObjects, 3, 24.0f);
}

bool NIsoMuon::MatchMuonToHighPtTriggerObject(
    const Muon &muon,
    const RVec<TrigObj> &triggerObjects
) const {
    // NanoAOD muon TrigObj quality bits: 10 = Mu50, 11 = Mu100.
    // In Run 2, bit 10 also covers the TkMu50 filter configured for 2016.
    // The event-level path OR is checked separately before this matching.
    return MatchMuonToTriggerBit(muon, triggerObjects, 10, 50.0f) ||
           MatchMuonToTriggerBit(muon, triggerObjects, 11, 100.0f);
}

float NIsoMuon::IsolatedTagPtThreshold() const {
    return DataEra == "2017" ? 29.0f : 26.0f;
}

TString NIsoMuon::EfficiencyBinTag(const Muon &probe) const {
    const std::array<double, 11> ptEdges = {
        10.0, 15.0, 20.0, 25.0, 30.0, 40.0,
        50.0, 60.0, 120.0, 200.0, 2000.0
    };
    const std::array<double, 5> etaEdges = {
        0.0, 0.9, 1.2, 2.1, 2.4
    };

    auto label = [](double value) {
        TString result = Form("%.1f", value);
        result.ReplaceAll(".", "p");
        result.ReplaceAll("-", "m");
        return result;
    };

    const double pt = probe.Pt();
    const double absoluteEta = std::fabs(probe.Eta());

    for (std::size_t etaIndex = 0;
         etaIndex + 1 < etaEdges.size();
         ++etaIndex) {
        if (!(etaEdges[etaIndex] <= absoluteEta &&
              absoluteEta < etaEdges[etaIndex + 1])) {
            continue;
        }
        for (std::size_t ptIndex = 0;
             ptIndex + 1 < ptEdges.size();
             ++ptIndex) {
            if (!(ptEdges[ptIndex] <= pt &&
                  pt < ptEdges[ptIndex + 1])) {
                continue;
            }
            return "Pt" + label(ptEdges[ptIndex]) +
                   "to" + label(ptEdges[ptIndex + 1]) +
                   "_AbsEta" + label(etaEdges[etaIndex]) +
                   "to" + label(etaEdges[etaIndex + 1]);
        }
    }
    return "";
}

void NIsoMuon::RunMuonIDEfficiency(
    const Event &ev,
    const RVec<Muon> &allMuons,
    const RVec<Jet> &allJets,
    const RVec<TrigObj> &triggerObjects
) {
    RVec<Muon> denominatorMuons;
    for (const auto &muon : allMuons) {
        if (!(muon.Pt() > 10.0)) continue;
        if (!(std::fabs(muon.Eta()) < kMuonEta)) continue;
        denominatorMuons.push_back(muon);
    }
    std::sort(denominatorMuons.begin(), denominatorMuons.end(), PtComparing);
    if (denominatorMuons.size() < 2) return;

    RVec<Jet> jets = SelectJets(
        allJets,
        Jet::JetID::TIGHT,
        kJetPt,
        kJetEta
    );
    std::sort(jets.begin(), jets.end(), PtComparing);
    if (jets.empty()) return;

    RVec<float> ptEdges = {
        10.0f, 15.0f, 20.0f, 25.0f, 30.0f, 40.0f,
        50.0f, 60.0f, 120.0f, 200.0f, 2000.0f
    };
    RVec<float> etaEdges = {0.0f, 0.9f, 1.2f, 2.1f, 2.4f};

    for (const auto &jet : jets) {
        for (std::size_t tagIndex = 0;
             tagIndex < denominatorMuons.size();
             ++tagIndex) {
            const Muon &tagMuon = denominatorMuons[tagIndex];
            if (!(tagMuon.Pt() > kLeadingMuonPt)) break;
            if (!tagMuon.PassID(Muon::MuonID::POG_MEDIUM)) continue;
            if (!(jet.DeltaR(tagMuon) < kMuonJetDR)) continue;
            if (!MatchMuonToHighPtTriggerObject(tagMuon, triggerObjects)) continue;

            for (std::size_t probeIndex = 0;
                 probeIndex < denominatorMuons.size();
                 ++probeIndex) {
                if (probeIndex == tagIndex) continue;
                const Muon &probe = denominatorMuons[probeIndex];

                if (!probe.isTracker()) continue;
                if (!(jet.DeltaR(probe) < kMuonJetDR)) continue;
                if (!(tagMuon.DeltaR(probe) > 0.05)) continue;
                if (!PassDileptonSign(DileptonSign::OS, tagMuon, probe)) continue;

                const float mass = (tagMuon + probe).M();
                if (!(2.0f < mass && mass < 5.0f)) continue;

                const TString binTag = EfficiencyBinTag(probe);
                if (binTag == "") continue;

                const float weight = CalculateEfficiencyWeight(ev);
                const bool passMediumID =
                    probe.PassID(Muon::MuonID::POG_MEDIUM);

                const TString denominatorRegion = "MuonIDEfficiency_DENOM";
                const TString numeratorRegion = "MuonIDEfficiency_NUM";

                FillHist(
                    denominatorRegion + "/Probe_absEta_Pt___" + denominatorRegion,
                    std::fabs(probe.Eta()),
                    probe.Pt(),
                    weight,
                    etaEdges,
                    ptEdges
                );
                if (passMediumID) {
                    FillHist(
                        numeratorRegion + "/Probe_absEta_Pt___" + numeratorRegion,
                        std::fabs(probe.Eta()),
                        probe.Pt(),
                        weight,
                        etaEdges,
                        ptEdges
                    );
                }

                const TString fitRegion =
                    "MuonIDEfficiency_" + binTag +
                    (passMediumID ? "_Pass" : "_Fail");
                FillHist(
                    fitRegion + "/DileptonJPsi_Mass___" + fitRegion,
                    mass,
                    weight,
                    300,
                    2.0,
                    5.0
                );
                return;
            }
        }
    }
}

void NIsoMuon::RunTriggerEfficiency(
    const Event &ev,
    const RVec<Muon> &allMuons,
    const RVec<Jet> &allJets,
    const RVec<TrigObj> &triggerObjects
) {
    RVec<Muon> muons;
    for (const auto &muon : allMuons) {
        if (!(muon.Pt() > 10.0)) continue;
        if (!(std::fabs(muon.Eta()) < kMuonEta)) continue;
        muons.push_back(muon);
    }
    std::sort(muons.begin(), muons.end(), PtComparing);
    if (muons.size() < 2) return;

    RVec<Jet> jets = SelectJets(
        allJets,
        Jet::JetID::TIGHT,
        kJetPt,
        kJetEta
    );
    std::sort(jets.begin(), jets.end(), PtComparing);
    if (jets.empty()) return;

    std::size_t tagIndex = muons.size();
    for (std::size_t index = 0; index < muons.size(); ++index) {
        const Muon &tagMuon = muons[index];
        if (!(tagMuon.Pt() > IsolatedTagPtThreshold())) continue;
        if (!tagMuon.PassID(Muon::MuonID::POG_TIGHT)) continue;
        if (!tagMuon.PassID(Muon::MuonID::POG_PFISO_TIGHT)) continue;
        if (!MatchMuonToIsolatedTriggerObject(tagMuon, triggerObjects)) continue;

        bool jetClean = true;
        for (const auto &jet : jets) {
            if (jet.DeltaR(tagMuon) < 0.5) {
                jetClean = false;
                break;
            }
        }
        if (!jetClean) continue;

        tagIndex = index;
        break;
    }
    if (tagIndex == muons.size()) return;

    RVec<float> ptEdges = {
        0.0f, 40.0f, 41.0f, 42.0f, 43.0f, 44.0f, 45.0f, 46.0f, 47.0f, 48.0f, 49.0f,
        50.0f, 51.0f, 52.0f, 53.0f, 54.0f, 55.0f, 56.0f, 57.0f, 58.0f, 59.0f, 60.0f, 70.0f,
        80.0f, 90.0f, 100.0f, 120.0f, 150.0f, 200.0f, 500.0f, 2000.0f
    };
    RVec<float> etaEdges = {0.0f, 0.9f, 1.2f, 2.1f, 2.4f};

    for (const auto &jet : jets) {
        for (std::size_t probeIndex = 0;
             probeIndex < muons.size();
             ++probeIndex) {
            if (probeIndex == tagIndex) continue;
            const Muon &probe = muons[probeIndex];

            if (!probe.PassID(Muon::MuonID::POG_MEDIUM)) continue;
            if (!(jet.DeltaR(probe) < kMuonJetDR)) continue;

            const TString binTag = EfficiencyBinTag(probe);
            if (binTag == "") continue;

            const float weight = CalculateEfficiencyWeight(ev);
            const bool passTarget =
                ev.PassTrigger(highPtEfficiencyTriggers) &&
                MatchMuonToHighPtTriggerObject(probe, triggerObjects);

            const TString denominatorRegion = "TriggerEfficiency_DENOM";
            const TString numeratorRegion = "TriggerEfficiency_NUM";

            FillHist(
                denominatorRegion + "/Probe_absEta_Pt___" + denominatorRegion,
                std::fabs(probe.Eta()),
                probe.Pt(),
                weight,
                etaEdges,
                ptEdges
            );
            if (passTarget) {
                FillHist(
                    numeratorRegion + "/Probe_absEta_Pt___" + numeratorRegion,
                    std::fabs(probe.Eta()),
                    probe.Pt(),
                    weight,
                    etaEdges,
                    ptEdges
                );
            }
            return;
        }
    }
}

float NIsoMuon::ReadTheoryWeight(const TString &kind, int index) const {
    if (kind == "PDFScale") {
        // NanoAOD stores the eight non-central scale weights in
        // LHEScaleWeight[0..7]; the central (muR,muF)=(1,1) point is omitted.
        //
        // Expose PDFScale0..8 in the Run-2 SKFlat ordering expected by the
        // downstream limit workflow:
        //   0: (1,1)
        //   1: (1,2)
        //   2: (1,0.5)
        //   3: (2,1)
        //   4: (2,2)
        //   5: (2,0.5)
        //   6: (0.5,1)
        //   7: (0.5,2)
        //   8: (0.5,0.5)
        if (index == 0) return 1.0f;

        if (nLHEScaleWeight < 8) {
            return std::numeric_limits<float>::quiet_NaN();
        }

        switch (index) {
            case 1: return LHEScaleWeight[4];
            case 2: return LHEScaleWeight[3];
            case 3: return LHEScaleWeight[6];
            case 4: return LHEScaleWeight[7];
            case 5: return LHEScaleWeight[5];
            case 6: return LHEScaleWeight[1];
            case 7: return LHEScaleWeight[2];
            case 8: return LHEScaleWeight[0];
            default:
                return std::numeric_limits<float>::quiet_NaN();
        }
    }

    if (kind == "PDFError") {
        int firstErrorIndex = 0;
        int availableErrors = nLHEPdfWeight;
        if (nLHEPdfWeight >= 103) {
            firstErrorIndex = 1;
            availableErrors = nLHEPdfWeight - 3;
        } else if (nLHEPdfWeight >= 102) {
            firstErrorIndex = 0;
            availableErrors = nLHEPdfWeight - 2;
        }

        if (index < 0 || index >= std::min(100, availableErrors)) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return LHEPdfWeight[firstErrorIndex + index];
    }

    if (kind == "PDFAlphaS") {
        int alphaSIndex = -1;
        if (nLHEPdfWeight >= 103) {
            alphaSIndex = 101 + index;
        } else if (nLHEPdfWeight >= 102) {
            alphaSIndex = 100 + index;
        }

        if (index < 0 || index > 1 ||
            alphaSIndex < 0 || alphaSIndex >= nLHEPdfWeight) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return LHEPdfWeight[alphaSIndex];
    }

    return std::numeric_limits<float>::quiet_NaN();
}

void NIsoMuon::executeEvent() {
    const Event ev = GetEvent();
    const RVec<Muon> allMuons = GetAllMuons();
    const RVec<Jet> allJets = GetAllJets();

    // Standard CMS event-quality cleaning.  In Run 3 PassNoiseFilter also
    // contains the dedicated early-Run-3 ECAL bad-calibration recipe.
    if (!PassNoiseFilter(allJets, ev)) return;

    // Apply the JME veto-map prescription before the physics selection.
    // PF-muon-overlap jets (DeltaR < 0.2) are deliberately excluded from the
    // veto-map probe, following the JME/SKNano prescription.  This is important
    // here because a muon-containing jet is an analysis object rather than a
    // jet that should be lepton-cleaned away.
    RVec<Jet> cleanedJets;
    cleanedJets.reserve(allJets.size());

    if (Run == 2) {
        // Run 2: remove only a jet that satisfies the loose veto-map probe
        // selection and lies in a vetoed detector region.
        for (const auto &jet : allJets) {
            bool probeVetoMap = jet.Pt() > 15.0f;
            probeVetoMap = probeVetoMap && jet.PassID(Jet::JetID::TIGHT);
            probeVetoMap = probeVetoMap &&
                (jet.Pt() > 50.0f || jet.PassID(Jet::JetID::PUID_LOOSE));
            probeVetoMap = probeVetoMap &&
                (jet.chEmEF() + jet.neEmEF() < 0.9f);

            for (const auto &muon : allMuons) {
                if (jet.DeltaR(muon) < 0.2f) {
                    probeVetoMap = false;
                    break;
                }
            }

            if (probeVetoMap &&
                myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), "jetvetomap")) {
                continue;
            }
            cleanedJets.push_back(jet);
        }
    } else {
        // Run 3: reject the complete event if any eligible non-muon-overlap
        // jet is in the veto map.  Otherwise retain the original jet collection.
        bool vetoEvent = false;
        for (const auto &jet : allJets) {
            bool probeVetoMap = jet.Pt() > 15.0f;
            probeVetoMap = probeVetoMap && std::fabs(jet.Eta()) < 5.0f;
            probeVetoMap = probeVetoMap && jet.PassID(Jet::JetID::TIGHT);
            probeVetoMap = probeVetoMap &&
                (jet.chEmEF() + jet.neEmEF() < 0.9f);

            for (const auto &muon : allMuons) {
                if (jet.DeltaR(muon) < 0.2f) {
                    probeVetoMap = false;
                    break;
                }
            }

            if (probeVetoMap &&
                myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), "jetvetomap")) {
                vetoEvent = true;
                break;
            }
        }
        if (vetoEvent) return;
        cleanedJets = allJets;
    }

    if (analysisMode == AnalysisMode::MuonIDEfficiency) {
        if (!ev.PassTrigger(highPtEfficiencyTriggers)) return;
        const RVec<TrigObj> triggerObjects = GetAllTrigObjs();
        RunMuonIDEfficiency(ev, allMuons, cleanedJets, triggerObjects);
        return;
    }

    if (analysisMode == AnalysisMode::TriggerEfficiency) {
        if (!ev.PassTrigger(isolatedMuonTriggers)) return;
        const RVec<TrigObj> triggerObjects = GetAllTrigObjs();
        RunTriggerEfficiency(ev, allMuons, cleanedJets, triggerObjects);
        return;
    }

    if (!ev.PassTrigger(highPtMuonTriggers)) return;

    RVec<GenJet> allGenJets;
    if (!IsDATA && runSyst) {
        allGenJets = GetAllGenJets();
    }

    gActiveBTagSource = "central";
    gActiveL1PrefireVariation = 0;

    RunVariation(
        ev,
        allMuons,
        cleanedJets,
        allGenJets,
        ExperimentalSyst::Central,
        ""
    );

    if (!IsDATA && runSyst) {
        // Non-b-tag experimental systematics remain unchanged.
        const std::array<std::pair<ExperimentalSyst, TString>, 12> variations = {{
            {ExperimentalSyst::JetResUp, "_Syst_JetResUp"},
            {ExperimentalSyst::JetResDown, "_Syst_JetResDown"},
            {ExperimentalSyst::JetEnUp, "_Syst_JetEnUp"},
            {ExperimentalSyst::JetEnDown, "_Syst_JetEnDown"},
            {ExperimentalSyst::MuonEnUp, "_Syst_MuonEnUp"},
            {ExperimentalSyst::MuonEnDown, "_Syst_MuonEnDown"},
            {ExperimentalSyst::MuonIDSFUp, "_Syst_MuonIDSFUp"},
            {ExperimentalSyst::MuonIDSFDown, "_Syst_MuonIDSFDown"},
            {ExperimentalSyst::MuonTriggerSFUp, "_Syst_MuonTriggerSFUp"},
            {ExperimentalSyst::MuonTriggerSFDown, "_Syst_MuonTriggerSFDown"},
            {ExperimentalSyst::PUUp, "_Syst_PUUp"},
            {ExperimentalSyst::PUDown, "_Syst_PUDown"}
        }};

        gActiveBTagSource = "central";
        gActiveL1PrefireVariation = 0;
        for (const auto &[variation, suffix] : variations) {
            RunVariation(
                ev,
                allMuons,
                cleanedJets,
                allGenJets,
                variation,
                suffix
            );
        }

        // L1 ECAL prefiring is relevant for 2016 and 2017.  The NanoAOD
        // branches already store the nominal/up/down event weights.  Keep this
        // as a separate era-specific nuisance in the statistical model.
        if (DataEra == "2016preVFP" || DataEra == "2016postVFP" ||
            DataEra == "2017") {
            gActiveL1PrefireVariation = +1;
            RunVariation(
                ev, allMuons, cleanedJets, allGenJets,
                ExperimentalSyst::Central, "_Syst_L1PrefireUp"
            );

            gActiveL1PrefireVariation = -1;
            RunVariation(
                ev, allMuons, cleanedJets, allGenJets,
                ExperimentalSyst::Central, "_Syst_L1PrefireDown"
            );
            gActiveL1PrefireVariation = 0;
        }

        // BTV fixed-WP prescription: retain heavy-/light-flavour and
        // correlated/uncorrelated components as independent templates.
        struct BTagVariationConfig {
            TString source;
            ExperimentalSyst direction;
            TString suffix;
        };

        const std::array<BTagVariationConfig, 8> bTagVariations = {{
            {"hf_corr",   ExperimentalSyst::BTagUp,   "_Syst_BTagHFCorrUp"},
            {"hf_corr",   ExperimentalSyst::BTagDown, "_Syst_BTagHFCorrDown"},
            {"hf_uncorr", ExperimentalSyst::BTagUp,   "_Syst_BTagHFUncorrUp"},
            {"hf_uncorr", ExperimentalSyst::BTagDown, "_Syst_BTagHFUncorrDown"},
            {"lf_corr",   ExperimentalSyst::BTagUp,   "_Syst_BTagLFCorrUp"},
            {"lf_corr",   ExperimentalSyst::BTagDown, "_Syst_BTagLFCorrDown"},
            {"lf_uncorr", ExperimentalSyst::BTagUp,   "_Syst_BTagLFUncorrUp"},
            {"lf_uncorr", ExperimentalSyst::BTagDown, "_Syst_BTagLFUncorrDown"}
        }};

        for (const auto &config : bTagVariations) {
            gActiveBTagSource = config.source;
            RunVariation(
                ev,
                allMuons,
                cleanedJets,
                allGenJets,
                config.direction,
                config.suffix
            );
        }
        gActiveBTagSource = "central";
    }

    if (!IsDATA && runXSecSyst) {
        gActiveBTagSource = "central";
        gActiveL1PrefireVariation = 0;
        for (int index = 0; index < 9; ++index) {
            const float factor = ReadTheoryWeight("PDFScale", index);
            if (!std::isfinite(factor)) continue;
            RunVariation(
                ev,
                allMuons,
                cleanedJets,
                allGenJets,
                ExperimentalSyst::Central,
                "_Syst_PDFScale" + TString::Itoa(index, 10),
                factor
            );
        }

        for (int index = 0; index < 100; ++index) {
            const float factor = ReadTheoryWeight("PDFError", index);
            if (!std::isfinite(factor)) continue;
            RunVariation(
                ev,
                allMuons,
                cleanedJets,
                allGenJets,
                ExperimentalSyst::Central,
                "_Syst_PDFError" + TString::Itoa(index, 10),
                factor
            );
        }

        for (int index = 0; index < 2; ++index) {
            const float factor = ReadTheoryWeight("PDFAlphaS", index);
            if (!std::isfinite(factor)) continue;
            RunVariation(
                ev,
                allMuons,
                cleanedJets,
                allGenJets,
                ExperimentalSyst::Central,
                "_Syst_PDFAlphaS" + TString::Itoa(index, 10),
                factor
            );
        }
    }
}

