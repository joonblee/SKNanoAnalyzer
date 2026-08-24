#ifndef NIsoMuon_h
#define NIsoMuon_h

#include "AnalyzerCore.h"
#include "correction.h"

#include <array>
#include <memory>

class NIsoMuon : public AnalyzerCore {
public:
    NIsoMuon();
    ~NIsoMuon();

    void initializeAnalyzer();
    void executeEvent();

private:
    // NIsoMuon efficiency modes patch 20260819_v2
    enum class AnalysisMode {
        NIsoDimuon,
        MuonIDEfficiency,
        TriggerEfficiency
    };

    enum class Category {
        BJet,
        LightJet
    };

    enum class DileptonSign {
        OS,
        SS
    };

    enum class ExperimentalSyst {
        Central,
        JetResUp,
        JetResDown,
        JetEnUp,
        JetEnDown,
        MuonEnUp,
        MuonEnDown,
        MuonIDSFUp,
        MuonIDSFDown,
        MuonTriggerSFUp,
        MuonTriggerSFDown,
        PUUp,
        PUDown,
        BTagUp,
        BTagDown
    };

    AnalysisMode analysisMode;
    bool runSyst;
    bool runXSecSyst;
    RVec<TString> highPtMuonTriggers;
    RVec<TString> triggerSFCoveredPaths;
    RVec<TString> highPtEfficiencyTriggers;
    RVec<TString> isolatedMuonTriggers;
    TString run2TriggerSFKey;
    std::unique_ptr<correction::CorrectionSet> run3MuonZCorrections; //!
    std::unique_ptr<correction::CorrectionSet> run3HighPtMuonCorrections; //!
    float deepJetMediumWP;
    float deepJetLooseWP;

    void ConfigureTriggers();
    void ConfigureMuonScaleFactors();
    void AddTriggerIfAvailable(const TString &name, bool coveredByTriggerSF);

    void RunVariation(
        const Event &ev,
        const RVec<Muon> &allMuons,
        const RVec<Jet> &allJets,
        const RVec<GenJet> &allGenJets,
        ExperimentalSyst syst,
        const TString &systSuffix,
        float theoryWeight = 1.0
    );

    bool FindCandidate(
        Category category,
        DileptonSign sign,
        const RVec<Muon> &muons,
        const RVec<Jet> &allJets,
        const RVec<Jet> &tagJets,
        Muon &leadingMuon,
        Muon &subleadingMuon,
        Jet &dimuonJet,
        Jet &tagJet
    ) const;

    bool PassDileptonSign(
        DileptonSign sign,
        const Muon &muon1,
        const Muon &muon2
    ) const;

    float CalculateWeight(
        const Event &ev,
        Category category,
        ExperimentalSyst syst,
        const RVec<Muon> &selectedMuons,
        const RVec<Jet> &eligibleJets,
        const Jet &tagJet,
        float theoryWeight
    );

    float CalculateMuonIDSF(
        const RVec<Muon> &muons,
        ExperimentalSyst syst
    ) const;

    float CalculateMuonTriggerSF(
        const Event &ev,
        const Muon &leadingMuon,
        ExperimentalSyst syst
    ) const;

    float CalculateBTagWeight(
        const RVec<Jet> &jets,
        JetTagging::JetFlavTaggerWP workingPoint,
        ExperimentalSyst syst
    );

    float CalculateEfficiencyWeight(const Event &ev);

    bool MatchMuonToTriggerBit(
        const Muon &muon,
        const RVec<TrigObj> &triggerObjects,
        int bit,
        float minimumObjectPt
    ) const;

    bool MatchMuonToIsolatedTriggerObject(
        const Muon &muon,
        const RVec<TrigObj> &triggerObjects
    ) const;

    bool MatchMuonToHighPtTriggerObject(
        const Muon &muon,
        const RVec<TrigObj> &triggerObjects
    ) const;

    float IsolatedTagPtThreshold() const;
    TString EfficiencyBinTag(const Muon &probe) const;

    void RunMuonIDEfficiency(
        const Event &ev,
        const RVec<Muon> &allMuons,
        const RVec<Jet> &allJets,
        const RVec<TrigObj> &triggerObjects
    );

    void RunTriggerEfficiency(
        const Event &ev,
        const RVec<Muon> &allMuons,
        const RVec<Jet> &allJets,
        const RVec<TrigObj> &triggerObjects
    );

    float ReadTheoryWeight(const TString &kind, int index) const;

    static TString CategoryName(Category category);
    static TString SignName(DileptonSign sign);
};

#endif
