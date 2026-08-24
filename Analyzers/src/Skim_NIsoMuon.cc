#include "Skim_NIsoMuon.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

Skim_NIsoMuon::Skim_NIsoMuon() : newtree(nullptr) {}
Skim_NIsoMuon::~Skim_NIsoMuon() {}

void Skim_NIsoMuon::addAvailableTrigger(const TString &name) {
    if (TriggerMap.find(name) != TriggerMap.end()) {
        triggerNames.push_back(name);
    }
}

void Skim_NIsoMuon::configureTriggers() {
    triggerNames.clear();

    RVec<TString> requested;
    if (DataYear == 2016) {
        requested = {
            "HLT_Mu50",
            "HLT_TkMu50",
            "HLT_IsoMu24",
            "HLT_IsoTkMu24"
        };
    } else if (DataYear == 2017) {
        requested = {
            "HLT_Mu50",
            "HLT_OldMu100",
            "HLT_TkMu100",
            "HLT_IsoMu27"
        };
    } else if (DataYear == 2018) {
        requested = {
            "HLT_Mu50",
            "HLT_OldMu100",
            "HLT_TkMu100",
            "HLT_IsoMu24"
        };
    } else if (DataYear == 2022 || DataYear == 2023) {
        requested = {
            "HLT_Mu50",
            "HLT_CascadeMu100",
            "HLT_HighPtTkMu100",
            "HLT_Mu50_L1SingleMuShower",
            "HLT_IsoMu24"
        };
    } else {
        throw std::runtime_error("[Skim_NIsoMuon] unsupported era: " + std::string(DataEra.Data()));
    }

    for (const auto &name : requested) {
        addAvailableTrigger(name);
    }

    if (triggerNames.empty()) {
        throw std::runtime_error("[Skim_NIsoMuon] none of the configured muon triggers is available");
    }

    std::cout << "[Skim_NIsoMuon] active triggers:";
    for (const auto &name : triggerNames) {
        std::cout << " " << name;
    }
    std::cout << std::endl;
}

void Skim_NIsoMuon::initializeAnalyzer() {
    GetOutfile()->cd();
    newtree = fChain->CloneTree(0);
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
    configureTriggers();
}

void Skim_NIsoMuon::executeEvent() {
    Event ev = GetEvent();
    if (!ev.PassTrigger(triggerNames)) {
        return;
    }

    RVec<Muon> muons = GetAllMuons();
    std::sort(muons.begin(), muons.end(), PtComparing);

    RVec<Muon> skimMuons;
    for (const auto &muon : muons) {
        if (muon.Pt() <= 8.) {
            continue;
        }
        if (std::fabs(muon.Eta()) >= 2.5) {
            continue;
        }
        skimMuons.push_back(muon);
    }

    if (skimMuons.size() < 2) {
        return;
    }
    if (skimMuons[0].Pt() <= 24.) {
        return;
    }

    const RVec<Jet> jets = GetAllJets();
    for (const auto &jet : jets) {
        if (jet.Pt() <= 20.) {
            continue;
        }
        if (std::fabs(jet.Eta()) >= 2.6) {
            continue;
        }

        for (const auto &muon : skimMuons) {
            if (jet.DeltaR(muon) < 0.5) {
                newtree->Fill();
                return;
            }
        }
    }
}

void Skim_NIsoMuon::WriteHist() {
    GetOutfile()->cd();
    newtree->Write();
}
