#ifndef Skim_NIsoMuon_h
#define Skim_NIsoMuon_h

#include "AnalyzerCore.h"

class Skim_NIsoMuon : public AnalyzerCore {
public:
    Skim_NIsoMuon();
    ~Skim_NIsoMuon();

    void initializeAnalyzer() override;
    void executeEvent() override;
    void WriteHist() override;

private:
    TTree *newtree;
    RVec<TString> triggerNames;

    void configureTriggers();
    void addAvailableTrigger(const TString &name);
};

#endif
