# SKNanoAnalyzer
**v1.0.1**

## Introduction
- For setting up the environment and starting the analysis, please refer [Getting Started](docs/GettingStarted.md).
- If you want to engage in the development and commit to the repository, please read [Development Guide](docs/DevelopmentGuide.md) first.

## To do
### Assigned
- [ ] FatJet and Tau: **Youngwan**
- [ ] Auto-Validation: Check the consistency of the output for a specific NanoAOD file: **Jin**

## Make Tutorial
- [ ] From make Training tree, Export ONNX, Use this model in the Analyzer : **Yeonjoon(Review: Youngwan)**
- [ ] Simple DY tutorial : **Eunsu(Review: Jin)**
- [ ] tt reconstruction using Kinematic Fit : **Minseok(Review: Byunghun)**
- [ ] tt reconstruction using KinFit Library : **Yeonjoon**
- [ ] NLO vs LO comparison in MadGraph Tutorial: **Jin**
- [ ] Tau+Fatjet Tutorial : **Youngwan**

## Useful links
- [MiniAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD)
- [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD)
- [NanoAOD-Tools (for postproc)](https://github.com/cms-sw/cmssw/tree/master/PhysicsTools/NanoAODTools)
- [SNU\_SKFlatProduction](https://docs.google.com/spreadsheets/d/12SfHxUR_1_-13Tzt6mI0slnenQeeCT40zDdggx_KPR8/edit?pli=1#gid=1739893496)
- [NanoAOD wiki@gitlab](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/home)
- [Run2 Summer20 UL](https://docs.google.com/spreadsheets/d/1ABl2p2uwr2EfEbolBEVNcKb_fIXigYY9sqCRT8XIi1Q/edit#gid=1318927481)
- [Run3 Common Bkg Sample Preparation](https://docs.google.com/spreadsheets/d/1xEbHtxzJVaWJpxDrJ7tSpxcVQMRIFbZeoeO7davNFMA/edit#gid=1396344675)
- [SKFlatMaker](https://github.com/CMSSNU/SKFlatMaker)
- [SKFlatAnalyzer](https://github.com/CMSSNU/SKFlatAnalyzer)
- [SKFlat Manual](https://jskim.web.cern.ch/jskim/SKFlat/Manual/Manual_SKFlat.pdf)
- [NanoAOD Branches](https://cms-nanoaod-integration.web.cern.ch/autoDoc/)
- [PdmV Run3 DATA full table](https://pdmv-pages.web.cern.ch/run_3_data/full_table.html)
- [CMS 13.6TeV SM cross sections](https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV)
- [LHC ttbar NNLO+NNLL cross sections](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO)
- [LHC Single Top NNLO+NNLL cross sections](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec)
- [LHC Higgs Cross Section](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/HiggsXSBR)
- [Rules for Run3 2024 dataset names](https://cms-pdmv.gitbook.io/project/mccontact/rules-for-run3-2024-dataset-names)
- [Rules for Run3 dataset names](https://cms-pdmv.gitbook.io/project/mccontact/rules-for-run3-dataset-names)


## NIsoMuon submission options

After setting up the environment, run from the repository root. `source run.sh`
defaults to **Run3, nominal, no top-pT reweighting**. `bash run.sh` accepts the
same options. Sourcing keeps script variables/functions and exits in a subshell.

```bash
source run.sh --help
source run.sh                              # Run3 nominal
source run.sh --era 2022EE --flags nominal,RunSyst
source run.sh --era Run3 --flags all --top-pt-reweight
source run.sh --era 2022EE --flags all --top-pt-reweight --dry-run
```

- `--era` (alias `--eras`) accepts a single era or a comma-separated list.
- `--flags nominal` is the default; `all` submits nominal, `RunSyst`, and
  `RunXSecSyst` independently. Existing efficiency modes remain selectable by name.
- `--top-pt-reweight` (optionally `true`/`false`) selects only POWHEG ttbar samples
  for the chosen modes. It does not submit data, other backgrounds, or signals.
- `--dry-run` prints the selected submission commands without building or submitting.

`TOP_PT_REWEIGHT=true bash run.sh` remains supported and now selects the default
**nominal** mode. Add `--flags all` to produce nominal and both systematic
collections. `--no-top-pt-reweight` overrides the environment setting.
Study outputs are kept below `NIsoMuon/TopPtReweight/<era>/`,
`NIsoMuon/TopPtReweight/RunSyst/<era>/`, and
`NIsoMuon/TopPtReweight/RunXSecSyst/<era>/`. Ordinary production paths are unchanged.

SKNano also accepts `Run2`, `Run3`, and `Run2+3` (alias `all`).
