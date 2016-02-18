```
mkdir NonResSignal
fggRunJobs.py --load RunJobs_Signal.json -H -D -P -n 500 -d NonResSignal -x cmsRun MakeTrees_FLASHgg.py  maxEvents=-1 -q 1nh --no-use-tarball
```   

Hadd all backgrounds:   
```
mkdir Hadd
hadd Hadd/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd Hadd/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd Hadd/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd Hadd/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd Hadd/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_*
hadd Hadd/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_*
```
