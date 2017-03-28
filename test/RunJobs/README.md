## Single Local Job Run
Before submitting batch jobs, try individual signal, background and data jobs.   
Examples:   

```
cmsRun MakeTrees_FLASHgg_MC.py maxEvents=-1 \
targetLumi=22.37e+3 doSelection=1 PURW=True \
puTarget=5.49e+04,4.58e+05,1.55e+06,2.72e+06,3.7e+06,5.61e+06,1.62e+07,3.89e+07,7.41e+07,1.64e+08,3.18e+08,5.21e+08,7.33e+08,9.76e+08,1.24e+09,1.43e+09,1.55e+09,1.59e+09,1.59e+09,1.56e+09,1.5e+09,1.42e+09,1.32e+09,1.19e+09,1.05e+09,9.08e+08,7.63e+08,6.27e+08,5.01e+08,3.87e+08,2.89e+08,2.08e+08,1.45e+08,9.66e+07,6.19e+07,3.8e+07,2.22e+07,1.23e+07,6.51e+06,3.26e+06,1.55e+06,7.01e+05,3.03e+05,1.27e+05,5.27e+04,2.34e+04,1.23e+04,8.3e+03,6.87e+03,6.29e+03 \
inputFiles=/store/group/phys_higgs/resonant_HH/RunII/MicroAOD/HHbbggSignal_Moriond17_v1/2_4_4/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/HHbbggSignal_Moriond17_v1-2_4_4-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170118_080928/0000/myMicroAODOutputFile_3.root \
outputFile=output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root doDoubleCountingMitigation=0

cmsRun MakeTrees_FLASHgg_MC.py maxEvents=100 campaign=RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2 targetLumi=22.37e+3 doSelection=1 PURW=True\ puTarget=5.49e+04,4.58e+05,1.55e+06,2.72e+06,3.7e+06,5.61e+06,1.62e+07,3.89e+07,7.41e+07,1.64e+08,3.18e+08,5.21e+08,7.33e+08,9.76e+08,1.24e+09,1.43e+09,1.55e+09,1.59e+09,1.59e+09,1.56e+09,1.5e+09,1.42e+09,1.32e+09,1.19e+09,1.05e+09,9.08e+08,7.63e+08,6.27e+08,5.01e+08,3.87e+08,2.89e+08,2.08e+08,1.45e+08,9.66e+07,6.19e+07,3.8e+07,2.22e+07,1.23e+07,6.51e+06,3.26e+06,1.55e+06,7.01e+05,3.03e+05,1.27e+05,5.27e+04,2.34e+04,1.23e+04,8.3e+03,6.87e+03,6.29e+03 \ inputFiles=/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/160707_151052/0000/myMicroAODOutputFile_1.root \ outputFile=output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root doDoubleCountingMitigation=0

cmsRun MakeTrees_FLASHgg_data.py maxEvents=5000 campaign=RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2 \
targetLumi=36.50e+3 doSelection=1 processId=Data processType=Data \  
inputFiles=/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/DoubleEG/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-Run2016B-23Sep2016-v2/161114_162452/0000/myMicroAODOutputFile_102.root \    
lumiMask=/afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/lumiJsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt \
outputFile=output_TestData.root
```

## Submission with flashgg tools

To submit a test job locally, try this:
```
fggRunJobs.py --load jsons/80X/RunJobs_Test.json -m 0 -n 1 -d NonResSignal -x cmsRun MakeTrees_FLASHgg_Test.py  maxEvents=-1 --no-use-tarball
```
You may find that it does not work because the correct ```dataset.json``` file does not exist. For that to work you need to create such dataset first, following instructions on the level above. Then, modify the ```json/RunJobs_Test.json``` file substituting the campaign name. Now it wshould work.

```
fggRunJobs.py --load <RunJobs_Signal.json> -H -D -P -n 500 -d <NonResSignal> -x cmsRun MakeTrees_FLASHgg_MC.py  maxEvents=-1 -q 1nh --no-use-tarball
fggRunJobs.py --load <RunJobs_Background.json> -H -D -P -n 500 -d <NonResSignal> -x cmsRun MakeTrees_FLASHgg_MC.py  maxEvents=-1 -q 1nh --no-use-tarball
fggRunJobs.py --load <RunJobs_Data.json> -H -D -P -n 500 -d <NonResSignal> -x cmsRun MakeTrees_FLASHgg_data.py  maxEvents=-1 -q 1nh --no-use-tarball
```

### Hadd all backgrounds:
```
mkdir Hadd
cd Hadd

#hadd all graviton samples
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-1000_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-1000_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-250_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-250_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-260_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-260_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-270_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-270_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-280_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-280_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-300_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-300_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-320_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-320_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-340_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-340_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-350_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-350_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-400_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-400_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-450_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-450_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-500_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-500_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-550_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-550_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-600_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-600_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-650_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-650_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-700_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-700_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-750_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-750_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-800_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-800_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-900_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-900_narrow_13TeV-madgraph_*

#hadd all non resonant samples
hadd output_GluGluToHHTo2B2G_node_10_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_10_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_11_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_11_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_12_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_12_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_13_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_13_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_2_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_2_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_3_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_3_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_4_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_4_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_5_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_5_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_6_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_6_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_7_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_7_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_8_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_8_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_9_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_9_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_box_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_box_13TeV-madgraph_*

#hadd all radion samples
hadd output_GluGluToRadionToHHTo2B2G_M-250_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-250_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-260_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-260_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-270_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-270_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-280_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-280_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-320_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-320_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-340_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-340_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-350_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-350_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-400_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-400_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-450_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-450_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-500_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-500_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-550_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-550_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-600_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-600_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-700_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-700_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-750_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-750_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-800_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-800_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph_*

#hadd all background samples
hadd output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root ../output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_*.root
hadd output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root
hadd output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root
hadd output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root
hadd output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root
hadd DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root ../output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ferrif-RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1-94d03017ce60e1991044b8807b85a209_USER_*.root
hadd output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root ../output_GluGluHToGG_M-125_13TeV_powheg_pythia8_*.root
hadd output_VBFHToGG_M-125_13TeV_powheg_pythia8.root ../output_VBFHToGG_M-125_13TeV_powheg_pythia8_*.root
hadd output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root ../output_ttHToGG_M125_13TeV_powheg_pythia8_v2_*.root
hadd output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root ../output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_*.root
hadd output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo.root ../output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo_*.root
hadd output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo.root ../output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo_*.root

#hadd data in separate datasets
hadd DoubleEG_Run2016B.root ../output_DoubleEG_*Run2016B*_USER_*.root
hadd DoubleEG_Run2016C.root ../output_DoubleEG_*Run2016C*_USER_*.root
hadd DoubleEG_Run2016D.root ../output_DoubleEG_*Run2016D*_USER_*.root
hadd DoubleEG_Run2016E.root ../output_DoubleEG_*Run2016E*_USER_*.root
hadd DoubleEG_Run2016F.root ../output_DoubleEG_*Run2016F*_USER_*.root
hadd DoubleEG_Run2016G.root ../output_DoubleEG_*Run2016G*_USER_*.root
hadd DoubleEG_Run2016H.root ../output_DoubleEG_*Run2016H*_USER_*.root
hadd DoubleEG.root DoubleEG_Run2016*

```
