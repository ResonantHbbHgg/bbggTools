### Run efficienci tree maker with fggRunJobs

```
fggRunJobs.py --load RunJobs_All_76X.json -H -D -P -n 500 -d <directory> -x cmsRun MakeEfficiencies_RunJob.py maxEvents=-1 --no-use-tarball
```

This will submit jobs to many lxplus machines, but it won't be on lxbatch. If you want to run on lxbatch (recommended), do:   

```
fggRunJobs.py --load RunJobs_All_76X.json -H -D -P -n 500 -d <directory> -x cmsRun MakeEfficiencies_RunJob.py maxEvents=-1 --no-use-tarball
```

After that, do hadd of all your files:   
```
cd <directory>
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
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-650_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-650_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-700_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-700_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-800_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-800_narrow_13TeV-madgraph_*
hadd output_GluGluToBulkGravitonToHHTo2B2G_M-900_narrow_13TeV-madgraph.root ../output_GluGluToBulkGravitonToHHTo2B2G_M-900_narrow_13TeV-madgraph_*

#hadd all non resonant samples
hadd output_GluGluToHHTo2B2G_node_10_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_10_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_11_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_11_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_12_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_12_13TeV-madgraph_*
hadd output_GluGluToHHTo2B2G_node_13_13TeV-madgraph.root ../output_GluGluToHHTo2B2G_node_13_13TeV-madgraph_*
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
hadd output_GluGluToRadionToHHTo2B2G_M-600_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-600_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-700_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-700_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-750_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-750_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-800_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-800_narrow_13TeV-madgraph_*
hadd output_GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph.root ../output_GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph_*

#hadd all background samples
hadd output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root ../output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_*
hadd output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root ../output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*
hadd output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root ../output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_*
hadd output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root ../output_GluGluHToGG_M-125_13TeV_powheg_pythia8_*
hadd output_VBFHToGG_M-125_13TeV_powheg_pythia8.root ../output_VBFHToGG_M-125_13TeV_powheg_pythia8_*
hadd output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root ../output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_*
hadd output_ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root ../output_ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_*

#hadd data in separate datasets
hadd output_DoubleEG_ferrif-RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015D-16Dec2015-v2-5edb826a5c37da64e53f7a0f39c62a72_USER.root ../output_DoubleEG_ferrif-RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015D-16Dec2015-v2-5edb826a5c37da64e53f7a0f39c62a72_USER_*
hadd output_DoubleEG_ferrif-RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015C_25ns-16Dec2015-v1-5edb826a5c37da64e53f7a0f39c62a72_USER.root ../output_DoubleEG_ferrif-RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015C_25ns-16Dec2015-v1-5edb826a5c37da64e53f7a0f39c62a72_USER_*

hadd DoubleEG.root ./output_DoubleEG_ferrif-RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015*

```
