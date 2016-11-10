#!/usr/bin/env bash
rm OutputAngles.root
python readLorentzVector.py  eosuser/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/DoubleEG.root DoubleEG
python readLorentzVector.py  eosuser/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root SMHH
python readLorentzVector.py  eosuser/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph.root Radion300
python readLorentzVector.py  eosuser/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToRadionToHHTo2B2G_M-700_narrow_13TeV-madgraph.root Radion700
python readLorentzVector.py  eosuser/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToBulkGravitonToHHTo2B2G_M-300_narrow_13TeV-madgraph.root Graviton300
python readLorentzVector.py  eosuser/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToBulkGravitonToHHTo2B2G_M-700_narrow_13TeV-madgraph.root Graviton700
