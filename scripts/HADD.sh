SIG='Signal_withTTHBDT'
BKG='BKG_withTTHBDT' 
DATA='Data_withTTHBDT'
TTBBLLNUNU='ttbbllnunu_withTTHBDT'

mkdir hadd_bkg

# ------------- Single Higgs -----------
#ggH
hadd hadd_bkg/output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root $BKG/output_GluGluHToGG_M-125_13TeV_powheg_pythia8_*.root

#VBF
hadd hadd_bkg/output_VBFHToGG_M-125_13TeV_powheg_pythia8.root $BKG/output_VBFHToGG_M-125_13TeV_powheg_pythia8_*.root

# VH
hadd hadd_bkg/output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root $BKG/output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_*.root

#ttH
hadd hadd_bkg/output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root $BKG/output_ttHToGG_M125_13TeV_powheg_pythia8_v2_*.root

#bbH
hadd hadd_bkg/output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo.root $BKG/output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo_*.root
hadd hadd_bkg/output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo.root  $BKG/output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo_*.root
hadd hadd_bkg/output_bbHToGG_M-125_4FS_all_13TeV_amcatnlo.root hadd_BKG/output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo.root hadd_BKG/output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo.root

# ------------- QCD backgrounds ---------
# ggbb
hadd hadd_bkg/output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root  $BKG/output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_*.root

# gjbb
hadd hadd_bkg/output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root $BKG/output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root
hadd hadd_bkg/output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root $BKG/output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root

#jjbb
hadd hadd_bkg/output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root $BKG/output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root
hadd hadd_bkg/output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root $BKG/output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_*.root

# ------------- ttbbllnunu ------------

hadd hadd_bkg/output_TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.root $TTBBLLNUNU/output_TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_*.root

# -------------- DATA -----------------

mkdir hadd_data

hadd hadd_data/DoubleEG_Run2016B.root $DATA/output_*Run2016B*.root

hadd hadd_data/DoubleEG_Run2016C.root $DATA/output_*Run2016C*.root

hadd hadd_data/DoubleEG_Run2016D.root $DATA/output_*Run2016D*.root

hadd hadd_data/DoubleEG_Run2016E.root $DATA/output_*Run2016E*.root

hadd hadd_data/DoubleEG_Run2016F.root $DATA/output_*Run2016F*.root

hadd hadd_data/DoubleEG_Run2016G.root $DATA/output_*Run2016G*.root

hadd hadd_data/DoubleEG_Run2016H.root $DATA/output_*Run2016H*.root

hadd hadd_data/DoubleEG.root hadd_data/DoubleEG_*.root

# -------------- 
NEWPROD='Mar292018_ForUpgrade_ttHBDT'
mkdir /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD
mkdir /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Background  
mkdir /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Background/Hadd  
mkdir /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Data  
mkdir /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Data/Hadd
mkdir /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Signal
mkdir /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Signal/Hadd

cp hadd_data/*.root /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Data/Hadd
#cp $DATA/*.sh /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Data/
#cp $DATA/*.log /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Data/
#cp $DATA/*.py /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Data/

cp hadd_bkg/*.root /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Background/Hadd
#cp $BKG/*.sh /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Background/
#cp $BKG/*.log /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Background/
#cp $BKG/*.py /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Background/

cp $SIG/*.root /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Signal/Hadd
#cp $SIG/*.sh /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Signal/
#cp $SIG/*.log /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Signal/
#cp $SIG/*.py /eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/$NEWPROD/Signal/
