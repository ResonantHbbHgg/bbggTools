import ROOT

# Get the data from the ROOT file
root_bkg = ROOT.TChain('bbggSelectionTree')
root_bkg.AddFile('/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar22_Moriond_Spring16Reg_v10/EGML_Data_v10/Hadd/DoubleEG.root')
print 'Signal events:',root_bkg.GetEntries()
root_sig = ROOT.TChain('bbggSelectionTree')
root_sig.AddFile('AllNonRes.root')
#root_sig.AddFile('AllRes.root')
#root_sig.AddFile('res_lowmass.root')
print 'Background events:',root_sig.GetEntries()

# Useful output information will be stored in a new root file:
f_out = ROOT.TFile("LearningOutput_BDT_LM400_GradBoost.root","RECREATE")

# Create the TMVA factory
ROOT.TMVA.Tools.Instance()
factory = ROOT.TMVA.Factory("TMVAClassification", f_out,"AnalysisType=Classification")

# Add the six variables to the TMVA factory as floats
mvaVars = [
##Main: 
'leadingJet_bDis',
'subleadingJet_bDis',
'diphotonCandidate.Pt()/(diHiggsCandidate.M())',
##Remove:
#          'leadingJet.Pt()/dijetCandidate.M()',
#	   'subleadingJet.Pt()/dijetCandidate.M()',
##Test:
          'fabs(CosThetaStar_CS)',
           'fabs(CosTheta_bb)',
           'fabs(CosTheta_gg)',
           'dijetCandidate.Pt()/(diHiggsCandidate.M())'
           ]

for x in mvaVars:
    factory.AddVariable(x,"F")

# Link the signal and background to the root_data ntuple
factory.AddSignalTree(root_sig)
factory.AddBackgroundTree(root_bkg)

# cuts defining the signal and background sample
sigCut = ROOT.TCut("(isSignal == 1 && (diHiggsCandidate.M() - dijetCandidate.M()-diphotonCandidate.M()+ 250) < 400) && leadingJet_bDis > -1 && subleadingJet_bDis > -1")
bkgCut = ROOT.TCut("(isSignal == 0 && (diHiggsCandidate.M() - dijetCandidate.M()-diphotonCandidate.M()+ 250) < 400) && leadingJet_bDis > -1 && subleadingJet_bDis > -1")
#sigCut = ROOT.TCut("(isSignal == 1 && (diHiggsCandidate.M() - dijetCandidate.M() + 125) > 350) ")
#bkgCut = ROOT.TCut("(isSignal == 0 && (diHiggsCandidate.M() - dijetCandidate.M() + 125) > 350) ")


# Prepare the training/testing signal/background
factory.PrepareTrainingAndTestTree(sigCut,bkgCut,"SplitMode=Random:NormMode=NumEvents:!V")

# Book the SVM method and train/test
#method = factory.BookMethod( ROOT.TMVA.Types.kSVM, "SVM", "C=1.0:Gamma=0.005:Tol=0.001:VarTransform=None" )
method = factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT", "UseRandomisedTrees=1:NTrees=1000:BoostType=Grad:NegWeightTreatment=IgnoreNegWeightsInTraining:MaxDepth=3:MinNodeSize=3:Shrinkage=0.1625")
#method = factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT", "NTrees=2000:BoostType=AdaBoost:AdaBoostBeta=0.6:UseRandomisedTrees=True:UseNVars=6:nCuts=2000:PruneMethod=CostComplexity:PruneStrength=-1")
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

