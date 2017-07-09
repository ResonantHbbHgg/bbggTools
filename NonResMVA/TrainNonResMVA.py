import ROOT

# Get the data from the ROOT file
root_bkg = ROOT.TChain('bbggSelectionTree')
root_bkg.AddFile('/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/April26_Mjj70to190/EGML_Data_Mjj70to190/Hadd/DoubleEG.root')
print 'Signal events:',root_bkg.GetEntries()
root_sig = ROOT.TChain('bbggSelectionTree')
#root_sig.AddFile('AllNonResMjj70.root')
root_sig.AddFile('LowMassResMjj70.root')
print 'Background events:',root_sig.GetEntries()

# Useful output information will be stored in a new root file:
f_out = ROOT.TFile("LearningOutput_BDT_LM500Res_GradBoost_Mjj70.root","RECREATE")

# Create the TMVA factory
ROOT.TMVA.Tools.Instance()
factory = ROOT.TMVA.Factory("TMVAClassification", f_out,"AnalysisType=Classification")

# Add the six variables to the TMVA factory as floats
mvaVars = [
##Main: 
'leadingJet_bDis',
'subleadingJet_bDis',
'diphotonCandidate.Pt()/(diHiggsCandidate.M())',
'diphotonCandidate.Pt()/dijetCandidate.Pt()',
#'fabs(dijetCandidate.Eta() - diphotonCandidate.Eta())',
#'diphotonCandidate.Pt()/dijetCandidate.Pt()',
#'fabs(diphotonCandidate.Pt() - dijetCandidate.Pt())/diphotonCandidate.Pt()',
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
sigCut = ROOT.TCut("(isSignal == 1 && (diHiggsCandidate.M() - dijetCandidate.M()-diphotonCandidate.M()+ 250) < 500) && leadingJet_bDis > -1 && subleadingJet_bDis > -1 && diphotonCandidate.Pt()/dijetCandidate.Pt() < 3 && diphotonCandidate.Pt()/(diHiggsCandidate.M()) < 1 && dijetCandidate.Pt()/(diHiggsCandidate.M()) < 1")
bkgCut = ROOT.TCut("(isSignal == 0 && (diHiggsCandidate.M() - dijetCandidate.M()-diphotonCandidate.M()+ 250) < 500) && leadingJet_bDis > -1 && subleadingJet_bDis > -1 && diphotonCandidate.Pt()/dijetCandidate.Pt() < 3 && diphotonCandidate.Pt()/(diHiggsCandidate.M()) < 1 && dijetCandidate.Pt()/(diHiggsCandidate.M()) < 1")
#sigCut = ROOT.TCut("(isSignal == 1 && (diHiggsCandidate.M() - dijetCandidate.M() + 125) > 350) ")
#bkgCut = ROOT.TCut("(isSignal == 0 && (diHiggsCandidate.M() - dijetCandidate.M() + 125) > 350) ")


# Prepare the training/testing signal/background
factory.PrepareTrainingAndTestTree(sigCut,bkgCut,"SplitMode=Random:NormMode=NumEvents:!V")

# Book the SVM method and train/test
#method = factory.BookMethod( ROOT.TMVA.Types.kSVM, "SVM", "C=1.0:Gamma=0.005:Tol=0.001:VarTransform=None" )
method = factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT", "UseRandomisedTrees=1:NTrees=1000:BoostType=Grad:NegWeightTreatment=IgnoreNegWeightsInTraining:MaxDepth=3:MinNodeSize=3:Shrinkage=0.1625")#:nCuts=200")
#method = factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT_2", "UseRandomisedTrees=1:NTrees=1000:BoostType=Grad:NegWeightTreatment=IgnoreNegWeightsInTraining:MaxDepth=3:MinNodeSize=3:Shrinkage=0.1625:nCuts=100")
#method = factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT_3", "UseRandomisedTrees=1:NTrees=1000:BoostType=Grad:NegWeightTreatment=IgnoreNegWeightsInTraining:MaxDepth=3:MinNodeSize=3:Shrinkage=0.1625:SeparationType=CrossEntropy")
#method = factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT", "NTrees=2000:BoostType=AdaBoost:AdaBoostBeta=0.6:UseRandomisedTrees=True:UseNVars=6:nCuts=2000:PruneMethod=CostComplexity:PruneStrength=-1")


factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

