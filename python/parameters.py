### Parameters to be used in the configuration files
import FWCore.ParameterSet.Config as cms
import flashgg.Taggers.flashggTags_cff as flashggTags

_benchmark  = cms.untracked.uint32(0)
_is2016			=	cms.untracked.uint32(1)
_doPhotonCR		=	cms.untracked.uint32(1)
_doSelectionTree	=	cms.untracked.uint32(1) 

##Tags and Objects
_triggerTag	=	cms.InputTag("TriggerResults", "", "HLT2")
_myTriggers	=	cms.untracked.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v")
#_DiPhotonTag	=	cms.untracked.InputTag('flashggDiPhotons')
_DiPhotonTag	=	cms.untracked.InputTag('flashggPreselectedDiPhotons')
_JetTag		=	cms.untracked.InputTag('flashggJets')
_metTag		=	cms.InputTag('flashggMets')
#_inputTagJets	=	flashggTags.UnpackedJetCollectionVInputTag,
#_inputTagJets	= 	flashggTags.UnpackedJetCollectionVInputTag
#_inputTagJets 	= 	flashggTags.UnpackedJetCollectionVInputTag
_GenTag		=	cms.untracked.InputTag('flashggGenPhotons')

_rhoFixedGridCollection		=	cms.InputTag('fixedGridRhoAll')

##Photon selection
_PhotonPtOverDiPhotonMass	=	cms.untracked.vdouble( 0.333, 0.25 ) #0: Pho1, 1: Pho2
_PhotonEta			=	cms.untracked.vdouble(1.479, 2.50) #0: EB boundary, 1: EE boundary (EE) (only (1) is used to cut on both)
_PhotonR9			=	cms.untracked.vdouble(0., 0.) #0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
_PhotonElectronVeto		=	cms.untracked.vint32(1, 1) #0: Pho1, 1: Pho2
_PhotonDoElectronVeto		=	cms.untracked.vint32(1 , 1) #0: Pho1, 1: Pho2
_PhotonDoID			=	cms.untracked.vint32(1 , 1) #0: Pho1, 1: Pho2
_PhotonDoISO			=	cms.untracked.vint32(1 , 1) #0: Pho1, 1: Pho2
_DiPhotonPt			=	cms.untracked.vdouble(0.) #0: lower boundary for dipho pt
_DiPhotonEta			=	cms.untracked.vdouble(100000) #0: upper boundary
_DiPhotonMassWindow		=	cms.untracked.vdouble(100., 180.) #0: DiPhoton mass window lower boundary, 1: upper boundary
_DiPhotonOnlyFirst		=	cms.untracked.uint32(0)#If you only want to look at first diphoton pair: 1; if all: 0
_PhotonWhichID			=	cms.untracked.vstring("loose", "loose") #Which photon ID working point? loose, medium or tight
_PhotonWhichISO			=	cms.untracked.vstring("loose", "loose") #Which photon ISO working point? loose, medium or tight
_DoMVAPhotonID			=	cms.untracked.uint32(1) #Do MVA photon ID instead of cut based
_MVAPhotonID			=	cms.untracked.vdouble(0.2, 0.2)#tight: (0.68, 0.60); loose: (0.2, 0.2) (2016 values)
#_PhotonMVAEstimator		=	cms.untracked.string("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values") #MVA user float
_PhotonMVAEstimator		=	cms.untracked.string("EGMPhotonMVA")
_doCustomPhotonMVA		=	cms.untracked.uint32(0)
if _doCustomPhotonMVA:
   _MVAPhotonID = cms.untracked.vdouble(0.07, -0.03)

##Jet selection
_JetPtOverDiJetMass	=	cms.untracked.vdouble(25., 0.0) #0: both jets have to pass [0] requirement. at least one jet has to pass [1] requirement
_JetEta			=	cms.untracked.vdouble(2.4, 2.4) #0: jet1, 1: jet2
_JetBDiscriminant	=	cms.untracked.vdouble(0., 0.8) #0: lowest b-tag requirement for any jet (default 0), standard b-tag cut (loose, medium, tight) 
_JetDoPUID		=	cms.untracked.vint32(1, 1) #0: jet1, 1: jet2
_n_bJets		=	cms.untracked.uint32(0) #Number of required jets passing requirements in JetBDiscriminant
_DiJetPt		=	cms.untracked.vdouble(0.) #0: lower boundary for dijet pt
_DiJetEta		=	cms.untracked.vdouble(20000000.) #0: upper boundary for dijet pt
_DiJetMassWindow	=	cms.untracked.vdouble(60., 180.) #0: DiJet mass window lower boundary, 1: upper boundary
_CandidateMassWindow	=	cms.untracked.vdouble(0., 25000.) #0: 4-candidate mass window lower boundary, 1: upper boundary
_CandidatePt		=	cms.untracked.vdouble(0.) #0 4-candidate pt lower bound
_CandidateEta		=	cms.untracked.vdouble(20000000.) #0 4-candidate eta upper bound
_bTagType		=	cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags') #string for btag algorithm
_JetDrPho		=	cms.untracked.vdouble(0.4, 0.4) #DeltaR between jets and photon
_CandidatesDeltaR	=	cms.untracked.vdouble(0.0) #DeltaR between dijet and diphoton candidates
_JetDoID		=	cms.untracked.vint32(1, 1) #Do jet ID? 0: jet1, 1: jet2
_doJetRegression	=	cms.untracked.uint32(1)
_bRegFileLeading      	=  	cms.untracked.FileInPath("flashgg/bbggTools/data/BRegression/2016/BDTG_15plus3_jetGenJet_nu_leading_summer16_3_24.weights.xml")
_bRegFileSubLeading     =  	cms.untracked.FileInPath("flashgg/bbggTools/data/BRegression/2016/BDTG_15plus3_jetGenJet_nu_trailing_summer16_3_24.weights.xml")


#H/E, Sigma_iEtaiEta, R9, electron veto
_phoIDlooseEB	=	cms.untracked.vdouble(0.05, 0.0102, 1.0)
_phoIDmediumEB	=	cms.untracked.vdouble(0.05, 0.0102, 1.0)
_phoIDtightEB	=	cms.untracked.vdouble(0.05, 0.0100, 1.0)
_phoIDlooseEE	=	cms.untracked.vdouble(0.05, 0.0274, 1.0)
_phoIDmediumEE	=	cms.untracked.vdouble(0.05, 0.0268, 1.0)
_phoIDtightEE	=	cms.untracked.vdouble(0.05, 0.0268, 1.0)

#charged, neutral, photon
_phoISOlooseEB	=	cms.untracked.vdouble(3.32, 1.92, 0.81)
_phoISOmediumEB	=	cms.untracked.vdouble(1.37, 1.06, 0.28)
_phoISOtightEB	=	cms.untracked.vdouble(0.76, 0.97, 0.08)
_phoISOlooseEE	=	cms.untracked.vdouble(1.97, 11.86, 0.83)
_phoISOmediumEE	=	cms.untracked.vdouble(1.10, 2.960, 0.39)
_phoISOtightEE	=	cms.untracked.vdouble(0.56, 2.090, 0.16)

#[0]*pho + [1]
_nhCorrEB	=	cms.untracked.vdouble(0.014, 0.000019)
_phCorrEB	=	cms.untracked.vdouble(0.0053, 0.)
_nhCorrEE	=	cms.untracked.vdouble(0.0139, 0.000025)
_phCorrEE	=	cms.untracked.vdouble(0.0034, 0.)

##Systematics
_jetSmear		=	cms.untracked.int32(0)
_randomLabel		=	cms.untracked.string("rnd_g_JER")
_jetScale		=	cms.untracked.int32(0)
##Photon corrections
_doPhotonScale		=	cms.untracked.int32(0) #-10: not applied | -1: -1sigma | 0: central | 1: 1sigma
_doPhotonExtraScale	=	cms.untracked.int32(0) #-10: not applied | -1: -1sigma | 0: central | 1: 1sigma
_doPhotonSmearing	=	cms.untracked.int32(0) #-10: not applied | -1: -1sigma | 0: central | 1: 1sigma

#_PhotonCorrectionFile	=	cms.untracked.string("EgammaAnalysis/ElectronTools/data/80X_ichepV2_2016_pho")
_PhotonCorrectionFile	=	cms.untracked.string("flashgg/Systematics/data/Moriond17_74x_pho")

##NonResCats
_addNonResMVA 			=	cms.untracked.uint32(1)
_NonResMVAWeights_LowMass 	=	cms.untracked.FileInPath("flashgg/bbggTools/data/NonResMVA/TMVAClassification_BDT.weights_LowMass_MX350_Mjj60.xml")
_NonResMVAWeights_HighMass 	= 	cms.untracked.FileInPath("flashgg/bbggTools/data/NonResMVA/TMVAClassification_BDT.weights_HighMass_MX350_Mjj60.xml")
_ResMVAWeights_LowMass 		=	cms.untracked.FileInPath("flashgg/bbggTools/data/NonResMVA/TMVAClassification_BDT.weights_ResLowMass_MX500_Mjj60.xml")
_ResMVAWeights_HighMass 	= 	cms.untracked.FileInPath("flashgg/bbggTools/data/NonResMVA/TMVAClassification_BDT.weights_ResHighMass_MX500_Mjj60.xml")
_NonResMVAVars 			=	cms.untracked.vstring('leadingJet_bDis','subleadingJet_bDis','diphotonCandidate.Pt()/(diHiggsCandidate.M())','fabs(CosThetaStar_CS)','fabs(CosTheta_bb)','fabs(CosTheta_gg)','dijetCandidate.Pt()/(diHiggsCandidate.M())')

