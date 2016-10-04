### Parameters to be used in the configuration files
import FWCore.ParameterSet.Config as cms
import flashgg.Taggers.flashggTags_cff as flashggTags


_doPhotonCR				=	cms.untracked.uint32(1)
_triggerTag				=	cms.InputTag("TriggerResults", "", "HLT")
_myTriggers				=	cms.untracked.vstring(
						"HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v",
						"HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v",
						"HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v"
													)
_DiPhotonTag				=	cms.untracked.InputTag('flashggDiPhotons')
#_DiPhotonTag				=	cms.untracked.InputTag('flashggPreselectedDiPhotons')
_JetTag						=	cms.untracked.InputTag('flashggJets')
_metTag					=	cms.InputTag('slimmedMETs')

#_inputTagJets					=	flashggTags.UnpackedJetCollectionVInputTag,
_rhoFixedGridCollection		=	cms.InputTag('fixedGridRhoAll')
_GenTag						=	cms.untracked.InputTag('flashggGenPhotons')
#0: Pho1, 1: Pho2
_PhotonPtOverDiPhotonMass	=	cms.untracked.vdouble( 0.333, 0.25 )
#0: First upper boundary (EB), 1: second upper boundary (EE) (only (1) is used to cut on both)
_PhotonEta					=	cms.untracked.vdouble(1.479, 2.50)
#0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
_PhotonR9					=	cms.untracked.vdouble(0., 0.)
#0: Pho1, 1: Pho2
_PhotonElectronVeto			=	cms.untracked.vint32(1, 1)
#0: Pho1, 1: Pho2
_PhotonDoElectronVeto		=	cms.untracked.vint32(1 , 1)
#0: Pho1, 1: Pho2
_PhotonDoID					=	cms.untracked.vint32(1 , 1)
#0: Pho1, 1: Pho2
_PhotonDoISO				=	cms.untracked.vint32(1 , 1)
#0: lower boundary for dipho pt
_DiPhotonPt					=	cms.untracked.vdouble(0.)
#0: upper boundary
_DiPhotonEta				=	cms.untracked.vdouble(100000)
#0: DiPhoton mass window lower boundary, 1: upper boundary
_DiPhotonMassWindow			=	cms.untracked.vdouble(100., 180.)
#If you only want to look at first diphoton pair: 1; if all: 0
_DiPhotonOnlyFirst			=	cms.untracked.uint32(0)

#0: both jets have to pass [0] requirement. at least one jet has to pass [1] requirement
_JetPtOverDiJetMass			=	cms.untracked.vdouble(25., 0.0)

#0: jet1, 1: jet2
_JetEta						=	cms.untracked.vdouble(2.4, 2.4)
#0: lowest b-tag requirement for any jet (default 0), standard b-tag cut (loose, medium, tight) 
_JetBDiscriminant			=	cms.untracked.vdouble(-50., 0.8)
#0: jet1, 1: jet2
_JetDoPUID					=	cms.untracked.vint32(1, 1)
#Number of required jets passing requirements in JetBDiscriminant
_n_bJets					=	cms.untracked.uint32(0)
#0: lower boundary for dijet pt
_DiJetPt					=	cms.untracked.vdouble(0.)
#0: upper boundary for dijet pt
_DiJetEta					=	cms.untracked.vdouble(20.)
#0: DiJet mass window lower boundary, 1: upper boundary
_DiJetMassWindow			=	cms.untracked.vdouble(80., 200.)
#0: 4-candidate mass window lower boundary, 1: upper boundary
_CandidateMassWindow		=	cms.untracked.vdouble(0., 25000.)
#0 4-candidate pt lower bound
_CandidatePt				=	cms.untracked.vdouble(0.)
#0 4-candidate eta upper bound
_CandidateEta				=	cms.untracked.vdouble(20.)
#string for btag algorithm
_bTagType					=	cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags')
#DeltaR between jets and photon
_JetDrPho					=	cms.untracked.vdouble(0.4, 0.4)
#DeltaR between dijet and diphoton candidates
_CandidatesDeltaR			=	cms.untracked.vdouble(0.0)

#_inputTagJets = flashggTags.UnpackedJetCollectionVInputTag

#H/E, Sigma_iEtaiEta, R9, electron veto
_phoIDlooseEB				=	cms.untracked.vdouble(0.05, 0.0102, 1.0)
_phoIDmediumEB				=	cms.untracked.vdouble(0.05, 0.0102, 1.0)
_phoIDtightEB				=	cms.untracked.vdouble(0.05, 0.0100, 1.0)

_phoIDlooseEE				=	cms.untracked.vdouble(0.05, 0.0274, 1.0)
_phoIDmediumEE				=	cms.untracked.vdouble(0.05, 0.0268, 1.0)
_phoIDtightEE				=	cms.untracked.vdouble(0.05, 0.0268, 1.0)

#charged, neutral, photon
_phoISOlooseEB				=	cms.untracked.vdouble(3.32, 1.92, 0.81)
_phoISOmediumEB				=	cms.untracked.vdouble(1.37, 1.06, 0.28)
_phoISOtightEB				=	cms.untracked.vdouble(0.76, 0.97, 0.08)

_phoISOlooseEE				=	cms.untracked.vdouble(1.97, 11.86, 0.83)
_phoISOmediumEE				=	cms.untracked.vdouble(1.10, 2.960, 0.39)
_phoISOtightEE				=	cms.untracked.vdouble(0.56, 2.090, 0.16)

#[0]*pho + [1]
_nhCorrEB					=	cms.untracked.vdouble(0.014, 0.000019)
_phCorrEB					=	cms.untracked.vdouble(0.0053, 0.)
_nhCorrEE					=	cms.untracked.vdouble(0.0139, 0.000025)
_phCorrEE					=	cms.untracked.vdouble(0.0034, 0.)

#if true, make selection tree, if false, make tree before selection
_doSelectionTree			= cms.untracked.uint32(1)

#Which photon ID working point? loose, medium or tight
_PhotonWhichID				=	cms.untracked.vstring("loose", "loose")

#Which photon ISO working point? loose, medium or tight
_PhotonWhichISO				=	cms.untracked.vstring("loose", "loose")

#Do jet ID? 0: jet1, 1: jet2
_JetDoID					=	cms.untracked.vint32(1, 1)

#Do MVA photon ID instead of cut based
_DoMVAPhotonID				=	cms.untracked.uint32(1)
#MVA cut EB/EE, values from Ming: (0.374, 0.472), EGM values: (0.374, 0.336)
_MVAPhotonID				=	cms.untracked.vdouble(0.374, 0.336)
#MVA user float
_PhotonMVAEstimator			=	cms.untracked.string("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values")

_doJetRegression			=	cms.untracked.uint32(1)


_bRegFile      =  cms.untracked.FileInPath("flashgg/bbggTools/data/BRegression/TMVARegression_BDTG.weights.xml")
 
#_bRegFile					=	cms.untracked.string("/afs/cern.ch/work/r/rateixei/work/DiHiggs/flashggJets/CMSSW_7_4_15/src/flashgg/bbggTools/Weights/BRegression/TMVARegression_BDTG.weights.xml")
#_bRegFile					=	 cms.untracked.string("/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/Weights/BRegression/BDTG_16plus2_jetGenJet_nu_7_6.weights.xml")

_jetSmear					=	cms.untracked.int32(0)

_randomLabel				=	cms.untracked.string("rnd_g_JER")

_jetScale					=	cms.untracked.int32(0)

