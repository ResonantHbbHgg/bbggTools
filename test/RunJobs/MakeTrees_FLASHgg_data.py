import FWCore.ParameterSet.Config as cms
from flashgg.bbggTools.pColors import *
import flashgg.Taggers.flashggTags_cff as flashggTags

process = cms.Process("bbggtree")
process.load("flashgg.bbggTools.bbggTree_cfi")
process.bbggtree.rho = cms.InputTag('fixedGridRhoAll')
process.bbggtree.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.bbggtree.puInfo=cms.InputTag("slimmedAddPileupInfo")
process.bbggtree.lumiWeight = cms.double(1.0)
process.bbggtree.intLumi = cms.double(1.0)
process.bbggtree.puReWeight=cms.bool(False)
#process.bbggtree.puBins=cms.vdouble()
#process.bbggtree.dataPu=cms.vdouble()
#process.bbggtree.mcPu=cms.vdouble()

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring("test.root")
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2000 )

# import flashgg customization
from flashgg.MetaData.JobConfig import customize
import FWCore.ParameterSet.VarParsing as VarParsing
# set default options if needed
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",2.58e+3)

customize.register('doSelection',
					False,
					VarParsing.VarParsing.multiplicity.singleton,
					VarParsing.VarParsing.varType.bool,
					'False: Make tree before selection; True: Make tree after selection')
customize.register('doDoubleCountingMitigation',
					False,
					VarParsing.VarParsing.multiplicity.singleton,
					VarParsing.VarParsing.varType.bool,
					'False: Do not remove double counted photons from QCD/GJet/DiPhotonJets; True: Remove double counted photons from QCD/GJet/DiPhotonJets')
customize.register('nPromptPhotons',
					0,
					VarParsing.VarParsing.multiplicity.singleton,
					VarParsing.VarParsing.varType.int,
					'Number of prompt photons to be selected - to use this, set doDoubleCountingMitigation=1')

# call the customization
customize(process)

maxEvents = 5
if customize.maxEvents:
        maxEvents = int(customize.maxEvents)

inputFile = "/store/user/rateixei/flashgg/RunIISpring15-50ns/Spring15BetaV2/GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150804_164203/0000/myMicroAODOutputFile_1.root" #RadFiles['650']
if customize.inputFiles:
        inputFile = customize.inputFiles

outputFile = "rest_rad700.root"
if customize.outputFile:
        outputFile = customize.outputFile

print customize.inputFiles, customize.outputFile, customize.maxEvents, customize.doSelection, customize.doDoubleCountingMitigation, customize.nPromptPhotons

#process.load("flashgg.bbggTools.bbggTree_cfi")
process.load("flashgg.Taggers.flashggTags_cff")
process.bbggtree.OutFileName = cms.untracked.string(outputFile)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
#### 2015 triggers
#process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v*",
#								"HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v*",
#								"HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v*") )
#### 2016 trigger
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
								) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#############   Geometry  ###############
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.load("Configuration.Geometry.GeometryECALHCAL_cff")
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)

process.dataRequirements = cms.Sequence()

process.dataRequirements += process.hltHighLevel
process.dataRequirements += process.eeBadScFilter


print bcolors.OKBLUE + "########################################################################" + bcolors.ENDC
print bcolors.OKBLUE + "#  _   _  _   _  _      _                     _____                    #" + bcolors.ENDC
print bcolors.OKBLUE + "# | | | || | | || |    | |                   |_   _|                   #" + bcolors.ENDC
print bcolors.OKBLUE + "# | |_| || |_| || |__  | |__    __ _   __ _    | |   _ __   ___   ___  #" + bcolors.ENDC
print bcolors.OKBLUE + "# |  _  ||  _  || '_ \ | '_ \  / _` | / _` |   | |  | '__| / _ \ / _ \ #" + bcolors.ENDC
print bcolors.OKBLUE + "# | | | || | | || |_) || |_) || (_| || (_| |   | |  | |   |  __/|  __/ #" + bcolors.ENDC
print bcolors.OKBLUE + "# \_| |_/\_| |_/|_.__/ |_.__/  \__, | \__, |   \_/  |_|    \___| \___| #" + bcolors.ENDC
print bcolors.OKBLUE + "#                               __/ |  __/ |                           #" + bcolors.ENDC
print bcolors.OKBLUE + "#                              |___/  |___/                            #" + bcolors.ENDC
print bcolors.OKBLUE + "########################################################################" + bcolors.ENDC
print bcolors.FAIL + "Is it performing the analysis selection?",customize.doSelection, bcolors.ENDC
print bcolors.FAIL + "Is it removing double counted photons from QCD/GJet/DiPhotonJets?",customize.doDoubleCountingMitigation, bcolors.ENDC
if customize.doSelection is True:
	process.bbggtree.doSelectionTree = cms.untracked.uint32(1)
if customize.doSelection is False:
	process.bbggtree.doSelectionTree = cms.untracked.uint32(0)
if customize.doDoubleCountingMitigation is True:
	process.bbggtree.doDoubleCountingMitigation = cms.untracked.uint32(1)
	process.bbggtree.nPromptPhotons = cms.untracked.uint32(customize.nPromptPhotons)
	print "Number of prompt photons in DiPhotonCandidate:",customize.nPromptPhotons
if customize.doDoubleCountingMitigation is False:
	process.bbggtree.doDoubleCountingMitigation = cms.untracked.uint32(0)

import flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi as flashggPhotonMVA
#flashggPhotonMVA.flashggUpdatedIdMVADiPhotons.reRunRegression = cms.bool(False)
process.load("flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi")

import flashgg.Taggers.flashggPreselectedDiPhotons_cfi as flashggPreSelection
process.load("flashgg.Taggers.flashggPreselectedDiPhotons_cfi")

process.p = cms.Path( process.dataRequirements
                      * flashggPhotonMVA.flashggUpdatedIdMVADiPhotons
                      * flashggPreSelection.flashggPreselectedDiPhotons
                      * flashggTags.flashggUnpackedJets
                      * process.bbggtree )

#process.p = cms.Path( process.dataRequirements
#                      * flashggPhotonMVA.flashggUpdatedIdMVADiPhotons
#                      * flashggTags.flashggUnpackedJets
#                      * process.bbggtree)
