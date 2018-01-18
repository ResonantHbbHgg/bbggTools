import FWCore.ParameterSet.Config as cms
from flashgg.bbggTools.pColors import *

process = cms.Process("bbggtree")
process.load("flashgg.bbggTools.bbggTree_cfi")
process.bbggtree.rho = cms.InputTag('fixedGridRhoAll')
process.bbggtree.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.bbggtree.puInfo=cms.InputTag("slimmedAddPileupInfo")
process.bbggtree.lumiWeight = cms.double(1.0)
process.bbggtree.intLumi = cms.double(1.0)
process.bbggtree.puReWeight=cms.bool(True)
process.bbggtree.puBins=cms.vdouble()
process.bbggtree.dataPu=cms.vdouble()
process.bbggtree.mcPu=cms.vdouble()
print "I'M HERE 1"

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
customize.register('PURW',
				1,
				VarParsing.VarParsing.multiplicity.singleton,
				VarParsing.VarParsing.varType.bool,
				"Do PU reweighting? Doesn't work on 76X")

customize.register('bench',
				0,
				VarParsing.VarParsing.multiplicity.singleton,
				VarParsing.VarParsing.varType.int,
				"Benchmark number for Non-Res weight. SHould be in range 1-12")



# call the customization
customize(process)

process.bbggtree.puReWeight=cms.bool( bool(customize.PURW) )
if customize.PURW == False:
	process.bbggtree.puTarget = cms.vdouble()
print "I'M HERE 2.0"


process.bbggtree.benchmark=cms.untracked.uint32(customize.bench)
if customize.bench>0:
  print "I'M HERE 2.1"
  process.bbggtree.getNonResGenInfo=cms.untracked.bool(True)

print "I'M HERE 2.2"

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

process.bbggtree.OutFileName = cms.untracked.string(outputFile)

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

process.bbggtree.doSigmaMdecorr = 1
process.bbggtree.sigmaMdecorrFile = cms.untracked.FileInPath("flashgg/Taggers/data/diphoMVA_sigmaMoMdecorr_split_Mgg40_180.root")
process.bbggtree.MVAResultTag = cms.InputTag('flashggDiPhotonMVA')

import flashgg.Taggers.flashggTags_cff as flashggTags
process.load("flashgg.Taggers.flashggTags_cff")

import flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi as flashggPhotonMVA
flashggPhotonMVA.flashggUpdatedIdMVADiPhotons.reRunRegression = cms.bool(False)
process.load("flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi")

import flashgg.Taggers.flashggPreselectedDiPhotons_cfi as flashggPreSelection
process.load("flashgg.Taggers.flashggPreselectedDiPhotons_cfi")

import flashgg.Taggers.flashggDiPhotonMVA_cfi as flashggDiPhotonMVA
process.load("flashgg.Taggers.flashggDiPhotonMVA_cfi")

process.p = cms.Path( flashggPhotonMVA.flashggUpdatedIdMVADiPhotons
                      * flashggPreSelection.flashggPreselectedDiPhotons
		      * flashggDiPhotonMVA.flashggDiPhotonMVA
                      * flashggTags.flashggUnpackedJets
                      * process.bbggtree )

#process.p = cms.Path( flashggPhotonMVA.flashggUpdatedIdMVADiPhotons
#                      * flashggTags.flashggUnpackedJets
#                      * process.bbggtree)
