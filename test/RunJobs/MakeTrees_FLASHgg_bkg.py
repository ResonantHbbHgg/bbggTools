import FWCore.ParameterSet.Config as cms
from flashgg.bbggTools.microAOD_RadFiles import *
from flashgg.bbggTools.microAOD_GravFiles import *
from flashgg.bbggTools.More_microAOD_DJet40Inf import *
from flashgg.bbggTools.pColors import *
import flashgg.Taggers.flashggTags_cff as flashggTags

##### Arguments
#import flashgg.bbggTools.VarParsing as opts
#options = opts.VarParsing('analysis')
#------- Add extra option
#options.register('doSelection',
#					False,
#					opts.VarParsing.multiplicity.singleton,
#					opts.VarParsing.varType.bool,
#					'False: Make tree before selection; True: Make tree after selection')
#options.register('doDoubleCountingMitigation',
#					False,
#					opts.VarParsing.multiplicity.singleton,
#					opts.VarParsing.varType.bool,
#					'False: Do not remove double counted photons from QCD/GJet/DiPhotonJets; True: Remove double counted photons from QCD/GJet/DiPhotonJets')
#options.register('nPromptPhotons',
#					0,
#					opts.VarParsing.multiplicity.singleton,
#					opts.VarParsing.varType.int,
#					'Number of prompt photons to be selected - to use this, set doDoubleCountingMitigation=1')
#
#-------
#
#options.parseArguments()

#maxEvents = 5
#if options.maxEvents:
#        maxEvents = int(options.maxEvents)

#inputFile = "/store/user/rateixei/flashgg/RunIISpring15-50ns/Spring15BetaV2/GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150804_164203/0000/myMicroAODOutputFile_1.root" #RadFiles['650']
#if options.inputFiles:
#        inputFile = options.inputFiles

#outputFile = "rest_rad700.root"
#if options.outputFile:
#        outputFile = options.outputFile
######
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


# call the customization
customize(process)

process.bbggtree.puReWeight=cms.bool( customize.PURW )
if customize.PURW == False:
	process.bbggtree.puTarget = cms.vdouble()
print "I'M HERE 2"

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


#process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2000 )

## Available mass points:
# RadFiles: 320, 340, 350, 400, 600, 650, 700
# GravFiles: 260, 270, 280, 320, 350, 500, 550
#NonResPhys14 = 'file:/afs/cern.ch/work/r/rateixei/work/DiHiggs/FLASHggPreSel/CMSSW_7_4_0_pre9/src/flashgg/MicroAOD/test/hhbbgg_hggVtx.root'

#process.source = cms.Source("PoolSource",
#    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#        customize.inputFiles
#    )
#)

print "I'M HERE 3"

#process.load("flashgg.bbggTools.bbggTree_cfi")
process.load("flashgg.Taggers.flashggTags_cff")
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

#process.p = cms.Path(process.bbggtree)

process.p = cms.Path(flashggTags.flashggUnpackedJets*process.bbggtree)
