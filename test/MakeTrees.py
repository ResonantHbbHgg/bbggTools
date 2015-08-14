import FWCore.ParameterSet.Config as cms
from flashgg.bbggTools.microAOD_RadFiles import *
from flashgg.bbggTools.microAOD_GravFiles import *
from flashgg.bbggTools.More_microAOD_DJet40Inf import *
from flashgg.bbggTools.pColors import *

##### Arguments
import flashgg.bbggTools.VarParsing as opts
options = opts.VarParsing('analysis')
#------- Add extra option
options.register('doSelection',
					False,
					opts.VarParsing.multiplicity.singleton,
					opts.VarParsing.varType.bool,
					'False: Make tree before selection; True: Make tree after selection')
#-------

options.parseArguments()

maxEvents = -1
if options.maxEvents:
        maxEvents = int(options.maxEvents)

inputFile = "/store/user/rateixei/flashgg/RunIISpring15-50ns/Spring15BetaV2/GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150804_164203/0000/myMicroAODOutputFile_1.root" #RadFiles['650']
if options.inputFiles:
        inputFile = options.inputFiles

outputFile = "rest_rad700.root"
if options.outputFile:
        outputFile = options.outputFile
######

process = cms.Process("bbggtree")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2000 )

## Available mass points:
# RadFiles: 320, 340, 350, 400, 600, 650, 700
# GravFiles: 260, 270, 280, 320, 350, 500, 550
NonResPhys14 = 'file:/afs/cern.ch/work/r/rateixei/work/DiHiggs/FLASHggPreSel/CMSSW_7_4_0_pre9/src/flashgg/MicroAOD/test/hhbbgg_hggVtx.root'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        inputFile
    )
)

process.load("flashgg.bbggTools.bbggTree_cfi")
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
print bcolors.FAIL + "Is it performing the analysis selection?",options.doSelection, bcolors.ENDC
if options.doSelection is True:
	print "HELOOOO"
	process.bbggtree.doSelectionTree = cms.untracked.uint32(1)
if options.doSelection is False:
	process.bbggtree.doSelectionTree = cms.untracked.uint32(0)

process.p = cms.Path(process.bbggtree)
