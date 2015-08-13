import FWCore.ParameterSet.Config as cms
#from mAOD_GravFiles import *
#from mAOD_RadFiles import *
from flashgg.bbggTools.microAOD_RadFiles import *

##### Arguments
import flashgg.bbggTools.VarParsing as opts
options = opts.VarParsing('analysis')
options.parseArguments()

maxEvents = -1
if options.maxEvents:
        maxEvents = int(options.maxEvents)

inputFile = RadFiles['700']
if options.inputFiles:
        inputFile = options.inputFiles

outputFile = "rest_rad700.root"
if options.outputFile:
        outputFile = options.outputFile
######

process = cms.Process("bbggtruthdiff")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2000 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( inputFile
#'file:/afs/cern.ch/work/r/rateixei/work/DiHiggs/FLASHggPreSel/CMSSW_7_4_0_pre9/src/flashgg/MicroAOD/test/hhbbgg_hggVtx.root' 
#'/store/user/rateixei/flashgg/RunIISpring15DR74/RunIISpring15MicroAODV1/GluGluToBulkGravitonToHHTo2B2G_M-260_narrow_13TeV-madgraph/RunIISpring15DR74-RunIISpring15MicroAODV1-v0-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150615_144924/0000/microAOD_1.root'
#        RadFiles['320']
#		GravFiles['270']
    )
)

process.load("flashgg.bbggTools.bbggTruthDiff_cfi")
process.bbggtruthdiff.OutFileName = cms.untracked.string(outputFile)


process.p = cms.Path(process.bbggtruthdiff)
