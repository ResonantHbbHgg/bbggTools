import FWCore.ParameterSet.Config as cms
from flashgg.bbggTools.microAOD_RadFiles import *
from flashgg.bbggTools.microAOD_GravFiles import *
from flashgg.bbggTools.More_microAOD_DJet40Inf import *

##### Arguments
import flashgg.bbggTools.VarParsing as opts
options = opts.VarParsing('analysis')
options.parseArguments()

maxEvents = -1
if options.maxEvents:
        maxEvents = int(options.maxEvents)

inputFile = XXXFiles['MASS']
if options.inputFiles:
        inputFile = options.inputFiles

outputFile = '/tmp/bbggTree_XXXMASS.root'
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


process.p = cms.Path(process.bbggtree)
