import FWCore.ParameterSet.Config as cms
from microAOD_RadFiles import *
from microAOD_GravFiles import *

process = cms.Process("bbggAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2000 )


## Available mass points:
# RadFiles: 320, 340, 350, 400, 600, 650, 700
# GravFiles: 260, 270, 280, 320, 350, 500, 550

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        XXXFiles['MASS']
    )
)

process.load("flashgg.bbggPlotter.CfiFile_cfi")
process.bbggAnalyzer.OutFileName = cms.untracked.string('XXXMASS.root')


process.p = cms.Path(process.bbggAnalyzer)
