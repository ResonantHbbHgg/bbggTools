import FWCore.ParameterSet.Config as cms
from microAOD_RadFiles import *
from microAOD_GravFiles import *
#from microAOD_GJet40Inf import *
from More_microAOD_DJet40Inf import *

process = cms.Process("bbggAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2000 )


## Available mass points:
# RadFiles: 320, 340, 350, 400, 600, 650, 700
# GravFiles: 260, 270, 280, 320, 350, 500, 550
NonResPhys14 = 'file:/afs/cern.ch/work/r/rateixei/work/DiHiggs/FLASHggPreSel/CMSSW_7_4_0_pre9/src/flashgg/MicroAOD/test/hhbbgg_hggVtx.root'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#	GJet40Inf
	RadFiles['320']
#        NonResPhys14
    )
)

process.load("flashgg.bbggPlotter.CfiFile_cfi")
#process.bbggAnalyzer.OutFileName = cms.untracked.string('test_GJet40Inf_INCLUSIVE.root')
process.bbggAnalyzer.OutFileName = cms.untracked.string('test_rad320.root')
#Number of required jets passing requirements in JetBDiscriminant (0: No cat, 1: cat0, 2: cat1)
#process.bbggAnalyzer.n_bJets=cms.untracked.uint32(2),

process.p = cms.Path(process.bbggAnalyzer)
