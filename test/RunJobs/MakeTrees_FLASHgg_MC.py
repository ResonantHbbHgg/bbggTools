import FWCore.ParameterSet.Config as cms
from flashgg.bbggTools.pColors import *
import flashgg.Taggers.flashggTags_cff as flashggTags

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

## Prepare photon ID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
process.egmPhotonIDs.physicsObjectSrc = cms.InputTag('flashggRandomizedPhotons')
process.photonMVAValueMapProducer.src = cms.InputTag('flashggRandomizedPhotons')
process.egmPhotonIsolation.srcToIsolate = cms.InputTag('flashggRandomizedPhotons')
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
#    getattr(process, idmod).physicsObjectSrc = cms.InputTag('flashggPhotons')

print "I'M HERE 3"

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

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

# Load the producer module to build full 5x5 cluster shapes and whatever 
# else is needed for IDs
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
#process.load("Configuration.Geometry.GeometryECALHCAL_cff")
#from RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi import *
# Load the producer for MVA IDs. Make sure it is also added to the sequence!
#from RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi import *
#process.photonMVAValueMapProducer.src = cms.InputTag('flashggRandomizedPhotons')
#process.photonIDValueMapProducer.src = cms.InputTag('flashggRandomizedPhotons')

from flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi import flashggUpdatedIdMVADiPhotons

process.p = cms.Path( flashggUpdatedIdMVADiPhotons
                      * flashggTags.flashggUnpackedJets
                      * process.bbggtree)
