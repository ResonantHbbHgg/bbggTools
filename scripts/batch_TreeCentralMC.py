import json, os
from flashgg.bbggTools.pColors import *


eosOutput = '/eos/cms/store/user/rateixei/HHbbgg/bbggTrees/bkg/'

GJets20to40 = ['GJets20to40',
'/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/mdonega-RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2-99bc95700bc743151c172a409f08df19/USER'
]

GJets40toInf = ['GJets40toInf',
'/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/mdonega-RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1-99bc95700bc743151c172a409f08df19/USER'
]

QCD30toInfMgg40to80 = ['QCD30toInfMgg40to80',
'/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/mdonega-RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1-99bc95700bc743151c172a409f08df19/USER'
]

QCD30to40 = ['QCD30to40',
'/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/mdonega-RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2-99bc95700bc743151c172a409f08df19/USER'
]
QCD40toInf = ['QCD40toInf',
'/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/mdonega-RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1-99bc95700bc743151c172a409f08df19/USER'
]

p_DYJets = ['DYJets',
'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/musella-EXOSpring15_v1-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2-aca3e57a154f8ee20c72275035ffe5ba/USER'
]
p_DiPhotonJetsSherpa = ['DiPhotonJetsSherpa',
'/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/crovelli-RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v3-99bc95700bc743151c172a409f08df19/USER'
]

bkgSample = GJets20to40

localDir = os.getcwd()

divisor = "myMicroAODOutputFile_"

data_file_location = localDir + '/../MetaData/microAODdatasets/Spring15BetaV2_MetaV3/datasets.json'
if 'crovelli' in bkgSample[1] or 'musella' in bkgSample[1]:
	data_file_location = localDir + '/../MetaData/microAODdatasets/Spring15BetaV2_MetaV3/bkgPasquale.json'
if 'musella' in bkgSample[1]:
	divisor = "diphotonsMicroAOD_"

data_file = open(data_file_location)
data = json.load(data_file)

dataFiles = data[bkgSample[1]]['files']

batchInit = '''#!/bin/bash
cd ''' + str(localDir) + '''
cd ../test/
eval \`scramv1 runtime -sh\`
'''

fileNameBase = "/tmp/bbggTree_" + bkgSample[0]

doSelection = True
if(doSelection):
	fileNameBase = "/tmp/bbggSelectionTree_" + bkgSample[0]
	eosOutput = '/eos/cms/store/user/rateixei/HHbbgg/bbggSelectionTrees/bkg/'


for files in dataFiles:
	if int(files['nevents']) == 0: continue
	print files['name'], files['nevents']
	fileNumber_temp = files['name'].split(divisor)[1]
	fileNumber = fileNumber_temp.split(".")[0]
	print fileNumber
	outputFile = fileNameBase + str(fileNumber) + ".root"
	extraCommand = ''
	if(doSelection):
		extraCommand = ' doSelection=1'
	command = "cmsRun MakeTrees.py inputFiles=" + str(files['name']) + " outputFile=" + outputFile + extraCommand
	print bcolors.OKBLUE + "Command to be issued in batch:" + bcolors.ENDC
	print bcolors.OKBLUE + bcolors.BOLD + "\t"+command + bcolors.ENDC
	print bcolors.FAIL + "Final output location: " + bcolors.BOLD + eosOutput + bcolors.ENDC
	batchFinish = "\n\n/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp " + outputFile + " " + eosOutput + "\nrm " + outputFile
	batchName = "TreeMaker_"+str(bkgSample[0]) + '_' +fileNumber+".sh"
	batchFile = open(batchName, "w")
	batchFile.write(batchInit)
	batchFile.write(command)
	batchFile.write(batchFinish)
	batchFile.close()
	os.system("chmod +x "+batchName)
	os.system("bsub -q 1nh -J t"+str(bkgSample[0]) + '_' +fileNumber+" < "+batchName)
	print ''


print "########################################################"
print "Done submitting all jobs to batch! Now wait patiently :)"
print "########################################################"