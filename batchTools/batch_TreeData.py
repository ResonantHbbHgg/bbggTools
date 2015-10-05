import json, os
from flashgg.bbggTools.pColors import *

localDir = os.getcwd()

batchInit = '''#!/bin/bash
cd ''' + str(localDir) + '''
cd ../test/
eval \`scramv1 runtime -sh\`
'''

data_file = open(localDir + '/../MetaData/microAODdatasets/Spring15BetaV2_MetaV3/datasets.json')
data = json.load(data_file)

eosOutput = '/eos/cms/store/user/rateixei/HHbbgg/bbggTrees/data/'

Run2015B = '/DoubleEG/mdonega-RunIISpring15-50ns-Spring15BetaV2-v0-Run2015B-PromptReco-v1-0433f895d0be63f5c271ad0870ad8023/USER'

dataFiles = data[Run2015B]['files']

fileNameBase = "/tmp/bbggTree_data_"

doSelection = True
if(doSelection):
	fileNameBase = "/tmp/bbggSelectionTree_data_"
	eosOutput = '/eos/cms/store/user/rateixei/HHbbgg/new_bbggSelectionTrees/data/'

for files in dataFiles:
	if int(files['nevents']) == 0: continue
	fileNumber_temp = files['name'].split("myMicroAODOutputFile_")[1]
	fileNumber = fileNumber_temp.split(".")[0]
	outputFile = fileNameBase + str(fileNumber) + ".root"
	extraCommand = ''
	if(doSelection):
		extraCommand = ' doSelection=1'
	command = "cmsRun MakeTrees.py inputFiles=" + str(files['name']) + " outputFile=" + outputFile + extraCommand
	print bcolors.OKBLUE + "Command to be issued in batch:" + bcolors.ENDC
	print bcolors.OKBLUE + bcolors.BOLD + "\t"+command + bcolors.ENDC
	print bcolors.FAIL + "Final output location: " + bcolors.BOLD + eosOutput + bcolors.ENDC
	batchFinish = "\n\n/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp " + outputFile + " " + eosOutput + "\nrm " + outputFile
	batchName = "TreeMaker_data_"+fileNumber+".sh"
	batchFile = open(batchName, "w")
	batchFile.write(batchInit)
	batchFile.write('\n' + command)
	batchFile.write(batchFinish)
	batchFile.close()
	os.system("chmod +x "+batchName)
	os.system("bsub -q 1nh -J tdata_"+fileNumber+" < "+batchName)
	print ''

print "########################################################"
print "Done submitting all jobs to batch! Now wait patiently :)"
print "########################################################"
