import json, os

data_file = open('datasets.json')
data = json.load(data_file)

eosOutput = '/eos/cms/store/user/rateixei/HHbbgg/bbggTrees/data/'

Run2015B = '/DoubleEG/mdonega-RunIISpring15-50ns-Spring15BetaV2-v0-Run2015B-PromptReco-v1-0433f895d0be63f5c271ad0870ad8023/USER'

dataFiles = data[Run2015B]['files']

localDir = os.getcwd()

batchInit = '''#!/bin/bash
cd ''' + str(localDir) + '''
cd ../test/
eval \`scramv1 runtime -sh\`
'''

for files in dataFiles:
	if int(files['nevents']) == 0: continue
	print files['name'], files['nevents']
	fileNumber_temp = files['name'].split("myMicroAODOutputFile_")[1]
	fileNumber = fileNumber_temp.split(".")[0]
	print fileNumber
	if int(fileNumber) == 10: continue
	outputFile = "/tmp/bbggTree_data_" + str(fileNumber) + ".root"
	command = "\ncmsRun MakeTrees.py inputFiles=" + str(files['name']) + " outputFile=" + outputFile
	batchFinish = "\n\n/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp " + outputFile + " " + eosOutput + "\nrm " + outputFile
	batchName = "TreeMaker_data_"+fileNumber+".sh"
	batchFile = open(batchName, "w")
	batchFile.write(batchInit)
	batchFile.write(command)
	batchFile.write(batchFinish)
	batchFile.close()
	os.system("chmod +x "+batchName)
	os.system("bsub -q 1nh -J tdata_"+fileNumber+" < "+batchName)
#	break
