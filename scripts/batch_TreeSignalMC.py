import json, os

localDir = os.getcwd()

batchInit = '''#!/bin/bash
cd ''' + str(localDir) + '''
cd ../test/
eval \`scramv1 runtime -sh\`
'''

data_file = open( localDir + '/../MetaData/microAODdatasets/Spring15BetaV2_MetaV3/signalMC.json')
data = json.load(data_file)

eosOutput = '/eos/cms/store/user/rateixei/HHbbgg/bbggSelectionTrees/signal/'

massPoints = ['Grav260', 'Grav270', 'Grav280', 'Grav320', 'Grav350',
				'Grav500', 'Grav550', 'Rad320', 'Rad340', 'Rad350', 
				'Rad400', 'Rad600', 'Rad650', 'Rad700']
				
for mPoint in massPoints:
	print mPoint
	dataFiles = data[mPoint]['files']
	rFiles = ''
	counter = 0
	for files in dataFiles:
		counter += 1
		if int(files['nevents']) == 0: continue
		rFiles += files['name']
		if counter < len(dataFiles): rFiles += ','
	print rFiles
	outputFile = "/tmp/bbggSelectionTree_" + str(mPoint) + ".root"
	command = "\ncmsRun MakeTrees.py inputFiles=" + str(rFiles) + " outputFile=" + outputFile + " doSelection=1"
	batchFinish = "\n\n/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp " + outputFile + " " + eosOutput + "\nrm " + outputFile
	batchName = "TreeMaker_"+str(mPoint)+".sh"
	batchFile = open(batchName, "w")
	batchFile.write(batchInit)
	batchFile.write(command)
	batchFile.write(batchFinish)
	batchFile.close()
	os.system("chmod +x "+batchName)
	os.system("bsub -q 1nh -J t"+str(mPoint)+" < "+batchName)
#	break
