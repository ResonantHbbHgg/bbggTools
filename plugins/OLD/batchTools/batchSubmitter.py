#!/usr/bin/env python

import json, os
from flashgg.bbggTools.pColors import *
from AvailableSamples import *
import getopt, sys

def main(argv):
	eosOutput = '/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/'
	campaign = ''
	sample = ''
	Type = '/bkg/'
	version = ""
	doSelection = 1
	jsonFile = ""
	prefix = ""
	try:
		opts, args = getopt.getopt(argv,"hc:s:o:div:nj:p:",["campaign=", "sample=", "outputDir=", "isData", "isSignal", "version=", "noSelection", "json=", "prefix="])
	except getopt.GetoptError:
		print 'batchSubmitter.py -c <campaign> -s <sample>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'batch_TreeMaker.py -c <campaign> -s <sample>'
			sys.exit()
		elif opt in ("-c", "--campaign"):
			campaign = arg
		elif opt in ("-s", "--sample"):
			sample = arg
		elif opt in ("-o", "--outputDir"):
			eosOutput = arg
		elif opt in ("-d", "--isData"):
			Type = '/data/'
		elif opt in ("-i", "--isSignal"):
			Type = '/sig/'
		elif opt in ("-v", "--version"):
			version = arg
		elif opt in ("-n", "--noSelection"):
			doSelection = 0
		elif opt in ("-j", "--json"):
			jsonFile = arg
		elif opt in ("-p", "--prefix"):
			prefix = arg

	if version == "":
		print "Please, specify version to be created (-v <version>)"
		sys.exit(2)

	eosOutput += version + Type
	
	if sample is '' and "sig" not in Type:
		print 'batchSubmitter.py -c <campaign> -s <sample>'
		sys.exit()
	
	localDir = os.getcwd()
	
	divisor = "myMicroAODOutputFile_"
	
	campaignFile = '../../MetaData/data/' + campaign + '/datasets.json'

	if "sig" in Type:
#		campaignFile = '../MetaData/RunIISpring15-ReMiniAOD-BetaV7/datasets.json'
		campaignFile = '../MetaData/Signal.json'

	data_file_location = campaignFile #localDir + '/../MetaData/microAODdatasets/Spring15BetaV2_MetaV3/datasets.json'

	if jsonFile != "":
		data_file_location = jsonFile
	
	data_file = open(data_file_location)
	data = json.load(data_file)
	
	gotIt = False

	outJson = ''
	if "data" in Type:
		outJson = open(prefix+"data_lumi.json", "w")
		outJson.write("{\n")
	jRuns = {}
	
	totEvents = 0
	for sName in data:
		if sample not in sName.split('/')[1]:
			continue

		gotIt = True
		dataFiles = data[sName]['files']

		print sName
		if "data" in Type:
			for f,dd in enumerate(dataFiles):
				runs = dd['lumis']
				for o,rr in enumerate(runs):
					lumis = runs[rr]
					if rr not in jRuns:
						jRuns[rr] = []
					for i,ll in enumerate(lumis):
						thisRun = '['+str(ll)+','+str(ll)+']'
						jRuns[rr].append(thisRun)
		
		batchInit = '''#!/bin/bash
		cd ''' + str(localDir) + '''
		cd ../test/
		eval `scramv1 runtime -sh`
		'''
		
		fileNameBase = "/tmp/bbggTree_" + prefix + sample
		
		extraCommand = ''
		if 'QCD' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=0 '
		if 'GJet' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=1 '
		if 'DiPhotonJet' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=2 '
		
		
		if(doSelection):
			fileNameBase = "/tmp/bbggSelectionTree_" + prefix + sample
			extraCommand += ' doSelection=1'
		
		if "sig" in Type:
			fileNameBase += sName.split("/")[1]

		for files in dataFiles:
			if int(files['nevents']) == 0: continue
			totEvents += int(files['nevents'])
			print files['name'], files['nevents']
			if "data" in Type:
				thisprefix = prefix + files['name'].split("/")[10] + "_"
				fileNameBase = "/tmp/bbggSelectionTree_" + thisprefix + sample
			fileNumber_temp = files['name'].split(divisor)[1]
			fileNumber = fileNumber_temp.split(".")[0]
			print fileNumber
			outputFile = fileNameBase + '_' + str(fileNumber) + ".root"
			command = "cmsRun MakeTrees.py inputFiles=" + str(files['name']) + " outputFile=" + outputFile + extraCommand
			print bcolors.OKBLUE + "Command to be issued in batch:" + bcolors.ENDC
			print bcolors.OKBLUE + bcolors.BOLD + "\t"+command + bcolors.ENDC
			print bcolors.FAIL + "Final output location: " + bcolors.BOLD + eosOutput + bcolors.ENDC
			batchFinish = "\n\n/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp " + outputFile + " " + eosOutput + "\nrm " + outputFile
			batchName = "TreeMaker_"+str(sample) + '_' +fileNumber+".sh"
			batchFile = open(batchName, "w")
			batchFile.write(batchInit)
			batchFile.write(command)
			batchFile.write(batchFinish)
			batchFile.close()
			os.system("chmod +x "+batchName)
			os.system("bsub -q 1nh -J t"+str(sample) + '_' +fileNumber+" < "+batchName)
			os.system("rm "+batchName)
			print ''
	
	if gotIt is False:
		print "Couldn't find your sample on the available samples. Please use showAvailableSamples.py to see what you can run in."
		sys.exit()
	
	print "########################################################"
	print "Done submitting all jobs to batch! Now wait patiently :)"
	print "Total number of events running on:", str(totEvents)
	print "########################################################"

	if "data" in Type:
		for i,run in enumerate(jRuns):
			if i > 0:
				outJson.write(',\n')
			outJson.write('"'+run+'":[')
			for j,r in enumerate(jRuns[run]):
				if j > 0:
					outJson.write(',')
				outJson.write(r)
			outJson.write(']')
		outJson.write('}')
	
if __name__ == "__main__":
	main(sys.argv[1:])
