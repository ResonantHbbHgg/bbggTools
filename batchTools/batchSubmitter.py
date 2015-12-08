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
	try:
		opts, args = getopt.getopt(argv,"hc:s:o:div:n",["campaign=", "sample=", "outputDir=", "isData", "isSignal", "version=", "noSelection"])
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
	
	data_file = open(data_file_location)
	data = json.load(data_file)
	
	gotIt = False
	
	for sName in data:
		if sample not in sName.split('/')[1]:
			continue

		gotIt = True
		dataFiles = data[sName]['files']
		
		batchInit = '''#!/bin/bash
		cd ''' + str(localDir) + '''
		cd ../test/
		eval `scramv1 runtime -sh`
		'''
		
		fileNameBase = "/tmp/bbggTree_" + sample
		
		extraCommand = ''
		if 'QCD' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=0 '
		if 'GJet' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=1 '
		if 'DiPhotonJet' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=2 '
		
		
		if(doSelection):
			fileNameBase = "/tmp/bbggSelectionTree_" + sample
			extraCommand += ' doSelection=1'
		
		if "sig" in Type:
			fileNameBase += sName.split("/")[1]

		for files in dataFiles:
			if int(files['nevents']) == 0: continue
			print files['name'], files['nevents']
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
	print "########################################################"
	
if __name__ == "__main__":
	main(sys.argv[1:])
