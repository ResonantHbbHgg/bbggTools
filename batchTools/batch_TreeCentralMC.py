import json, os
from flashgg.bbggTools.pColors import *
from AvailableSamples import *
import getopt, sys

def main(argv):
	eosOutput = '/eos/cms/store/user/rateixei/HHbbgg/RunIISpring15-50ns-Spring15BetaV4/'
	campaign = ''
	sample = ''
	Type = 'bkg/'
	try:
		opts, args = getopt.getopt(argv,"hc:s:o:d",["campaign=", "sample=", "outputDir=", "isData"])
	except getopt.GetoptError:
		print 'batch_TreeMaker.py -c <campaign> -s <sample>'
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
			Type = 'data/'

	eosOutput += Type
	
	if campaign is '' or sample is '':
		print 'batch_TreeMaker.py -c <campaign> -s <sample>'
		sys.exit()
	
	localDir = os.getcwd()
	
	divisor = "myMicroAODOutputFile_"
	
	campaignFile = '../../MetaData/data/' + campaign + '/datasets.json'
	data_file_location = campaignFile #localDir + '/../MetaData/microAODdatasets/Spring15BetaV2_MetaV3/datasets.json'
	# if 'crovelli' in bkgSample[1] or 'musella' in bkgSample[1]:
	# 	data_file_location = localDir + '/../MetaData/microAODdatasets/Spring15BetaV2_MetaV3/bkgPasquale.json'
	# if 'musella' in bkgSample[1]:
	# 	divisor = "diphotonsMicroAOD_"
	
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
		eval \`scramv1 runtime -sh\`
		'''
		
		fileNameBase = "/tmp/bbggTree_" + sample
		
		extraCommand = ''
		if 'QCD' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=0 '
		if 'GJet' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=1 '
		if 'DiPhotonJet' in sample:
			extraCommand += ' doDoubleCountingMitigation=1 nPromptPhotons=2 '
		
		
		doSelection = True
		if(doSelection):
			fileNameBase = "/tmp/bbggSelectionTree_" + sample
			extraCommand += ' doSelection=1'
		
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
#			os.system("bsub -q 1nh -J t"+str(sample) + '_' +fileNumber+" < "+batchName)
			print ''
	
	if gotIt is False:
		print "Couldn't find your sample on the available samples. Please use showAvailableSamples.py to see what you can run in."
		sys.exit()
	
	print "########################################################"
	print "Done submitting all jobs to batch! Now wait patiently :)"
	print "########################################################"
	
if __name__ == "__main__":
	main(sys.argv[1:])
