import json, os
from flashgg.bbggTools.pColors import *
import getopt, sys

def main(argv):
	campaign = ''
	try:
		opts, args = getopt.getopt(argv,"hc:",["campaign="])
	except getopt.GetoptError:
		print 'showAvailableDatasets.py -c <campaign>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'showAvailableDatasets.py -c <campaign>'
			sys.exit()
		elif opt in ("-c", "--campaign"):
			campaign = arg

	if campaign is '':
		print "You haven't specified the campaign!"
		print 'showAvailableDatasets.py -c <campaign>'
		return

	campaignFile = '../../MetaData/data/' + campaign + '/datasets.json'
	data_file = open(campaignFile)
	data = json.load(data_file)

	for x in data:
		dsName = x.split('/')[1]
		dsType = ''
		if 'Reco' in x.split('/')[2]:
			dsType = 'DATA'
		else:
			dsType = 'MC'
		print dsType,'\t',dsName
		

if __name__ == "__main__":
	main(sys.argv[1:])