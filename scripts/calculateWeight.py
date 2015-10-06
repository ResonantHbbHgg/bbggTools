import json, os
from flashgg.bbggTools.pColors import *
import getopt, sys

def main(argv):
	campaign = ''
	sample = ''
	try:
		opts, args = getopt.getopt(argv,"hc:s:",["campaign=", "sample="])
	except getopt.GetoptError:
		print 'showAvailableDatasets.py -c <campaign> -s <sample>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'showAvailableDatasets.py -c <campaign> -s <sample>'
			sys.exit()
		elif opt in ("-c", "--campaign"):
			campaign = arg
		elif opt in ("-s", "--sample"):
			sample = arg

	if campaign is '':
		print "You haven't specified the campaign!"
		print 'showAvailableDatasets.py -c <campaign>'
		return

	campaignFile = '../../MetaData/data/' + campaign + '/datasets.json'
	data_file = open(campaignFile)
	data = json.load(data_file)

	sWeights = {}

	for x in data:
		if sample not in x: continue
#		print data[x]
		sWeights[x] = {}
		sWeights[x]['weight'] = 0
		sWeights[x]['events'] = 0
		for w in data[x]['files']:
#			print sWeights[x]['weight']
			sWeights[x]['weight'] += w['weights']
			sWeights[x]['events'] += w['totEvents']

	print sWeights

if __name__ == "__main__":
	main(sys.argv[1:])
