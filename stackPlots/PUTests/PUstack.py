from pullClass import *
from ROOT import *
import json, os
import shutil

doBlind = False

data_file = open("datasets_NoQCD_weights_Martina.json")
data = json.load(data_file)

version = "S15V7_PhotonCR"

fLocation = "/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/" + version + "/Hadd/"

prefix = "PhoIDloose_"
dirSuffix = "25ns_2600pb_photonCR_pureweight_PU_Martina"
dirPrefix = "/afs/cern.ch/user/r/rateixei/www/HHBBGG/"
dirName = dirPrefix + dirSuffix
isPhoCR = 0

if not os.path.exists(dirName):
	print dirName, "doesn't exist, creating it..."
	os.makedirs(dirName)
	shutil.copy2(dirPrefix + "index.php", dirName+"/index.php")
	if os.path.exists(dirName):
		print dirName, "now exists!"


lumi = 2580 #in pb

datasets = []
signals = []

for bkg in data['background']:
	fLocation = "/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/" + version + "/Hadd/"
#	if "QCD" in bkg['name']:
#		continue
	if "DiPho" in bkg['name']:
		fLocation = ''
	if "DiPho" not in bkg['name']:
		continue
	tempfile = TFile(fLocation+bkg['file'], "READ")
	temptree = tempfile.Get('bbggSelectionTree')
	normalization = (lumi*bkg['xsec']*bkg['sfactor'])/(bkg['weight'])
	arr = [bkg['name'], bkg['legend'], bkg['file'], normalization, bkg['color']]
	dataset = [tempfile, temptree, arr]
	datasets.append(dataset)

'''
for i,bkg in enumerate(data['signal']):
	tempfile = TFile.Open(fLocation+bkg['file'])
	temptree = tempfile.Get('bbggSelectionTree')
	normalization = 800
	arr = ["signal_"+str(i), bkg['legend'], bkg['file'], normalization, bkg['color']]
	dataset = [tempfile, temptree, arr]
	signals.append(dataset)
'''

print data['data']
f = TFile.Open( data['data'] )
t = f.Get('bbggSelectionTree')

bhists = {}


Cut = " diphotonCandidate.M() > 100 && diphotonCandidate.M() < 180"
Cut += " && dijetCandidate.M() > 60 && dijetCandidate.M() < 200"
Cut += " && diHiggsCandidate.M() > 225 && diHiggsCandidate.M() < 350"
#Cut += " && (leadingJet.pt()/dijetCandidate.M()) > 0.3333"
#Cut += " && (subleadingJet.pt()/dijetCandidate.M()) > 0.25"
#Cut += " && leadingPhoton.pt() > 35 && subleadingPhoton.pt() > 35 "
#Cut += " && leadingJet.pt() > 35 && subleadingJet.pt() > 35 "
#Cut += " && leadingJet.Eta() < 2. && leadingJet.Eta() > -2"
#Cut += " && leadingPhotonID[0] == 1 "
#Cut += " && subleadingPhotonID[0] == 1 "
#Cut += " && leadingPhotonEVeto == 0 "
#Cut += " && subleadingPhotonEVeto == 0 "
#Cut += " && leadingPhotonEVeto == 1 "
#Cut += " && subleadingPhotonEVeto == 1 "
#Cut += " && leadingPhotonISO[0] == 1 "
#Cut += " && subleadingPhotonISO[0] == 1 "
#Cut += " && leadingJet.pt() > 45 && subleadingJet.pt() > 45 "
#Cut += " && leadingJet_bDis > 0.9 && subleadingJet_bDis > 0.9"
if doBlind == True:
	Cut += " && !((diphotonCandidate.M() > 120 && diphotonCandidate.M() < 130))"
if isPhoCR == 1:
	Cut += " && (isPhotonCR == 1)"
if isPhoCR == 0:
	Cut += " && (isSignal == 1)"
weight = ""
weight += "( genTotalWeight )*"
#weight += "( genTotalWeight/fabs(genTotalWeight) )*"

#cut = TCut(weight+"("+Cut+")")
cut_data = TCut(Cut)

#variable = "diHiggsCandidate.Pt()"
#varName = "DiHiggs Candidate p_{T} (GeV)"

plots = []

nbin = 50

dr = "sqrt( (leadingPhoton.Eta() - subleadingPhoton.Eta())*(leadingPhoton.Eta() - subleadingPhoton.Eta()) + (leadingPhoton.Phi() - subleadingPhoton.Phi())*(leadingPhoton.Phi() - subleadingPhoton.Phi()) )"

plots.append(["nvtx", "nvtx", "Number of vertices", 52, 0, 52])
plots.append(["mgg", "diphotonCandidate.M()", "diphomass", 50, 100, 180])
plots.append(["mjj", "dijetCandidate.M()", "M(jj) (GeV)", 40, 100, 200])

#variable = "diphotonCandidate.M()"
#varName = "DiPhoton Candidate Mass (GeV)"
for plot in plots:
	print plot
	h1 = TH1F("h1", "Histograms", plot[3], plot[4], plot[5])
	h1.Reset()
	variable = plot[1]
	varName = plot[2]
	stack = myStack('test'+plot[0], varName, varName, dirName, lumi)
	if isPhoCR == 1:
		stack.makePhoCR()

	for bkg in datasets:
		cut = TCut(weight+"("+Cut+")")
		if 'DiPho' in bkg[2][0]:
			cut = TCut("puweight*"+weight+"("+Cut+")")
		print "Adding", bkg[0], 'to stack!'
	#	bfiles[bkg[0]] = TFile(bkg[2])
	#	btrees[bkg[0]] = bfiles[bkg[0]].Get('bbggSelectionTree')
		bhists[bkg[2][0]] = h1.Clone(bkg[2][0])
		bhists[bkg[2][0]].Reset()
		bkg[1].Draw(variable+'>>'+bkg[2][0], cut, 'hist')
		bhists[bkg[2][0]].Scale(bkg[2][3])
		bhists[bkg[2][0]].SetFillColor(bkg[2][4])
#		bhists[bkg[2][0]].SetLineColor(bkg[2][4])
		bhists[bkg[2][0]].SetLineColor(1)
		bhists[bkg[2][0]].SetLineWidth(1)
#		bhists[bkg[2][0]].SetFillStyle(3001)
		stack.addHist(bhists[bkg[2][0]], bkg[2][1], bkg[2][3])

	for bkg in signals:
		print "Adding", bkg[0], 'to stack!'
		bhists[bkg[2][0]] = h1.Clone(bkg[2][0])
		bhists[bkg[2][0]].Reset()
		bkg[1].Draw(variable+'>>'+bkg[2][0], cut, 'hist')
		total = bhists[bkg[2][0]].Integral()
		scale = bkg[2][3]/total
		bhists[bkg[2][0]].Scale(scale)
		bhists[bkg[2][0]].SetLineColor(bkg[2][4])
		bhists[bkg[2][0]].SetLineWidth(2)
#		stack.addSignal(bhists[bkg[2][0]], bkg[2][1], bkg[2][3])

	print 'Adding data to stack!'
	#h1.Reset()
#	f = TFile('data.root')
#	t = f.Get('bbggSelectionTree')
	data = h1.Clone("dat")
	data.Reset()
	t.Draw(variable+'>>dat', cut_data)
	data.SetMarkerStyle(20)
	data.SetMarkerSize(0.8)
	data.Sumw2()
	data.SetMarkerColor(1)
	data.SetLineColor(1)
	data.SetLineWidth(2)
	stack.addData(data, "Data")

	stack.drawStack(prefix + plot[0])
#	stack.drawStack("DCM_" + plot[0])
	del stack

