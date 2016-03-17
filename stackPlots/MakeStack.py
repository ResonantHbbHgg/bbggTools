from pullClass import *
from ROOT import *
from colors import *
import json, os
import shutil
from configs import *

data = json.load(data_file)

if not os.path.exists(dirName):
	print dirName, "doesn't exist, creating it..."
	os.makedirs(dirName)
	shutil.copy2(dirPrefix + "index.php", dirName+"/index.php")
	if os.path.exists(dirName):
		print dirName, "now exists!"

yieldsFile = open(dirPrefix+dirSuffix+"/yields.txt", "w")

datasets = []
signals = []

if addHiggs:
	for bkg in data['higgs']:
		fLocation = higgsLocation
		if addbbH == False and "bbH" in bkg['name']:
			continue
		tempfile = TFile(fLocation+bkg['file'], "READ")
		temptree = tempfile.Get('bbggSelectionTree')
		normalization = (lumi*bkg['xsec']*bkg['sfactor'])/(bkg['weight'])
		arr = [bkg['name'], bkg['legend'], bkg['file'], normalization, bkg['color']]
		dataset = [tempfile, temptree, arr]
		datasets.append(dataset)


for bkg in data['background']:
	fLocation = bkgLocation
#	if "QCD" in bkg['name']:
#		continue
#	if "DiPho" in bkg['name']:
#		fLocation = ''
#	if "DiPho" not in bkg['name']:
#		continue
	tempfile = TFile(fLocation+bkg['file'], "READ")
	temptree = tempfile.Get('bbggSelectionTree')
	normalization = (lumi*bkg['xsec']*bkg['sfactor'])/(bkg['weight'])
	arr = [bkg['name'], bkg['legend'], bkg['file'], normalization, bkg['color']]

	dataset = [tempfile, temptree, arr]
	datasets.append(dataset)


for i,bkg in enumerate(data['signal']):
	fLocation = signalLocation
	tempfile = TFile.Open(fLocation+bkg['file'])
	temptree = tempfile.Get('bbggSelectionTree')
	normalization = 50
	arr = ["signal_"+str(i), bkg['legend'], bkg['file'], normalization, bkg['color']]
	dataset = [tempfile, temptree, arr]
	signals.append(dataset)


print data['data']
fLocation = dataLocation
f = TFile.Open( fLocation+data['data'] )
t = f.Get('bbggSelectionTree')

bhists = {}

CutSignal = Cut
if doBlind == True:
	Cut += " && !((diphotonCandidate.M() > 120 && diphotonCandidate.M() < 130))"
if isPhoCR == True:
	Cut += " && (isPhotonCR == 1)"
	CutSignal += " && (isPhotonCR == 1)"
if isPhoCR == False:
	Cut += " && (isSignal == 1)"
	CutSignal += " && (isSignal == 1)"
if doSignalRegion == True:
	CutSignal += " && ( leadingJet_bDis > 0.8 || subleadingJet_bDis > 0.8 ) "
	Cut += " && ( leadingJet_bDis > 0.8 || subleadingJet_bDis > 0.8 ) "
if doJetCR == True:
	Cut += " && leadingJet_bDis < 0.80 && subleadingJet_bDis < 0.8 "
weight = ""
weight += "( genTotalWeight )*"
cut_data = TCut(Cut)

yieldsFile.write("Extra cuts on selection tree:\n")
yieldsFile.write(Cut+"\nYields:\n")

for i,plot in enumerate(plots):
	print plot
	h1 = TH1F("h1"+str(i), "Histograms", plot[3], plot[4], plot[5])
	h1.Reset()
	variable = plot[1]
	varName = plot[2]
	stack = myStack('test'+plot[0], varName, varName, dirName, lumi)
	if hideData == True:
		stack.hideData()
	if isPhoCR == 1:
		stack.makePhoCR()
	if doJetCR == 1:
		stack.makeJetCR()

	#addingbbh
	
	if(addbbH):
		bbhist = 0
		bbN2 = 0
		bbN3 = 0
		for bkg in datasets:
			cut = TCut(weight+"("+Cut+")")
			if 'bbH' not in bkg[2][0]:
				continue
			print "Adding", bkg[0], 'to stack!'
		#	bfiles[bkg[0]] = TFile(bkg[2])
		#	btrees[bkg[0]] = bfiles[bkg[0]].Get('bbggSelectionTree')
			bhists[bkg[2][0]] = h1.Clone(bkg[2][0])
			bhists[bkg[2][0]].Reset()
			bkg[1].Draw(variable+'>>'+bkg[2][0], cut, 'hist')
			bhists[bkg[2][0]].Sumw2()
			bhists[bkg[2][0]].Scale(bkg[2][3])
			bhists[bkg[2][0]].SetFillColor(bkg[2][4])
	#		bhists[bkg[2][0]].SetLineColor(bkg[2][4])
			bhists[bkg[2][0]].SetLineColor(1)
			print bkg[2][1]
			bhists[bkg[2][0]].SetLineWidth(0)
			if(bbhist == 0):
				bbhist = bhists[bkg[2][0]].Clone()
				bbN2 = bkg[2][1]
				bbN3 = bkg[2][3]
			else:
				bbhist.Sumw2()
				bhists[bkg[2][0]].Sumw2()
				bbhist.Add(bhists[bkg[2][0]])
	#		bhists[bkg[2][0]].SetFillStyle(3001)
		stack.addHist( bbhist, bbN2, bbN3 )


	allbkg = []
	for BK in bkgsToAdd:
		nInThis = 0
		for bkg in datasets:
			if BK[0] not in bkg[2][0]:
				continue
			if BK[0] in bkg[2][0]:
				nInThis += 1

			cut = TCut(weight+"("+Cut+")")
			if 'bbH' in bkg[2][0]:
				continue
			print "Adding", bkg[0], 'to stack!'
			if(nInThis == 1):
				bhists[ BK[0] ] = h1.Clone(BK[0])
				bhists[ BK[0] ].Reset()
				bhists[ BK[0] ].Sumw2()
				bhists[ BK[0] ].SetFillColor(bkg[2][4])
				bhists[ BK[0] ].SetLineColor(1)
				bhists[ BK[0] ].SetLineWidth(0)
			tempHist = bhists[BK[0]].Clone(BK[0] + "_Temp")
			tempHist.Reset()
			tempHist.Sumw2()

			bkg[1].Draw(variable+'>>'+BK[0] + "_Temp", cut, 'hist')
			tempHist.Sumw2()
			tempHist.Scale(bkg[2][3])

			bhists[ BK[0] ].Add(tempHist)

			if nInThis == BK[2]:
				stack.addHist(bhists[ BK[0] ], bkg[2][1], bkg[2][3])

	for bkg in signals:
		print "Adding", bkg[0], 'to stack!'
		bhists[bkg[2][0]] = h1.Clone(bkg[2][0])
		bhists[bkg[2][0]].Reset()
		bkg[1].Draw(variable+'>>'+bkg[2][0], TCut(CutSignal), 'hist')
		total = bhists[bkg[2][0]].Integral()
		nEvtsBeforeScale = bhists[bkg[2][0]].Integral()
		scale = 1.
		if bkg[2][3] > 0 and total > 0:
			scale = bkg[2][3]/total
		bhists[bkg[2][0]].Scale(scale)
		bhists[bkg[2][0]].SetLineColor(bkg[2][4])
		bhists[bkg[2][0]].SetLineWidth(2)
		stack.addSignal(bhists[bkg[2][0]], bkg[2][1], bkg[2][3])
		if(i == 0):
			Yield = nEvtsBeforeScale #bhists[bkg[2][0]].Integral()
			err = sqrt(nEvtsBeforeScale)
			yieldsFile.write(bkg[2][1] + "\t&\t" + str(Yield)+ "\t \pm \t" + str(err) + "\n")

	print 'Adding data to stack!'
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
	del stack

yieldsFile.close()
