from pullClass import *
from ROOT import *
from colors import *
import json, os
import shutil
from configs import *

doBlind = True
doSignalRegion = True
doJetCR = False
isPhoCR = False
addHiggs = True
hideData = False
addbbH = False

data_file = open("datasets/datasets76X.json")
data = json.load(data_file)

bkgsToAdd = [['vbf', "VBF H(#gamma#gamma)", 1], ['vh', "VH(#gamma#gamma)", 1], ['ggf', "ggH(#gamma#gamma)", 1], ['tth', "ttH(#gamma#gamma)", 1], ["DYJet", "Z/#gamma*+Jets", 1], ["GJets", "#gamma+Jets", 2], ["DiPhoJet", "#gamma#gamma+Jets", 1], ["QCD", "QCD", 2] ]

version = "correctcuts_71"

fLocation = "/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/" + version + "/Hadd/"

prefix = ""
dirSuffix = "1_3_0_SR_PhoMVAID_76X"
dirPrefix = "/afs/cern.ch/user/r/rateixei/www/HHBBGG/"
dirName = dirPrefix + dirSuffix

#yieldsFile = open(dirPrefix+dirSuffix+"/yields.txt", "w")

if not os.path.exists(dirName):
	print dirName, "doesn't exist, creating it..."
	os.makedirs(dirName)
	shutil.copy2(dirPrefix + "index.php", dirName+"/index.php")
	if os.path.exists(dirName):
		print dirName, "now exists!"

yieldsFile = open(dirPrefix+dirSuffix+"/yields.txt", "w")

lumi = 2700 #in pb

datasets = []
signals = []

if addHiggs:
	for bkg in data['higgs']:
		fLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_hig/Hadd/"
		if addbbH == False and "bbH" in bkg['name']:
			continue
		tempfile = TFile(fLocation+bkg['file'], "READ")
		temptree = tempfile.Get('bbggSelectionTree')
		normalization = (lumi*bkg['xsec']*bkg['sfactor'])/(bkg['weight'])
		arr = [bkg['name'], bkg['legend'], bkg['file'], normalization, bkg['color']]
		dataset = [tempfile, temptree, arr]
		datasets.append(dataset)


for bkg in data['background']:
	fLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_bkg/Hadd/"
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
	fLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_sig/Hadd/"
	tempfile = TFile.Open(fLocation+bkg['file'])
	temptree = tempfile.Get('bbggSelectionTree')
	normalization = 50
	arr = ["signal_"+str(i), bkg['legend'], bkg['file'], normalization, bkg['color']]
	dataset = [tempfile, temptree, arr]
	signals.append(dataset)


print data['data']
fLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_data/Hadd/"
f = TFile.Open( fLocation+data['data'] )
t = f.Get('bbggSelectionTree')

bhists = {}


Cut = " diphotonCandidate.M() > 100 && diphotonCandidate.M() < 180"
Cut += " && dijetCandidate.M() > 60 && dijetCandidate.M() < 180"
#Cut += " && diHiggsCandidate.M() > 280 && diHiggsCandidate.M() < 320"
#Cut += " && (leadingJet.pt()/dijetCandidate.M()) > 0.3333"
#Cut += " && (subleadingJet.pt()/dijetCandidate.M()) > 0.25"
#Cut += " && leadingPhoton.pt() > 35 && subleadingPhoton.pt() > 35 "
#Cut += " && leadingJet.pt() > 35 && subleadingJet.pt() > 35 "
#Cut += " && leadingJet.Eta() < 2. && leadingJet.Eta() > -2"
#Cut += " && leadingPhotonID[1] == 1 "
#Cut += " && subleadingPhotonID[1] == 1 "
#Cut += " && leadingPhotonISO[1] == 1 "
#Cut += " && subleadingPhotonISO[1] == 1 "
#Cut += " && leadingPhotonEVeto == 0 "
#Cut += " && subleadingPhotonEVeto == 0 "
#Cut += " && leadingPhotonEVeto == 1 "
#Cut += " && subleadingPhotonEVeto == 1 "
#Cut += " && leadingPhotonISO[0] == 1 "
#Cut += " && subleadingPhotonISO[0] == 1 "
#Cut += " && leadingJet.pt() > 45 && subleadingJet.pt() > 45 "
#Cut += " && leadingJet_bDis > 0.9 && subleadingJet_bDis > 0.9"
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
#	Cut += " && (( leadingJet_bDis > 0.8 && subleadingJet_bDis < 0.8 ) || (leadingJet_bDis < 0.8 && subleadingJet_bDis > 0.8) ) && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) > 290 && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) < 310 "
#	Cut += " && (( leadingJet_bDis < 0.8 && subleadingJet_bDis > 0.8 )) && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) > 290 && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) < 310 "
#	CutSignal += " && (( leadingJet_bDis > 0.8 && subleadingJet_bDis < 0.8 ) || (leadingJet_bDis < 0.8 && subleadingJet_bDis > 0.8) ) && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) > 290 && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) < 310"
#	CutSignal += " && (( leadingJet_bDis < 0.8 && subleadingJet_bDis > 0.8 )) && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) > 290 && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) < 310"
	CutSignal += " && ( leadingJet_bDis > 0.8 || subleadingJet_bDis > 0.8 ) "
	Cut += " && ( leadingJet_bDis > 0.8 || subleadingJet_bDis > 0.8 ) "
if doJetCR == True:
	Cut += " && leadingJet_bDis < 0.80 && subleadingJet_bDis < 0.8 "
weight = ""
weight += "( genTotalWeight )*"
#weight += "( genTotalWeight/fabs(genTotalWeight) )*"

#cut = TCut(weight+"("+Cut+")")
cut_data = TCut(Cut)

#variable = "diHiggsCandidate.Pt()"
#varName = "DiHiggs Candidate p_{T} (GeV)"

plots = []

yieldsFile.write("Extra cuts on selection tree:\n")
yieldsFile.write(Cut+"\nYields:\n")

nbin = 40

dr = "sqrt( (leadingPhoton.Eta() - subleadingPhoton.Eta())*(leadingPhoton.Eta() - subleadingPhoton.Eta()) + (leadingPhoton.Phi() - subleadingPhoton.Phi())*(leadingPhoton.Phi() - subleadingPhoton.Phi()) )"

plots.append(["j1ratio_dijet", "leadingJet.Pt()/dijetCandidate.M()", "p_{T}(j_{1})/M(jj)", nbin, 0.1, 1.5])
plots.append(["dijet_deta", "fabs(leadingJet.Eta() - subleadingJet.Eta())", "#Delta#eta between jets", nbin, 0, 5])
plots.append(["diPho_Mass", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 100, 180])
#plots.append(["costhetastar", "fabs(CosThetaStar)", "|cos#theta*|", nbin, 0, 1])
plots.append(["leadingPhoton_pt", "leadingPhoton.pt()", "Leading Photon Pt (GeV)", nbin, 30, 300])
plots.append(["nvtx", "nvtx", "Number of vertices", 30, 0, 30])
plots.append(["diPho_Mass", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 100, 180])
plots.append(["diPho_Mass_HM", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 80, 2000])
plots.append(["diJet_Mass", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 60, 180])
plots.append(["diJet_Mass_Limit", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 60, 180])
plots.append(["diJet_Mass_HM", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 80, 2000])
plots.append(["leadingJet_pt", "leadingJet.pt()", "Leading Jet Pt (GeV)", nbin, 15, 200] )
plots.append(["subleadingJet_pt", "subleadingJet.pt()", "Subleading Jet Pt (GeV)", nbin, 15, 80] )
plots.append(["leadingJet_eta", "leadingJet.Eta()", "Leading Jet Eta", nbin, -3, 3] )
plots.append(["leadingPhoton_eta", "leadingPhoton.Eta()", "Leading Photon Eta", nbin, -3, 3] )
plots.append(["dicandidate_Mass", "diHiggsCandidate.M()", "DiHiggs Candidate Mass (GeV)", nbin, 100, 1000])
plots.append(["MX", "diHiggsCandidate.M() - dijetCandidate.M() + 125.", "M_{X} (GeV)", nbin, 100, 1000])
plots.append(["MKF", "diHiggsCandidate_KF.M()", "M_{KinFit}(bb#gamma#gamma) (GeV)", nbin, 100, 1000])
plots.append(["dicandidate_Mass_Limit", "diHiggsCandidate.M()", "DiHiggs Candidate Mass (GeV)", nbin, 225, 350])
plots.append(["dicandidate_Mass_HM", "diHiggsCandidate.M()", "DiHiggs Candidate Mass (GeV)", nbin, 250, 5000])
plots.append(["btagSum", "leadingJet_bDis+subleadingJet_bDis", "Sum of b-tag of jet pair", nbin, 0, 2])
plots.append(["subleadingPhoton_pt", "subleadingPhoton.pt()", "Subleading Photon Pt (GeV)", nbin, 30, 150])
plots.append(["subleadingPhoton_eta", "subleadingPhoton.eta()", "Subleading Photon Eta)", nbin, -3, 3])
plots.append(["dr_photons", dr, "#DeltaR between photons", nbin, 0, 10])
plots.append(["leadingPho_MVA", "leadingPhotonIDMVA", "Leading Photon #gammaMVA discriminant", nbin, 0, 1])
plots.append(["subleadingPho_MVA", "subleadingPhotonIDMVA", "SubLeading Photon #gammaMVA discriminant", nbin, 0, 1])


#variable = "diphotonCandidate.M()"
#varName = "DiPhoton Candidate Mass (GeV)"
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

#			if(i == 0):
#				Yield = bhists[bkg[2][0]].Integral()
#				err = sqrt(nEvtsBeforeScale)*bkg[2][3]
#				yieldsFile.write(bkg[2][1] + "\t&\t" + str(Yield)+ "\t \pm \t" + str(err) + "\n")


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

yieldsFile.close()
