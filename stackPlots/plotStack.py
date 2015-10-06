from pullClass import *
from ROOT import *
import json, os

data_file = open("datasets.json")
data = json.load(data_file)

'''
bkgs = {}

#bkgs[0] = ['diphoJets', '#gamma#gamma+Jets', 'diphoJets.root', 2.43798295134949343e-04, kGreen]
#bkgs[0] = ['diphoJets', '#gamma#gamma+Jets', 'diphoJets.root', 1*1.49648458328955920e-03, kGreen+1] #weights
bkgs[0] = ['diphoJets', '#gamma#gamma+Jets', 'diphoJets.root', 1*2.517683e-03, kGreen+1] #weights
#bkgs[1] = ['gj_20to40', '#gamma+Jets (20to40)', 'gjets20to40.root', 1.58317937176184038e-03, kRed]
bkgs[1] = ['gj_20to40', '#gamma+Jets (20to40)', 'gjets20to40.root', 0.00621686274241627, kRed+1]
#bkgs[2] = ['gj_40toIn', '#gamma+Jets (40toInf)', 'gjets40Inf.root', 2.03614370860702052e-03, kRed+2]
bkgs[2] = ['gj_40toIn', '#gamma+Jets (40toInf)', 'gjets40Inf.root', 1*0.008181706230436655, kRed+3]
#bkgs[3] = ['qcd_30to40', 'QCD (30to40)', 'qcd30to40.root', 2.22683548783303153e-01, kBlue]
bkgs[3] = ['qcd_30to40', 'QCD (30to40-Mgg>80)', 'qcd30to40.root', 0.82925311506336, kBlue+1]
#bkgs[4] = ['qcd_40toIn', 'QCD (40toInf)', 'qcd40toinf.root', 1*8.70970120699850359e-01, kCyan]
bkgs[4] = ['qcd_40toIn', 'QCD (40toInf-Mgg>80)', 'qcd40toinf.root', 3.6937595545044473, kCyan+1]
#bkgs[5] = ['qcd_30toIn', 'QCD (30toInf-40<Mgg<80)', 'QCD_30toInf.root', 4.318834176, kCyan+2]
#bkgs[6] = ['dyjets', 'DYJets', 'DYJets.root', 8.70383969100297673e-03, kYellow]
bkgs[5] = ['dyjets', 'DYJets', 'DYJets.root', 1*6.32980151931769395e-07, kYellow+1] #weights
'''

lumi = 150.0 #in pb

datasets = []

for bkg in data['background']:
	tempfile = TFile(bkg['file'])
	temptree = tempfile.Get('bbggSelectionTree')
	normalization = (lumi*bkg['xsec']*bkg['sfactor'])/(bkg['weight'])
	arr = [bkg['name'], bkg['legend'], bkg['file'], normalization, bkg['color']]
	dataset = [tempfile, temptree, arr]
	datasets.append(dataset)

'''
datasets = []

for bkg in bkgs:
	tempfile = TFile(bkgs[bkg][2])
	temptree = tempfile.Get("bbggSelectionTree")
	dataset = [tempfile, temptree, bkgs[bkg]]
	datasets.append(dataset)
'''
print data['data']
f = TFile( data['data'] )
t = f.Get('bbggSelectionTree')

bhists = {}


Cut = " diphotonCandidate.M() > 80 "
Cut += " && dijetCandidate.M() > 80 "
Cut += " && leadingPhoton.pt() > 35 && subleadingPhoton.pt() > 35 "
Cut += " && leadingJet.pt() > 35 && subleadingJet.pt() > 35 "
#Cut += " && leadingJet.Eta() < 2. && leadingJet.Eta() > -2"
Cut += " && leadingPhotonID[0] == 1 "
Cut += " && subleadingPhotonID[0] == 1 "
#Cut += " && leadingPhotonEVeto == 1 "
#Cut += " && subleadingPhotonEVeto == 1 "
#Cut += " && leadingPhotonEVeto == 1 "
#Cut += " && subleadingPhotonEVeto == 1 "
#Cut += " && leadingPhotonISO[0] == 1 "
#Cut += " && subleadingPhotonISO[0] == 1 "
#Cut += " && leadingJet.pt() > 45 && subleadingJet.pt() > 45 "
#Cut += " && leadingJet_bDis > 0.9 && subleadingJet_bDis > 0.9"
weight = ""
weight += "( genTotalWeight )*"
#weight += "( genTotalWeight/fabs(genTotalWeight) )*"

cut = TCut(weight+"("+Cut+")")
cut_data = TCut(Cut)

#variable = "diHiggsCandidate.Pt()"
#varName = "DiHiggs Candidate p_{T} (GeV)"

plots = {}

nbin = 50

plots[0] = ["diPho_Mass", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 80, 250]
plots[1] = ["diJet_Mass", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 80, 500]
plots[2] = ["leadingJet_pt", "leadingJet.pt()", "Leading Jet Pt (GeV)", nbin, 30, 400] 
plots[3] = ["leadingJet_eta", "leadingJet.Eta()", "Leading Jet Eta", nbin, -4, 4] 
plots[4] = ["leadingPhoton_eta", "leadingPhoton.Eta()", "Leading Photon Eta", nbin, -4, 4] 
plots[5] = ["dicandidate_Mass", "diHiggsCandidate.M()", "DiHiggs Candidate Mass (GeV)", nbin, 250, 2000]
plots[6] = ["leadingPhoton_pt", "leadingPhoton.pt()", "Leading Photon Pt (GeV)", nbin, 30, 300]
plots[7] = ["btagSum", "leadingJet_bDis+subleadingJet_bDis", "Sum of b-tag of jet pair", nbin, 0, 2]
#plots[5] = ["dr_photons", "leadingPhoton.Dot(subleadingPhoton)", "Leading Photon Eta", 20, 0, 5] 

#variable = "diphotonCandidate.M()"
#varName = "DiPhoton Candidate Mass (GeV)"
for plot in plots:
	print plots[plot]
	h1 = TH1F("h1", "Histograms", plots[plot][3], plots[plot][4], plots[plot][5])
	h1.Reset()
	variable = plots[plot][1]
	varName = plots[plot][2]
	stack = myStack('test'+plots[plot][0], varName, varName)

	for bkg in datasets:
		print "Adding", bkg[0], 'to stack!'
	#	bfiles[bkg[0]] = TFile(bkg[2])
	#	btrees[bkg[0]] = bfiles[bkg[0]].Get('bbggSelectionTree')
		bhists[bkg[2][0]] = h1.Clone(bkg[2][0])
		bhists[bkg[2][0]].Reset()
		bkg[1].Draw(variable+'>>'+bkg[2][0], cut, 'hist')
		bhists[bkg[2][0]].Scale(bkg[2][3])
		bhists[bkg[2][0]].SetFillColor(bkg[2][4])
		bhists[bkg[2][0]].SetLineColor(bkg[2][4])
		bhists[bkg[2][0]].SetLineWidth(0)
		bhists[bkg[2][0]].SetFillStyle(3001)
		stack.addHist(bhists[bkg[2][0]], bkg[2][1], bkg[2][3])

	print 'Adding data to stack!'
	#h1.Reset()
#	f = TFile('data.root')
#	t = f.Get('bbggSelectionTree')
	data = h1.Clone("dat")
	data.Reset()
	t.Draw(variable+'>>dat', cut_data)
	data.SetMarkerStyle(8)
	data.Sumw2()
	data.SetMarkerColor(1)
	data.SetLineColor(1)
	data.SetLineWidth(3)
	stack.addData(data, "Data")

	stack.drawStack("DCM_phoIDLoose_25ns_" + plots[plot][0])
#	stack.drawStack("DCM_" + plots[plot][0])
	del stack

