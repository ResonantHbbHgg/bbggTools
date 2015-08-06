##
##	Workspace builder
##

from ROOT import *
from ROOT import RooFit, RooRealVar, RooGaussian, RooDataSet, RooArgList, RooTreeData
from ROOT import RooCmdArg, RooArgSet, kFALSE, RooLinkedList, RooArgusBG, RooAddPdf
from ROOT import RooAbsPdf, RooFormulaVar, RooCategory, RooRealConstant, RooSuperCategory
from ROOT import RooMappedCategory, RooThresholdCategory, RooTruthModel, RooDecay, RooGaussModel
from ROOT import RooAddModel, RooWorkspace

#def BuildWorkspace(file):
file = "../trees/Tree_Grav500.root"
File = TFile(file, "READ")
if(File.IsZombie()):
	print "[BuildWorkspace] File not found!"

bbggTree = File.Get("bbggTree")
if(bbggTree.IsZombie()):
	print "[BuildWorkspace] Tree not found!"

mgg = TH1F("mgg", "Invariant mass of #gamma#gamma candidate", 20, 115, 135)
bbggTree.Draw("diphotonCandidate.mass()>>mgg")

x = RooRealVar("x", "M(#gamma#gamma)", 115, 135)
data = RooDataHist("data", "Dataset with x", RooArgList(x), mgg)

W = RooWorkspace()
W.factory("RooCBShape::CBall(x[115,135], mean[120,130], sigma[0.01,10], alpha[0,10000],n[0,100000])")
W.import(data)
W.Print()

frame = x.frame()
RooAbsData.plotOn(data, frame)
c1 = TCanvas("b", "b", 900, 800)
frame.Draw()
c1.SaveAs("roofit.pdf")

c0 = TCanvas("a", "a", 900, 800)
mgg.Draw()
c0.SaveAs("test.pdf")
