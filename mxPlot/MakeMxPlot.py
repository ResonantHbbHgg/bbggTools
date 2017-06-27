from ROOT import *
from flashgg.bbggTools.MyCMSStyle import *

cNiceRed = TColor.GetColor('#FA4912')
cNicePurple = TColor.GetColor('#885BB2')

#gROOT.SetBatch()

masses = [280, 350, 450, 600, 800]

filesLoc = '/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Signal_Mjj70_NewMVA/Hadd/output_GluGluToBulkGravitonToHHTo2B2G_M-MASS_narrow_13TeV-madgraph.root'

axes = [
[211, 349],
[261, 439],
[341, 559],
[451, 749],
[601, 999]
]

hists = []
hists_org = []
files = []
trees = []

for mm in masses:
  tf = TFile(filesLoc.replace("MASS", str(mm)))
  files.append(tf)
  tr = tf.Get("bbggSelectionTree")
  trees.append(tr)

ofile = TFile("out.root", "RECREATE")

for imm,mm in enumerate(masses):
  hh = TH1F('h_'+str(mm), ';;Normalized Yields [A.U.]', 40, axes[imm][0], axes[imm][1])
  hh.SetLineWidth(2)
  hh.SetLineColor(cNiceRed)
  hh_org = TH1F('o_'+str(mm), ';;Normalized Yields [A.U.]', 40, axes[imm][0], axes[imm][1])
#  hh_org.SetLineStyle(kDotted)
  hh_org.SetLineWidth(2)
  hh_org.SetLineColor(cNicePurple)
  tf = files[imm]
  tree = trees[imm]
  print tree.GetEntries("isSignal")
  tree.Draw("(diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)>>h_"+str(mm), "isSignal")
  tree.Draw("(diHiggsCandidate.M())>>o_"+str(mm), "isSignal")
  hh_org.Print()
  print hh.Integral()
  hists.append(hh)
  hists_org.append(hh_org)

gStyle.SetOptStat(0)
ofile.cd()
cc = TCanvas("cc", "cc", 1000, 350)
#cc.SetLeftMargin(0.15)
margin = 0.02
cc.SetBottomMargin(0)
cc.SetRightMargin(0.012)
cc.SetLeftMargin(0.005)
#cc.Divide(len(masses),1,0,0)
legs = []
for mm,M in enumerate(masses):
  hists[mm].GetXaxis().SetNdivisions(504,1)
  hists[mm].GetXaxis().SetLabelFont(43)
  hists[mm].GetXaxis().SetLabelSize(15)
#  if mm == 0:
#    hists[mm].GetYaxis().SetLabelFont(43)
#    hists[mm].GetYaxis().SetLabelSize(15)  
#    hists[mm].GetYaxis().SetTitleOffset()
  #  SetAxisTextSizes(hists[mm])
  #  pad = cc.cd(mm+1)
  pad = TPad(str(M), str(M), margin-0.01+float(mm)*(1-margin)/float(len(masses)), 0.05, margin-0.01+float(mm+1)*(1-margin)/float(len(masses)), 1)
  pad.SetLeftMargin(0)
  pad.SetRightMargin(0)
  pad.cd()
  hnorm = hists[mm].Integral()
  hists[mm].Scale(1./hnorm)
  onorm = hists_org[mm].Integral()
  hists_org[mm].Scale(1./onorm)
  hists[mm].SetMaximum(hists[mm].GetMaximum()*1.2)
  hists[mm].Draw("")
  hists_org[mm].Draw("same")
  tlt = TLatex()
  tlt.SetTextFont(43)
  tlt.SetTextSize(18)
  tlt.DrawLatexNDC(0.04, 0.84, "M = " + str(M) + ' GeV')
#  else: tlt.DrawLatexNDC(0.04+0.035, 0.85, "M = " + str(M) + ' GeV')
  if mm == 0:
    leg = TLegend(0.68, 0.7, 0.95, 0.85)
    leg.SetTextFont(43)
    leg.SetTextSize(15)
    leg.AddEntry(hists[mm], '#tilde{M}_{X}', 'l')
    leg.AddEntry(hists_org[mm], 'M(#gamma#gammajj)', 'l')
    leg.SetBorderSize(0)
    legs.append(leg)
    legs[mm].Draw("same")
  hists[mm].Write()
  hists_org[mm].Write()
  cc.cd()
  pad.Draw()
cc.cd()
tl = TLatex()
tl.SetTextFont(43)
tl.SetTextSize(18)
tl.DrawLatexNDC(0.78, 0.03, "Reconstructed Mass [GeV]")
DrawCMSLabels(cc, '', 1, 0, 0.055)
cc.SaveAs("/afs/cern.ch/user/r/rateixei/www/HHBBGG/MXAnalyzer/MX.pdf")
