from ROOT import *
import math
from flashgg.bbggTools.MyCMSStyle import *

cNiceRed = TColor.GetColor('#FA4912')
cNicePurple = TColor.GetColor('#885BB2')

#gROOT.SetBatch()

masses = [280, 450, 800]

filesLoc = '/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Signal_Mjj70_NewMVA/Hadd/output_GluGluToBulkGravitonToHHTo2B2G_M-MASS_narrow_13TeV-madgraph.root'

axes = []

#axes = [
#[211, 349],
#[341, 559],
#[601, 999]
#]


hists = []
hists_org = []
files = []
trees = []

W = 140

for mm in masses:
  tf = TFile(filesLoc.replace("MASS", str(mm)))
  files.append(tf)
  tr = tf.Get("bbggSelectionTree")
  trees.append(tr)
  axes.append([mm-W/2, mm+W/2])

ofile = TFile("out.root", "RECREATE")

gStyle.SetOptStat(0)
for imm,mm in enumerate(masses):
  hh = TH1F('h_'+str(mm), ';;', 40, axes[imm][0], axes[imm][1])
  hh.SetLineWidth(3)
  hh.SetLineColor(cNiceRed)
  hh_org = TH1F('o_'+str(mm), ';;', 40, axes[imm][0], axes[imm][1])
#  hh_org.SetLineStyle(kDotted)
  hh_org.SetLineWidth(3)
  hh_org.SetLineColor(cNicePurple)
  hh_org.SetLineStyle(kDashed)
  tf = files[imm]
  tree = trees[imm]
  print tree.GetEntries("isSignal")
  tree.Draw("(diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)>>h_"+str(mm), "isSignal")
  tree.Draw("(diHiggsCandidate.M())>>o_"+str(mm), "isSignal")
  hh_org.Print()
  print hh.Integral()
  hists.append(hh)
  hists_org.append(hh_org)


ofile.cd()
cc = TCanvas("cc", "cc", 600, 350)
#SetPadStyle(cc)
cc.SetLeftMargin(0.08)
cc.SetRightMargin(0.01)
margin = 0.06
extramargin = 0.09
#cc.SetBottomMargin(0)
#cc.SetRightMargin(0.012)
#cc.SetLeftMargin(0.005)
#cc.Divide(len(masses),1,0,0)
legs = []

Max = 0

for mm,M in enumerate(masses):
  if mm > 10:
    continue
  hists[mm].GetXaxis().SetNdivisions(504,1)
  hists[mm].GetXaxis().SetLabelFont(43)
  hists[mm].GetXaxis().SetLabelSize(15)
#  if mm == 0:
#    hists[mm].GetYaxis().SetLabelFont(43)
  hists[mm].GetXaxis().SetLabelSize(15)  
#    hists[mm].GetYaxis().SetTitleOffset()
  #  SetAxisTextSizes(hists[mm])
  #  pad = cc.cd(mm+1)

  lm = len(masses)
  padW = (1-margin)/(lm+extramargin)

  # padW is defind in such a way to keep the same aspect off all 3 plots

  print ' pad Width ', padW

  if mm == 0:
    pad = TPad(str(M), str(M), margin+float(mm)*padW-0.01, 0.05, margin+extramargin*padW+float(mm+1)*padW-0.01, 1-0.02)
    pad.SetLeftMargin(extramargin)
    hists[mm].GetYaxis().SetLabelSize(0.06)
  else :
    pad = TPad(str(M), str(M), margin+extramargin*padW+float(mm)*padW-0.01, 0.05, margin+extramargin*padW+float(mm+1)*padW-0.01, 1-0.02)
    pad.SetLeftMargin(0)
  
  pad.SetRightMargin(0) 
  hists[mm].GetYaxis().SetNdivisions(305) 

  SetPadStyle(pad)
  pad.cd()
  hnorm = hists[mm].Integral()
  hists[mm].Scale(1./hnorm)
  onorm = hists_org[mm].Integral()
  hists_org[mm].Scale(1./onorm)
  if Max < hists[mm].GetMaximum():
    Max = hists[mm].GetMaximum()
  hists[mm].SetMaximum(hists[mm].GetMaximum()*1.2)
  hists[mm].Draw("")
  hists_org[mm].Draw("same")
  tlt = TLatex()
  tlt.SetTextFont(43)
  tlt.SetTextSize(15)
  if mm == 0:
    tlt.DrawLatexNDC(0.13, 0.82, "m_{X} = " + str(M) + ' GeV')
  else :
    tlt.DrawLatexNDC(0.04, 0.82, "m_{X} = " + str(M) + ' GeV')
#  else: tlt.DrawLatexNDC(0.04+0.035, 0.85, "M = " + str(M) + ' GeV')
  if mm == 0:
    leg = TLegend(0.65, 0.65, 0.95, 0.80)
    leg.SetTextFont(43)
    leg.SetTextSize(15)
    leg.AddEntry(hists[mm], '#tilde{M}_{X}', 'l')
    leg.AddEntry(hists_org[mm], 'm_{#gamma#gammajj}', 'l')
    leg.SetBorderSize(0)
    legs.append(leg)
    legs[mm].Draw("same")
  hists[mm].Write()
  hists_org[mm].Write()
  cc.cd()
  pad.Draw()

for mm,M in enumerate(masses):
  hists[mm].SetMaximum(Max*1.2)

cc.Update()

cc.cd()
tl = TLatex()
tl.SetTextFont(43)
tl.SetTextSize(16)
tl.DrawLatexNDC(0.68, 0.05, "Reconstructed mass [GeV]")
tl.SetTextAngle(90)
tl.DrawLatexNDC(0.04, 0.39, "Events fraction /(7 GeV)")
#DrawCMSLabels(cc, '', 1)
DrawCMSLabels(cc, '', 1, 0, 0.055)
cc.SaveAs("MX.pdf")
