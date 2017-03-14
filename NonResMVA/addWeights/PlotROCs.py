from ROOT import *

files = [
#["roc_SM3sig_HM.root", "SM+Node3 HM", kGreen+2],
#["roc_SM3sig_Together_HM.root", "SM+Node3 HM+LM", kBlue],
#["roc_SMsig_HM.root", "SM HM", kRed],
#["roc_SMsig_Together_HM.root", "SM HM+LM", kViolet],
#["roc_HHTagger_HM.root", "All Nodes HM", kOrange]
#["roc_HHTagger_Rad300.root", "All Nodes HM", kOrange],
#["roc_ResTrain_Rad300.root", "Low Mass Res Training", kViolet]
["roc_HHTagger_HM_BSR.root", "Control Region Training", kOrange],
["roc_HBCRTrain_HM_BSR.root", "Blinded CR Training", kViolet],
["roc_HBSRTrain_HM_BSR.root", "Blinded SR Training", kGreen+2]
]

dummy = TFile("dummy.root", "RECREATE")
grs = []

for f in files:
  ff = TFile(f[0])
  gr = ff.Get("ROC")
  dummy.cd()
  gr.SetLineColor(f[2])
  gr.Write()
  grs.append([gr, f[1]])

c = TCanvas("c", "c", 800, 600)
leg = TLegend(0.2, 0.2, 0.5, 0.6)
leg.SetLineWidth(0)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
for ig,g in enumerate(grs):
  if ig == 0:
    g[0].Draw("AL")
    g[0].GetXaxis().SetLimits(0.7, 1)
    g[0].GetXaxis().SetTitle("Signal Efficiency")
    g[0].GetYaxis().SetTitle("Background Rejection (1-Eff)")
    c.Update()
  else:
    g[0].Draw("same")
  g[0].SetLineWidth(2)
  c.Update()
  leg.AddEntry(g[0], g[1], "l")
leg.Draw("same")
c.SaveAs("ROC_BCR.pdf")
