from ROOT import *

files = [
#["EffCR_SM3sig_HM.root", "SM+Node3 HM", kGreen+2],
#["EffCR_SM3sig_HM_AllHH.root", "SM+Node3 HM+LM", kBlue],
#["EffCR_SMsig_HM.root", "SM HM", kRed],
#["EffCR_SMsig_HM_AllHH.root", "SM HM+LM", kViolet],
#["EffCR_HHTagger_HM.root", "All Nodes HM", kOrange]
["EffCR_HHTagger_LM.root", "All Nodes LM", kOrange],
["EffCR_ResLM.root", "Low Mass Res Training", kViolet]
#["EffBSR_HHTagger_HM.root", "Control Region Training", kOrange],
#["EffBSR_BCR_HM_AllHH.root", "Blinded CR Training", kViolet],
#["EffBSR_BSR_HM_AllHH.root", "Blinded SR Training", kGreen+2]
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
c.SetLogy()
leg = TLegend(0.2, 0.6, 0.5, 0.89)
leg.SetLineWidth(0)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
f1 = TF1("f1", "pol1", 0,1)
f1.SetParameters(0,1)
f1.SetLineColor(kBlack)
f1.SetLineStyle(kDashed)
f1.SetLineWidth(2)
for ig,g in enumerate(grs):
  if ig == 0:
    g[0].Draw("AL")
    g[0].GetYaxis().SetLimits(0.015, 0.2)
    g[0].GetYaxis().SetRangeUser(0.015, 0.2)
    g[0].GetXaxis().SetLimits(0.7,1)
    g[0].GetXaxis().SetTitle("Signal Efficiency")
    g[0].GetYaxis().SetTitle("Background Efficiency")
    c.Update()
  else:
    g[0].Draw("same")
  g[0].SetLineWidth(2)
#  f1.Draw("same")
  c.Update()
  leg.AddEntry(g[0], g[1], "l")
leg.Draw("same")
c.SaveAs("reff_res.pdf")
