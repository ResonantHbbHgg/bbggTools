from ROOT import *
from configs import *
from array import array 

gROOT.ProcessLine("gErrorIgnoreLevel = 5000;")

legendSize = 0.03

grs = []
effs = {}
for st in steps:
    effs[st] = []

for i in extraStepsN:
    effs[i] = []

for n,i in enumerate(xaxis):
    print "analyzing file:", ffolder+files[n]
    thisFile = TFile(ffolder+files[n])
    thisEffHist = thisFile.Get("h_Efficiencies")
    totalEvs = thisEffHist.GetBinContent(1)
    for st in steps:
        effs[st].append(float(thisEffHist.GetBinContent(st+1))/float(totalEvs))
    thisTree = thisFile.Get("bbggSelectionTree")
    for ss,st in enumerate(extraSteps):
        thisName = "thisHisto"
        thisHisto = TH1F(thisName, "thisHisto", 1000, -1000, 1000)
        thisTree.Draw("leadingJet_bDis>>"+thisName, "gen_NRW*("+st+")")
        thisIndex = len(steps)+ss
        effs[extraStepsN[ss]].append(thisHisto.Integral()/float(totalEvs))

leg = TLegend(0.15, 0.7, 0.85, 0.89)
leg.SetHeader(header)
leg.SetFillStyle(0)
leg.SetLineWidth(0)
leg.SetNColumns(3)
leg.SetTextSize(legendSize)
if(DEBUG): print effs
for n,i in enumerate(steps+extraStepsN):
    if(DEBUG): print xaxis
    if(DEBUG): print effs[i]
    gr = TGraph(len(xaxis), array('d', xaxis), array('d', effs[i]))
    gr.SetLineColor(kBlack)
    gr.SetMarkerColor(kBlack)
    gr.SetMarkerStyle(24)
    gr.SetMarkerSize(1)
    grs.append(gr)
    leg.AddEntry(grs[n], stepLegs[i], "lp")

c0 = TCanvas("c", "c", 800, 800)
c0.SetGridy()
c0.SetGridx()
for i,gr in enumerate(grs):
    if i == 0:
        gr.Draw("A"+drawOpt)
        gr.SetMaximum(1.3)
        gr.SetMinimum(0.0)
        gr.SetTitle("")
        gr.GetXaxis().SetTitle(xtitle)
        gr.GetYaxis().SetTitle("Acceptance x Efficiency")
        gr.GetYaxis().SetTitleOffset(1.25)
        gr.GetXaxis().SetLimits(xmin,xmax)
        c0.Update()
    else:
        gr.Draw(drawOpt+"same")

grs2 = []
leg2 = TLegend(0.15, 0.7, 0.85, 0.89)
leg2.SetHeader(header)
leg2.SetFillStyle(0)
leg2.SetLineWidth(0)
leg2.SetNColumns(3)
leg2.SetTextSize(legendSize)
for n,i in enumerate(steps+extraStepsN):
    if(DEBUG): print xaxis
    if(DEBUG): print effs[i]
    gr = TGraph(len(xaxis), array('d', xaxis), array('d', effs[i]))
    gr.SetMarkerColorAlpha(TColor.GetColor(colors[n]),0.8)
#    gr.SetMarkerColor(n+1)
    gr.SetMarkerStyle(9)
    gr.SetMarkerSize(1)
    grs2.append(gr)
    leg2.AddEntry(grs2[n], stepLegs[i], "p")

for gr in grs2:
   gr.Draw(drawOpt+"same")

leg.Draw("same")
leg2.Draw("same")
c0.SaveAs(outLoc+outname)
c0.SaveAs(outLoc+outname.replace("pdf", "png"))
