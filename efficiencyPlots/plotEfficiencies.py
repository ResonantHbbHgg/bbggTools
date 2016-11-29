from ROOT import *
from configs import *
from array import array 

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
        thisTree.Draw("leadingJet_bDis>>"+thisName, "("+st+")")
        thisIndex = len(steps)+ss
        effs[extraStepsN[ss]].append(thisHisto.Integral()/float(totalEvs))

leg = TLegend(0.15, 0.78, 0.9, 0.89)
leg.SetHeader("pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma")
leg.SetFillStyle(0)
leg.SetLineWidth(0)
leg.SetNColumns(3)
print effs
for n,i in enumerate(steps+extraStepsN):
    print xaxis
    print effs[i]
    gr = TGraph(len(xaxis), array('d', xaxis), array('d', effs[i]))
    gr.SetLineColor(n+1)
    gr.SetMarkerColor(n+1)
    grs.append(gr)
    leg.AddEntry(grs[n], stepLegs[i], "lp")

c0 = TCanvas("c", "c", 800, 800)
for i,gr in enumerate(grs):
    if i == 0:
        gr.Draw("APL*")
        gr.SetMaximum(1.2)
        gr.SetMinimum(0.0)
        gr.SetTitle("")
        gr.GetXaxis().SetTitle(xtitle)
        gr.GetYaxis().SetTitle("Acceptance x Efficiency")
        gr.GetYaxis().SetTitleOffset(1.25)
        c0.Update()
    else:
        gr.Draw("PL*same")
leg.Draw("same")
c0.SaveAs(outname)

