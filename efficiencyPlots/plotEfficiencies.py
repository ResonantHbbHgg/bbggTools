from ROOT import *
from configs import *
from array import array 
from flashgg.bbggTools.MyCMSStyle import *

gROOT.ProcessLine("gErrorIgnoreLevel = 5000;")

legendSize = 0.03

markerStyles = [ 24, 25, 26, 27, 28, 30, 32 ]
markerStylesFull = [ 20, 21, 22, 33, 34, 29, 23 ]

axes = {}

grs = []
effs = {}
for st in steps:
    effs[st] = []
    axes[st] = []

for i in extraStepsN:
    effs[i] = []
    axes[i] = []

for n,i in enumerate(xaxis):
    print "analyzing file:", ffolder+files[n]
    thisFile = TFile(ffolder+files[n])
    thisEffHist = thisFile.Get("h_Efficiencies")
    totalEvs = 1
    if thisEffHist != None:
      totalEvs = 0.5*thisEffHist.GetBinContent(1)
      for st in steps:
          effs[st].append(float(thisEffHist.GetBinContent(st+1))/float(totalEvs))
          axes[st].append(int(i))
    thisTree = thisFile.Get("bbggSelectionTree")
    for ss,st in enumerate(extraSteps):
        if 'HM' in st and 'Res' in st and int(i) < 600: continue
        if 'LM' in st and 'Res' in st and int(i) > 600: continue
        thisName = "thisHisto"
        thisHisto = TH1F(thisName, "thisHisto", 1000, -1000, 1000)
        thisTree.Draw("leadingJet_bDis>>"+thisName, "1*("+st+")")
        print thisHisto.Integral(), st
        thisIndex = len(steps)+ss
        effs[extraStepsN[ss]].append(thisHisto.Integral()/float(totalEvs))
        axes[extraStepsN[ss]].append(int(i))

leg = TLegend(0.13, legymin, 0.87, 0.89)
leg.SetFillColorAlpha(kWhite, 0.8)
leg.SetHeader(header)
Header = leg.GetListOfPrimitives().First()
Header.SetTextSize(.035)
leg.SetLineWidth(0)
leg.SetNColumns(3)
leg.SetTextSize(legendSize)
if(DEBUG): print effs
for n,i in enumerate(steps+extraStepsN):
    if(DEBUG): print xaxis
    if(DEBUG): print effs[i]
    gr = TGraph(len(axes[i]), array('d', axes[i]), array('d', effs[i]))
    gr.SetLineColor(kBlack)
    gr.SetMarkerColor(kBlack)
    gr.SetMarkerStyle(markerStyles[n])
    gr.SetMarkerSize(2)
    grs.append(gr)
    leg.AddEntry(grs[n], stepLegs[i], "lp")

c0 = TCanvas("c", "c", 800, 750)
c0.SetGridy()
c0.SetGridx()
for i,gr in enumerate(grs):
    if i == 0:
        gr.Draw("A"+drawOpt)
        gr.SetMaximum(ymax)
        gr.SetMinimum(ymin)
        gr.SetTitle("")
        gr.GetXaxis().SetTitle(xtitle)
        gr.GetYaxis().SetTitle("Acceptance x Efficiency")
        gr.GetYaxis().SetTitleOffset(1.25)
        SetAxisTextSizes(gr, yoff)
        gr.GetXaxis().SetLimits(xmin,xmax)
        c0.Update()
    else:
        gr.Draw(drawOpt+"same")

grs2 = []
leg2 = leg.Clone("leg2")
leg2.Clear()
leg2.SetHeader(header)
leg2.SetFillStyle(0)
Header = leg2.GetListOfPrimitives().First()
Header.SetTextSize(.035)
leg2.SetTextSize(legendSize)
for n,i in enumerate(steps+extraStepsN):
    if(DEBUG): print xaxis
    if(DEBUG): print effs[i]
    gr = TGraph(len(axes[i]), array('d', axes[i]), array('d', effs[i]))
#    gr = TGraph(len(xaxis), array('d', xaxis), array('d', effs[i]))
    gr.SetMarkerColorAlpha(TColor.GetColor(colors[n]),0.8)
#    gr.SetMarkerColor(n+1)
    gr.SetMarkerStyle(markerStylesFull[n])
    gr.SetMarkerSize(2)
    grs2.append(gr)
    leg2.AddEntry(grs2[n], stepLegs[i], "p")

for gr in grs2:
   gr.Draw(drawOpt+"same")

leg.Draw("same")
leg2.Draw("same")
DrawCMSLabels(c0, '', 1)
c0.SaveAs(outLoc+outname)
c0.SaveAs(outLoc+outname.replace("pdf", "png"))
