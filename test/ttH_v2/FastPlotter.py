#!/usr/bin/env python
from ROOT import *
#gROOT.ProcessLine(".L ../python/MyCMSStyle.py")
#SetAxisTextSizes(obj
#gROOT.ProcessLine("setTDRStyle();")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

import os

def buildLegend(x1, y1, x2, y2):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetLineWidth(0)
    leg.SetFillStyle(0)
    return leg

os.system("python ../../scripts/readLorentzVector_2017Signal_ttH.py Signal_v2/output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph_0.root bbggSelectionTree 1 1 1")

os.system("python ../../scripts/readLorentzVector_2017Signal_ttH.py SingleHiggsBackground/ttH_v2/output_ttHToGG_M125_13TeV_powheg_pythia8_v2_0.root bbggSelectionTree 1 1 1")


fB = TFile("SingleHiggsBackground_ttH_v2_output_ttHToGG_M125_13TeV_powheg_pythia8_v2_0.root")
fS = TFile("Signal_v2_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph_0.root")

listMET = [["hMET_low_mp;1", 1], ["hMET_low_hp;1", 2], ["hMET_high_mp;1", 3], ["hMET_high_hp;1", 4]]

listsumEt = [["hsumEt_low_mp;1", 1], ["hsumEt_low_hp;1", 2], ["hsumEt_high_mp;1", 3], ["hsumEt_high_hp;1", 4]]

listMET_over_sumEt = [["hMET_over_sumEt_low_mp;1", 1], ["hMET_over_sumEt_low_hp;1", 2], ["hMET_over_sumEt_high_mp;1", 3], ["hMET_over_sumEt_high_hp;1", 4]]

listMET_njets4 = [["hMET_njets4_low_mp;1", 1], ["hMET_njets4_low_hp;1", 2], ["hMET_njets4_high_mp;1", 3], ["hMET_njets4_high_hp;1", 4]]

listNjets = [["hNjets_low_mp;1", 1], ["hNjets_low_hp;1", 2], ["hNjets_high_mp;1", 3], ["hNjets_high_hp;1", 4]]




listMET_nleptonsg1 = [["hMET_nleptonsg1_low_mp;1", 1], ["hMET_nleptonsg1_low_hp;1", 2], ["hMET_nleptonsg1_high_mp;1", 3], ["hMET_nleptonsg1_high_hp;1", 4]]

listNleptons = [["hNleptons_low_mp;1", 1], ["hNleptons_low_hp;1", 2], ["hNleptons_high_mp;1", 3], ["hNleptons_high_hp;1", 4]]



listXtt0 = [["hXtt0_low_mp;1", 1], ["hXtt0_low_hp;1", 2], ["hXtt0_high_mp;1", 3], ["hXtt0_high_hp;1", 4]]
listXtt1 = [["hXtt1_low_mp;1", 1], ["hXtt1_low_hp;1", 2], ["hXtt1_high_mp;1", 3], ["hXtt1_high_hp;1", 4]]

listMjjW0 = [["hMjjW0_low_mp;1", 1], ["hMjjW0_low_hp;1", 2], ["hMjjW0_high_mp;1", 3], ["hMjjW0_high_hp;1", 4]]
listMjjW1 = [["hMjjW1_low_mp;1", 1], ["hMjjW1_low_hp;1", 2], ["hMjjW1_high_mp;1", 3], ["hMjjW1_high_hp;1", 4]]

listMjjbt0 = [["hMjjbt0_low_mp;1", 1], ["hMjjbt0_low_hp;1", 2], ["hMjjbt0_high_mp;1", 3], ["hMjjbt0_high_hp;1", 4]]
listMjjbt1 = [["hMjjbt1_low_mp;1", 1], ["hMjjbt1_low_hp;1", 2], ["hMjjbt1_high_mp;1", 3], ["hMjjbt1_high_hp;1", 4]]







enum = [["MET.png", listMET, "S", "", ""], ["MET_njets4.png", listMET_njets4, "S","", "Nj <= 4"], 
        ["MET_over_sumEt.png", listMET_over_sumEt, "B", "", ""],
        ["sumEt.png", listsumEt, "S", "", ""],
        ["MET_nleptonsg1.png", listMET_nleptonsg1, "S","", "Nlept => 1"], 
        ["Njets.png", listNjets, "S","", ""],  
        ["Nleptons.png", listNleptons, "S","LOG", ""], 
        ["Nleptons.png", listNleptons, "S","", ""], 
        ["Xtt0.png", listXtt0, "B","", "NJ >=4"],  ["MjjW0.png", listMjjW0, "B","", "NJ >=4"],  ["Mjjbt0.png", listMjjbt0, "B","", "NJ >=4"],
        ["Xtt1.png", listXtt1, "B","", "NJ >=6"],  ["MjjW1.png", listMjjW1, "B","", "NJ >=6"],  ["Mjjbt1.png", listMjjbt1, "B","", "NJ >=6"]]

###
canv = TCanvas("canv", "canv", 500, 500)

title = TText()
title.SetNDC()

for plotTOplot in enum:
    canv.Clear()
    canv.Divide(2,2)
    for hist in plotTOplot[1]:
        canv.cd(hist[1])
        print plotTOplot[0], " ", hist[1], " ", hist[0]
        hist_to_draw_B = fB.Get(hist[0])
        hist_to_draw_S = fS.Get(hist[0])

        hist_to_draw_B.SetLineColor(kBlue)
        hist_to_draw_B.SetLineStyle(1)
        hist_to_draw_S.SetLineColor(kRed)    
        hist_to_draw_S.SetLineStyle(2)    

        intB = hist_to_draw_B.GetEntries()
        intS = hist_to_draw_S.GetEntries()
        
        if intS > 0:
            hist_to_draw_S.Scale(1./hist_to_draw_S.GetEntries())
        if intB > 0:
            hist_to_draw_B.Scale(1./hist_to_draw_B.GetEntries())

        if "S" in plotTOplot[2]:
            hist_to_draw_S.Draw()
            hist_to_draw_B.Draw("SAME")
        else:
            hist_to_draw_B.Draw()
            hist_to_draw_S.Draw("SAME")

      


        title.DrawText(0.5, 0.85, hist[0])

        leg = buildLegend(0.60, 0.5, 0.90, 0.84)
        leg.SetHeader(plotTOplot[4])
        leg.AddEntry(hist_to_draw_S, "HH", "l")
        leg.AddEntry(hist_to_draw_B, "ttH", "l")
        leg.Draw()
        
        if "LOG" in plotTOplot[3]:
            gPad.SetLogy()

    if "LOG" in plotTOplot[3]:
        canv.SaveAs("LOG_"+plotTOplot[0])
    else: 
        canv.SaveAs(plotTOplot[0])


###
