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

os.system("python ../scripts/readLorentzVector_2017Signal.py Signal/output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root bbggSelectionTree 1 1 1")

os.system("python ../scripts/readLorentzVector_2017Signal.py SingleHiggsBackground/ttH/output_ttHToGG_M125_13TeV_powheg_pythia8_v2_0.root bbggSelectionTree 1 1 1")


fB = TFile("SingleHiggsBackground_ttH_output_ttHToGG_M125_13TeV_powheg_pythia8_v2_0.root")
fS = TFile("Signal_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root")

listMET = [["hMET_low_mp;1", 1], ["hMET_low_hp;1", 2], ["hMET_high_mp;1", 3], ["hMET_high_hp;1", 4]]

listMET_njets4 = [["hMET_njets4_low_mp;1", 1], ["hMET_njets4_low_hp;1", 2], ["hMET_njets4_high_mp;1", 3], ["hMET_njets4_high_hp;1", 4]]

listNjets = [["hNjets_low_mp;1", 1], ["hNjets_low_hp;1", 2], ["hNjets_high_mp;1", 3], ["hNjets_high_hp;1", 4]]



enum = [["MET.png", listMET], ["MET_njets4.png", listMET_njets4], ["Njets.png", listNjets]]

###
canv = TCanvas("canv", "canv", 500, 500)

title = TText()
title.SetNDC()

for pair in enum:
    canv.Clear()
    canv.Divide(2,2)
    for hist in pair[1]:
        canv.cd(hist[1])
        print pair[0], " ", hist[1], " ", hist[0]
        hist_to_draw_B = fB.Get(hist[0])
        hist_to_draw_S = fS.Get(hist[0])

        hist_to_draw_B.SetLineColor(kBlue)
        hist_to_draw_B.SetLineStyle(1)
        hist_to_draw_S.SetLineColor(kRed)    
        hist_to_draw_S.SetLineStyle(2)    

        hist_to_draw_B.Scale(1./hist_to_draw_B.GetEntries())
        hist_to_draw_S.Scale(1./hist_to_draw_S.GetEntries())
        
        hist_to_draw_S.Draw()
        hist_to_draw_B.Draw("SAME")


        title.DrawText(0.5, 0.85, hist[0])

        leg = buildLegend(0.60, 0.5, 0.80, 0.84)
        leg.AddEntry(hist_to_draw_S, "HH", "l")
        leg.AddEntry(hist_to_draw_B, "ttH", "l")
        leg.Draw()


    canv.SaveAs(pair[0])


###
