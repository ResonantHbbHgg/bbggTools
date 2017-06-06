#!/usr/bin/env python

import os, sys, ROOT, math
from array import array

rndm = ROOT.TRandom(1)

hggbr = 2.27e-3
genDiPhoFilterFactor = 1./(1 - 0.06)
lumi = 36.5*1000

#L, M, T = 0.5426, 0.8484, 0.9535

def readLorentzVector():



  f = ROOT.TFile.Open('/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar30_ForApproval/'+sys.argv[1], "READ")

  print sys.argv[1]

  treeName = sys.argv[2]

  sigma = float(sys.argv[3])
  nEvt = float(sys.argv[4])
  bWeightSigma = int(sys.argv[5])

  print sigma, ' ', nEvt

  if bWeightSigma: 
    weightSigma = 1./nEvt*hggbr*genDiPhoFilterFactor*sigma*lumi
  else :
    weightSigma = 1.

  hDijet_low_2b = ROOT.TH1F("hDijet_low_2b", "; mjj; nevents", 15, 50, 200);
  hDijet_high_2b = ROOT.TH1F("hDijet_high_2b", "; mjj; nevents", 15, 50, 200);   
  hDijet_low_1b = ROOT.TH1F("hDijet_low_1b", "; mjj; nevents", 15, 50, 200);
  hDijet_high_1b = ROOT.TH1F("hDijet_high_1b", "; mjj; nevents", 15, 50, 200);

  hDijet_2b = ROOT.TH2F("hDijet_2b", "; mjj; mjjgg; nevents", 15, 50, 200, 30, 250, 1000);  


  mychain = f.Get("bbggSelectionTree")
  entries = mychain.GetEntriesFast()


#  t.Branch( 'HHTagger_LM', HHTagger_LM, 'HHTagger_LM/F' )
#  t.Branch( 'HHTagger_HM', HHTagger_HM, 'HHTagger_HM/F' )
#  t.Branch( 'MX', MX, 'MX/F' )

  LMMPC = 0
  HMMPC = 0
  LMHPC = 0
  HMHPC = 0
 
  LMMbtag = 0
  HMMbtag = 0
  LMHbtag = 0
  HMHbtag = 0


  for event in mychain:

    dijet = event.dijetCandidate
    diphoton = event.diphotonCandidate

    mass_jj = dijet.M()
    mass_gg = diphoton.M()

    MX = event.MX
    HHTagger_HM = event.HHTagger_HM
    HHTagger_LM = event.HHTagger_LM

    leadingJet_bDis = event.leadingJet_bDis
    subleadingJet_bDis = event.subleadingJet_bDis
    genTotalWeight = event.genTotalWeight


    bCuts = mass_jj > 80 and mass_jj < 180 and  mass_gg > 100 and  mass_gg < 180
    bLM = MX < 350
    bHM = MX > 350


    bLMMPC = bCuts and bLM and (HHTagger_LM > 0.850 and HHTagger_LM < 0.985 and subleadingJet_bDis > 0.5426 and leadingJet_bDis > 0.5426)
    bHMMPC = bCuts and bHM and (MX > 350 and HHTagger_HM > 0.675 and HHTagger_HM < 0.970)

    bLMHPC = bCuts and bLM and (MX < 350 and HHTagger_LM > 0.985 and subleadingJet_bDis > 0.5426 and leadingJet_bDis > 0.5426)
    bHMHPC = bCuts and bHM and (MX > 350 and HHTagger_HM > 0.970)

    if bLMMPC:
      LMMPC += genTotalWeight
    elif bHMMPC:
      HMMPC += genTotalWeight
    elif bLMHPC:
      LMHPC += genTotalWeight
    elif bHMHPC:
      HMHPC += genTotalWeight

    bHMHbtag = bCuts and bHM and ( subleadingJet_bDis > 0.8484 and leadingJet_bDis > 0.8484 )
    bLMHbtag = bCuts and bLM and ( subleadingJet_bDis > 0.8484 and leadingJet_bDis > 0.8484 )

    bHMMbtag = bCuts and bHM and ( ( leadingJet_bDis > 0.8484 and subleadingJet_bDis > 0.5426 and subleadingJet_bDis < 0.8484 ) or ( subleadingJet_bDis > 0.8484 and leadingJet_bDis > 0.5426 and leadingJet_bDis < 0.8484 ) )
    bLMMbtag = bCuts and bLM and ( ( leadingJet_bDis > 0.8484 and subleadingJet_bDis > 0.5426 and subleadingJet_bDis < 0.8484 ) or ( subleadingJet_bDis > 0.8484 and leadingJet_bDis > 0.5426 and leadingJet_bDis < 0.8484 ) )

    if bLMMbtag:
       LMMbtag += genTotalWeight
    elif bHMMbtag:
       HMMbtag += genTotalWeight
    elif bLMHbtag:
       LMHbtag += genTotalWeight
    elif bHMHbtag:
       HMHbtag += genTotalWeight



  print 'HMHPC ',  HMHPC
  print 'HMMPC ',  HMMPC
  print 'LMHPC ',  LMHPC
  print 'LMMPC ',  LMMPC

  print ''
  print 'HMHPC \t HMMPC'
  print '%.2f' % (HMHPC*weightSigma), '\t', '%.2f' % (HMMPC*weightSigma)
  print 'LMHPC \t LMMPC'
  print '%.2f' % (LMHPC*weightSigma), '\t', '%.2f' % (LMMPC*weightSigma)
  print ''

  print 'HMHbtag ',  HMHbtag
  print 'HMMbtag ',  HMMbtag
  print 'LMHbtag ',  LMHbtag
  print 'LMMbtag ',  LMMbtag

  print ''
  print 'HMHbtag \t HMMbtag'
  print '%.2f' % (HMHbtag*weightSigma), '\t', '%.2f' % (HMMbtag*weightSigma)
  print 'LMHbtag \t LMMbtag'
  print '%.2f' % (LMHbtag*weightSigma), '\t', '%.2f' % (LMMbtag*weightSigma)
  print ''


#    b2btag = ( and subleadingJet_bDis > 0.80);
###    b0btag = (leadingJet_bDis < 0.80 and subleadingJet_bDis < 0.80);
#    b1btag = (not b2btag) and (not b0btag);
###    b1btag = (leadingJet_bDis > 0.80 and subleadingJet_bDis > 0.46 and subleadingJet_bDis < 0.80) or (leadingJet_bDis > 0.46 and subleadingJet_bDis > 0.80 and leadingJet_bDis < 0.80);

#    if mass_hh < 350 and b2btag:
#      hDijet_low_2b.Fill(mass_jj)
#    elif mass_hh > 350 and b2btag:
#      hDijet_high_2b.Fill(mass_jj)	
#    elif mass_hh < 350 and b1btag:
#      hDijet_low_1b.Fill(mass_jj)
#      hDijet_flavour1_1b_low.Fill(mass_jj, leadingJet_flavour)	
#      hDijet_flavour2_1b_low.Fill(mass_jj, subleadingJet_flavour)
#      hpT_flavour1_1b_low.Fill(pT1, leadingJet_flavour)	
#      hpT_flavour2_1b_low.Fill(pT2, subleadingJet_flavour)
#    elif (mass_hh > 350 and b1btag):
#      hDijet_high_1b.Fill(mass_jj)
#      hDijet_flavour1_1b_high.Fill(mass_jj, leadingJet_flavour)	
#      hDijet_flavour2_1b_high.Fill(mass_jj, subleadingJet_flavour)
#      hpT_flavour1_1b_high.Fill(pT1, leadingJet_flavour)	
#      hpT_flavour2_1b_high.Fill(pT2, subleadingJet_flavour)



#  hcosThetaStar_vs_hcosThetaStarHgg.Write("cosThetaStar_vs_cosThetaStarHgg_"+treeName)
#  hcosThetaStar_vs_hcosThetaStarHbb.Write("cosThetaStar_vs_cosThetaStarHbb_"+treeName)
#  hcosThetaStarHgg_vs_hcosThetaStarHbb.Write("cosThetaStarHgg_vs_cosThetaStarHbb_"+treeName)

#  hcosThetaStarHgg_vs_hmgg.Write("cosThetaStarHgg_vs_hmgg_"+treeName)

#  hcosThetaStarHgg_vs_hpTg2.Write("cosThetaStarHgg_vs_hpTg2_"+treeName)

#  corrHH_Hgg =  hcosThetaStar_vs_hcosThetaStarHgg.GetCorrelationFactor();
#  corrHH_Hbb =  hcosThetaStar_vs_hcosThetaStarHbb.GetCorrelationFactor();
#  corrHgg_Hbb =  hcosThetaStarHgg_vs_hcosThetaStarHbb.GetCorrelationFactor();


readLorentzVector()
