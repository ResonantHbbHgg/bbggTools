#!/usr/bin/env python

import os, sys, ROOT, math
from array import array

rndm = ROOT.TRandom(1)

hggbr = 2.27e-3
genDiPhoFilterFactor = 1./(1 - 0.06)
lumi = 36.5*1000

#L, M, T = 0.5426, 0.8484, 0.9535

def readLorentzVector():



#  f = ROOT.TFile.Open('/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar30_ForApproval/'+sys.argv[1], "READ")
  f = ROOT.TFile.Open('/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Dec13_ForTTH_rejection/'+sys.argv[1], "READ")

  print sys.argv[1]

  treeName = sys.argv[2]

  outfilename = sys.argv[1].replace('/','_')

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


  hNjets_low_mp = ROOT.TH1F("hNjets_low_mp", "MX < 350 && MP; njets; nevents", 14, -0., 13.5);
  hNjets_low_hp = ROOT.TH1F("hNjets_low_hp", "MX < 350 && HP; njets; nevents", 14, -0., 13.5);   
  hNjets_high_mp = ROOT.TH1F("hNjets_high_mp", "MX > 350 && MP; njets; nevents", 14, -0.5, 13.5);
  hNjets_high_hp = ROOT.TH1F("hNjets_high_hp", "MX > 350 && HP; njets; nevents", 14, -0.5, 13.5);


  hMET_low_mp = ROOT.TH1F("hMET_low_mp", "; MET; nevents", 20, 0, 200);
  hMET_low_hp = ROOT.TH1F("hMET_low_hp", "; MET; nevents", 20, 0, 200);   
  hMET_high_mp = ROOT.TH1F("hMET_high_mp", "; MET; nevents", 20, 0, 200);
  hMET_high_hp = ROOT.TH1F("hMET_high_hp", "; MET; nevents", 20, 0, 200);


  hMET_njets4_low_mp = ROOT.TH1F("hMET_njets4_low_mp", "; MET; nevents", 20, 0, 200);
  hMET_njets4_low_hp = ROOT.TH1F("hMET_njets4_low_hp", "; MET; nevents", 20, 0, 200);   
  hMET_njets4_high_mp = ROOT.TH1F("hMET_njets4_high_mp", "; MET; nevents", 20, 0, 200);
  hMET_njets4_high_hp = ROOT.TH1F("hMET_njets4_high_hp", "; MET; nevents", 20, 0, 200);



  mychain = f.Get("bbggSelectionTree")
  entries = mychain.GetEntriesFast()


#  t.Branch( 'HHTagger_LM', HHTagger_LM, 'HHTagger_LM/F' )
#  t.Branch( 'HHTagger_HM', HHTagger_HM, 'HHTagger_HM/F' )
#  t.Branch( 'MX', MX, 'MX/F' )

  LMMPC = 0
  HMMPC = 0
  LMHPC = 0
  HMHPC = 0


  nLMMPC = 0
  nHMMPC = 0
  nLMHPC = 0
  nHMHPC = 0
 

  nLMMPC_clean = 0
  nHMMPC_clean = 0
  nLMHPC_clean = 0
  nHMHPC_clean = 0

  LMMbtag = 0
  HMMbtag = 0
  LMHbtag = 0
  HMHbtag = 0


  for event in mychain:

    dijet = event.dijetCandidate
    diphoton = event.diphotonCandidate
    pMET = event.MET

    MET = pMET.Pt()

    mass_jj = dijet.M()
    mass_gg = diphoton.M()

    MX = event.MX
    HHTagger_HM = event.HHTagger_HM
    HHTagger_LM = event.HHTagger_LM

    leadingJet_bDis = event.leadingJet_bDis
    subleadingJet_bDis = event.subleadingJet_bDis
    genTotalWeight = event.genTotalWeight

    njets = event.njets

    bCuts = mass_jj > 80 and mass_jj < 180 and  mass_gg > 100 and  mass_gg < 180
    bLM = MX < 350
    bHM = MX > 350


    bLMMPC = bCuts and bLM and (HHTagger_LM > 0.850 and HHTagger_LM < 0.985 and subleadingJet_bDis > 0.5426 and leadingJet_bDis > 0.5426)
    bHMMPC = bCuts and bHM and (MX > 350 and HHTagger_HM > 0.675 and HHTagger_HM < 0.970)

    bLMHPC = bCuts and bLM and (MX < 350 and HHTagger_LM > 0.985 and subleadingJet_bDis > 0.5426 and leadingJet_bDis > 0.5426)
    bHMHPC = bCuts and bHM and (MX > 350 and HHTagger_HM > 0.970)

    if bLMMPC:
      LMMPC += genTotalWeight
      nLMMPC += 1
      if njets < 7 and MET < 90:
        nLMMPC_clean += 1
      hNjets_low_mp.Fill(njets)
      hMET_low_mp.Fill(MET)
      if njets < 5:
        hMET_njets4_low_mp.Fill(MET)
    elif bHMMPC:
      HMMPC += genTotalWeight
      nHMMPC += 1
      if njets < 7 and MET < 90:
        nHMMPC_clean += 1
      hNjets_high_mp.Fill(njets)
      hMET_high_mp.Fill(MET)
      if njets < 5:
        hMET_njets4_high_mp.Fill(MET)
    elif bLMHPC:
      LMHPC += genTotalWeight
      nLMHPC += 1 
      if njets < 7 and MET < 90:
        nLMHPC_clean += 1
      hNjets_low_hp.Fill(njets)
      hMET_low_hp.Fill(MET)
      if njets < 5:
        hMET_njets4_low_hp.Fill(MET)
    elif bHMHPC:
      HMHPC += genTotalWeight
      nHMHPC += 1
      if njets < 7 and MET < 90:
        nHMHPC_clean += 1
      hNjets_high_hp.Fill(njets)
      hMET_high_hp.Fill(MET)
      if njets < 5:
        hMET_njets4_high_hp.Fill(MET)




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

  print 'nHMHPC ',  nHMHPC, ' nHMHPC_clean ',  nHMHPC_clean, " rate ", nHMHPC_clean/(nHMHPC+1.0)
  print 'nHMMPC ',  nHMMPC, ' nHMMPC_clean ',  nHMMPC_clean, " rate ", nHMMPC_clean/(nHMMPC+1.0)
  print 'nLMHPC ',  nLMHPC, ' nLMHPC_clean ',  nLMHPC_clean, " rate ", nLMHPC_clean/(nLMHPC+1.0)
  print 'nLMMPC ',  nLMMPC, ' nLMMPC_clean ',  nLMMPC_clean, " rate ", nLMMPC_clean/(nLMMPC+1.0)


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
  outfile = ROOT.TFile(outfilename, "RECREATE")
  hMET_low_mp.Write()
  hMET_low_hp.Write()
  hMET_high_mp.Write()
  hMET_high_hp.Write()

  hMET_njets4_low_mp.Write()
  hMET_njets4_low_hp.Write()
  hMET_njets4_high_mp.Write()
  hMET_njets4_high_hp.Write()

  hNjets_low_mp.Write()
  hNjets_low_hp.Write()
  hNjets_high_mp.Write()
  hNjets_high_hp.Write()
  outfile.Close()

readLorentzVector()
