#!/usr/bin/env python

import os, sys, ROOT, math
from array import array
import numpy as n

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
  outfile = ROOT.TFile(outfilename, "RECREATE")

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


  hNjets_low_mp = ROOT.TH1F("hNjets_low_mp", "MX < 350 && MP; njets; nevents", 14, -0.5, 13.5);
  hNjets_low_hp = ROOT.TH1F("hNjets_low_hp", "MX < 350 && HP; njets; nevents", 14, -0.5, 13.5);   
  hNjets_high_mp = ROOT.TH1F("hNjets_high_mp", "MX > 350 && MP; njets; nevents", 14, -0.5, 13.5);
  hNjets_high_hp = ROOT.TH1F("hNjets_high_hp", "MX > 350 && HP; njets; nevents", 14, -0.5, 13.5);

  hNleptons_low_mp = ROOT.TH1F("hNleptons_low_mp", "MX < 350 && MP; nleptons; nevents", 5, -0.5, 4.5);
  hNleptons_low_hp = ROOT.TH1F("hNleptons_low_hp", "MX < 350 && HP; nleptons; nevents", 5, -0.5, 4.5);   
  hNleptons_high_mp = ROOT.TH1F("hNleptons_high_mp", "MX > 350 && MP; nleptons; nevents", 5, -0.5, 4.5);
  hNleptons_high_hp = ROOT.TH1F("hNleptons_high_hp", "MX > 350 && HP; nleptons; nevents", 5, -0.5, 4.5);

  hMET_low_mp = ROOT.TH1F("hMET_low_mp", "; MET; nevents", 20, 0, 200);
  hMET_low_hp = ROOT.TH1F("hMET_low_hp", "; MET; nevents", 20, 0, 200);   
  hMET_high_mp = ROOT.TH1F("hMET_high_mp", "; MET; nevents", 20, 0, 200);
  hMET_high_hp = ROOT.TH1F("hMET_high_hp", "; MET; nevents", 20, 0, 200);

  hsumEt_low_mp = ROOT.TH1F("hsumEt_low_mp", "; sumEt; nevents", 20, 0, 1000);
  hsumEt_low_hp = ROOT.TH1F("hsumEt_low_hp", "; sumEt; nevents", 20, 0, 1000);   
  hsumEt_high_mp = ROOT.TH1F("hsumEt_high_mp", "; sumEt; nevents", 20, 0, 1000);
  hsumEt_high_hp = ROOT.TH1F("hsumEt_high_hp", "; sumEt; nevents", 20, 0, 1000);

  hMET_over_sumEt_low_mp = ROOT.TH1F("hMET_over_sumEt_low_mp", "; MET_over_sumEt; nevents",   40, 0, 5);
  hMET_over_sumEt_low_hp = ROOT.TH1F("hMET_over_sumEt_low_hp", "; MET_over_sumEt; nevents",   40, 0, 5);   
  hMET_over_sumEt_high_mp = ROOT.TH1F("hMET_over_sumEt_high_mp", "; MET_over_sumEt; nevents", 40, 0, 5);
  hMET_over_sumEt_high_hp = ROOT.TH1F("hMET_over_sumEt_high_hp", "; MET_over_sumEt; nevents", 40, 0, 5);

  hMET_nleptonsg1_low_mp = ROOT.TH1F("hMET_nleptonsg1_low_mp", "; MET; nevents", 20, 0, 200);
  hMET_nleptonsg1_low_hp = ROOT.TH1F("hMET_nleptonsg1_low_hp", "; MET; nevents", 20, 0, 200);   
  hMET_nleptonsg1_high_mp = ROOT.TH1F("hMET_nleptonsg1_high_mp", "; MET; nevents", 20, 0, 200);
  hMET_nleptonsg1_high_hp = ROOT.TH1F("hMET_nleptonsg1_high_hp", "; MET; nevents", 20, 0, 200);

  hMET_njets4_low_mp = ROOT.TH1F("hMET_njets4_low_mp", "; MET; nevents", 20, 0, 200);
  hMET_njets4_low_hp = ROOT.TH1F("hMET_njets4_low_hp", "; MET; nevents", 20, 0, 200);   
  hMET_njets4_high_mp = ROOT.TH1F("hMET_njets4_high_mp", "; MET; nevents", 20, 0, 200);
  hMET_njets4_high_hp = ROOT.TH1F("hMET_njets4_high_hp", "; MET; nevents", 20, 0, 200);

  hXtt0_low_mp = ROOT.TH1F("hXtt0_low_mp", "; Xtt0; nevents", 60, 0, 60);
  hXtt0_low_hp = ROOT.TH1F("hXtt0_low_hp", "; Xtt0; nevents", 60, 0, 60);   
  hXtt0_high_mp = ROOT.TH1F("hXtt0_high_mp", "; Xtt0; nevents", 60, 0, 60);
  hXtt0_high_hp = ROOT.TH1F("hXtt0_high_hp", "; Xtt0; nevents", 60, 0, 60);

  hXtt1_low_mp = ROOT.TH1F("hXtt1_low_mp", "; Xtt1; nevents", 60, 0, 60);
  hXtt1_low_hp = ROOT.TH1F("hXtt1_low_hp", "; Xtt1; nevents", 60, 0, 60);   
  hXtt1_high_mp = ROOT.TH1F("hXtt1_high_mp", "; Xtt1; nevents", 60, 0, 60);
  hXtt1_high_hp = ROOT.TH1F("hXtt1_high_hp", "; Xtt1; nevents", 60, 0, 60);

  hMjjW0_low_mp = ROOT.TH1F("hMjjW0_low_mp", "; MjjW0; nevents", 36, 20, 200);
  hMjjW0_low_hp = ROOT.TH1F("hMjjW0_low_hp", "; MjjW0; nevents", 36, 20, 200);   
  hMjjW0_high_mp = ROOT.TH1F("hMjjW0_high_mp", "; MjjW0; nevents", 36, 20, 200);
  hMjjW0_high_hp = ROOT.TH1F("hMjjW0_high_hp", "; MjjW0; nevents", 36, 20, 200);

  hMjjW1_low_mp = ROOT.TH1F("hMjjW1_low_mp", "; MjjW1; nevents", 36, 20, 200);
  hMjjW1_low_hp = ROOT.TH1F("hMjjW1_low_hp", "; MjjW1; nevents", 36, 20, 200);   
  hMjjW1_high_mp = ROOT.TH1F("hMjjW1_high_mp", "; MjjW1; nevents", 36, 20, 200);
  hMjjW1_high_hp = ROOT.TH1F("hMjjW1_high_hp", "; MjjW1; nevents", 36, 20, 200);

  hMjjbt0_low_mp = ROOT.TH1F("hMjjbt0_low_mp", "; Mjjbt0; nevents", 40, 100, 300);
  hMjjbt0_low_hp = ROOT.TH1F("hMjjbt0_low_hp", "; Mjjbt0; nevents", 40, 100, 300);   
  hMjjbt0_high_mp = ROOT.TH1F("hMjjbt0_high_mp", "; Mjjbt0; nevents", 40, 100, 300);
  hMjjbt0_high_hp = ROOT.TH1F("hMjjbt0_high_hp", "; Mjjbt0; nevents", 40, 100, 300);

  hMjjbt1_low_mp = ROOT.TH1F("hMjjbt1_low_mp", "; Mjjbt1; nevents", 40, 100, 300);
  hMjjbt1_low_hp = ROOT.TH1F("hMjjbt1_low_hp", "; Mjjbt1; nevents", 40, 100, 300);   
  hMjjbt1_high_mp = ROOT.TH1F("hMjjbt1_high_mp", "; Mjjbt1; nevents", 40, 100, 300);
  hMjjbt1_high_hp = ROOT.TH1F("hMjjbt1_high_hp", "; Mjjbt1; nevents", 40, 100, 300);



  mychain = f.Get("bbggSelectionTree")
  entries = mychain.GetEntriesFast()

  nsumEt = n.zeros(1, dtype=float)
  nMET = n.zeros(1, dtype=float)
  nXtt0 = n.zeros(1, dtype=float)
  nXtt1 = n.zeros(1, dtype=float)
  nnjets = n.zeros(1, dtype=int)
  nnmus = n.zeros(1, dtype=int)
  nnelecs = n.zeros(1, dtype=int)

  nm_over_ptjj = n.zeros(1, dtype=float)
  nm_over_ptgg = n.zeros(1, dtype=float)

  ndPhi1 = n.zeros(1, dtype=float)
  ndPhi2 = n.zeros(1, dtype=float)

  nCosThetaStar_CS = n.zeros(1, dtype=float)
  nPhoJetMinDr = n.zeros(1, dtype=float)
  nCosTheta_bb = n.zeros(1, dtype=float)


# --------------- BDT tree definition ---------------

  tBDT = ROOT.TTree("BDT", "Tree for ggHH vs ttH BDT")

  tBDT.Branch("sumEt", nsumEt, "sumEt/D")
  tBDT.Branch("MET", nMET, "MET/D")

  tBDT.Branch("Xtt0", nXtt0, "Xtt0/D")
  tBDT.Branch("Xtt1", nXtt1, "Xtt1/D")

  tBDT.Branch("njets", nnjets, "njets/I")
  tBDT.Branch("nmus", nnmus, "nmus/I")
  tBDT.Branch("nelecs", nnelecs, "nelecs/I")

  tBDT.Branch("m_over_ptjj", nm_over_ptjj, "m_over_ptjj/D")
  tBDT.Branch("m_over_ptgg", nm_over_ptgg, "m_over_ptgg/D")

  tBDT.Branch("dPhi1", ndPhi1, "dPhi1/D")
  tBDT.Branch("dPhi2", ndPhi2, "dPhi2/D")

  tBDT.Branch("CosTheta_bb", nCosTheta_bb, "CosTheta_bb/D");
  tBDT.Branch("CosThetaStar_CS", nCosThetaStar_CS, "CosThetaStar_CS/D");

  tBDT.Branch("PhoJetMinDr", nPhoJetMinDr, "PhoJetMinDr/D");
# --------------- --------------- ---------------

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

  nHMHPC_filter_L = 0
  nHMHPC_filter_LMET = 0
  nHMHPC_filter_LMET_JET = 0
  nHMHPC_filter_LMET_JET_H = 0
  nHMHPC_filter_LMET_JET_HFH = 0

  nHMHPC_filter_LMET_NJ3 = 0

  nHMHPC_filter_NJ3 = 0
  nHMHPC_filter_NJ4to5 = 0
  nHMHPC_filter_NJ6to7 = 0
  nHMHPC_filter_NJ8 = 0

  nHMHPC_filter_LMET_H = 0
  nHMHPC_filter_LMET_FH = 0


  LMMbtag = 0
  HMMbtag = 0
  LMHbtag = 0
  HMHbtag = 0


  for event in mychain:

    dijet = event.dijetCandidate
    diphoton = event.diphotonCandidate

    leadingJet = event.leadingJet
    subleadingJet = event.subleadingJet

    CosThetaStar_CS = event.CosThetaStar_CS 
    CosTheta_bb = event.CosTheta_bb

    pMET = event.MET

    MET = pMET.Pt()

    PhoJetMinDr = event.PhoJetMinDr

    v3MET = ROOT.TVector3(pMET.Px(), pMET.Py(), 0)
    v3j1  = ROOT.TVector3(leadingJet.Px(), leadingJet.Py(), 0)
    v3j2  = ROOT.TVector3(subleadingJet.Px(), subleadingJet.Py(), 0)

    dPhi1 = ROOT.TMath.ACos(v3MET.Dot(v3j1)/(MET*v3j1.Mag()+1e-10))
    dPhi2 = ROOT.TMath.ACos(v3MET.Dot(v3j2)/(MET*v3j2.Mag()+1e-10))


    nelecs = event.nelecs
    nmus = event.nmus

    nleptons = nelecs + nmus
    
    Xtt0 = event.Xtt0
    Xtt1 = event.Xtt1

    MjjW0 = event.MjjW0
    MjjW1 = event.MjjW1

    Mjjbt0 = event.Mjjbt0
    Mjjbt1 = event.Mjjbt1


    mass_jj = dijet.M()
    mass_gg = diphoton.M()

    MX = event.MX
    HHTagger_HM = event.HHTagger_HM
    HHTagger_LM = event.HHTagger_LM

    leadingJet_bDis = event.leadingJet_bDis
    subleadingJet_bDis = event.subleadingJet_bDis
    genTotalWeight = event.genTotalWeight

    njets = event.njets

    sumEt = event.sumEt
    sumBjets = leadingJet.pt() + subleadingJet.pt()

    m_over_ptjj = dijet.Pt()/dijet.M()
    m_over_ptgg = diphoton.Pt()/diphoton.M()

    nsumEt[0] = sumEt
    nMET[0] = MET
    nXtt0[0] = Xtt0
    nXtt1[0] = Xtt1
    nnjets[0] = njets
    nnmus[0] = nmus
    nnelecs[0] = nelecs
    nm_over_ptjj[0] = m_over_ptjj
    nm_over_ptgg[0] = m_over_ptgg

    ndPhi1[0] = dPhi1
    ndPhi2[0] = dPhi2

    nCosTheta_bb[0] = CosTheta_bb
    nCosThetaStar_CS[0] = CosThetaStar_CS 
    
    nPhoJetMinDr[0] = PhoJetMinDr



    tBDT.Fill()
    

    # ============================================================================


    bCuts = mass_jj > 80 and mass_jj < 180 and  mass_gg > 100 and  mass_gg < 180
    bLM = MX < 350
    bHM = MX > 350


    bLMMPC = bCuts and bLM and (HHTagger_LM > 0.850 and HHTagger_LM < 0.985 and subleadingJet_bDis > 0.5426 and leadingJet_bDis > 0.5426)
    bHMMPC = bCuts and bHM and (MX > 350 and HHTagger_HM > 0.675 and HHTagger_HM < 0.970)

    bLMHPC = bCuts and bLM and (MX < 350 and HHTagger_LM > 0.985 and subleadingJet_bDis > 0.5426 and leadingJet_bDis > 0.5426)
    bHMHPC = bCuts and bHM and (MX > 350 and HHTagger_HM > 0.970)

# ---------------------------------------------------------

    filterSL_L = nleptons > 0
    filterSL_MET = MET > 130
    filterSL_H = njets > 3 and njets < 6 and ((Xtt0 < 5 and MET > 60) or (MET > 90 and Xtt0 < 10) or MET > 110)
    filterFH = njets > 5 and njets < 8 and ((Xtt0 < 5 and Xtt1 < 10) or (Xtt0 < 10 and Xtt1 < 20 and MET > 70) or MET > 100)
    filterFH_JETS = njets > 7

# ---------------------------------------------------------

    if bLMMPC:
      LMMPC += genTotalWeight
      nLMMPC += 1
      if not filterSL_L and not filterSL_MET and not filterSL_H and not filterFH and not filterFH_JETS:
        nLMMPC_clean += 1
      hNjets_low_mp.Fill(njets)
      hMET_low_mp.Fill(MET)
      hNleptons_low_mp.Fill(nleptons)

      hsumEt_low_mp.Fill(sumEt)
      hMET_over_sumEt_low_mp.Fill(MET/(sumBjets+0.01))

      if nleptons > 0:
        hMET_nleptonsg1_low_mp.Fill(MET)
      if njets < 5:
        hMET_njets4_low_mp.Fill(MET)
      if njets > 3:
        hXtt0_low_mp.Fill(Xtt0)
        hMjjW0_low_mp.Fill(MjjW0)
        hMjjbt0_low_mp.Fill(Mjjbt0)
      if njets > 5:
        hXtt1_low_mp.Fill(Xtt1)
        hMjjW1_low_mp.Fill(MjjW1)
        hMjjbt1_low_mp.Fill(Mjjbt1)


    elif bHMMPC:
      HMMPC += genTotalWeight
      nHMMPC += 1
      if not filterSL_L and not filterSL_MET  and not filterSL_H and not filterFH and not filterFH_JETS:
        nHMMPC_clean += 1
      hNjets_high_mp.Fill(njets)
      hMET_high_mp.Fill(MET)
      hNleptons_high_mp.Fill(nleptons)

      hsumEt_high_mp.Fill(sumEt)
      hMET_over_sumEt_high_mp.Fill(MET/(sumBjets+0.01))

      if nleptons > 0:
        hMET_nleptonsg1_high_mp.Fill(MET)
      if njets < 5:
        hMET_njets4_high_mp.Fill(MET)
      if njets > 3:
        hXtt0_high_mp.Fill(Xtt0)
        hMjjW0_high_mp.Fill(MjjW0)
        hMjjbt0_high_mp.Fill(Mjjbt0)
      if njets > 5:
        hXtt1_high_mp.Fill(Xtt1)
        hMjjW1_high_mp.Fill(MjjW1)
        hMjjbt1_high_mp.Fill(Mjjbt1)


    elif bLMHPC:
      LMHPC += genTotalWeight
      nLMHPC += 1 
      if not filterSL_L and not filterSL_MET  and not filterSL_H and not filterFH and not filterFH_JETS:
        nLMHPC_clean += 1
      hNjets_low_hp.Fill(njets)
      hMET_low_hp.Fill(MET)
      hNleptons_low_hp.Fill(nleptons)

      hsumEt_low_hp.Fill(sumEt)
      hMET_over_sumEt_low_hp.Fill(MET/(sumBjets+0.01))

      if nleptons > 0:
        hMET_nleptonsg1_low_hp.Fill(MET)
      if njets < 5:
        hMET_njets4_low_hp.Fill(MET)
      if njets > 3:
        hXtt0_low_hp.Fill(Xtt0)
        hMjjW0_low_hp.Fill(MjjW0)
        hMjjbt0_low_hp.Fill(Mjjbt0)
      if njets > 5:
        hXtt1_low_hp.Fill(Xtt1)
        hMjjW1_low_hp.Fill(MjjW1)
        hMjjbt1_low_hp.Fill(Mjjbt1)


    elif bHMHPC:
      HMHPC += genTotalWeight
      nHMHPC += 1
      if not filterSL_L and not filterSL_MET  and not filterSL_H and not filterFH and not filterFH_JETS:
        nHMHPC_clean += 1

      if not filterSL_L:
        nHMHPC_filter_L += 1
      if not filterSL_L and not filterSL_MET:
        nHMHPC_filter_LMET += 1
      if not filterSL_L and not filterSL_MET and not filterFH_JETS:
        nHMHPC_filter_LMET_JET += 1
      if not filterSL_L and not filterSL_MET and not filterFH_JETS and not filterSL_H:
        nHMHPC_filter_LMET_JET_H += 1
      if not filterSL_L and not filterSL_MET and not filterFH_JETS and not filterSL_H and not filterFH:
        nHMHPC_filter_LMET_JET_HFH += 1

      if not filterSL_L and not filterSL_MET and njets < 4:
        nHMHPC_filter_LMET_NJ3 += 1

      if njets < 4:
        nHMHPC_filter_NJ3 += 1

      if njets > 3 and njets < 6:
        nHMHPC_filter_NJ4to5 += 1

      if njets > 5 and njets < 8:
        nHMHPC_filter_NJ6to7 += 1

      if njets > 7:
        nHMHPC_filter_NJ8 += 1



      if not filterSL_L and not filterSL_MET and not filterSL_H and  njets > 3 and njets < 6:
        nHMHPC_filter_LMET_H += 1

      if not filterSL_L and not filterSL_MET and not filterFH and  njets > 5 and njets < 8:
        nHMHPC_filter_LMET_FH += 1



      hNjets_high_hp.Fill(njets)
      hMET_high_hp.Fill(MET)
      hNleptons_high_hp.Fill(nleptons)

      hsumEt_high_hp.Fill(sumEt)
      hMET_over_sumEt_high_hp.Fill(MET/(sumBjets+0.01))

      if nleptons > 0:
        hMET_nleptonsg1_high_hp.Fill(MET)
      if njets < 5:
        hMET_njets4_high_hp.Fill(MET)
      if njets > 3:
        hXtt0_high_hp.Fill(Xtt0)
        hMjjW0_high_hp.Fill(MjjW0)
        hMjjbt0_high_hp.Fill(Mjjbt0)
      if njets > 5:
        hXtt1_high_hp.Fill(Xtt1)
        hMjjW1_high_hp.Fill(MjjW1)
        hMjjbt1_high_hp.Fill(Mjjbt1)



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

  print 'High Mass High Purity'
  print '---------------------'
  print 'lept > 0 ', '%.3f' % (nHMHPC_filter_L/(nHMHPC+1.0))
  print 'lept > 0 && MET > 130 ', '%.3f' %(nHMHPC_filter_LMET/(nHMHPC+1.0))
  print 'lept > 0 && MET > 130 && Njets > 8',  '%.3f' %(nHMHPC_filter_LMET_JET/(nHMHPC+1.0))
  print 'lept > 0 && MET > 130 && Njets > 8 && Wjj SL cut',  '%.3f' %(nHMHPC_filter_LMET_JET_H/(nHMHPC+1.0))
  print 'lept > 0 && MET > 130 && Njets > 8 && Wjj SL cut && Full H cut ',  '%.3f' %(nHMHPC_filter_LMET_JET_HFH/(nHMHPC+1.0))
  print '---------------------'
  print 'ttH hard to remove: Njets < 4',  '\t \t %.3f' %(nHMHPC_filter_NJ3/(nHMHPC+1.0))
  print '    lept > 0 && MET > 130 && Njets < 4',      '\t\t  relative %.3f' %(nHMHPC_filter_LMET_NJ3/(nHMHPC_filter_NJ3+1.0))
  print '---------------------'
  print 'ttH hard to remove: Njets > 3 and Njets < 6',  '\t \t %.3f' %(nHMHPC_filter_NJ4to5/(nHMHPC+1.0))
  print '    lept > 0 && MET > 130 && Wjj SL cut',     '\t \t relative %.3f' %(nHMHPC_filter_LMET_H/(nHMHPC_filter_NJ4to5+1.0))
  print '---------------------'
  print 'Njets > 6 and Njets < 8',  '\t \t %.3f' %(nHMHPC_filter_NJ6to7/(nHMHPC+1.0))
  print '    lept > 0 && MET > 130 && Fully Hadronic', '\t \t relative  %.3f' %(nHMHPC_filter_LMET_FH/(nHMHPC_filter_NJ6to7+1.0))
  print '---------------------'
  print 'Njets > 7',  '\t \t %.3f' %(nHMHPC_filter_NJ8/(nHMHPC+1.0))
  print '    removed Njets > 7: \t \t 0' 
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


  hNleptons_low_mp.Write()
  hNleptons_low_hp.Write()
  hNleptons_high_mp.Write()
  hNleptons_high_hp.Write()

  hMET_low_mp.Write()
  hMET_low_hp.Write()
  hMET_high_mp.Write()
  hMET_high_hp.Write()

  hMET_nleptonsg1_low_mp.Write()
  hMET_nleptonsg1_low_hp.Write()
  hMET_nleptonsg1_high_mp.Write()
  hMET_nleptonsg1_high_hp.Write()

  hMET_njets4_low_mp.Write()
  hMET_njets4_low_hp.Write()
  hMET_njets4_high_mp.Write()
  hMET_njets4_high_hp.Write()

  hXtt0_low_mp.Write()
  hXtt0_low_hp.Write()
  hXtt0_high_mp.Write()
  hXtt0_high_hp.Write()

  hXtt1_low_mp.Write()
  hXtt1_low_hp.Write()
  hXtt1_high_mp.Write()
  hXtt1_high_hp.Write()

  hMjjW0_low_mp.Write()
  hMjjW0_low_hp.Write()
  hMjjW0_high_mp.Write()
  hMjjW0_high_hp.Write()

  hMjjW1_low_mp.Write()
  hMjjW1_low_hp.Write()
  hMjjW1_high_mp.Write()
  hMjjW1_high_hp.Write()

  hMjjbt0_low_mp.Write()
  hMjjbt0_low_hp.Write()
  hMjjbt0_high_mp.Write()
  hMjjbt0_high_hp.Write()

  hMjjbt1_low_mp.Write()
  hMjjbt1_low_hp.Write()
  hMjjbt1_high_mp.Write()
  hMjjbt1_high_hp.Write()

  hsumEt_low_mp.Write()
  hsumEt_low_hp.Write()
  hsumEt_high_mp.Write()
  hsumEt_high_hp.Write()

  hMET_over_sumEt_low_mp.Write()
  hMET_over_sumEt_low_hp.Write()
  hMET_over_sumEt_high_mp.Write()
  hMET_over_sumEt_high_hp.Write()

  tBDT.Write()

  outfile.Close()

readLorentzVector()
