#!/usr/bin/env python

import os, sys, ROOT, math
from array import array

def readLorentzVector():

  f = ROOT.TFile.Open(sys.argv[1], "READ")

  treeName = sys.argv[2]

  hDijet_low_2b = ROOT.TH1F("hDijet_low_2b", "; mjj; nevents", 15, 50, 200);
  hDijet_high_2b = ROOT.TH1F("hDijet_high_2b", "; mjj; nevents", 15, 50, 200);   
  hDijet_low_1b = ROOT.TH1F("hDijet_low_1b", "; mjj; nevents", 15, 50, 200);
  hDijet_high_1b = ROOT.TH1F("hDijet_high_1b", "; mjj; nevents", 15, 50, 200);

# =========================================== #

  hDijet_low_1b_flavour5 = ROOT.TH1F("hDijet_low_1b_flavour5", "; mjj; nevents", 15, 50, 200);
  hDijet_high_1b_flavour5 = ROOT.TH1F("hDijet_high_1b_flavour5", "; mjj; nevents", 15, 50, 200); 
  hDijet_low_1b_flavournon5 = ROOT.TH1F("hDijet_low_1b_flavournon5", "; mjj; nevents", 15, 50, 200);
  hDijet_high_1b_flavournon5 = ROOT.TH1F("hDijet_high_1b_flavournon5", "; mjj; nevents", 15, 50, 200); 

  hDijet_low_2b_flavour5 = ROOT.TH1F("hDijet_low_2b_flavour5", "; mjj; nevents", 15, 50, 200);
  hDijet_high_2b_flavour5 = ROOT.TH1F("hDijet_high_2b_flavour5", "; mjj; nevents", 15, 50, 200); 
  hDijet_low_2b_flavournon5 = ROOT.TH1F("hDijet_low_2b_flavournon5", "; mjj; nevents", 15, 50, 200);
  hDijet_high_2b_flavournon5 = ROOT.TH1F("hDijet_high_2b_flavournon5", "; mjj; nevents", 15, 50, 200); 

# =========================================== #
# =========================================== #
# =========================================== #

  hDijet_2b = ROOT.TH2F("hDijet_2b", "; mjj; mjjgg; nevents", 15, 50, 200, 30, 250, 1000);  
  hDijet_1b = ROOT.TH2F("hDijet_1b", "; mjj; mjjgg; nevents", 15, 50, 200, 30, 250, 1000);  

  hDijet_2b_flavour5 = ROOT.TH2F("hDijet_2b_flavour5", "; mjj; mjjgg; nevents", 15, 50, 200, 30, 250, 1000);  
  hDijet_1b_flavour5 = ROOT.TH2F("hDijet_1b_flavour5", "; mjj; mjjgg; nevents", 15, 50, 200, 30, 250, 1000);  

  hDijet_2b_flavournon5 = ROOT.TH2F("hDijet_2b_flavournon5", "; mjj; mjjgg; nevents", 15, 50, 200, 30, 250, 1000);  
  hDijet_1b_flavournon5 = ROOT.TH2F("hDijet_1b_flavournon5", "; mjj; mjjgg; nevents", 15, 50, 200, 30, 250, 1000);  


  hDijet_btag = ROOT.TH2F("hDijet_btag", "; mjj; btag1+btag2; nevents", 15, 50, 200, 30, -0.1, 2.1);  
  hDijet_flavour1_1b_low = ROOT.TH2F("hDijet_flavour1_1b_low", "; mjj; flavour1; nevents", 15, 50, 200, 50, -25, 25);  
  hDijet_flavour2_1b_low = ROOT.TH2F("hDijet_flavour2_1b_low", "; mjj; flavour2; nevents", 15, 50, 200, 50, -25, 25);  
  hDijet_flavour1_1b_high = ROOT.TH2F("hDijet_flavour1_1b_high", "; mjj; flavour1; nevents", 15, 50, 200, 50, -25, 25);  
  hDijet_flavour2_1b_high = ROOT.TH2F("hDijet_flavour2_1b_high", "; mjj; flavour2; nevents", 15, 50, 200, 50, -25, 25);  

  hpT_flavour1_1b_low = ROOT.TH2F("hpT_flavour1_1b_low", "; p_{T}; flavour1; nevents", 36, 10, 100, 50, -25, 25);  
  hpT_flavour2_1b_low = ROOT.TH2F("hpT_flavour2_1b_low", "; p_{T}; flavour2; nevents", 36, 10, 100, 50, -25, 25); 
  hpT_flavour1_1b_high = ROOT.TH2F("hpT_flavour1_1b_high", "; p_{T}; flavour1; nevents", 36, 10, 100, 50, -25, 25);  
  hpT_flavour2_1b_high = ROOT.TH2F("hpT_flavour2_1b_high", "; p_{T}; flavour2; nevents", 36, 10, 100, 50, -25, 25); 

  hDefcosThetaStar = ROOT.TH1F("hDefcosThetaStar", "; cos (#theta^{*}_{H}); nevents", 40, -1, 1);  
  hcosThetaStar = ROOT.TH1F("hcosThetaStar", "; cos (#theta^{*}_{H}); nevents", 40, -1, 1);  
  hDeltaCosThetaStar = ROOT.TH1F("hDeltaCosThetaStar", "; cos (#theta^{*}_{H, CS}) - cos (#theta^{*}_{H}); nevents", 20, -0.5, 0.5);  
  hcosThetaStarCS = ROOT.TH1F("hcosThetaStarCS", "; cos (#theta^{*}_{H}) in CS frame; nevents", 40, -1, 1);  
  hcosThetaStarHgg = ROOT.TH1F("hcosThetaStarHgg", "; cos (#theta^{*}_{#gamma}); nevents", 40, -1, 1);  
  hcosThetaStarHbb = ROOT.TH1F("hcosThetaStarHbb", "; cos (#theta^{*}_{b}); nevents", 40, -1, 1);  

  hPhi =  ROOT.TH1F("hPhi", "; #phi; nevents", 16, -4, 4);  
  hPhi1 =  ROOT.TH1F("hPhi1", "; #phi_{1}; nevents", 16, -4, 4);  

  hcosThetaStar_vs_hcosThetaStarHgg = ROOT.TH2F("hcosThetaStar_vs_hcosThetaStarHgg", "; cos (#theta^{*}_{H}); cos (#theta^{*}_{#gamma})", 20, -1, 1, 20, -1, 1); 
  hcosThetaStar_vs_hcosThetaStarHbb = ROOT.TH2F("hcosThetaStar_vs_hcosThetaStarHbb", "; cos (#theta^{*}_{H}); cos (#theta^{*}_{b})", 20, -1, 1, 20, -1, 1); 
  hcosThetaStarHgg_vs_hcosThetaStarHbb = ROOT.TH2F("hcosThetaStarHgg_vs_hcosThetaStarHbb", "; cos (#theta^{*}_{H}); cos (#theta^{*}_{b})", 20, -1, 1, 20, -1, 1); 

  mychain = f.Get("bbggSelectionTree")
  entries = mychain.GetEntriesFast()


  f = ROOT.TFile( 'Angles_'+treeName+'.root', 'recreate' )
  t = ROOT.TTree( treeName, 'Tree with Angles' )

  acosThetaStar = array( 'f', [ 0 ] )
  acosThetaStarHgg = array( 'f', [ 0 ] )
  acosThetaStarHbb = array( 'f', [ 0 ] )
  aPhi = array( 'f', [ 0 ] )
  aPhi1 = array( 'f', [ 0 ] )
  amgg = array( 'f', [ 0 ])
  ambb = array( 'f', [ 0 ])
  amhh = array( 'f', [ 0 ])


  t.Branch( 'cosThetaStar', acosThetaStar, 'cosThetaStar/F' )
  t.Branch( 'cosThetaStarHgg', acosThetaStarHgg, 'cosThetaStarHgg/F' )
  t.Branch( 'cosThetaStarHbb', acosThetaStarHbb, 'cosThetaStarHbb/F' )
  t.Branch( 'Phi', aPhi, 'Phi/F' )
  t.Branch( 'Phi1', aPhi1, 'Phi1/F' )
  t.Branch( 'mgg', amgg, 'mgg/F' )
  t.Branch( 'mbb', ambb, 'mbb/F' )
  t.Branch( 'mhh', amhh, 'mhh/F' )






  for event in mychain:
    dijet = event.dijetCandidate
   # print " dijetx %f dijety %f dijetz %f e %f m %f" % (dijet.Px(), dijet.Py(), dijet.Pz(), dijet.E(), dijet.M())
    diphoton = event.diphotonCandidate
#    print " diphotonx %f diphotony %f diphotonz %f e %f m %f" % (diphoton.Px(), diphoton.Py(), diphoton.Pz(), diphoton.E(), diphoton.M())
    dihiggs = event.diHiggsCandidate
#    print " dihiggsx %f dihiggsy %f dihiggsz %f e %f m %f" % (dihiggs.Px(), dihiggs.Py(), dihiggs.Pz(), dihiggs.E(), dihiggs.M())


    mass_jj = dijet.M()
    mass_gg = diphoton.M()
    mass_hh = dihiggs.M()

    leadingJet_bDis = event.leadingJet_bDis
    subleadingJet_bDis = event.subleadingJet_bDis

    leadingJet_flavour = event.leadingJet_flavour
    subleadingJet_flavour = event.subleadingJet_flavour

    b2btag = (leadingJet_bDis > 0.80 and subleadingJet_bDis > 0.80);
    b0btag = (leadingJet_bDis < 0.80 and subleadingJet_bDis < 0.80);
#    b1btag = (not b2btag) and (not b0btag);
    b1btag = (leadingJet_bDis > 0.80 and subleadingJet_bDis > 0.46 and subleadingJet_bDis < 0.80) or (leadingJet_bDis > 0.46 and subleadingJet_bDis > 0.80 and leadingJet_bDis < 0.80);


    leadingJet = event.leadingJet
    pT1 = leadingJet.pt();

    subleadingJet = event.subleadingJet
    pT2 = subleadingJet.pt();

    leadingPhoton = event.leadingPhoton
    subleadingPhoton = event.subleadingPhoton

    
    defCosThetaStar = event.CosThetaStar

    bCS = False

    fCosThetaStar = CosThetaStar(diphoton, dihiggs, False)
    fCosThetaStarCS = CosThetaStar(diphoton, dihiggs, True)
    fCosThetaStarHgg = CosThetaStarh1(leadingPhoton,diphoton)
    fCosThetaStarHbb = CosThetaStarh1(leadingJet,dijet)

    vPhi = {} 
    vPhi = Phi(leadingPhoton, subleadingPhoton, leadingJet, subleadingJet, diphoton, dijet, dihiggs)

    acosThetaStar[0] = fCosThetaStar
    acosThetaStarHgg[0] = fCosThetaStarHgg
    acosThetaStarHbb[0] = fCosThetaStarHbb
    aPhi[0] = vPhi[0]
    aPhi1[0] = vPhi[1]

    ambb[0] = mass_jj
    amgg[0] = mass_gg
    amhh[0] = mass_hh

    t.Fill()


#    print " dihiggs Cos Theta* ", fCosThetaStar, " hgg Cos Theta* ", fCosThetaStarHgg, " hbb Cos Theta* ", fCosThetaStarHbb, " phi ", fPhi, " phi1 ", fPhi1

    hDefcosThetaStar.Fill(defCosThetaStar)
    hcosThetaStarCS.Fill(fCosThetaStarCS)
    hDeltaCosThetaStar.Fill(fCosThetaStar-fCosThetaStarCS)
    hcosThetaStar.Fill(fCosThetaStar)
    hcosThetaStarHgg.Fill(fCosThetaStarHgg)
    hcosThetaStarHbb.Fill(fCosThetaStarHbb)
    
    hPhi.Fill(vPhi[0])
    hPhi1.Fill(vPhi[1])

    hcosThetaStar_vs_hcosThetaStarHgg.Fill(fCosThetaStar, fCosThetaStarHgg)
    hcosThetaStar_vs_hcosThetaStarHbb.Fill(fCosThetaStar, fCosThetaStarHbb)
    hcosThetaStarHgg_vs_hcosThetaStarHbb.Fill(fCosThetaStarHgg, fCosThetaStarHbb)

    if mass_hh < 350 and b2btag:
      hDijet_low_2b.Fill(mass_jj)
    elif mass_hh > 350 and b2btag:
      hDijet_high_2b.Fill(mass_jj)	
    elif mass_hh < 350 and b1btag:
      hDijet_low_1b.Fill(mass_jj)
      hDijet_flavour1_1b_low.Fill(mass_jj, leadingJet_flavour)	
      hDijet_flavour2_1b_low.Fill(mass_jj, subleadingJet_flavour)
      hpT_flavour1_1b_low.Fill(pT1, leadingJet_flavour)	
      hpT_flavour2_1b_low.Fill(pT2, subleadingJet_flavour)
    elif (mass_hh > 350 and b1btag):
      hDijet_high_1b.Fill(mass_jj)
      hDijet_flavour1_1b_high.Fill(mass_jj, leadingJet_flavour)	
      hDijet_flavour2_1b_high.Fill(mass_jj, subleadingJet_flavour)
      hpT_flavour1_1b_high.Fill(pT1, leadingJet_flavour)	
      hpT_flavour2_1b_high.Fill(pT2, subleadingJet_flavour)


    if (mass_hh < 350 and b1btag and abs(subleadingJet_flavour) == 5 and abs(leadingJet_flavour) == 5):
      hDijet_low_1b_flavour5.Fill(mass_jj)
    elif (mass_hh < 350 and b1btag and (abs(subleadingJet_flavour) != 5 or abs(leadingJet_flavour) != 5)):
      hDijet_low_1b_flavournon5.Fill(mass_jj)
    elif (mass_hh > 350 and b1btag and abs(subleadingJet_flavour) == 5 and abs(leadingJet_flavour) == 5):
      hDijet_high_1b_flavour5.Fill(mass_jj)
    elif (mass_hh > 350 and b1btag and (abs(subleadingJet_flavour) != 5 or abs(leadingJet_flavour) != 5)):
      hDijet_high_1b_flavournon5.Fill(mass_jj)


    if (mass_hh < 350 and b2btag and abs(subleadingJet_flavour) == 5 and abs(leadingJet_flavour) == 5):
      hDijet_low_2b_flavour5.Fill(mass_jj)
    elif (mass_hh < 350 and b2btag and (abs(subleadingJet_flavour) != 5 or abs(leadingJet_flavour) != 5)):
      hDijet_low_2b_flavournon5.Fill(mass_jj)
    elif (mass_hh > 350 and b2btag and abs(subleadingJet_flavour) == 5 and abs(leadingJet_flavour) == 5):
      hDijet_high_2b_flavour5.Fill(mass_jj)
    elif (mass_hh > 350 and b2btag and (abs(subleadingJet_flavour) != 5 or abs(leadingJet_flavour) != 5)):
      hDijet_high_2b_flavournon5.Fill(mass_jj)


    if b2btag:
      hDijet_2b.Fill(mass_jj, mass_hh)	
      if (abs(subleadingJet_flavour) == 5 and abs(leadingJet_flavour) == 5):
        hDijet_2b_flavour5.Fill(mass_jj, mass_hh)
      else :
        hDijet_2b_flavournon5.Fill(mass_jj, mass_hh)
    elif b1btag:
      hDijet_1b.Fill(mass_jj, mass_hh)	
      if (abs(subleadingJet_flavour) == 5 and abs(leadingJet_flavour) == 5):
        hDijet_1b_flavour5.Fill(mass_jj, mass_hh)
      else :
        hDijet_1b_flavournon5.Fill(mass_jj, mass_hh)
  
    if mass_hh < 350:
       hDijet_btag.Fill(mass_jj, leadingJet_bDis+subleadingJet_bDis);


  output = ROOT.TFile("OutputAngles.root", "UPDATE")



#  hDijet_low_2b.Write()
#  hDijet_high_2b.Write()
#  hDijet_low_1b.Write()
#  hDijet_high_1b.Write()
#  hDijet_2b.Write()
#  hDijet_1b.Write()
#  hDijet_btag.Write()
#  hDijet_flavour1_1b_high.Write()
#  hDijet_flavour2_1b_high.Write()
#  hDijet_flavour1_1b_low.Write()
#  hDijet_flavour2_1b_low.Write()
#  hDijet_low_1b_flavour5.Write()
#  hDijet_high_1b_flavour5.Write()
#  hDijet_low_1b_flavournon5.Write()
#  hDijet_high_1b_flavournon5.Write()
#  hDijet_low_2b_flavour5.Write()
#  hDijet_high_2b_flavour5.Write()
#  hDijet_low_2b_flavournon5.Write()
#  hDijet_high_2b_flavournon5.Write()

#  hDijet_2b_flavournon5.Write()
#  hDijet_1b_flavournon5.Write()

#  hDijet_2b_flavour5.Write()
#  hDijet_1b_flavour5.Write()

#  hpT_flavour1_1b_low.Write()
#  hpT_flavour2_1b_low.Write()
#  hpT_flavour1_1b_high.Write()
#  hpT_flavour2_1b_high.Write()

  hcosThetaStar_vs_hcosThetaStarHgg.Write("cosThetaStar_vs_cosThetaStarHgg_"+treeName)
  hcosThetaStar_vs_hcosThetaStarHbb.Write("cosThetaStar_vs_cosThetaStarHbb_"+treeName)
  hcosThetaStarHgg_vs_hcosThetaStarHbb.Write("cosThetaStarHgg_vs_cosThetaStarHbb_"+treeName)

  corrHH_Hgg =  hcosThetaStar_vs_hcosThetaStarHgg.GetCorrelationFactor();
  corrHH_Hbb =  hcosThetaStar_vs_hcosThetaStarHbb.GetCorrelationFactor();
  corrHgg_Hbb =  hcosThetaStarHgg_vs_hcosThetaStarHbb.GetCorrelationFactor();

  print "Correlation factor ", treeName, " HH_Hgg", corrHH_Hgg, " HH_Hbb", corrHH_Hbb, " Hgg_Hbb", corrHgg_Hbb

  hcosThetaStarCS.Write("cosThetaStarCS_"+treeName)
  hcosThetaStar.Write("cosThetaStar_"+treeName)
  hDeltaCosThetaStar.Write("DeltaCosThetaStar_"+treeName)
  hcosThetaStarHgg.Write("cosThetaStarHgg_"+treeName)
  hcosThetaStarHbb.Write("cosThetaStarHbb_"+treeName)
    
  hPhi.Write("Phi_"+treeName)
  hPhi1.Write("Phi1_"+treeName)


  output.Close();
  f.Write()
  f.Close()

def CSaxis(H):

  p1 = ROOT.TLorentzVector();
  p2 = ROOT.TLorentzVector();
  p1.SetPxPyPzE(0, 0,  6500, 6500);
  p2.SetPxPyPzE(0, 0,  -6500, 6500);

  boost_H = -H.BoostVector();
  p1.Boost(boost_H)
  p2.Boost(boost_H)


  p1v3 = p1.BoostVector().Unit(); 
  p2v3 = p2.BoostVector().Unit();
  CSaxisv3 = (p1v3 - p2v3).Unit();

  return CSaxisv3;


def CosThetaStar(diphoton, dihiggs, bCS):

  H = ROOT.TLorentzVector()
  H.SetPxPyPzE(dihiggs.px(), dihiggs.py(),  dihiggs.pz(), dihiggs.e());
  boost_H = -H.BoostVector()

  h1 = ROOT.TLorentzVector()
  h1.SetPxPyPzE(diphoton.px(), diphoton.py(),  diphoton.pz(), diphoton.e());

  h1.Boost(boost_H)
  h1_vect = h1.Vect().Unit();

  fCSaxis = CSaxis(H)


  cosThetaStar = h1.CosTheta()
  if bCS: 
    cosThetaStar = fCSaxis.Dot(h1_vect);


  return cosThetaStar;


def CosThetaStarh1(parton, higgs):

  part1 = ROOT.TLorentzVector()
  part1.SetPxPyPzE(parton.px(), parton.py(),  parton.pz(), parton.e())

#  print " partonx %f partony %f partonz %f e %f m %f" % (parton.Px(), parton.Py(), parton.Pz(), parton.E(), parton.M())
#  print " part1x %f part1y %f part1z %f e %f m %f" % (part1.Px(), part1.Py(), part1.Pz(), part1.E(), part1.M())

  h1 = ROOT.TLorentzVector()
  h1.SetPxPyPzE(higgs.px(), higgs.py(),  higgs.pz(), higgs.e());
  boost_h1 = -h1.BoostVector()

#  print " higgsx %f higgsy %f higgsz %f e %f m %f" % (higgs.Px(), higgs.Py(), higgs.Pz(), higgs.E(), higgs.M())
#  print " h1x %f h1y %f h1z %f e %f m %f" % (h1.Px(), h1.Py(), h1.Pz(), h1.E(), h1.M())

  part1.Boost(boost_h1)



  cosThetaStarh1 = part1.CosTheta()


  return cosThetaStarh1;


def norm_planes_hi(leadingPhoton, subleadingPhoton, leadingJet, subleadingJet, diphoton, dijet, dihiggs):

  g1 = ROOT.TLorentzVector();
  g1.SetPxPyPzE(leadingPhoton.px(), leadingPhoton.py(),  leadingPhoton.pz(), leadingPhoton.e());

  g2 = ROOT.TLorentzVector();
  g2.SetPxPyPzE(subleadingPhoton.px(), subleadingPhoton.py(),  subleadingPhoton.pz(), subleadingPhoton.e());

  b1 = ROOT.TLorentzVector();
  b1.SetPxPyPzE(leadingJet.px(), leadingJet.py(),  leadingJet.pz(), leadingJet.e());

  b2 = ROOT.TLorentzVector();
  b2.SetPxPyPzE(subleadingJet.px(), subleadingJet.py(),  subleadingJet.pz(), subleadingJet.e());

  hgg = ROOT.TLorentzVector();
  hgg.SetPxPyPzE(diphoton.px(), diphoton.py(),  diphoton.pz(), diphoton.e());

  hbb = ROOT.TLorentzVector();
  hbb.SetPxPyPzE(dijet.px(), dijet.py(),  dijet.pz(), dijet.e());

  H = ROOT.TLorentzVector();
  H.SetPxPyPzE(dihiggs.px(), dihiggs.py(),  dihiggs.pz(), dihiggs.e());

  boost_H = -H.BoostVector();

  partons = {}
  partons3v = {}
  

  partons[0] = g1
  partons[1] = g2
  partons[2] = b1
  partons[3] = b2
  
  

  list = {0,1,2,3};

  for i in list: 
    partons[i].Boost(boost_H)
    partons3v[i] = partons[i].Vect().Unit();

  vnorm = {}

  list = {0,1};
  for i in list: 
    vnorm[i] = (partons3v[i*2].Cross(partons3v[i*2+1])).Unit();
  
  
  return vnorm;


def Phi(leadingPhoton, subleadingPhoton, leadingJet, subleadingJet, diphoton, dijet, dihiggs):

  vPhi = {}
# Define hgg direction

  hgg = ROOT.TLorentzVector();
  hgg.SetPxPyPzE(diphoton.px(), diphoton.py(),  diphoton.pz(), diphoton.e());

  H = ROOT.TLorentzVector();
  H.SetPxPyPzE(dihiggs.px(), dihiggs.py(),  dihiggs.pz(), dihiggs.e());

  boost_H = - H.BoostVector();
  hgg.Boost(boost_H)
  hgg_vect = hgg.Vect();


#  Calculate the normal to Hgg and Hbb decay plane
  vnorm = norm_planes_hi(leadingPhoton, subleadingPhoton, leadingJet, subleadingJet, diphoton, dijet, dihiggs);


# calculate Phi

  dsignh1 = hgg_vect.Dot(vnorm[1].Cross(vnorm[0]))/(abs(hgg_vect.Dot(vnorm[1].Cross(vnorm[0]))));
  vPhi[0] = dsignh1*(-1)*math.acos(vnorm[0].Dot(vnorm[1]));

#======================================================================

# Define z direction
  p1 = ROOT.TLorentzVector();
  p1.SetPxPyPzE(0, 0,  6500, 6500);
  z_vect = p1.Vect();

#  Calcuate the normal to the zz' plane
  zzprime = (z_vect.Cross(hgg_vect)).Unit();

# calculate Phi1
  dsignh1 = hgg_vect.Dot(zzprime.Cross(vnorm[0]))/(abs(hgg_vect.Dot(zzprime.Cross(vnorm[0]))));
  vPhi[1] = dsignh1*math.acos(zzprime.Dot(vnorm[0]));

  return vPhi;


readLorentzVector()
