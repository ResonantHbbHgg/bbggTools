from ROOT import *
from array import array

import argparse
parser =  argparse.ArgumentParser(description='Add Classification BDT weights')
parser.add_argument('-F', '--file', dest='File', required=True, type=str)
parser.add_argument('-L', '--low', dest='Low', required=True, type=str)
parser.add_argument('-H', '--high', dest='High', required=True, type=str)
parser.add_argument('-O', '--output', dest='Out', required=True, type=str)

opt = parser.parse_args()


reader_HM = TMVA.Reader()
reader_LM = TMVA.Reader()

ljbt = array('f', [0])
sjbt = array('f', [0])
ctscs = array('f', [0])
ctsbb = array('f', [0])
ctsgg = array('f', [0])
jjratio = array('f', [0])
ggratio = array('f', [0])

reader_HM.AddVariable('leadingJet_bDis', ljbt)
reader_HM.AddVariable('subleadingJet_bDis', sjbt)
reader_HM.AddVariable('diphotonCandidate.Pt()/(diHiggsCandidate.M())', ggratio)
reader_HM.AddVariable("fabs(CosThetaStar_CS)", ctscs)
reader_HM.AddVariable('fabs(CosTheta_bb)', ctsbb)
reader_HM.AddVariable('fabs(CosTheta_gg)', ctsgg)
reader_HM.AddVariable('dijetCandidate.Pt()/(diHiggsCandidate.M())', jjratio)

reader_LM.AddVariable('leadingJet_bDis', ljbt)
reader_LM.AddVariable('subleadingJet_bDis', sjbt)
reader_LM.AddVariable('diphotonCandidate.Pt()/(diHiggsCandidate.M())', ggratio)
reader_LM.AddVariable("fabs(CosThetaStar_CS)", ctscs)
reader_LM.AddVariable('fabs(CosTheta_bb)', ctsbb)
reader_LM.AddVariable('fabs(CosTheta_gg)', ctsgg)
reader_LM.AddVariable('dijetCandidate.Pt()/(diHiggsCandidate.M())', jjratio)

reader_HM.BookMVA("BDT",opt.High)
reader_LM.BookMVA("BDT",opt.Low)

FilesToRedo = opt.File.split(',')

for f in FilesToRedo:
  print f
#  print "ReReco/"+f.split("/")[12]
  infilename = f
  infile = TFile(infilename)
  intree = infile.Get("bbggSelectionTree")

  outfile = TFile(opt.Out+f.split("/")[len(f.split("/"))-1], "RECREATE")
  outtree = intree.CloneTree(0)
  bdt_lm = array('f', [0])
  bdt_hm = array('f', [0])
  bdt = array('f', [0])
  _bdt_lm = outtree.Branch('bdt_lm', bdt_lm, 'bdt_lm/F')
  _bdt_hm = outtree.Branch('bdt_hm', bdt_hm, 'bdt_hm/F')
  _bdt = outtree.Branch('bdt', bdt, 'bdt/F')
  nentries = intree.GetEntries()

  for i in range(0, nentries):
     if i%1000 == 0: print i
     intree.GetEntry(i)

     ctscs[0] = abs(intree.CosThetaStar_CS)
     ctsbb[0] = abs(intree.CosTheta_bb)
     ctsgg[0] = abs(intree.CosTheta_gg)
     ljbt[0] = intree.leadingJet_bDis
     sjbt[0] = intree.subleadingJet_bDis
     jjratio[0] = intree.dijetCandidate.Pt()/intree.diHiggsCandidate.M()
     ggratio[0] = intree.diphotonCandidate.Pt()/intree.diHiggsCandidate.M()

     bdt_hm[0] = reader_HM.EvaluateMVA("BDT")
     bdt_lm[0] = reader_LM.EvaluateMVA("BDT")
     bdt[0] = bdt_hm[0]
     if (intree.diHiggsCandidate.M() - intree.dijetCandidate.M() + 125) < 400: bdt[0] = bdt_lm[0]
     outtree.Fill()

  outfile.cd()
  outtree.Write()
  outfile.Close()
  
