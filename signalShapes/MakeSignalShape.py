from ROOT import *
gROOT.SetBatch(True)

def DSCrystalBall(xx,par):
  x = float(xx[0])
  N = float(par[0]) #normalization
  mu = float(par[1]) #mean
  sig = float(par[2]) #sigma
  ##tail parameters
  a1 = float(par[3])
  p1 = float(par[4])
  a2 = float(par[5])
  p2 = float(par[6])

  u = (x-mu)/sig
  A1 = TMath.Power(p1/TMath.Abs(a1), p1)*TMath.Exp(-a1*a1/2)
  A2 = TMath.Power(p2/TMath.Abs(a2), p2)*TMath.Exp(-a2*a2/2)
  B1 = p1/TMath.Abs(a1) - TMath.Abs(a1)
  B2 = p2/TMath.Abs(a2) - TMath.Abs(a2)

  Res = -99
  if u < -a1: Res = A1*TMath.Power(B1-u,-p1)
  elif u < a2: Res = TMath.Exp(-u*u/2)
  else: Res = A2*TMath.Power(B2+u,-p2)

  Result = Res*N

  return Result


import argparse
parser =  argparse.ArgumentParser(description='Add Classification BDT weights')
parser.add_argument('-F', '--file', dest='File', required=True, type=str)
parser.add_argument('-L', '--label', dest='Label', required=True, type=str)
parser.add_argument('-C', '--cut', dest='cut', type=str, default='')

opt = parser.parse_args()

dscb_mjj = TF1("dscb_mjj", DSCrystalBall, 60, 180, 7)
#dscb_mjj = TF1("dscb_mjj", DSCrystalBall, 80, 200, 7)
dscb_mjj.SetParameters(1, 120, 20, 1.5, 2, 1.5, 2)
dscb_mjj.SetParLimits(2, 0.00001, 50.)
dscb_mjj.SetParLimits(3, 0.1, 5.)
dscb_mjj.SetParLimits(5, 1.0001, 5.)

dscb_mgg = TF1("dscb_mgg", DSCrystalBall, 115, 135, 7)
dscb_mgg.SetParameters(1, 125, 1.5, 1.5, 2, 1.5, 2)
dscb_mgg.SetParLimits(2, 0.00001, 10.)
dscb_mgg.SetParLimits(3, 1.0001, 5.)
dscb_mgg.SetParLimits(5, 1.0001, 5.)

tf = TFile(opt.File)
tree = tf.Get("bbggSelectionTree")

h_mjj = TH1F("mjj", ";M(jj) [GeV]; Normalized Yields", 40, 60, 180)
#h_mjj = TH1F("mjj", ";M(jj) [GeV]; Normalized Yields", 40, 80, 200)
h_mgg = TH1F("mgg", ";M(#gamma#gamma) [GeV]; Normalized Yields", 20, 115, 135)

h_mjj.Sumw2()
h_mgg.Sumw2()

tree.Draw("dijetCandidate.M()>>mjj", "isSignal "+opt.cut,'goff') 
tree.Draw("diphotonCandidate.M()>>mgg", "isSignal "+opt.cut,'goff') 

mjj_n = h_mjj.Integral()
h_mjj.Scale(1./mjj_n)
mgg_n = h_mgg.Integral()
h_mgg.Scale(1./mgg_n)

mjjfit = h_mjj.Fit("dscb_mjj", "S")
mggfit = h_mgg.Fit("dscb_mgg", "S")

c = TCanvas("c", "c", 800, 600)
gStyle.SetOptStat(0)
gStyle.SetOptFit(1)
h_mjj.Draw("E")
dscb_mjj.Draw("same")
c.SaveAs("mjj"+opt.Label+".pdf")
c.SaveAs("mjj"+opt.Label+".png")

c.Clear()
#h_mgg.GetXaxis().SetRangeUser(115, 135)
#h_mgg.GetXaxis().SetLimits(115, 135)
h_mgg.Draw("E")
dscb_mgg.Draw("same")
c.SaveAs("mgg"+opt.Label+".pdf")
c.SaveAs("mgg"+opt.Label+".png")

oFile = TFile(opt.Label+".root", "RECREATE")
oFile.cd()
h_mjj.Write()
h_mgg.Write()
mjjfit.Write()
mggfit.Write()
dscb_mjj.Write()
dscb_mgg.Write()
oFile.Close()


