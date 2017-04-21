from ROOT import *
import sys,math
gROOT.SetBatch(True)

def PoissonianError(Exp, Obs):
  alpha = 1 - 0.6827
  if Obs < 1E-5: Obs = 0
  if Exp < 1E-5: Exp = 0
  if Obs == 0: L = 0
  else: L = Math.gamma_quantile(alpha/2., Obs, 1.)
  U = Math.gamma_quantile_c(alpha/2., Obs+1, 1.)
  diff = Obs - Exp
  if diff > 0: err = Obs-L
  else: err = U-Obs
#  print err, Obs, Exp, L, U, Obs-Exp, (Obs-Exp)/err
  return err

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

def DSCrystalBall2D(xx,par):
  x=float(xx[0])
  y=float(xx[1])
  pararr = [ par[i] for i in range(0,14) ]
#  print pararr
  parx=pararr[:len(pararr)/2]
  pary=pararr[len(pararr)/2:]
#  print parx
  fx=DSCrystalBall([x],parx)
  fy=DSCrystalBall([y],pary)
  return fx*fy


import argparse
parser =  argparse.ArgumentParser(description='Add Classification BDT weights')
parser.add_argument('-F', '--file', dest='File', required=True, type=str)
parser.add_argument('-L', '--label', dest='Label', required=True, type=str)
parser.add_argument('-C', '--cut', dest='cut', type=str, default=' (1>0) ')
parser.add_argument('--fixPars', dest='fpars', type=str, default='')
parser.add_argument('--ggH', dest='ggh', action='store_true', default=False)
opt = parser.parse_args()

myvars = ['mgg', 'mjj']
Pars = {'mgg':{}, 'mjj':{}}
for vv in myvars:
  if vv not in opt.fpars: continue
  fp1 = opt.fpars.split(';')
  for ff in fp1:
    if vv in ff:
      print vv
      mggParsAll = ff.split(":")[1].split(',')
      for mggv in mggParsAll:
        mggp = mggv.split("=")
        Pars[vv][mggp[0]] = mggp[1]
    

if not opt.ggh:
  dscb_mjj = TF1("dscb_mjj", DSCrystalBall, 60, 180, 7)
  #dscb_mjj = TF1("dscb_mjj", DSCrystalBall, 80, 200, 7)
  dscb_mjj.SetParameters(1, 120, 15, 1.5, 2, 1.5, 2)
  dscb_mjj.SetParLimits(2, 5., 500.)
  dscb_mjj.SetParLimits(3, 0.5, 5.)
  dscb_mjj.SetParLimits(5, 1.0001, 5.)
if opt.ggh:
  dscb_mjj = TF1("dscb_mjj", "[0] + [1]*x + [2]*x*x", 60, 180)
  dscb_mjj.SetParameters(8.0, -0.08, 0.0001)

dscb_mgg = TF1("dscb_mgg", DSCrystalBall, 115, 135, 7)
dscb_mgg.SetParameters(1, 125, 1.5, 1.5, 2, 1.5, 2)
dscb_mgg.SetParLimits(2, 0.00001, 10.)
dscb_mgg.SetParLimits(3, 1.0001, 5.)
dscb_mgg.SetParLimits(5, 1.0001, 5.)

if not opt.ggh:
  dscb_mjj_mgg = TF2("dscb_mjj_mgg", DSCrystalBall2D, 60,180,115,135,14)
  dscb_mjj_mgg.SetParameter(0,0.5)
  dscb_mjj_mgg.SetParameter(1, 120)
  dscb_mjj_mgg.SetParameter(2, 15)
  dscb_mjj_mgg.SetParameter(3, 1.5)
  dscb_mjj_mgg.SetParameter(4, 2)
  dscb_mjj_mgg.SetParameter(5, 1.5)
  dscb_mjj_mgg.SetParameter(6, 2)
  #dscb_mjj_mgg.SetParameter(7, 0.5)
  dscb_mjj_mgg.FixParameter(7, 1.0)
  dscb_mjj_mgg.SetParameter(8, 125)
  dscb_mjj_mgg.SetParameter(9, 8.98766e-01)
  dscb_mjj_mgg.SetParameter(10, 1.5)
  dscb_mjj_mgg.SetParameter(11, 2)
  dscb_mjj_mgg.SetParameter(12, 1.5)
  dscb_mjj_mgg.SetParameter(13, 2)
  dscb_mjj_mgg.SetParLimits(2, 5., 50.)
  dscb_mjj_mgg.SetParLimits(3, 0.1, 5.)
  dscb_mjj_mgg.SetParLimits(5, 1.0001, 5.)
  dscb_mjj_mgg.SetParLimits(9, 8.58766e-01, 9.166e-01)
  dscb_mjj_mgg.SetParLimits(10, 1.0001, 5.)
  dscb_mjj_mgg.SetParLimits(12, 1.0001, 5.)

for par in Pars['mjj']:
  print "Fixing mgg parameter",par,"to value",Pars['mjj'][par]
  dscb_mjj.FixParameter(int(par), float(Pars['mjj'][par]))
  if not opt.ghh: dscb_mjj_mgg.FixParameter(int(par), float(Pars['mjj'][par]))

for par in Pars['mgg']:
  print "Fixing mgg parameter",par,"to value",Pars['mgg'][par]
  dscb_mgg.FixParameter(int(par), float(Pars['mgg'][par]))
  if not opt.ghh: dscb_mjj_mgg.FixParameter(int(par)+7, float(Pars['mjj'][par]))


tf = TFile(opt.File)
tree = tf.Get("bbggSelectionTree")

h_mjj = TH1F("mjj", ";M(jj) [GeV]; Normalized Yields", 40, 60, 180)
#h_mjj = TH1F("mjj", ";M(jj) [GeV]; Normalized Yields", 40, 80, 200)
h_mgg = TH1F("mgg", ";M(#gamma#gamma) [GeV]; Normalized Yields", 20, 115, 135)
h_mjj_mgg = TH2F("mjj_mgg", ";M(jj) [GeV]; M(#gamma#gamma) [GeV]", 40,60,180,20,115,135)
r_mjj_mgg = TH2F("r_mjj_mgg", "Residues;M(jj) [GeV]; M(#gamma#gamma) [GeV]", 40,60,180,20,115,135)

h_mjj.Sumw2()
h_mgg.Sumw2()
h_mjj_mgg.Sumw2()

tree.Draw("dijetCandidate.M()>>mjj", "isSignal && "+opt.cut,'goff') 
tree.Draw("diphotonCandidate.M()>>mgg", "isSignal && "+opt.cut,'goff') 
tree.Draw("diphotonCandidate.M():dijetCandidate.M()>>mjj_mgg", "isSignal && "+opt.cut, 'goff')

mjj_n = h_mjj.Integral()
h_mjj.Scale(1./mjj_n)
mgg_n = h_mgg.Integral()
h_mgg.Scale(1./mgg_n)
mjj_mgg_n = h_mjj_mgg.Integral()
h_mjj_mgg.Scale(1./mjj_mgg_n)

mjjfit = h_mjj.Fit("dscb_mjj", "S")
mggfit = h_mgg.Fit("dscb_mgg", "S")
mjjmggfit = TObject()#h_mjj_mgg.Fit("dscb_mjj_mgg", "S")
if not opt.ggh:
  for ip in range(0,7):
    dscb_mjj_mgg.SetParameter(ip, dscb_mjj.GetParameter(ip))
    dscb_mjj_mgg.SetParameter(ip+7, dscb_mgg.GetParameter(ip))

for xx in range(1, 41):
 for yy in range(1,21):
  mc_c = h_mjj_mgg.GetBinContent(xx,yy)
  xval = h_mjj_mgg.GetXaxis().GetBinCenter(xx)
  yval = h_mjj_mgg.GetYaxis().GetBinCenter(yy)
#  print xval,yval
  if not opt.ggh:
    pd_c = dscb_mjj_mgg.Eval(xval,yval)
    err = abs(PoissonianError(pd_c,mc_c))
    resi = abs(pd_c - mc_c)/err
    r_mjj_mgg.SetBinContent(xx,yy,resi)

c = TCanvas("c", "c", 800, 600)
gStyle.SetOptStat(0)
gStyle.SetOptFit(1)
h_mjj.SetMinimum(0)
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

if not opt.ggh:
  dscb_mjj_mgg_projy = TF12("dscb_mjj_mgg_projy",dscb_mjj_mgg,125,"y")
  dscb_mjj_mgg_projx = TF12("dscb_mjj_mgg_projy",dscb_mjj_mgg,125,"x")

oFile = TFile(opt.Label+".root", "RECREATE")
oFile.cd()
h_mjj.Write()
h_mgg.Write()
h_mjj_mgg.Write()
r_mjj_mgg.Write()
mjjfit.Write()
mggfit.Write()
mjjmggfit.Write()
dscb_mjj.Write()
dscb_mgg.Write()
if not opt.ggh:
  dscb_mjj_mgg.Write()
  dscb_mjj_mgg_projy.Write()
  dscb_mjj_mgg_projx.Write()
oFile.Close()


