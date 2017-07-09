from ROOT import TCanvas
from array import array
from flashgg.bbggTools.MyCMSStyle import *

def DrawNoPull(data, bkg, legend, fileName, varName, dirName, lumi, signals, SUM, ControlRegion, hideData, year):
  frame = data.Clone(data.GetName()+'_frame')
  frame.Reset()
  frame.GetXaxis().SetTitle(varName)
  myMax = max(data.GetMaximum(), SUM.GetMaximum())
  frame.SetMaximum(myMax*1.7)
  frame.SetMinimum(0.05)
#  frame.GetYaxis().SetMaxDigits(2)

  frame.GetYaxis().SetTitle("Events")
  if "GeV" in varName:
    nbins = frame.GetXaxis().GetNbins()
    binslow = frame.GetXaxis().GetBinLowEdge(1)
    binsup = frame.GetXaxis().GetBinUpEdge(nbins)
    perbin = (float(binsup) - float(binslow))/float(nbins)
    thisLabel = "Events/("+str(perbin)+" GeV)"
    frame.GetYaxis().SetTitle(thisLabel)
    if binsup > 999:
      frame.GetXaxis().SetLimits(binslow, 999)

  SetAxisTextSizes(frame)
  SetGeneralStyle()
  tc = TCanvas('tc', 'tc', 800, 700)
  SetPadStyle(tc)
  frame.Draw()
  bkg.Draw("hist same")
  SUM.Draw("E2 same")
  if hideData==False: data.Draw("E same")

  for si in signals:
#    print si
    si[0].Draw("hist  same")

  tc.Update()
  tc.RedrawAxis()

  for ll in legend:
    ll.Draw("same")

  Lumi = str(lumi/1000.)

  l1 = DrawCMSLabels(tc, Lumi, 0, 0.08)
  tc.Update()

  tc.SaveAs(dirName+"/NP_" + fileName + ".pdf")
  tc.SaveAs(dirName+"/NP_" + fileName + ".png")

  for ll in l1:
    ll.Delete()

  frame.SetMaximum( myMax*100000 )
  tc.SetLogy()
  tc.Update()

  l2 = DrawCMSLabels(tc, Lumi, 0, 0)
  
  tc.SaveAs(dirName+"/NP_LOG_" + fileName + ".pdf")
  tc.SaveAs(dirName+"/NP_LOG_" + fileName + ".png")
