#!/bin/python
from ROOT import *
from copy import deepcopy
from math import *
from array import array

#
# The following function was copied (with permission) from the Latino's GitHub:
# https://github.com/latinos/LatinoAnalysis/blob/master/ShapeAnalysis/scripts/mkPlot.py#L437
# Thanks Andrea!
#

# _____________________________________________________________________________
# --- poissonian error bayesian 1sigma band
#                                     1/0   1/0
def GetPoissError(numberEvents, down, up):
	alpha = (1-0.6827)
	L = 0
	if numberEvents!=0 : 
		L = Math.gamma_quantile (alpha/2,numberEvents,1.)
	U = 0
	if numberEvents==0 :
		U = Math.gamma_quantile_c (alpha,numberEvents+1,1.) 
	else :
		U = Math.gamma_quantile_c (alpha/2,numberEvents+1,1.)
         
       # the error
	L = numberEvents - L
	if numberEvents > 0 :
		U = U - numberEvents
       #else :
         #U = 1.14 # --> bayesian interval Poisson with 0 events observed
         #1.14790758039 from 10 lines above
         
	if up and not down :
		return U
	if down and not up :
		return L
	if up and down :
		return [L,U]

def doPull(bkg, data, stack):
	data_nbins = data.GetNbinsX()
	bkg_nbins = bkg.GetNbinsX()
	if data_nbins != bkg_nbins:
		print "Your data and background histograms don't match!"
		return 0
	px = []
	pxh = []
	pxl = []
	pyh = []
	pyl = []
	py = []
	pe = []
	pex = []
	pey = []
	largest = 0
	UpEdge = 0
	LowEdge = 0
	bWidth = 0
	for x in xrange(0, data_nbins):
		LowEdge = data.GetXaxis().GetBinLowEdge(data_nbins - x)
		UpEdge = data.GetXaxis().GetBinUpEdge(x+1)
		n_data = data.GetBinContent(x+1)
		n_bkg = bkg.GetBinContent(x+1)
		bin_center = data.GetBinCenter(x+1)
		bWidth = data.GetBinWidth(x+1)*0.5
		px.append(bin_center)
		pey.append(0)
		pex.append(bWidth)
		pxh.append(0)
		pxl.append(0)
		yErrors = GetPoissError(n_data, 1, 1)
		pyl.append(yErrors[0])
		pyh.append(yErrors[1])
		if n_bkg == 0 :
			if n_data == 0:
				py.append(1)
				pe.append(0)
			else:
				py.append(5)
				pe.append(5)
		else:
			sigma = (float(n_data))/float((n_bkg))
#			erro = sigma*sqrt(1./float(n_data) + 1./float(n_bkg))
#			erro = sqrt( float(n_bkg) )
			erro = bkg.GetBinError(x+1)
			py.append(sigma)
			pe.append(erro)
			print "bin center:", bin_center, "bin width",bWidth, "n data:",n_data, "nbkg:",n_bkg, "sigma:",sigma, "erro:",erro
			if abs(sigma) > largest:
				largest = abs(sigma)
			if abs(erro) > largest:
				largest = abs(erro)
	pX = array('d', px)
	pXH = array('d', pxh)
	pXL = array('d', pxl)
	pY = array('d', py)
	pYH = array('d', pyh)
	pYL = array('d', pyl)
	pE = array('d', pe)
	pEY = array('d', pey)
	pEX = array('d', pex)
	pullData = TGraphAsymmErrors(len(pX), pX, pY, pXL, pXH, pYL, pYH) #points
	pullError = TGraphErrors(len(pX), pX, pEY, pEX, pE) #bars
	pullError.SetMaximum(largest*1.2)
	pullError.SetMinimum(largest*(-1.2))
	pullData.SetMaximum(largest*1.2)
	pullData.SetMinimum(largest*(-1.2))
	if largest == 0:
		pullError.SetMaximum(1.)
		pullError.SetMinimum(-1.)	
		pullData.SetMaximum(1.)
		pullData.SetMinimum(-1.)
	pullError.SetTitle("")
	pullData.SetTitle("")
	pullError.SetFillColor(kGray)
	pullError.SetFillStyle(3001)
	pullError.SetLineWidth(0)
	pullError.SetLineColor(0)
	pullData.SetMarkerStyle(8)
	pullAll = [pullData, pullError, LowEdge - bWidth*0.18, UpEdge+bWidth*0.18]
	print LowEdge, UpEdge
	return pullAll

def SaveNoPull(data, bkg, fileName):
	c0 = TCanvas("c0", "c0", 1000, 800)
	c0.cd()
	data.Draw()
	bkg.Draw('same')
	c0.SaveAs('/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/noPull_'+str(fileName) + ".pdf")
	c0.SaveAs('/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/noPull_'+str(fileName) + ".png")

def SaveWithPull(data, bkg, legend, pullH, pullE, fileName, varName):
	data.SetStats(0)
	Font = 43
	labelSize = 20
	titleSize = 24
	
	bkg.GetXaxis().SetTitleFont(Font)
	bkg.GetXaxis().SetTitleSize(0)
	bkg.GetXaxis().SetLabelFont(Font)
	bkg.GetXaxis().SetLabelSize(0)

	bkg.GetYaxis().SetTitleFont(Font)
	bkg.GetYaxis().SetTitleSize(titleSize)
	bkg.GetYaxis().SetLabelFont(Font)
	bkg.GetYaxis().SetLabelSize(labelSize)
	bkg.GetYaxis().SetTitle("Events")
	bkg.GetYaxis().SetTitleOffset(1.8)

#	gStyle.SetPadTickX(1)
#	gStyle.SetPadTickY(1)

	pullE.GetXaxis().SetLabelFont(Font)
	pullE.GetXaxis().SetLabelSize(labelSize)
	pullE.GetXaxis().SetLabelOffset(0.01)

	pullE.GetXaxis().SetTitle(varName)
	pullE.GetXaxis().SetTitleFont(Font)
	pullE.GetXaxis().SetTitleSize(titleSize)
	pullE.GetXaxis().SetTitleOffset(5.0)

	pullE.GetYaxis().SetTitle("Data/Background")
	pullE.GetYaxis().SetTitleFont(Font)
	pullE.GetYaxis().SetTitleSize(10)
	pullE.GetYaxis().SetTitleOffset(1.8)
	pullE.GetYaxis().CenterTitle()
	pullE.GetYaxis().SetLabelFont(Font)
	pullE.GetYaxis().SetLabelSize(8)

	pullE.GetYaxis().SetNdivisions(6)


	ratio = 0.2
	epsilon = 0.00
	c1 = TCanvas("c1", "c1", 900, 800)
	SetOwnership(c1,False) #If I don't put this, I get memory leak problems...
	p1 = TPad("pad1","pad1", 0, float(ratio - epsilon), 1, 1)
	SetOwnership(p1,False)
	p1.SetBottomMargin(epsilon)
	p2 = TPad("pad2","pad2",0,0,1,float(ratio*(1-epsilon)) )
	SetOwnership(p2,False)
	p2.SetFillColor(0)
	p2.SetFillStyle(0)
	p2.SetTopMargin(0.06)
	p2.SetBottomMargin(0.35)
	p2.SetGridy()
	p1.cd()
	bkg.Draw('hist')
	data.Draw('Esame')
	tlatex = TLatex()
	baseSize = 25
	tlatex.SetNDC()
	tlatex.SetTextAngle(0)
	tlatex.SetTextColor(kBlack)
	tlatex.SetTextFont(63)
	tlatex.SetTextAlign(31)
	tlatex.SetTextSize(baseSize)
	tlatex.DrawLatex(0.215, 0.91, "CMS")
	tlatex.SetTextFont(53)
	tlatex.SetTextSize(baseSize - 5)
	tlatex.DrawLatex(0.325, 0.91, "Preliminary")
	tlatex.SetTextFont(43)
	tlatex.SetTextSize(baseSize - 5)
	tlatex.DrawLatex(0.88, 0.91, "L = 150 pb^{-1} (13 TeV, 25 ns)")

	for leg in legend:
		leg.Draw('same')
	p2.cd()
	pullE.Draw("A5")
	pullH.Draw("Psame")
	c1.cd()
	p2.Draw()
	p1.Draw()
	c1.Update()
	c1.SaveAs("/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/" + fileName + ".pdf")
	c1.SaveAs("/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/" + fileName + ".png")

	p1.cd()
	p1.SetLogy()
	p1.Update()
	c1.cd()
	c1.Update()
	c1.SaveAs("/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/LOG_" + fileName + ".pdf")
	c1.SaveAs("/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/LOG_" + fileName + ".png")

def SavePull(pullH, pullE, LowEdge, UpEdge):
	ca = TCanvas("ca", "ca", 1000, 800)
	ca.cd()
	pullE.Draw("A2")
	pullE.GetXaxis().SetLimits(LowEdge,UpEdge)
	pullE.GetXaxis().SetTitle
	ca.Update()
	pullH.Draw("Psame")
#		pullH.Draw("A*")
#		pullH.GetXaxis().SetRangeUser(LowEdge,UpEdge)
#		pullE.Draw("same2")
	ca.SaveAs("/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/pull.pdf")
	ca.SaveAs("/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/pull.png")
	

def MakeLegend(HistList, DataHist):
	nBoxes = 3
	nLegPerBoxMin = 4
	legends = []
	y0 = 0.65
	y1 = 0.89
	ystep = (y1 - y0)/4
	xi = 0.08
	xf = 0.81
	xstep = (xf - xi)/nBoxes
	xgap = xstep/2 + 0.03

	totalHists = len(HistList) + 1

	nPerGroup = (totalHists)//nBoxes
	if (totalHists)%nBoxes > 0:
		nPerGroup += 1

	if nPerGroup < nLegPerBoxMin:
		nPerGroup = nLegPerBoxMin

	nLeg3 = min( totalHists, nPerGroup)
	nLeg2 = min( totalHists - nLeg3, nPerGroup)
	print "nLeg2: ", totalHists, nLeg3
	nLeg1 = totalHists - (nLeg2 + nLeg3)

	yLeg3 = max( y1 - nLeg3*ystep, 0.65)
	yLeg2 = max( y1 - nLeg2*ystep, 0.65)
	print "yLeg2 = max( 0.9 - ", nLeg2*ystep, ", 0.65)"
	yLeg1 = max( y1 - nLeg1*ystep, 0.65)
	print yLeg3, yLeg2, yLeg1

        leg1 = TLegend( xi + xstep*0 + xgap*1, yLeg1, xi + xstep*1 - xgap*1, y1)
        leg2 = TLegend( xi + xstep*1 + xgap*1, yLeg2, xi + xstep*2 - xgap*1, y1)
        leg3 = TLegend( xi + xstep*2 + xgap*1, yLeg3, xi + xstep*3 - xgap*1, y1)
#	textFont = 63
#	textSize = 18
#	leg1.SetTextFont(textFont)
#	leg1.SetTextSize(textSize)
#	leg2.SetTextFont(textFont)
#	leg2.SetTextSize(textSize)
#	leg3.SetTextFont(textFont)
#	leg3.SetTextSize(textSize)

	leg3.AddEntry(DataHist, " Data (40/pb)", "lep")
	
	for i,hist in enumerate(HistList):
		if i+1 < nPerGroup:
			leg3.AddEntry(hist[0], " " +hist[1], "f")
		if i+1 >= nPerGroup and i+1 < nPerGroup*2:
			leg2.AddEntry(hist[0], " " +hist[1], "f")
		if i+1 >= nPerGroup*2 and i+1 < nPerGroup*3:
			leg1.AddEntry(hist[0], " " +hist[1], "f")
		
	legends.append(leg3)
	legends.append(leg2)
	legends.append(leg1)
        textFont = 63
        textSize = 18
	for leg in legends:
		leg.SetFillStyle(0)
		leg.SetLineWidth(0)
		leg.SetBorderSize(0)
		leg.SetTextFont(textFont)
		leg.SetTextSize(textSize)

	return legends
