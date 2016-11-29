#!/bin/python
from ROOT import *
from copy import deepcopy
from math import *
from array import array

def getRatio(hist1, hist2):
	graph = TGraphAsymmErrors(hist1)
	npoint = 0
	for i in xrange(0, hist1.GetNbinsX()):
		Bin = i+1
		b1 = hist1.GetBinContent(Bin)
		b2 = hist2.GetBinContent(Bin)

		if b1 == 0 or b2 == 0:
			continue

		ratio = b1/b2

		b1sq = b1*b1
		b2sq = b2*b2

		e1sq_up = hist1.GetBinErrorUp(Bin)*hist1.GetBinErrorUp(Bin)
		e2sq_up = hist2.GetBinErrorUp(Bin)*hist2.GetBinErrorUp(Bin)

		e1sq_low = hist1.GetBinErrorLow(Bin)*hist1.GetBinErrorLow(Bin)
		e2sq_low = hist2.GetBinErrorLow(Bin)*hist2.GetBinErrorLow(Bin)

		error_up = sqrt((e1sq_up * b2sq + e2sq_up * b1sq) / (b2sq * b2sq))
		error_low = sqrt((e1sq_low * b2sq + e2sq_low * b1sq) / (b2sq * b2sq))

		graph.SetPoint(npoint, hist1.GetBinCenter(Bin), ratio)
		graph.SetPointError(npoint, 0, 0, error_low, error_up)
		npoint += 1
	graph.Set(npoint)
	return graph

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
	thisStack = stack.Clone("thisStackClone")
	for x in xrange(0, data_nbins):
		LowEdge = data.GetXaxis().GetBinLowEdge(data_nbins - x)
		UpEdge = data.GetXaxis().GetBinUpEdge(x+1)
		n_data = data.GetBinContent(x+1)
		n_bkg = bkg.GetBinContent(x+1)
		bin_center = data.GetBinCenter(x+1)
		bWidth = data.GetBinWidth(x+1)*0.5
#		bWidth = data.GetBinWidth(x+1)*1
		px.append(bin_center)
		pey.append(1)
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
#			erro = (bkg.GetBinContent(x+1)+bkg.GetBinError(x+1))/bkg.GetBinContent(x+1)
			erro = bkg.GetBinError(x+1)
			err = thisStack.GetStack().Last().GetBinError(x+1)
			cont = thisStack.GetStack().Last().GetBinContent(x+1)
			erro = (cont+err)/cont - 1
#			erro = stack.GetStack().Last().GetBinError(x+1)
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
#	pullData = TGraphAsymmErrors(len(pX), pX, pY, pXL, pXH, pYL, pYH) #points
	pullData = getRatio(data, bkg) #TGraphAsymmErrors(len(pX), pX, pY, pXL, pXH, pYL, pYH) #points
	pullError = TGraphErrors(len(pX), pX, pEY, pEX, pE) #bars
	Largest = ceil(largest*1.2)
	pullError.SetMaximum(1.99)
	pullError.SetMinimum(0.01)	
	pullData.SetMaximum(1.99)
	pullData.SetMinimum(0.01)
	pullError.SetTitle("")
	pullData.SetTitle("")
#	pullError.SetFillColorAlpha(kRed, 0.5)
	pullError.SetFillColorAlpha(kGray+2, 0.5)
	pullError.SetFillStyle(1001)
#	pullError.SetLineWidth(0)
#	pullError.SetLineColor(0)
#	pullData.SetMarkerStyle(8)
#	pullData.SetMarkerColorAlpha(0,0)
	pullAll = [pullData, pullError, LowEdge - bWidth*0.18, UpEdge+bWidth*0.18]
	print LowEdge, UpEdge
	return pullAll

def SaveNoPull(data, bkg, fileName):
	c0 = TCanvas("c0", "c0", 1000, 800)
	c0.cd()
	data.Draw()
	bkg.Draw('histsame')
#	c0.SaveAs('/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/noPull_'+str(fileName) + ".pdf")
#	c0.SaveAs('/afs/cern.ch/user/r/rateixei/www/HHBBGG/TestBench/noPull_'+str(fileName) + ".png")
#	c0.Delete()

def SaveWithPull(data, bkg, legend, pullH, pullE, fileName, varName, dirName, lumi, signals, SUM, ControlRegion, hideData, year):
	gStyle.SetHatchesLineWidth(1)
	data.SetStats(0)
	Font = 43
	labelSize = 20
	titleSize = 24

	bkg.Draw('hist')
	bkg.SetTitle("")
	bkg.SetMinimum(0.001)
	print bkg.GetNhists()

	#### Configure thstack
	bkg.GetHistogram().GetXaxis().SetNdivisions(515)


	ratio = 0.2
	epsilon = 0.0
	c1 = TCanvas("c2", "c2", 900, 800)
	SetOwnership(c1,False) #If I don't put this, I get memory leak problems...
	p1 = TPad("pad1","pad1", 0, float(ratio - epsilon), 1, 1)
	SetOwnership(p1,False)
	p1.SetBottomMargin(epsilon)
	p2 = TPad("pad2","pad2",0,0,1,float(ratio*(1-epsilon)) )
	SetOwnership(p2,False)
	p2.SetFillColor(0)
	p2.SetFillStyle(0)
	p2.SetTopMargin(0.0)
	p2.SetBottomMargin(0.35)
	p2.SetGridy()
	p1.cd()
	bkg.Draw('hist')
	bkg.SetTitle("")
	bkg.SetMinimum(0.001)
	print bkg.GetNhists()

	#### Configure thstack
	bkg.GetHistogram().GetXaxis().SetNdivisions(515)
	if(hideData == False):
		bkg.GetHistogram().GetXaxis().SetTitleFont(Font)
		bkg.GetHistogram().GetXaxis().SetTitleSize(0)
		bkg.GetHistogram().GetXaxis().SetLabelFont(Font)
		bkg.GetHistogram().GetXaxis().SetLabelSize(0)
	if(hideData):
		bkg.GetHistogram().GetXaxis().SetTitleFont(Font)
		bkg.GetHistogram().GetXaxis().SetTitleSize(titleSize)
		bkg.GetHistogram().GetXaxis().SetLabelFont(Font)
		bkg.GetHistogram().GetXaxis().SetLabelSize(labelSize)
		bkg.GetHistogram().GetXaxis().SetTitle(varName)
		

	bkg.GetYaxis().SetTitleFont(Font)
	bkg.GetYaxis().SetTitleSize(titleSize)
	bkg.GetYaxis().SetLabelFont(Font)
	bkg.GetYaxis().SetLabelSize(labelSize)
	bkg.GetYaxis().SetTitle("Events")
	if "GeV" in varName:
		nbins = bkg.GetXaxis().GetNbins()
		binslow = bkg.GetXaxis().GetBinLowEdge(1)
		binsup = bkg.GetXaxis().GetBinUpEdge(nbins)
		perbin = (float(binsup) - float(binslow))/float(nbins)
		thisLabel = "Events/("+str(perbin)+" GeV)"
		bkg.GetYaxis().SetTitle(thisLabel)

	bkg.GetYaxis().SetTitleOffset(1.75)

	pullE.GetXaxis().SetLabelFont(Font)
	pullE.GetXaxis().SetLabelSize(labelSize)
	pullE.GetXaxis().SetLabelOffset(0.01)

	pullE.GetXaxis().SetTitle(varName)
	pullE.GetXaxis().SetTitleFont(Font)
	pullE.GetXaxis().SetTitleSize(titleSize)
	pullE.GetXaxis().SetTitleOffset(5)
	pullE.GetXaxis().SetLabelOffset(0.05)
	pullE.GetXaxis().SetNdivisions(515)
	pullE.GetXaxis().SetTickLength(0.15)

	pullE.GetYaxis().SetTitle("Data/MC")
	pullE.GetYaxis().SetTitleFont(Font)
	pullE.GetYaxis().SetTitleSize(18)
	pullE.GetYaxis().SetTitleOffset(1.5)
	pullE.GetYaxis().CenterTitle()
	pullE.GetYaxis().SetLabelFont(Font)
	pullE.GetYaxis().SetLabelSize(15)

	pullE.GetYaxis().SetNdivisions(504)
	pullE.GetXaxis().SetRangeUser(bkg.GetHistogram().GetXaxis().GetXmin(), bkg.GetHistogram().GetXaxis().GetXmax())
	### Fin

	p1.Update()
	gStyle.SetHatchesLineWidth(2)
	SUM.Draw("E2 same")
#	SUM.Draw("a3 same")
	if(hideData == False):
		data.Draw('E same')
	tlatex = TLatex()
	baseSize = 25
	tlatex.SetNDC()
	tlatex.SetTextAngle(0)
	tlatex.SetTextColor(kBlack)
	tlatex.SetTextFont(63)
	tlatex.SetTextAlign(11)
	tlatex.SetTextSize(25)
	tlatex.DrawLatex(0.11, 0.91, "CMS")
	tlatex.SetTextFont(53)
	tlatex.DrawLatex(0.18, 0.91, "Preliminary")
	tlatex.SetTextFont(43)
	tlatex.SetTextSize(23)
#	tlatex.DrawLatex(0.65, 0.91,"L = 2.70 fb^{-1} (13 TeV)")
	Lumi = "L = " + str(lumi) + " pb^{-1} (13 TeV)"#, "+ year + ")"
	if lumi > 1000:
		llumi = float(lumi)/1000.
		Lumi = "L = " + str(llumi) + " fb^{-1} (13 TeV)"#, "+year+")"
#	tlatex.DrawLatex(0.14, 0.82, Lumi)
	tlatex.SetTextAlign(31)
	tlatex.DrawLatex(0.9, 0.91, Lumi)
	tlatex.SetTextAlign(11)

	if ControlRegion != "":
		tlatex.SetTextFont(63)
		tlatex.DrawLatex(0.14, 0.78, ControlRegion)

	for leg in legend:
		leg.Draw('same')

	for h in signals:
		h[0].Draw("same hist")

	p2.cd()
#	pullE.GetYaxis().SetNdivisions(4, False)
	if(hideData==False):
		pullE.Draw("AF2")
	if(hideData):
		pullE.Draw("A")
	Line = TLine(bkg.GetHistogram().GetXaxis().GetXmin(), 1., bkg.GetHistogram().GetXaxis().GetXmax(), 1.)
	Line.SetLineColor(kRed)
	if(hideData==False):
		Line.Draw()
		pullH.Draw("Psame")
	c1.cd()
	p2.Draw()
	p1.Draw()
	c1.Update()
	c1.SaveAs(dirName+"/" + fileName + ".pdf")
	c1.SaveAs(dirName+"/" + fileName + ".png")

	p1.cd()
	p1.SetLogy()
	GeneralMaximus = max(bkg.GetMaximum(), data.GetMaximum())
	if(hideData):
		GeneralMaximus = bkg.GetMaximum()*10
	GeneralMaximus = max(GeneralMaximus, 1E-5)
	print GeneralMaximus, log10(abs(GeneralMaximus))
	GenMax = pow(10, log10(abs(GeneralMaximus))*3.)
	bkg.SetMaximum(GenMax)
	GenMin = bkg.GetMinimum()
	if GenMin == 0:
		GenMin = 1
	bkg.SetMinimum(0.01)
	if(hideData):
		bkg.SetMinimum(0.001)
	p1.Update()
	c1.cd()
	c1.Update()
	c1.SaveAs(dirName+"/LOG_" + fileName + ".pdf")
	c1.SaveAs(dirName+"/LOG_" + fileName + ".png")
#	c1.Delete()

	print "Expected number of events (MC):", SUM.Integral()
	print "Observed number of events (DATA):", data.Integral()

def SavePull(pullH, pullE, LowEdge, UpEdge, dirName):
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
#	ca.SaveAs(dirName + "/pull.pdf")
#	ca.SaveAs(dirName + "/pull.png")
#	ca.Delete()

def MakeLegend(HistList, DataHist, lumi, Signals, SUM):
	newList = []
	newLegs = []
	for h in HistList:
		if h[1] in newLegs:
			continue
		else:
			newList.append(h)
			newLegs.append(h[1])

	nMaxPerBox = (len(newList)+len(Signals)+2)/3
	if (3*nMaxPerBox < len(newList)+len(Signals)+2):
		nMaxPerBox += 1
	lenPerHist = 0.27/float(nMaxPerBox)

	legends = []
#	leg1 = TLegend(0.65, 0.65, 0.71, 0.65+lenPerHist*float(nMaxPerBox))
	leg1 = TLegend(0.68, 0.85-lenPerHist*float(nMaxPerBox), 0.74, 0.89)
#	if (3*nMaxPerBox > len(newList)+len(Signals)+2):
#		nMaxPerBox -= 1
	leg2 = TLegend(0.43, 0.85-lenPerHist*float(nMaxPerBox), 0.49, 0.89)

	leg3 = TLegend(0.13, 0.85-lenPerHist*float(nMaxPerBox), 0.19, 0.89)
	
	legends.append(leg1)
	legends.append(leg2)
	legends.append(leg3)

	data_lumi = " Data (" + str(lumi) + " pb^{-1})"
	if lumi > 1000:
		llumi = float(lumi)/1000.
		data_lumi = " Data (" + str(llumi) + " fb^{-1})"
#	leg1.AddEntry(DataHist, data_lumi, "lep")
#	leg1.AddEntry(SUM, " Stat. Uncertainty", "lf")
	
	allLegs = []
	allLegs.append([DataHist, "Data (" + str(llumi) + " fb^{-1})"])
	allLegs.append([SUM,  "Stat. Uncertainty"])
	allLegs += newList + Signals

	nMaxPerBox = (len(newList)+len(Signals)+2)/3
	if (3*nMaxPerBox < len(newList)+len(Signals)+2):
		nMaxPerBox += 1

	for i,l in enumerate(allLegs):
		Type = 'f'
		if 'Data' in l[1]:
			Type = 'lep'
		if 'Stat' in l[1]:
			Type = 'f'
		iLeg = i//nMaxPerBox
		legends[iLeg].AddEntry(l[0], ' '+l[1], Type)

        textFont = 43
        textSize = 20
	for leg in legends:
		leg.SetFillStyle(0)
		leg.SetLineWidth(0)
		leg.SetBorderSize(0)
		leg.SetTextFont(textFont)
		leg.SetTextSize(textSize)

	return legends
