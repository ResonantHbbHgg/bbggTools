#!/bin/python
from ROOT import *
from copy import deepcopy
from math import *
from array import array
from pullUtils import *


class myStack:
	myHistograms = []
	tStack = ''
	mySumHists = ''
	myData = ''
	dataLegend = ''
	varName = ''
	dirName = ''
	lumi = ''
	isPhoCR = 0
	year = ""
	def __init__(self, name, title, varName, dirName, lumi):
		self.myHistograms = []
		self.mySignals = []
		self.tStack = ''
		self.mySumHists = ''
		self.myData = ''
		self.dataLegend = ''
		self.varName = ''
		self.tStack = THStack(str(name), str(title))
		SetOwnership(self.tStack, False)
		self.varName = varName
		self.dirName = dirName
		self.lumi = lumi
		self.SUM = ''
		self.isPhoCR = 0
		self.isJetCR = 0
		self.legends = []
		self.hideData_ = 0
	def hideData(self):
		self.hideData_ = 1
	def setYear(self, Year):
		self.year = Year
	def makeJetCR(self):
		self.isJetCR = 1
	def makePhoCR(self):
		self.isPhoCR = 1
	def addHist(self, hist, legend, norm):
		self.myHistograms.append([deepcopy(hist), str(legend), norm])
		SetOwnership( self.myHistograms[ len(self.myHistograms)-1 ][0], False )
		self.legends.append(legend)
	def addSignal(self, hist, legend, norm):
		self.mySignals.append([deepcopy(hist), str(legend), norm])
		SetOwnership( self.mySignals[len(self.mySignals)-1][0], False)
	def addData(self, hist, legend):
		self.myData = deepcopy(hist)
		SetOwnership( self.myData, False)
		self.dataLegend = legend
	def alreadyHas(self, legend):
		if legend in self.legends:
			return 1
		else:
			return 0
	def drawStack(self, fileName):
		if len(self.myHistograms) < 1:
			print 'Your list of histograms is empty!'
			return 0
		if self.myData == '':
			print "You haven't added a data histogram!"
			return 0
		mySumHists = self.myHistograms[0][0].Clone("mySumHists")
		mySumHists.Reset()
		for Hist in self.myHistograms:
			hist = Hist[0]
			self.tStack.Add(hist, "hist")
			mySumHists.Add(hist,1)
		generalMaximus = max( self.myData.GetMaximum(), self.tStack.GetMaximum() )
		self.tStack.SetMaximum(generalMaximus*1.5)
		pullHandE = doPull(mySumHists, self.myData, self.tStack)
		pullH = pullHandE[0]
		pullE = pullHandE[1]
		LowEdge = pullHandE[2]
		UpEdge = pullHandE[3]
		self.SUM = self.tStack.GetStack().Last().Clone("SUM")
		self.SUM.SetLineWidth(0)
		self.SUM.SetFillColorAlpha(kGray+2, 0.5)
		self.SUM.SetMarkerColorAlpha(0,0)
#		self.SUM.SetLineColorAlpha(kGray+2,0.5)
		self.SUM.SetFillStyle(1001)
		legend = MakeLegend(self.myHistograms, self.myData, self.lumi, self.mySignals, self.SUM)
		ControlRegion = ""
		if self.isPhoCR == 1:
			ControlRegion = "Fake Photon CR"
		if self.isJetCR == 1:
			ControlRegion = "Light Jets CR"
		print self.tStack.GetNhists()
		SaveWithPull(self.myData, self.tStack, legend, pullH, pullE, fileName, self.varName, self.dirName, self.lumi, self.mySignals, self.SUM, ControlRegion, self.hideData_, self.year)
#		gROOT.EndOfProcessCleanups()

