from pullClass import *
from ROOT import *
import json, os
import shutil
from configs import *
import resource

gStyle.SetOptStat(0)
TGaxis.SetMaxDigits(3)

#yellow 4aa6db 641555

#SetMemoryPolicy( kMemoryStrict )

dummyTFile = TFile("dummy.root", "RECREATE")

NiceBlue = '#51A7F9'
NiceBlueDark = '#2175E0'
NiceGreen = '#6FBF41'
NiceGreen2 = '#35DC3D'
NiceYellow = '#FBE12A'
NiceYellow2 = '#FEEA01'
NiceOrange = '#F0951A'
NiceRed = '#FA4912'
NicePurple = '#885BB2'
NicePaleYellow = '#FFFF66'
NiceMidnight = '#000080'
NiceTangerine = '#FF8000'
cNiceBlue = TColor.GetColor('#51A7F9')
cNiceBlueDark = TColor.GetColor('#2175E0')
cNiceGreen = TColor.GetColor('#6FBF41')
cNiceGreen2 = TColor.GetColor(NiceGreen2)
cNiceGreenDark = TColor.GetColor('#008040')
cNiceYellow = TColor.GetColor('#FBE12A')
cNiceYellow2 = TColor.GetColor(NiceYellow2)
cNiceOrange = TColor.GetColor('#F0951A')
cNiceRed = TColor.GetColor('#FA4912')
cNicePurple = TColor.GetColor('#885BB2')
cNicePaleYellow = TColor.GetColor('#FFFF66')
cNiceMidnight = TColor.GetColor('#000080')
cNiceTangerine = TColor.GetColor('#FF8000')

myCols = [cNiceMidnight, cNiceGreenDark, kRed, kMagenta+1]

if not os.path.exists(dirName):
        print dirName, "doesn't exist, creating it..."
        os.makedirs(dirName)
#        shutil.copy2(dirPrefix + "index.php", dirName+"/index.php")
        if os.path.exists(dirName):
                print dirName, "now exists!"

datasets = json.load(data_file)

Trees = {}

CutSignal = Cut
if doBlind == True:
	Cut += " && !((diphotonCandidate.M() > 115 && diphotonCandidate.M() < 135))"
if isPhoCR == True:
	Cut += " && (isPhotonCR == 1)"
	CutSignal += " && (isPhotonCR == 1)"
if isPhoCR == False:
	Cut += " && (isSignal == 1)"
	CutSignal += " && (isSignal == 1)"
if doSignalRegion == True:
	CutSignal += " && ( leadingJet_bDis > "+ str(BTAG) + " || subleadingJet_bDis > "+str(BTAG)+" ) "
	Cut += " && ( leadingJet_bDis > "+ str(BTAG) + " || subleadingJet_bDis > "+ str(BTAG) + " ) "
if doJetCR == True:
	Cut += " && leadingJet_bDis < "+ str(BTAG) + " && subleadingJet_bDis < "+ str(BTAG) + " "
weightedcut = ""
weightedcut += "( genTotalWeight )*"
if (doPUweight):
	weightedcut += "( puweight )*("+Cut+")"
	weightedcut_vbf =  weightedcut.replace("puweight", "1>0")

weightedCut = TCut(weightedcut)
cut_data = TCut(Cut)
cut_signal = TCut(weightedcut.replace("!((diphotonCandidate.M() > 115 && diphotonCandidate.M() < 135))", "(1>0)"))
cut_signal_vbf = TCut(weightedcut_vbf.replace("!((diphotonCandidate.M() > 115 && diphotonCandidate.M() < 135))", "(1>0)"))

for plot in plots:
    Histos = []
    variable = plot[1]
    varName = plot[2]
    thisStack = 0
    thisHist = 0
    thisStack = myStack('test'+plot[0], varName, varName, dirName, lumi)
    if hideData == True:
        thisStack.hideData()
    if isPhoCR == 1:
        thisStack.makePhoCR()
    if doJetCR == 1:
        thisStack.makeJetCR()

    thisStack.setYear(year)

    modelHist = TH1F(plot[0]+"_hist", "", plot[3], plot[4], plot[5])

    backgroundHists = []
    for background in datasets["background"]:
        if not addbbH and 'bbH' in background: continue
        if not dyjets and "DYJ" in background: continue
        if "QCD" in background: continue
        print background
        thisName = plot[0]+"_hist"+"_"+background
        thisHist = modelHist.Clone(thisName)
#        thisHist.Clear()
#        thisHist.Sumw2()
        thisHist.SetLineColor(TColor.GetColor(datasets["background"][background]["color"]))
        thisHist.SetFillColor(TColor.GetColor(datasets["background"][background]["color"]))
        for i,fi in enumerate(datasets["background"][background]["files"]):
            print fi
            thisTreeLoc = fi["file"]
            if thisTreeLoc not in Trees:
                Trees[thisTreeLoc] = TChain("bbggSelectionTree")
                Trees[thisTreeLoc].AddFile(bkgLocation+thisTreeLoc)
                SetOwnership( Trees[thisTreeLoc], True )
            locName = thisName+str(i)
            locHist = thisHist.Clone(locName)
            thisWeightedCut = weightedCut
            if "QCD" in thisTreeLoc:
                thisWeightedCut = TCut(weightedcut.replace("isSignal == 1", "isSignal == 0"))
            Trees[thisTreeLoc].Draw(plot[1]+">>"+locName, thisWeightedCut)
            locHist.Scale(MCSF*lumi*fi["xsec"]*fi["sfactor"]/fi["weight"])
            thisHist.Add(locHist)
            Histos.append(locHist)
	    print thisName, " ",  thisHist.Integral()
#            thisFile = TFile(plot[0]+"_"+fi["file"], "RECREATE")
#            thisFile.cd()
#            thisHist.Write()
#            thisFile.Close()
#            dummyTFile.cd()
            
        backgroundHists.append([thisHist, datasets["background"][background]["legend"], datasets["background"][background]["position"]])
        del thisHist
    OrderedBackgrounds = sorted(backgroundHists, key=lambda x: x[2], reverse=True)
    print OrderedBackgrounds
    for background in OrderedBackgrounds: thisStack.addHist(background[0], background[1], background[2])
    
    for isi, signal in enumerate(datasets["signal"]):
        thisName = plot[0]+"_Signal_"+"_"+signal['name']
#        modelHist.Clear()
#        thisHist = 0
 #       thisHist = modelHist
#        thisHist.SetName(thisName)
#        thisHist.Sumw2()
        thisHist = modelHist.Clone(thisName)
#        thisHist.SetLineColor(signal["color"])
        thisHist.SetLineColor(myCols[isi])
        thisHist.SetLineWidth(3)
        thisHist.SetLineStyle(signal["style"])
        thisTreeLoc = signal["file"]
        if thisTreeLoc not in Trees:
            Trees[thisTreeLoc] = TChain("bbggSelectionTree")
            Trees[thisTreeLoc].AddFile(signalLocation+thisTreeLoc)
            SetOwnership( Trees[thisTreeLoc], True )
	if "VBF" in thisName:
		Trees[thisTreeLoc].Draw(plot[1]+">>"+thisName, cut_signal_vbf)
	else:
		Trees[thisTreeLoc].Draw(plot[1]+">>"+thisName, cut_signal)

	print thisName, " ",  thisHist.Integral()
        thisHist.Scale(signalFactor*lumi*signal["xsec"]*signal["sfactor"]/signal["weight"])
	print thisName, " ",  thisHist.Integral()
        Histos.append(thisHist)
        thisStack.addSignal(thisHist, signal["legend"], lumi*signal["xsec"]*signal["sfactor"]/signal["weight"])
        del thisHist
    
    
    print datasets['data']
    dataName = plot[0]+"_hist"+"_data"
    modelHist.Clear()
    dataHist = 0
#    dataHist = modelHist
    dataHist = modelHist.Clone(dataName)
#    dataHist.Sumw2()
    if datasets['data'] not in Trees:
        Trees[datasets['data']] = TChain("bbggSelectionTree")
        Trees[datasets['data']].AddFile(dataLocation+datasets['data'])
        SetOwnership( Trees[datasets['data']], True )
    Trees[datasets['data']].Draw(plot[1]+">>"+dataName, cut_data)
    dataHist.SetMarkerStyle(20)
    dataHist.SetMarkerSize(0.8)
    dataHist.SetMarkerColor(1)
    dataHist.SetLineColor(1)
    dataHist.SetLineWidth(2)
    dataHist.SetBinErrorOption(TH1.kPoisson)
    thisStack.addData(dataHist, "Data")
#    thisFile = TFile(plot[0]+"_data.root", "RECREATE")
#    thisFile.cd()
#    dataHist.Write()
#    thisFile.Close()
#    dummyTFile.cd()

    thisStack.drawStack(prefix + plot[0])
    
    del thisStack 


dummyTFile.Close()
os.system("rm dummy.root")

#for tt in Trees:
#    Trees[tt].IsA().Destructor( Trees[tt] )

#    gROOT.GetListOfCanvases().Delete()
#    print Trees
#    print 'Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
#    gROOT.GetListOfCleanups().Delete()
    

