#####################################
## input: aabb minitrees to 13 benchmarks
## output: one validation sample SM + 12 benchmarks + 72 outliers + 52 lambda-only scan
###################################
#!/usr/bin/env python
import os, sys, time,math
import ROOT
from ROOT import TH1D
from array import array
from ROOT import TLatex,TPad,TList,TH1,TH1F,TH2F,TH1D,TH2D,TFile,TTree,TCanvas,TLegend,SetOwnership,gDirectory,TObject,gStyle,gROOT,TLorentzVector,TGraph,TMultiGraph,TColor,TAttMarker,TLine,TDatime,TGaxis,TF1,THStack,TAxis,TStyle,TPaveText,TAttFill
import math
import bisect
from optparse import OptionParser
import numpy as np
from ROOT import TTree, TFile, AddressOf, gROOT



#inputLM = "NRDir_LM_350/"
#inputHM = "NRDir_HM_350/"
inputLM = "/afs/cern.ch/work/z/zghiche/public/ForXanda/NRDir_LM_350/"
inputHM = "/afs/cern.ch/work/z/zghiche/public/ForXanda/NRDir_HM_350/"
os.system('mkdir V3benchmarks')
os.system('mkdir V3outliers')
os.system('mkdir lambdaonly')

file = "LT_output_GluGluToHHTo2B2G_node_"
fileH="_13TeV-madgraph_HIGHMASS.root"
fileL="_13TeV-madgraph_LOWMASS.root"

nclu=12
nout=6
Nlambda=52

Nvarbinmass=91 
Nvarbincosth=11
#a = np.ones((3,2))        # a 2D array with 3 rows, 2 columns, filled with ones
#Matrix1 = [[[None]*Nvarbincosth]*Nvarbinmass]*nclu # of the benchmark
#Matrix2high = [[None]*Nvarbincosth]*Nvarbinmass # of the fullsim set
#Matrix2low = [[None]*Nvarbincosth]*Nvarbinmass # of the fullsim set

Matrix1 = np.zeros((nclu,Nvarbinmass,Nvarbincosth))  # of the benchmark
Matrix1out = np.zeros((nclu,nout,Nvarbinmass,Nvarbincosth))  # of the ouliers
Matrix2out = np.zeros((Nlambda,Nvarbinmass,Nvarbincosth))  # of the lambda-only

Matrix2high = np.zeros((Nvarbinmass,Nvarbincosth)) # of the fullsim set
Matrix2low = np.zeros((Nvarbinmass,Nvarbincosth))# of the fullsim set

Matrix1validationhigh = np.zeros((Nvarbinmass,Nvarbincosth))# of the fullsim set
Matrix1validationlow = np.zeros((Nvarbinmass,Nvarbincosth))# of the fullsim set
      
######################################
### read the reweiting maps to the benchmarks
lambdaonly = "Distros_lambdaOnly_5p_50000ev_52sam_13TeV_JHEP_50K.root"
envelope = "Distros_envelope_5p_20000ev_6sam_13TeV.root"
benchmarks = "Distros_benchmarks_5p_500000ev_12sam_13TeV_JHEPv3.root" # "Distros_5p_500000ev_12sam_13TeV_JHEP_500K.root"

######################################
# benchmarks
######################################
Filebenchmarks = ROOT.TFile.Open(benchmarks, 'read')
histo2D = [None]*nclu 
counter=0
# the file already contains the 2D histos
for cluster in range(0,nclu) :
    # just take the histo
    # bin1.SetBins(90,0.,1800.,10,-1,1.);
    #Filebenchmarks->cd("clu"+str(cluster))
    histoname = str(cluster)+"_bin1"
    histo2D[cluster] =  Filebenchmarks.Get(histoname)
    #c=ROOT.TCanvas()
    #histo2D[cluster].Draw()
    #c.Print("teste.png")
    #print histo2D[cluster].GetSize()
    Xaxis = histo2D[cluster].GetXaxis()
    Yaxis = histo2D[cluster].GetYaxis()
    for ibin in range(0,int(Nvarbinmass)) :
       for ibin2 in range(0,int(Nvarbincosth)) :
           mass = float(1810- 20*(ibin)) 
           cost = 1.1- 0.2*(ibin2)
           #print str(mass) + " " + str(cost)
           binspot= histo2D[cluster].FindBin(mass,cost)
           nevBin = histo2D[cluster].GetBinContent(int(binspot))
           Matrix1[cluster][ibin][ibin2] =  nevBin
           counter = counter + histo2D[cluster].GetBinContent(int(binspot))
           #print (ibin,ibin2,Matrix1[cluster][ibin][ibin2],histo2D[cluster].GetBinContent(int(binspot)))
    print str(cluster) + " " +str(counter) + " events to reweinght to benchmark"
    counter=0
Filebenchmarks.Close()
print "read benchmarks"
######################################
# envelope - the sample 0 is the benchmarks
######################################
Fileenvelope = ROOT.TFile.Open(envelope, 'read')
histo2Dout = [[None]*nout]*nclu 
counter=0 
for cluster in range(0,nclu) :
  histoPath= "clu"+str(cluster+1)
  gDirectory.cd(histoPath)
  for out in range(0,nout) :
    histoname = str(out)+"_bin1"
    histo2Dout[cluster][out] =  gDirectory.Get(histoname)
    Xaxis = histo2Dout[cluster][out].GetXaxis()
    Yaxis = histo2Dout[cluster][out].GetYaxis()
    for ibin in range(0,int(Nvarbinmass)) :
        for ibin2 in range(0,int(Nvarbincosth)) :
            mass = float(1810- 20*(ibin)) 
            cost = 1.1- 0.2*(ibin2)
            #print str(mass) + " " + str(cost)
            binspot= histo2Dout[cluster][out].FindBin(mass,cost)
            nevBin = histo2Dout[cluster][out].GetBinContent(int(binspot))
            Matrix1out[cluster][out][ibin][ibin2] =  nevBin
            counter = counter + histo2Dout[cluster][out].GetBinContent(int(binspot))
    #print (ibin,ibin2,Matrix1[cluster][ibin][ibin2],histo2D[cluster].GetBinContent(int(binspot)))
    #print str(cluster+1) +" "+ str(out) +" "+str(counter) + " events to reweinght to benchmark"
    print str(cluster) + " "+ str(out) + " " +str(counter) + " events to reweinght to benchmark"
    counter=0
    #print histoPath + " " + histoname 
  gDirectory.Close()
Fileenvelope.Close()
print "read outliers"
######################################
# lambda-only
######################################
Filelambdaonly = ROOT.TFile.Open(lambdaonly, 'read')
histo2Dlam = [None]*Nlambda
counter=0
print lambdaonly
# the file already contains the 2D histos
gDirectory.cd("lamdaOnly")
for cluster in range(0,Nlambda) :
    histoname = str(cluster)+"_bin1"
    histo2Dlam[cluster] =  gDirectory.Get(histoname)
    c=ROOT.TCanvas()
    histo2Dlam[cluster].Draw()
    c.Print("teste_lam.png")
    Xaxis = histo2Dlam[cluster].GetXaxis()
    Yaxis = histo2Dlam[cluster].GetYaxis()
    for ibin in range(0,int(Nvarbinmass)) :
        for ibin2 in range(0,int(Nvarbincosth)) :
            mass = float(1810- 20*(ibin)) 
            cost = 1.1- 0.2*(ibin2)
            #print str(mass) + " " + str(cost)
            binspot= histo2Dlam[cluster].FindBin(mass,cost)
            nevBin = histo2Dlam[cluster].GetBinContent(int(binspot))
            Matrix2out[cluster][ibin][ibin2] =  nevBin
            counter = counter + histo2Dlam[cluster].GetBinContent(int(binspot))
    #print (ibin,ibin2,Matrix1[cluster][ibin][ibin2],histo2D[cluster].GetBinContent(int(binspot)))
    print str(cluster) + " " +str(counter) + " events to reweinght to benchmark"
    counter=0
gDirectory.Close()
Filelambdaonly.Close()
######################################
'''
for cluster in range(0,1) :
   for ibin in range(0,int(Nvarbinmass)) :
      for ibin2 in range(0,int(Nvarbincosth)) :
          print (ibin,ibin2,Matrix1[cluster][ibin][ibin2])
'''          


######################################
## pass over all the events, make a matrix  weights
Sum2D = ROOT.TH2D("Sum2D","Sum2D",90,0.,1800.,10,-1,1.)
countereventHs=0
countereventLs=0
countereventHa=0
countereventLa=0
onebin=0
for files in range(0,13) :
   fileLow = inputLM + file + str(files) + fileL 
   fileHigh = inputHM + file + str(files) + fileH
   print fileLow 
   FileLow = ROOT.TFile.Open(fileLow, 'read')
   FileHigh = ROOT.TFile.Open(fileHigh, 'read')
   # ####################
   # high
   countereventHs = countereventHs+ FileHigh.TCVARS.GetEntries()
   for event in FileHigh.TCVARS:
       Mtot =  event.mtot # change to the gen mtot
       Cost = abs(event.cut_based_ct)-0.1 # change to the gen cost
       #print Cost
       Sum2D.Fill(Mtot,Cost)
       # find the bin
       imass=-1 
       icosth=-1
       for ibin in range(0,int(Nvarbinmass)) :
           #print (float(1800- 20*(ibin)),Mtot)
           if (Mtot > float(1800- 20*(ibin)) and imass < 0): 
               #print "bingo"
              imass=ibin
              continue
       for ibin in range(0,int(Nvarbincosth)) :
           #print (float(1.0- 0.2*(ibin)),Cost)
           if (Cost >= float(1.0- 0.2*(ibin)) and icosth < 0 ): 
               #print "bingo"
              icosth = ibin
              continue
       #print str(imass) + " " + str(icosth)
       if (imass < 0 or icosth < 0): print "no bin found - making fullsimed events matrix" 
       elif (imass > 0 and icosth > 0) :
           #old = Matrix2high[imass][icosth]
           Matrix2high[imass][icosth] +=  1
           countereventHa += 1
           continue
       if(files==0) : Matrix1validationhigh[imass][icosth] +=  1
   # ####################
   # low
   countereventLs = countereventLs+ FileLow.TCVARS.GetEntries()
   for event in FileLow.TCVARS:
        Mtot =  event.mtot # change to the gen mtot
        Cost = abs(event.cut_based_ct)-0.1 # change to the gen cost
        #Sum2D.Fill(Mtot,Cost)
        # find the bin
        imass=-1 
        icosth=-1
        for ibin in range(0,int(Nvarbinmass)) :
            #print (float(1500)- 30*(ibin),massV[i][S1],float(imass))
            if (Mtot > float(1800- 20*(ibin)) and float(imass) < 0): 
                imass=ibin
                continue
        for ibin in range(0,int(Nvarbincosth)) :
            if (Cost >= 1.0- 0.2*(ibin) and icosth < 0 ): 
                icosth = ibin
                continue
        #print str(imass) + " " + str(icosth)
        if (imass < 0 or icosth < 0): print "no bin found - making fullsimed events matrix" 
        else : 
          #print Matrix2low[imass][icosth]
          Matrix2low[imass][icosth] += 1
          countereventLa += 1
        if(files==0) : Matrix1validationlow[imass][icosth] +=  1
   FileLow.Close()
   FileHigh.Close()
#######################################################
c=ROOT.TCanvas()
Sum2D.Draw()
c.Print("Sum2D.png")
#######################################
#print onebin
#print Matrix2high[69][1]
countereventH=0
countereventL=0
countereventbench=0
for ibin in range(0,int(Nvarbinmass)) :
    for ibin2 in range(0,int(Nvarbincosth)) :
        #print Matrix2[ibin][ibin2]
        #print str(ibin) + " , " +str(ibin2)+" , "+ str(Matrix1[0][ibin][ibin2]) + " " + str(Matrix2low[ibin][ibin2]) + " " + str(Matrix2high[ibin][ibin2])
        countereventH +=  Matrix2high[ibin][ibin2]
        countereventL += Matrix2low[ibin][ibin2]
        countereventbench = countereventbench + Matrix1[0][ibin][ibin2]
print "we have "+str(countereventL)+" events to reweight in the low mass category ("+str(countereventLs)+", "+str(countereventLa)+"), and"
print "we have "+str(countereventH)+" events to reweight in the high mass category ("+str(countereventHs)+", "+str(countereventHa)+"), and"
print "we have "+str(countereventbench)+" events to reweight in each benchmark"
#######################################
# make the minitrees 
# Create a struct
gROOT.ProcessLine(\
                  "struct MyStruct{\
                  Int_t cut_based_ct;\
                  Double_t mtotgen;\
                  Double_t costgen;\
                  Double_t mtot;\
                  Double_t mtotTHweighted;\
                  Double_t mgg;\
                  Double_t mjj;\
                  Double_t evWeight;\
                  Double_t evWeightTh;\
                  Double_t evWeightTot;\
                  };")
from ROOT import MyStruct
################################
# benchmarks
benchmarksV3 = "V3benchmarks/V3_LT_output_GluGluToHHTo2B2G_"
fcluH = [None]*int(nclu+1) 
fcluL = [None]*int(nclu+1) 
################################
for cluster in range(0,nclu+1) : # cluster = -1 (file =0) is for the validation sample (0 of V1)
    name = None
    if (cluster ==0) : name = "box_validation"
    elif (cluster <13): name = "JHEPv3_benchmark_" # benchmarks
    fcluL[cluster] = ROOT.TFile(benchmarksV3+name+str(cluster)+fileL,"RECREATE")   
    t = TTree('TCVARS','My test tree')
    s = MyStruct()
    t.Branch('cut_based_ct',AddressOf(s,'cut_based_ct'),'cut_based_ct/I')
    t.Branch('mtotgen',AddressOf(s,'mtotgen'),'mtotgen/D')
    t.Branch('costgen',AddressOf(s,'costgen'),'costgen/D')
    t.Branch('mtot',AddressOf(s,'mtot'),'mtot/D')
    t.Branch('mgg',AddressOf(s,'mgg'),'mgg/D')
    t.Branch('mjj',AddressOf(s,'mjj'),'mjj/D')
    t.Branch('evWeight',AddressOf(s,'evWeight'),'evWeight/D')
    t.Branch('evWeightTh',AddressOf(s,'evWeightTh'),'evWeightTh/D')
    t.Branch('evWeightTot',AddressOf(s,'evWeightTot'),'evWeightTot/D')
    #
    for files in range(0,13) :
       fileLow = inputLM + file + str(files) + fileL 
       FileLow = ROOT.TFile.Open(fileLow, 'read')
       ################################
       # low
       for event in FileLow.TCVARS:
           #
           ###
           Mtotgen =  event.mtot # change to the gen mtot
           Costgen = event.cut_based_ct # change to the gen cost
           Mtot =  event.mtot 
           Mgg =  event.mgg 
           Cat = event.cut_based_ct 
           Mjj =  event.mjj 
           evWei = event.evWeight
           # I calculate weight by event and save to a tree
           weight = 0 #np.zeros((nclu))
           weightvalidation=0
           ###############
           # find the bin
           imass=-1 
           icosth=-1
           for ibin in range(0,int(Nvarbinmass)) :
               #print (float(1500)- 30*(ibin),massV[i][S1],float(imass))
               if (Mtotgen > float(1800- 20*(ibin)) and float(imass) < 0): 
                  imass=ibin
                  continue
           for ibin in range(0,int(Nvarbincosth)) :
               if (Costgen > 1.0- 0.2*(ibin) and icosth < 0 ): 
                  icosth = ibin
                  continue
               if (Costgen ==-1.0 and icosth < 0 ): icosth = 10
           #print str(Matrix1[imass][icosth]) + " " + str(Matrix2low[imass][icosth])
           if (imass < 0 or icosth < 0): print "no bin found - weight event by event " +str(Mtotgen)+" "+str(Costgen)
           elif (Matrix2low[imass][icosth] > 0 ) :
               if (cluster ==0) : weightvalidation = float(Matrix1validationlow[imass][icosth]/Matrix2low[imass][icosth])
               else : weight = float(Matrix1[cluster-1][imass][icosth]/Matrix2low[imass][icosth])                
           s.cut_based_ct = int(Cat)
           s.mtotgen = float(Mtotgen)
           s.costgen = float(Costgen)
           s.mtot = float(Mtot)
           s.mgg = float(Mgg)
           s.mjj = float(Mjj)
           s.evWeight = float(evWei)
           s.evWeightTh = float(weight)
           s.evWeightTot = float(weight*evWei)
           t.Fill()
    fcluL[cluster].Write()
    fcluL[cluster].Close()
    print "saved "+str(benchmarksV3+name+str(cluster)+fileL)
    ################################
    #
    ###############################
    fcluH[cluster] = ROOT.TFile(benchmarksV3+name+str(cluster+1)+fileH,"RECREATE")   
    tH = TTree('TCVARS','My test tree')
    tH.Branch('cut_based_ct',AddressOf(s,'cut_based_ct'),'cut_based_ct/I')
    tH.Branch('mtotgen',AddressOf(s,'mtotgen'),'mtotgen/D')
    tH.Branch('costgen',AddressOf(s,'costgen'),'costgen/D')
    tH.Branch('mtot',AddressOf(s,'mtot'),'mtot/D')
    tH.Branch('mgg',AddressOf(s,'mgg'),'mgg/D')
    tH.Branch('mjj',AddressOf(s,'mjj'),'mjj/D')
    tH.Branch('evWeight',AddressOf(s,'evWeight'),'evWeight/D')
    tH.Branch('evWeightTh',AddressOf(s,'evWeightTh'),'evWeightTh/D')
    tH.Branch('evWeightTot',AddressOf(s,'evWeightTot'),'evWeightTot/D')
    for files in range(0,13) :
       #################################
       # high
       fileHigh = inputHM + file + str(files) + fileH        
       FileHigh = ROOT.TFile.Open(fileHigh, 'read') 
       for event in FileHigh.TCVARS:
           #
           ###
           Mtotgen =  event.mtot # change to the gen mtot
           Costgen = event.cut_based_ct # change to the gen cost
           Mtot =  event.mtot 
           Mgg =  event.mgg 
           Cat = event.cut_based_ct 
           Mjj =  event.mjj 
           evWei = event.evWeight
           # I calculate weight by event and save to a tree
           weight = 0 #np.zeros((nclu))
           weightvalidation=0
           ###############
           # find the bin
           imass=-1 
           icosth=-1
           for ibin in range(0,int(Nvarbinmass)) :
               #print (float(1500)- 30*(ibin),massV[i][S1],float(imass))
               if (Mtotgen > float(1800- 20*(ibin)) and float(imass) < 0): 
                   imass=ibin
                   continue
           for ibin in range(0,int(Nvarbincosth)) :
               if (Costgen > 1.0- 0.2*(ibin) and icosth < 0 ): 
                   icosth = ibin
                   continue
               if (Costgen ==-1.0 and icosth < 0 ): icosth = 10
           #print str(Matrix1[imass][icosth]) + " " + str(Matrix2low[imass][icosth])
           if (imass < 0 or icosth < 0): print "no bin found - weight event by event " +str(Mtotgen)+" "+str(Costgen)
           elif (Matrix2low[imass][icosth] > 0 ) :
              if (cluster ==0) : weightvalidation = float(Matrix1validationhigh[imass][icosth]/Matrix2high[imass][icosth])
              elif (cluster <13): weight = float(Matrix1[cluster-1][imass][icosth]/Matrix2high[imass][icosth]) # benchmarks
           s.cut_based_ct = int(Cat)
           s.mtotgen = float(Mtotgen)
           s.costgen = float(Costgen)
           s.mtot = float(Mtot)
           s.mgg = float(Mgg)
           s.mjj = float(Mjj)
           s.evWeight = float(evWei)
           s.evWeightTh = float(weight)
           s.evWeightTot = float(weight*evWei)
           tH.Fill()
    fcluH[cluster].Write()
    fcluH[cluster].Close()
    print "saved "+str(benchmarksV3+name+str(cluster)+fileH)

###################################################
# lambda only
##################################################
lambdaonlyFolder = "lambdaonly/V3_LT_output_GluGluToHHTo2B2G_"
fcluH = [None]*Nlambda
fcluL = [None]*Nlambda 
name = "lambda_only"
for cluster in range(0,52) : 
    fcluL[cluster] = ROOT.TFile(lambdaonlyFolder+name+str(cluster+1)+fileL,"RECREATE")   
    t = TTree('TCVARS','My test tree')
    s = MyStruct()
    t.Branch('cut_based_ct',AddressOf(s,'cut_based_ct'),'cut_based_ct/I')
    t.Branch('mtotgen',AddressOf(s,'mtotgen'),'mtotgen/D')
    t.Branch('costgen',AddressOf(s,'costgen'),'costgen/D')
    t.Branch('mtot',AddressOf(s,'mtot'),'mtot/D')
    t.Branch('mgg',AddressOf(s,'mgg'),'mgg/D')
    t.Branch('mjj',AddressOf(s,'mjj'),'mjj/D')
    t.Branch('evWeight',AddressOf(s,'evWeight'),'evWeight/D')
    t.Branch('evWeightTh',AddressOf(s,'evWeightTh'),'evWeightTh/D')
    t.Branch('evWeightTot',AddressOf(s,'evWeightTot'),'evWeightTot/D')
    #
    for files in range(0,13) :
       fileLow = inputLM + file + str(files) + fileL 
       FileLow = ROOT.TFile.Open(fileLow, 'read')
       ################################
       # low
       for event in FileLow.TCVARS:
           #
           ###
           Mtotgen =  event.mtot # change to the gen mtot
           Costgen = event.cut_based_ct # change to the gen cost
           Mtot =  event.mtot 
           Mgg =  event.mgg 
           Cat = event.cut_based_ct 
           Mjj =  event.mjj 
           evWei = event.evWeight
           # I calculate weight by event and save to a tree
           weight = 0 #np.zeros((nclu))
           weightvalidation=0
           ###############
           # find the bin
           imass=-1 
           icosth=-1
           for ibin in range(0,int(Nvarbinmass)) :
               #print (float(1500)- 30*(ibin),massV[i][S1],float(imass))
               if (Mtotgen > float(1800- 20*(ibin)) and float(imass) < 0): 
                  imass=ibin
                  continue
           for ibin in range(0,int(Nvarbincosth)) :
               if (Costgen > 1.0- 0.2*(ibin) and icosth < 0 ): 
                  icosth = ibin
                  continue
               if (Costgen ==-1.0 and icosth < 0 ): icosth = 10
           #print str(Matrix1[imass][icosth]) + " " + str(Matrix2low[imass][icosth])
           if (imass < 0 or icosth < 0): print "no bin found - weight event by event " +str(Mtotgen)+" "+str(Costgen)
           elif (Matrix2low[imass][icosth] > 0 ) :
               weight = float(Matrix2out[cluster][imass][icosth]/Matrix2low[imass][icosth])                
           s.cut_based_ct = int(Cat)
           s.mtotgen = float(Mtotgen)
           s.costgen = float(Costgen)
           s.mtot = float(Mtot)
           s.mgg = float(Mgg)
           s.mjj = float(Mjj)
           s.evWeight = float(evWei)
           s.evWeightTh = float(weight)
           s.evWeightTot = float(weight*evWei)
           t.Fill()
    fcluL[cluster].Write()
    fcluL[cluster].Close()
    print "saved "+str(lambdaonlyFolder+name+str(cluster+1)+fileL)
    ################################
    #
    ###############################
    fcluH[cluster] = ROOT.TFile(lambdaonlyFolder+str(cluster+1)+fileH,"RECREATE")   
    tH = TTree('TCVARS','My test tree')
    tH.Branch('cut_based_ct',AddressOf(s,'cut_based_ct'),'cut_based_ct/I')
    tH.Branch('mtotgen',AddressOf(s,'mtotgen'),'mtotgen/D')
    tH.Branch('costgen',AddressOf(s,'costgen'),'costgen/D')
    tH.Branch('mtot',AddressOf(s,'mtot'),'mtot/D')
    tH.Branch('mgg',AddressOf(s,'mgg'),'mgg/D')
    tH.Branch('mjj',AddressOf(s,'mjj'),'mjj/D')
    tH.Branch('evWeight',AddressOf(s,'evWeight'),'evWeight/D')
    tH.Branch('evWeightTh',AddressOf(s,'evWeightTh'),'evWeightTh/D')
    tH.Branch('evWeightTot',AddressOf(s,'evWeightTot'),'evWeightTot/D')
    for files in range(0,13) :
       #################################
       # high
       fileHigh = inputHM + file + str(files) + fileH        
       FileHigh = ROOT.TFile.Open(fileHigh, 'read') 
       for event in FileHigh.TCVARS:
           #
           ###
           Mtotgen =  event.mtot # change to the gen mtot
           Costgen = event.cut_based_ct # change to the gen cost
           Mtot =  event.mtot 
           Mgg =  event.mgg 
           Cat = event.cut_based_ct 
           Mjj =  event.mjj 
           evWei = event.evWeight
           # I calculate weight by event and save to a tree
           weight = 0 #np.zeros((nclu))
           weightvalidation=0
           ###############
           # find the bin
           imass=-1 
           icosth=-1
           for ibin in range(0,int(Nvarbinmass)) :
               #print (float(1500)- 30*(ibin),massV[i][S1],float(imass))
               if (Mtotgen > float(1800- 20*(ibin)) and float(imass) < 0): 
                   imass=ibin
                   continue
           for ibin in range(0,int(Nvarbincosth)) :
               if (Costgen > 1.0- 0.2*(ibin) and icosth < 0 ): 
                   icosth = ibin
                   continue
               if (Costgen ==-1.0 and icosth < 0 ): icosth = 10
           #print str(Matrix1[imass][icosth]) + " " + str(Matrix2low[imass][icosth])
           if (imass < 0 or icosth < 0): print "no bin found - weight event by event " +str(Mtotgen)+" "+str(Costgen)
           elif (Matrix2low[imass][icosth] > 0 ) :
              weight = float(Matrix2out[cluster][imass][icosth]/Matrix2high[imass][icosth]) # benchmarks
           s.cut_based_ct = int(Cat)
           s.mtotgen = float(Mtotgen)
           s.costgen = float(Costgen)
           s.mtot = float(Mtot)
           s.mgg = float(Mgg)
           s.mjj = float(Mjj)
           s.evWeight = float(evWei)
           s.evWeightTh = float(weight)
           s.evWeightTot = float(weight*evWei)
           tH.Fill()
    fcluH[cluster].Write()
    fcluH[cluster].Close()
    print "saved "+ lambdaonlyFolder+name+str(cluster+1)+fileH

###################################################
# outliers
##################################################
outliersV3 = "V3outliers/V3_LT_output_GluGluToHHTo2B2G_"
fcluH = [[None]*nout]*nclu 
fcluL = [[None]*nout]*nclu 
name = "JHEPv3_cluster_"
name2="_outlier_"
# low mass
for cluster in range(0,nclu) : 
  for out in range(0,nout) :
    fcluL[cluster][out] = ROOT.TFile(outliersV3+name+str(cluster+1)+name2+ str(out) +fileL,"RECREATE")   
    t = TTree('TCVARS','My test tree')
    s = MyStruct()
    t.Branch('cut_based_ct',AddressOf(s,'cut_based_ct'),'cut_based_ct/I')
    t.Branch('mtotgen',AddressOf(s,'mtotgen'),'mtotgen/D')
    t.Branch('costgen',AddressOf(s,'costgen'),'costgen/D')
    t.Branch('mtot',AddressOf(s,'mtot'),'mtot/D')
    t.Branch('mgg',AddressOf(s,'mgg'),'mgg/D')
    t.Branch('mjj',AddressOf(s,'mjj'),'mjj/D')
    t.Branch('evWeight',AddressOf(s,'evWeight'),'evWeight/D')
    t.Branch('evWeightTh',AddressOf(s,'evWeightTh'),'evWeightTh/D')
    t.Branch('evWeightTot',AddressOf(s,'evWeightTot'),'evWeightTot/D')
    #t.Branch('mtotTHweighted',AddressOf(s,'mtotTHweighted'),'mtotTHweighted/D')
    for files in range(0,13) :
       fileLow = inputLM + file + str(files) + fileL 
       FileLow = ROOT.TFile.Open(fileLow, 'read')
       ################################
       # low
       for event in FileLow.TCVARS:
           #
           ###
           Mtotgen =  event.mtot # change to the gen mtot
           Costgen = event.cut_based_ct # change to the gen cost
           Mtot =  event.mtot 
           Mgg =  event.mgg 
           Cat = event.cut_based_ct 
           Mjj =  event.mjj 
           evWei = event.evWeight
           # I calculate weight by event and save to a tree
           weight = 0 #np.zeros((nclu))
           weightvalidation=0
           ###############
           # find the bin
           imass=-1 
           icosth=-1
           for ibin in range(0,int(Nvarbinmass)) :
               #print (float(1500)- 30*(ibin),massV[i][S1],float(imass))
               if (Mtotgen > float(1800- 20*(ibin)) and float(imass) < 0): 
                  imass=ibin
                  continue
           for ibin in range(0,int(Nvarbincosth)) :
               if (Costgen > 1.0- 0.2*(ibin) and icosth < 0 ): 
                  icosth = ibin
                  continue
               if (Costgen ==-1.0 and icosth < 0 ): icosth = 10
           #print str(Matrix1[imass][icosth]) + " " + str(Matrix2low[imass][icosth])
           if (imass < 0 or icosth < 0): print "no bin found - weight event by event " +str(Mtotgen)+" "+str(Costgen)
           elif (Matrix2low[imass][icosth] > 0 ) :
               weight = float(Matrix1out[cluster][out][imass][icosth]/Matrix2low[imass][icosth])                
           s.cut_based_ct = int(Cat)
           s.mtotgen = float(Mtotgen)
           s.costgen = float(Costgen)
           s.mtot = float(Mtot)
           s.mgg = float(Mgg)
           s.mjj = float(Mjj)
           s.evWeight = float(evWei)
           s.evWeightTh = float(weight)
           s.evWeightTot = float(weight*evWei)
           #s.mtotTHweighted= float(Mtot)
           t.Fill()
    fcluL[cluster][out].Write()
    fcluL[cluster][out].Close()
    print "saved "+str(outliersV3+name+str(cluster+1)+name2+ str(out) +fileL)
    ################################
    #
    ###############################
# high mass
for cluster in range(0,nclu) : 
    for out in range(0,nout) :    
      fcluH[cluster][out] = ROOT.TFile(outliersV3+name+str(cluster+1)+name2+ str(out) +fileH,"RECREATE")   
      tHH = TTree('TCVARS','My test tree')
      s = MyStruct()
      tHH.Branch('cut_based_ct',AddressOf(s,'cut_based_ct'),'cut_based_ct/I')
      tHH.Branch('mtotgen',AddressOf(s,'mtotgen'),'mtotgen/D')
      tHH.Branch('costgen',AddressOf(s,'costgen'),'costgen/D')
      tHH.Branch('mtot',AddressOf(s,'mtot'),'mtot/D')
      tHH.Branch('mgg',AddressOf(s,'mgg'),'mgg/D')
      tHH.Branch('mjj',AddressOf(s,'mjj'),'mjj/D')
      tHH.Branch('evWeight',AddressOf(s,'evWeight'),'evWeight/D')
      tHH.Branch('evWeightTh',AddressOf(s,'evWeightTh'),'evWeightTh/D')
      tHH.Branch('evWeightTot',AddressOf(s,'evWeightTot'),'evWeightTot/D')
      for files in range(0,13) :
         #################################
         # high
         fileHigh = inputHM + file + str(files) + fileH        
         FileHigh = ROOT.TFile.Open(fileHigh, 'read') 
         for event in FileHigh.TCVARS:
             ###
             Mtotgen =  event.mtot # change to the gen mtot
             Costgen = event.cut_based_ct # change to the gen cost
             Mtot =  event.mtot 
             Mgg =  event.mgg 
             Cat = event.cut_based_ct 
             Mjj =  event.mjj 
             evWei = event.evWeight
             # I calculate weight by event and save to a tree
             weight = 0 #np.zeros((nclu))
             weightvalidation=0
             ###############
             # find the bin
             imass=-1 
             icosth=-1
             for ibin in range(0,int(Nvarbinmass)) :
                 #print (float(1500)- 30*(ibin),massV[i][S1],float(imass))
                 if (Mtotgen > float(1800- 20*(ibin)) and float(imass) < 0): 
                     imass=ibin
                     continue
             for ibin in range(0,int(Nvarbincosth)) :
                 if (Costgen > 1.0- 0.2*(ibin) and icosth < 0 ): 
                     icosth = ibin
                     continue
                 if (Costgen ==-1.0 and icosth < 0 ): icosth = 10
             #print str(Matrix1[imass][icosth]) + " " + str(Matrix2low[imass][icosth])
             if (imass < 0 or icosth < 0): print "no bin found - weight event by event " +str(Mtotgen)+" "+str(Costgen)
             elif (Matrix2low[imass][icosth] > 0 ) : weight = float(Matrix1out[cluster][out][imass][icosth]/Matrix2high[imass][icosth]) # benchmarks
             s.cut_based_ct = int(Cat)
             s.mtotgen = float(Mtotgen)
             s.costgen = float(Costgen)
             s.mtot = float(Mtot)
             s.mgg = float(Mgg)
             s.mjj = float(Mjj)
             s.evWeight = float(evWei)
             #s.evWeightTh = float(weight)
             s.evWeightTot = float(weight*evWei)
             tHH.Fill()
      fcluH[cluster][out].Write()
      fcluH[cluster][out].Close()
      print "saved "+str(outliersV3+name+str(cluster+1)+name2+ str(out) +fileH)
       


