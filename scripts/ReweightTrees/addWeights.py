#!/usr/bin/env python
from ROOT import TFile, TTree, TH1D
from array import array

# This file should have been created by another script:
wFile = TFile("fout.root","read")

h = {}
h['wBox'] = wFile.Get('e_wBox')

# Here we get the Hists that contain Evnt-to-weight data
# also, we create a leaves for a new tree
leavesBench = "wBox/D"
leavesOut   = ""
for b in xrange(1,13):
    name = 'bench'+str(b)
    h[name] = wFile.Get('e_'+name)
    leavesBench+=':%s/D'%name
    for o in xrange(1,7):
        name = 'bench'+str(b)+'_out'+str(o)
        h[name] = wFile.Get('e_'+name)
        if len(leavesOut):
            leavesOut+=':%s/D'%name
        else:
            leavesOut+='%s/D'%name


# Now lets loop over all Limit trees (one file per Node)
# and add the weigts leaves for each benchmark/outlier

DirWithLimitTrees = '/afs/cern.ch/user/a/andrey/public/HH/HHnonRes_LimitTrees/'

#for n in ['SM','box','2','3','4','5','6','7','8','9','10','11','12','13']:
for n in ['SM','3']: # <- use for tests
  print '\t * Doing Node', n

  fName = 'LT_output_GluGluToHHTo2B2G_node_'+n+'_13TeV-madgraph_0.root'
  fSig = TFile(DirWithLimitTrees+fName,'read')

  # fSig.Print()
  tSig = fSig.Get('TCVARS')
  
  #tSig.Print()


  newFile = TFile("Wei_"+fName,"RECREATE")
  newTree = tSig.CloneTree(0)


  wArrBench = array('d',[0]*13)
  wArrOut   = array('d',[0]*(12*6))

  newBranch1 = newTree.Branch("benchmarks" , wArrBench, leavesBench)
  newBranch2 = newTree.Branch("outliers" ,   wArrOut,   leavesOut)

  alist, i = [], 0

  if n=='SM':
    fNum = 0
  elif n=='box':
    fNum = 1
  else:
    fNum = int(n)

  while tSig.GetEntry(i):
    i += 1

    eNum = tSig.GetLeaf( "o_evt" ).GetValue()
    #if eNum < 10:
    #  print fNum, eNum

    for b in xrange(1,13):
      name = 'bench'+str(b)
      wArrBench[b-1] = h[name].GetBinContent(int(fNum*50000+eNum))
    for b in xrange(1,13):
      for o in xrange(1,7):
        name = 'bench'+str(b)+'_out'+str(o)
        wArrOut[(b-1)*6+o-1] = h[name].GetBinContent(int(fNum*50000+eNum))

    newTree.Fill()

  newTree.Write()

