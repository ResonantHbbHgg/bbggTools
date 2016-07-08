#!/usr/bin/env python
from ROOT import TFile, TTree, TH1D
from array import array

fSig = TFile('/afs/cern.ch/user/a/andrey/public/HH/NonResAll.root','read')

tSig = fSig.Get('fsDir/TCVARS')
tSig.Print()


wFile = TFile("fout.root","read")
#wFile = TFile("teste_weight.root","read")
#wTree = wFile.Get("TCVARS")
h = {}
h['wBox'] = wFile.Get('e_wBox')

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

newFile = TFile("newFilename.root","RECREATE")
newTree = tSig.CloneTree(0)


wArrBench = array('d',[0]*13)
wArrOut   = array('d',[0]*(12*6))

newBranch1 = newTree.Branch("benchmarks" , wArrBench, leavesBench)
newBranch2 = newTree.Branch("outliers" ,   wArrOut,   leavesOut)

alist, i = [], 0
while tSig.GetEntry(i):
    i += 1

    fNum = tSig.GetLeaf( "file" ).GetValue()
    eNum = tSig.GetLeaf( "o_evt" ).GetValue()
    if eNum < 10:
        print fNum, eNum

    for b in xrange(1,13):
        name = 'bench'+str(b)
        wArrBench[b-1] = h[name].GetBinContent(int(fNum*50000+eNum))
    for b in xrange(1,13):
        for o in xrange(1,7):
            name = 'bench'+str(b)+'_out'+str(o)
            wArrOut[(b-1)*6+o-1] = h[name].GetBinContent(int(fNum*50000+eNum))

    newTree.Fill()

newTree.Write()
