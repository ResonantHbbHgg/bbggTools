#!/usr/bin/env python
from ROOT import *
gROOT.SetBatch()


f = TFile('/afs/cern.ch/user/a/andrey/public/HH/NonRes_weights_from_Xandra_v0.root','OPEN')
t = f.Get('TCVARS')


t.Draw('mHH>>h(200,0,2500)','file==0')
c1.SaveAs('mHH_Box_file.png')
t.Draw('mHH>>h(200,0,2500)','wBox')
c1.SaveAs('mHH_Box_wei.png')

'''
for n in range(1,13):
    print n
    t.Draw('mHH>>h(200,0,2500)','file=='+str(n))
    c1.SaveAs('mHH_'+str(n)+'_file.png')
   
    t.Draw('mHH>>h(200,0,2500)','bench'+str(n))
    c1.SaveAs('mHH_'+str(n)+'_weight.png')

'''

#t.Draw('(file*50000+evt)>>e_Box(650000,1,650001)','wBox')
#c1.SaveAs('evt.png')

fout = TFile('fout.root','recreate')

name = 'wBox'
t.Draw('(file*50000+evt)>>e_'+name+'(650000,1,650001)',name)
gPad.GetPrimitive('e_'+name).Write()
for b in xrange(1,13):
    name = 'bench'+str(b)
    t.Draw('(file*50000+evt)>>e_'+name+'(650000,1,650001)',name)
    gPad.GetPrimitive('e_'+name).Write()
    for o in xrange(1,7):
        name = 'bench'+str(b)+'_out'+str(o)
        t.Draw('(file*50000+evt)>>e_'+name+'(650000,1,650001)',name)
        gPad.GetPrimitive('e_'+name).Write()
        


# Manually cross check some weights from few events
for a in xrange(46000, 46020):
    print a, e_wBox.GetBinContent(12*50000+a)
t.Scan('file:evt:wBox','file==12&&evt>46000&&evt<46020')

