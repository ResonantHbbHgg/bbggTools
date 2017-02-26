#!/usr/bin/env python
import os
import sys

print '###########################'
print '## bbggTools After Build ##'
print '###########################'

env = os.environ['CMSSW_BASE']+'/src/flashgg/'
bbggMetaData = env+'/bbggTools/MetaData/'
flashggMetaData = env+'/MetaData/data/'
dirs = next(os.walk(bbggMetaData))[1]

print 'Copying directories in', bbggMetaData, 'to', flashggMetaData
for d in dirs:
  print '\t', d
  td = bbggMetaData+d
  command = 'cp -r '+td+' '+flashggMetaData
  os.system(command)


print 'Adding signal cross sections to flashgg cross sections list'
mycross = file(bbggMetaData+'HHSignalCrossSections.txt', 'r')
flashggcross = file(flashggMetaData+'cross_sections.json')
outcross = ''
if flashggcross.read().find('GluGluToHHTo2B2G'):
  print '\t HH signal cross sections already in flashgg cross sections, only need to run after build once!'
else:
  for l in flashggcross:
   outcross += l
   if '{' in l and '}' not in l:
     for c in mycross:
        outcross += c
  os.system('mv ' + flashggMetaData+'cross_sections.json ' + flashggMetaData+'old_cross_sections.json')
  newcross = file(flashggMetaData+'cross_sections.json', 'w+')
  newcross.write(outcross)
  newcross.close()

print '###########'
print '## Done! ##'
print '###########'
