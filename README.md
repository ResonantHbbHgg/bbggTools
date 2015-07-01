bbggTools
=========

Welcome to the HH->bbgg RunII analysis!

This package is a set of tools to be used to perform the resonant and non-resonant HH->bbgg analysis for RunII.
It should be used under the FLASHgg framework and not as a standalone package.

Before cloning this package, get the FLASHgg framework: https://github.com/cms-analysis/flashgg  
Also check their Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/FLASHggFramework#New_simple_instructions  
Current CMSSW version to be used (as of 01July2015): CMSSW_7_4_6_patch2

Follow their instructions (including compilation).

Then, go to the flashgg directory and clone bbggTools, and compile again:
```
cd $CMSSW_BASE/src/flashgg
git clone git@github.com:ResonantHbbHgg/bbggTools.git bbggTools
cd bbggTools
scramv1 b -j 10
```

#### MetaData
On the MetaData directory, you will find the JSON files that are used to produce microAOD with FLASHgg from the miniAOD datasets. We only need to produce the signal microAOD for ourselves, the rest should be inherited from Hgg. For more information about how to create MicroAOD, see the MetaData README.md.

#### Test
On the test directory, you will find the different scripts to run our analysis, such as the flat tree producer, plots producer, etc. See directory README for more information.
