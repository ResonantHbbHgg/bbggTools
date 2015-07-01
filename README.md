bbggTools
=========

Welcome to the HH->bbgg RunII analysis!

This package is a set of tools to be used to perform the resonant and non-resonant HH->bbgg analysis for RunII.
It should be used under the FLASHgg framework and not as a standalone package.

Before cloning this package, get the FLASHgg framework:
https://github.com/cms-analysis/flashgg

Follow their instructions (including compilation).

Then, go to the flashgg directory and clone bbggTools, and compile again:
```
cd $CMSSW_BASE/src/flashgg
git clone git@github.com:ResonantHbbHgg/bbggTools.git bbggTools
cd bbggTools
scramv1 b -j 10
```

#### MetaData
On the MetaData directory, you will find the JSON files that are used to produce microAOD with FLASHgg from the miniAOD datasets. We only need to produce the signal microAOD for ourselves, the rest should be inherited from Hgg.

#### Test
On the test directory, you will find the different scripts to run our analysis, such as the flat tree producer, plots producer, etc. See directory README for more information.
