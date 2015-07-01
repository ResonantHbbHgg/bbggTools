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
git clone git@github.com:ResonantHbbHgg/bbggTools.git
cd bbggTools
scramv1 b -j 10
```
