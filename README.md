bbggTools
=========

Welcome to the HH->bbgg RunII analysis!

This package is a set of tools to be used to perform the resonant and non-resonant HH->bbgg analysis for RunII.
It should be used under the FLASHgg framework and not as a standalone package.

Before cloning this package, get the FLASHgg framework: https://github.com/cms-analysis/flashgg  
Also check their Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/FLASHggFramework#Instructions_for_users  
Current CMSSW version to be used (as of 12October2015): CMSSW_7_4_12

Follow their instructions (including compilation).

Then, go to the flashgg directory and clone bbggTools, and compile again:
```
cd $CMSSW_BASE/src/flashgg
git clone git@github.com:ResonantHbbHgg/bbggTools.git bbggTools
cd bbggTools
./Compile.sh
```

Due to some problems with unused variables in limit codes (memory handling in RooFit), we need to set the propert compilation variables to ignore unused variables. That is done in Compile.sh. For cleaning the area, you can still do simply scramv1 b clean.


## Important scripts

### Produce Selection Trees

Selection trees are flat trees containing the selected diphoton and the selected dijet objects (along with the photons and the jets). It also contains the event weight, some simple ID information on each object and the b-tagging for the jets. It's used to make stack plots and to generate limit trees.

To generate one single tree, do:
```
cd test
cmsRun MakeTrees.py doSelection=1 inputFiles=<input microAOD> outputFile=<output flat tree>
```

Due to the prompt photon double counting, there are two extra options:
1) doDoubleCountingMitigation: =1 for QCD/GJet/DiPhotonJets.
2) nPromptPhotons: =2 for DiPhotonJets, =1 for GJet, =0 for QCD

To generate multiple trees (by dataset) on batch, do:
```
cd batchTools
python batch_TreeCentralMC.py -c <Campaign> -s <Sample>
```

Campaign: the current microAOD campaign. It's the folder where the datasets.json file is stored: flashgg/MetaData/data/<Campaign>/datasets.json
Sample: name of the sample you want to process. List of important samples is listed on AvailableSamples.py. 

Example:
```
python batch_TreeCentralMC.py -c RunIISpring15-25ns -s GJet_Pt-20to40 
python batch_TreeCentralMC.py -c RunIISpring15-25ns -s DoubleEG --isData
```

The script will automatically set doDoubleCountingMitigation and nPromptPhotons for the correct samples.

In batch_TreeCentralMC.py, don't forget to modify the output location to your folder!

### Produce Limit Trees

The Limit Tree is a created from Selection Trees and store only the needed variables to produce the final limits (diphoton invariant mass, dijet invariant mass, 4-body invariant mass, event weight, category).

Example:
```
LimitTreeMaker <input text file> <output location>
```

The input text file should be a list of Selection Trees root files. One limit tree file will be generated for each line in the text file. You can hadd things after, if necessary. The output location is the folder where you want the limit tree files to be stored.

#### MetaData
On the MetaData directory, you will find the JSON files that are used to produce microAOD with FLASHgg from the miniAOD datasets. We only need to produce the signal microAOD for ourselves, the rest should be inherited from Hgg. For more information about how to create MicroAOD, see the MetaData README.md.

#### Test
On the test directory, you will find the different scripts to run our analysis, such as the flat tree producer, plots producer, etc. See directory README for more information.
