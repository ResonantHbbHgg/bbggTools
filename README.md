bbggTools
=========

Welcome to the HH->bbgg RunII analysis!

This package is a set of tools to be used to perform the resonant and non-resonant HH->bbgg analysis for RunII.
It should be used under the FLASHgg framework and not as a standalone package.

Before cloning this package, get the FLASHgg framework: https://github.com/cms-analysis/flashgg  
Also check their Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/FLASHggFramework#Instructions_for_users  
ALWAYS FOLLOW THE INSTRUCTIONS LISTED ON "INSTRUCTIONS FOR USERS".  

Follow their instructions (including compilation).

Then, go to the flashgg directory and clone bbggTools, and compile again:
```
cd $CMSSW_BASE/src/flashgg
git clone git@github.com:ResonantHbbHgg/bbggTools.git bbggTools
cd bbggTools
scramv1 b -j 10
cmsenv
bbggAfterBuild.py
```   
   
Attention: run bbggAfterBuild.py only once. It will copy the directories under bbggTools/MetaData to flashgg/MetaData/data (signal datasets) and will add the signal cross sections to the flashgg list of cross sections.   


## Important steps

### Produce MicroAOD

1) Go to flashgg/MetaData/work

2) Create a json file that lists the datasets you want to process, for example:  
```
{
    "data" : ["/DoubleElectron/CMSSW_7_0_6_patch1-GR_70_V2_AN1_RelVal_zEl2012D-v1/MINIAOD"
              ],
    "sig"  : ["/GluGluToHToGG_M-125_13TeV-powheg-pythia6/Spring14miniaod-PU20bx25_POSTLS170_V5-v2/MINIAODSIM",
              ],
    "bkg"  : ["/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM",
              ]
}
```  
You can find such a json file, with all of our signal samples, under bbggTools/MetaData/SignalForMicroAOD.json. All needed background is processed centrally, so no need to worry about that.  

It is important to list the signal samples under "sig". This will trigger the signal customizations from FLASHgg, which includes the calculation of the PDF weights.  

Because of some differences regarding the storage of PDF weights in the LHE headers created by different generators, running on our signal samples **is not out of the box**. There are available unstructions for two different versions of the code:  

A) For **Spring15BetaV7**:  
You need to change the input pdfset from "PDF_variation" to "NNPDF30_lo_as_0130_nf_4.LHgrid" (the PDF set of our signal samples) in https://github.com/cms-analysis/flashgg/blob/Spring15BetaV7/MicroAOD/python/flashggPDFWeightObject_cfi.py#L5   

B) For the new versions (**1_1_X**):  
The functionality to pass the PDF set name was removed, so running our signal samples with the signal customization proved painful. Unfortunately, in this tag, we need to run without the PDF weights. So, list the signal samples under the "bkg" section of the json file.  

C) For upcoming versions:   
The FLASHgg guys have been kind enough to restore that functionality for us. Once the pull request with the needed modifications is merged, I'll update the instructions.   

3) After these modifications, you need to prepare the crab jobs:   
```
./prepareCrabJobs.py -C <Campaign> -U 5 -L 25 -s <signal json> -V <flashggVersion> -p ${CMSSW_BASE}/src/flashgg/MicroAOD/test/microAODstd.py --outputPath /store/group/phys_higgs/resonant_HH/RunII/MicroAOD
```   
Try to use a Campaign name that is different from the standard FLASHgg campaigns, so that a different catalog is created (below). For more information on these commands, see flashgg/MetaData.   

To submit everything, do:   
```
cd <microAODCampaignName>
echo crabConfig_*.py | xargs -n 1 crab sub
```   
There is the possibility that Crab won't like the configuration files because the request name is too long. In case that happens, do something similar to:   
```
sed -i 's/1_1_1_GluGluToBulkGraviton/1_1_1_Grav/' *
sed -i 's/1_1_1_GluGluToRadion/1_1_1_Rad/' *
sed -i 's/narrow_13TeV-madgraph_RunIISpring/13TeV-madgraph_RunIISpring/' *
```   
Where "1_1_1" is the FLASHgg version. This has no impact on the crab jobs themselves, only on the name of the job.   

#### Update Catalog With Produced Samples
(These instructions are explained more here: https://github.com/cms-analysis/flashgg/tree/master/MetaData)   

In order to run the analyzer in the newly created MicroAODs, you need to create a catalog with these files after your crab jobs are done. To do that, you need to run:   
```
fggManageSamples.py -C <Campaign> -V <flashggVersion> import
```  
This will create a datasets.json file under flashgg/MetaData/data/Campaign/. If you think there might be duplicated samples on your catalog, do:   
```
fggManageSamples.py -C <Campaign> review
```   
Finally, in order to have the weights on the samples catalog, do:
```
fggManageSamples.py -C <Campaign> check
```   
The datasets catalog file created will then be used to produce selection trees.   

**It is also good to include the newly created samples in the cross-section database under flashgg/MetaData/data/cross_sections.json. A file like that one, with the available samples today (Jan 14) can be found under bbggTools/MetaData/cross_sections.json. Just copy the file from our repo on flashgg!**

### Produce Selection Trees

#### Produce Selection Trees Using FLASHgg Functionalities
```
cd bbggTools/test/RunJobs/
```   
Create a json file for the datasets you want to process <ToProcess.json>. Those already exist for our available signal samples, background samples and data under bbggTools/test/RunJobs/. Create a directory output for your jobs <ToProcessDir>. Then do:   
```
fggRunJobs.py --load ToProcess.json -H -D -P -n 500 -d ToProcessDir -x cmsRun MakeTrees_FLASHgg.py  maxEvents=-1 -q 1nh --no-use-tarball
```   
This will submit the jobs to LSF batch. You can quit the job watcher if you want, and reload or get a summary of the jobs it with:
```
fggRunJobs.py --load ToProcessDir/config.json --cont
fggRunJobs.py --load ToProcessDir/config.json --summary
```

#### Old way (not maintained)

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
