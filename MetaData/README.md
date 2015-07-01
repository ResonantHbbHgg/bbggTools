Producing MicroAOD on FLASHgg with MetaData
===========================================

To produce microAOD from datasets in a JSON file, take the JSON file to flashgg/MetaData/work and do:
```
./prepareCrabJobs.py -s HHbbgg.json -o /store/user/rateixei/flashgg/ -C RunIISpring15DR74 -V RunIISpring15MicroAODV1
```
Where:  
-C = Campaign  
-V = Version  
-o = output folder on EOS

This will create a folder on flashgg/MetaData/work with all the crab submition files. To submit all of them at the same time, do:
```
echo crabConfig_*.py | xargs -n 1 crab sub
```

To check the status of all of them at the same time, do:
```
echo crab_* | xargs -n 1 crab status
```
