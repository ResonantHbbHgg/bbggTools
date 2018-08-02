PROD='withTTHBDT'

fggRunJobs.py --load jsons/80X/RunJobs_NonResSignal_HHbbggSignal_Moriond17_v1.json -m 0 -n 1 -d Signal_$PROD -x cmsRun MakeTrees_FLASHgg_MC_2017.py maxEvents=-1 --no-use-tarball -q 1nd

fggRunJobs.py --load jsons/80X/RunJobs_Data_ReMiniAOD2016.json -m 0 -n 1 -d Data_$PROD -x cmsRun MakeTrees_FLASHgg_data_2017.py maxEvents=-1 --no-use-tarball -q 1nw --hadd --njobs 50

fggRunJobs.py --load jsons/80X/RunJobs_Background_ReMiniAOD2016.json -m 0 -n 10 -d BKG_$PROD -x cmsRun MakeTrees_FLASHgg_MC_2017.py maxEvents=-1 --no-use-tarball -q 1nd

fggRunJobs.py --load jsons/80X/RunJobs_tteebbnunu.json -m 0 -n 100 -d ttbbllnunu_$PROD -x cmsRun MakeTrees_FLASHgg_MC_2017.py maxEvents=-1 --no-use-tarball -q 1nd



