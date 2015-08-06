#!/bin/bash

sed "s/XXX/Rad/g" ${CMSSW_BASE}/src/flashgg/bbggTools/test/TEMPLATE_MakeTrees.py > TEMPLATERAD_MakeTrees.py
for RadMass in 320 340 350 400 600 650 700
do
	echo "Running on Radion sample with M = "$RadMass
	sed "s/MASS/${RadMass}/g" TEMPLATERAD_MakeTrees.py > RadMass_${RadMass}_MakeTrees.py
	cmsRun RadMass_${RadMass}_MakeTrees.py >& LOG_TREES_RadMass_${RadMass}.out &
done
rm TEMPLATERAD_MakeTrees.py
#rm RadMass*MakeTrees.py

sed "s/XXX/Grav/g" ${CMSSW_BASE}/src/flashgg/bbggTools/test/TEMPLATE_MakeTrees.py > TEMPLATEGRAV_MakeTrees.py
for GravMass in 260 270 280 320 350 500 550
do
        echo "Running on Graviton sample with M = "$GravMass
        sed "s/MASS/${GravMass}/g" TEMPLATEGRAV_MakeTrees.py > GravMass_${GravMass}_MakeTrees.py
	cmsRun GravMass_${GravMass}_MakeTrees.py >& LOG_TREES_GravMass_${GravMass}.out &
done
rm TEMPLATEGRAV_MakeTrees.py
#rm GravMass*MakeTrees.py
