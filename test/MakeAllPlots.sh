#!/bin/bash

sed "s/XXX/Rad/g" TEMPLATE_MakePlots.py > TEMPLATERAD_MakePlots.py
for RadMass in 320 340 350 400 600 650 700
do
	echo "Running on Radion sample with M = "$RadMass
	sed "s/MASS/${RadMass}/g" TEMPLATERAD_MakePlots.py > RadMass_${RadMass}_MakePlots.py
	cmsRun RadMass_${RadMass}_MakePlots.py >& LOG_RadMass_${RadMass}.out &
done
rm TEMPLATERAD_MakePlots.py
#rm RadMass*MakePlots.py

sed "s/XXX/Grav/g" TEMPLATE_MakePlots.py > TEMPLATEGRAV_MakePlots.py
for GravMass in 260 270 280 320 350 500 550
do
        echo "Running on Graviton sample with M = "$GravMass
        sed "s/MASS/${GravMass}/g" TEMPLATEGRAV_MakePlots.py > GravMass_${GravMass}_MakePlots.py
	cmsRun GravMass_${GravMass}_MakePlots.py >& LOG_GravMass_${GravMass}.out &
done
rm TEMPLATEGRAV_MakePlots.py
#rm GravMass*MakePlots.py
