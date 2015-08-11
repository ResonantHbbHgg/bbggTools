#!/bin/bash

sed "s/XXX/Rad/g" ${CMSSW_BASE}/src/flashgg/bbggTools/test/TEMPLATE_MakePlots.py > TEMPLATERAD_MakePlots.py
for RadMass in 320 340 350 400 600 650 700
do
	echo "Running on Radion sample with M = "$RadMass
	sed "s/MASS/${RadMass}/g" TEMPLATERAD_MakePlots.py > RadMass_${RadMass}_MakePlots.py
	cat > RadMass_${RadMass}_MakePlots.sh << EOF
	#!/bin/bash
	cd /afs/cern.ch/work/r/rateixei/work/DiHiggs/Jun15/CMSSW_7_4_1/src/flashgg/bbggTools/test
	eval \`scramv1 runtime -sh\`
	cmsRun RadMass_${RadMass}_MakePlots.py
EOF
	chmod a+rwx RadMass_${RadMass}_MakePlots.sh
	bsub -J Plots_R${RadMass} -q 1nh < RadMass_${RadMass}_MakePlots.sh 
done
rm TEMPLATERAD_MakePlots.py
#rm RadMass*MakePlots.py

sed "s/XXX/Grav/g" ${CMSSW_BASE}/src/flashgg/bbggTools/test/TEMPLATE_MakePlots.py > TEMPLATEGRAV_MakePlots.py
for GravMass in 260 270 280 320 350 500 550
do
        echo "Running on Graviton sample with M = "$GravMass
        sed "s/MASS/${GravMass}/g" TEMPLATEGRAV_MakePlots.py > GravMass_${GravMass}_MakePlots.py
	cat > GravMass_${GravMass}_MakePlots.sh << EOF
        #!/bin/bash
        cd /afs/cern.ch/work/r/rateixei/work/DiHiggs/Jun15/CMSSW_7_4_1/src/flashgg/bbggTools/test
        eval \`scramv1 runtime -sh\`
	cmsRun GravMass_${GravMass}_MakePlots.py
EOF
	chmod a+rwx GravMass_${GravMass}_MakePlots.sh
	bsub -J Plots_G${GravMass} -q 1nh < GravMass_${GravMass}_MakePlots.sh 
done
rm TEMPLATEGRAV_MakePlots.py
#rm GravMass*MakePlots.py
