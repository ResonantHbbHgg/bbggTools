#!/bin/bash

sed "s/XXX/Rad/g" TEMPLATE_MakeTrees.py > TEMPLATERAD_MakeTrees.py
for RadMass in 320 340 350 400 600 650 700
do
	echo "Running on Radion sample with M = "$RadMass
	sed "s/MASS/${RadMass}/g" TEMPLATERAD_MakeTrees.py > RadMass_${RadMass}_MakeTrees.py
	cat > RadMass_${RadMass}_MakeTrees.sh << EOF
	#!/bin/bash
	cd /afs/cern.ch/work/r/rateixei/work/DiHiggs/Jun15/CMSSW_7_4_1/src/flashgg/bbggTools/test
	eval \`scramv1 runtime -sh\`
	cmsRun RadMass_${RadMass}_MakeTrees.py
EOF
	chmod a+rwx RadMass_${RadMass}_MakeTrees.sh
	bsub -J Tree_R${RadMass} -q 1nh < RadMass_${RadMass}_MakeTrees.sh 
done
rm TEMPLATERAD_MakeTrees.py
#rm RadMass*MakeTrees.py

sed "s/XXX/Grav/g" TEMPLATE_MakeTrees.py > TEMPLATEGRAV_MakeTrees.py
for GravMass in 260 270 280 320 350 500 550
do
        echo "Running on Graviton sample with M = "$GravMass
        sed "s/MASS/${GravMass}/g" TEMPLATEGRAV_MakeTrees.py > GravMass_${GravMass}_MakeTrees.py
	cat > GravMass_${GravMass}_MakeTrees.sh << EOF
        #!/bin/bash
        cd /afs/cern.ch/work/r/rateixei/work/DiHiggs/Jun15/CMSSW_7_4_1/src/flashgg/bbggTools/test
        eval \`scramv1 runtime -sh\`
	cmsRun GravMass_${GravMass}_MakeTrees.py
EOF
	chmod a+rwx GravMass_${GravMass}_MakeTrees.sh
	bsub -J Tree_G${GravMass} -q 1nh < GravMass_${GravMass}_MakeTrees.sh 
done
rm TEMPLATEGRAV_MakeTrees.py
#rm GravMass*MakeTrees.py
