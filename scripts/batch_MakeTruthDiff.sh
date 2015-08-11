#!/bin/bash

queue="1nh"
dir=${PWD}

sed "s/XXX/Rad/g" ${CMSSW_BASE}/src/flashgg/bbggTools/test/TEMPLATE_MakeTruthDiff.py > TEMPLATERAD_MakeTruthDiff.py
for RadMass in 320 340 350 400 600 650 700
do
	echo "Running on Radion sample with M = "$RadMass
	sed "s/MASS/${RadMass}/g" TEMPLATERAD_MakeTruthDiff.py > RadMass_${RadMass}_MakeTruthDiff.py
	cat > RadMass_${RadMass}_MakeTruthDiff.sh << EOF
	#!/bin/bash
	cd ${dir}
	eval \`scramv1 runtime -sh\`
	cmsRun RadMass_${RadMass}_MakeTruthDiff.py
EOF
	chmod a+rwx RadMass_${RadMass}_MakeTruthDiff.sh
	bsub -J TrDif_R${RadMass} -q ${queue} < RadMass_${RadMass}_MakeTruthDiff.sh 
done
rm TEMPLATERAD_MakeTruthDiff.py
#rm RadMass*MakeTruthDiff.py

sed "s/XXX/Grav/g" ${CMSSW_BASE}/src/flashgg/bbggTools/test/TEMPLATE_MakeTruthDiff.py > TEMPLATEGRAV_MakeTruthDiff.py
for GravMass in 260 270 280 320 350 500 550
do
        echo "Running on Graviton sample with M = "$GravMass
        sed "s/MASS/${GravMass}/g" TEMPLATEGRAV_MakeTruthDiff.py > GravMass_${GravMass}_MakeTruthDiff.py
	cat > GravMass_${GravMass}_MakeTruthDiff.sh << EOF
        #!/bin/bash
        cd ${dir}
        eval \`scramv1 runtime -sh\`
	cmsRun GravMass_${GravMass}_MakeTruthDiff.py
EOF
	chmod a+rwx GravMass_${GravMass}_MakeTruthDiff.sh
	bsub -J TrDif_G${GravMass} -q ${queue} < GravMass_${GravMass}_MakeTruthDiff.sh 
done
rm TEMPLATEGRAV_MakeTruthDiff.py
#rm GravMass*MakeTruthDiff.py
