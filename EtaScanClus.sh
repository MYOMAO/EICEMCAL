#!/bin/sh
i=11
N=40

etainit=1.25
step=0.25

while [ $i -lt $N ]
do

    eta=$(echo "${etainit} + ${step} * $i" | bc) ;

	echo "eta now " $eta	

	root -b -l -q Fun4All_G4_EICDetector.C'('5000,1,${eta}')'

	mv G4EICDetector.root_g4femc_eval.root GammaPosAna/G4EICDetector.root_g4femc_eval_${i}.root
	mv ShowerInfo.root GammaPosAna/ShowerInfo_${i}.root



	i=$(($i+1))

done


