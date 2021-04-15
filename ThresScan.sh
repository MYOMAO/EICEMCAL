#!/bin/sh
i=17
N=18

ethresinit=0.0
step=0.001

while [ $i -lt $N ]
do

    ethres=$(echo "${ethresinit} + ${step} * $i" | bc) ;

	echo "Index = " $i	
	echo "Energy Threshold Now " $ethres	"   GeV"

	echo "${ethres}" > EThresholdScan/EThres.dat


	root -b -l -q Fun4All_G4_EICDetector.C'('20000,1,50')'

	mv G4EICDetector.root_g4femc_eval.root EThresholdScan/G4EICDetector.root_g4femc_eval_${i}.root


	i=$(($i+1))

done


