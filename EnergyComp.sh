#!/bin/sh
i=0
N=4
NEvent=4000
config="Quarter"

e=(1 2 4 8)



while [ $i -lt $N ]
do


	echo "E now " ${e[i]}	



	root -b -l -q Fun4All_G4_EICDetector.C'('${NEvent},${e[i]},2')'

#	mv ShowerInfo.root EnergyScanLay/Unity/ShowerInfo_${e[i]}.root
#	mv G4EICDetector.root_DSTReader.root EnergyScanLay/Unity/G4EICDetector.root_DSTReader_${e[i]}.root

	mv ShowerInfo.root EnergyScanLay/${config}/ShowerInfo_${e[i]}.root
	mv G4EICDetector.root_DSTReader.root EnergyScanLay/${config}/G4EICDetector.root_DSTReader_${e[i]}.root
	mv G4EICDetector.root_g4femc_eval.root EnergyScanLay/${config}/G4EICDetector.root_g4femc_eval_${e[i]}.root

	
	i=$(($i+1))


done




