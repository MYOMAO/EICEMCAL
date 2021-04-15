#!/bin/sh
i=0
Ni=4
j=0
Nj=3

e=(1 2 4 8)
thick=(10 14 15)


while [ $i -lt $Ni ]
do

	echo "Thickness now " ${thick[j]}	

	cp EnergyTestConfig/PHG4ForwardEcalDetector_${thick[i]}.cc  g4detectors/PHG4ForwardEcalDetector.cc
	

	cd g4detectors/
	autogen.sh --prefix=$MYINSTALL
	make
	make install
	cd ..  

	clear


	while [ $j -lt $Nj ]
	do


		echo "E now " ${e[j]}	



		root -b -l -q Fun4All_G4_EICDetector.C'('6000,${e[j]},2')'

		mv ShowerInfo.root CompLay/${e[j]}/ShowerInfo_35to${thick[i]}.root

		#eta=$(echo "${etainit} + ${step} * $i" | bc) ;

		j=$(($j+1))

	done

	i=$(($i+1))


done


