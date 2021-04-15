echo "-------------- CHANGING SOFTWARE NOW BRO --------------------"

source /opt/sphenix/core/bin/sphenix_setup.csh -n
mkdir install
setenv MYINSTALL $PWD/install/
setenv LD_LIBRARY_PATH $MYINSTALL/lib:$LD_LIBRARY_PATH
set path = ( $MYINSTALL/bin $path )

cd g4detectors/
#cd coresoftware/simulation/g4simulation/g4detectors/
autogen.sh --prefix=$MYINSTALL
make
make install
cd ../

echo "-------------- BUILD g4eval ------------------------"


#cd g4eval
#autogen.sh --prefix=$MYINSTALL
#make
#make install
#cd ..


#echo "-------------- DONE CHANGING SOFTWARE NOW BRO --------------------"
