#!/bin/bash

filename='./data/blocks.dat'
numfiles=0
while read -r dir option
do
    if [ ! -z "$dir" ]; then
	numfiles=$(($numfiles+1))
	echo "taking data from directory ../$dir"
	echo "option is $option"
	if [[ $option == "var" ]]; then
	    cp -f ../$dir/me_blocks/variational.me ./ts1.me
	    sed -ie "s/Variational/ /g" ts1.me
	else
	    cp -f ../$dir/me_blocks/ts1.me .
	fi
#	exit
	sed -i '1d' ts1.me
	sed -i '$ d' ts1.me
	sed -i '$ d' ts1.me
#	sed -i '$ d' ts1.me
	cp -f ts1.me var_$numfiles.me
#	exit
	if [ $numfiles == 1 ]; then
	    cp -f ../$dir/me_blocks/ts1_en.me ./me_files
	    echo "0." > ./me_files/de1_TSvar.me
	fi
    fi
done <  $filename

echo "numfiles is $numfiles"

echo " Barrier TS REACS PRODS" > variational.me
echo " Variational" >> variational.me

for nstru in $(seq "$numfiles" )
do
    echo " !***** block $nstru" >> variational.me
    cat variational.me  var_$nstru.me > variational1.me
    cp -f variational1.me variational.me
done
#echo  " End" >> variational.me 
cp -f  variational.me ./me_files
rm -f var_*.me
rm -f  variational1.me
exit 0


