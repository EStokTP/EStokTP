#!/bin/bash

if [ "$#" != "3" ] 
then
    echo " usage  ./get_files.com  #from #to #directory"
    exit
fi

FILE1=$2
FILE2=$1
DIR_FROM=$3

echo "this will copy all the files for  $FILE2"
echo "taking the data from files of type $FILE1"
echo "located in the directory  $DIR_FROM"
echo "do you want to proceed? (type yes)"

read answer

if [ "$answer" = "yes" ]
then
    echo "dir_from is  $DIR_FROM"
    ls $DIR_FROM/output/${FILE1}_opt_* > temp.log
    grep -c ^ temp.log > num.log
    NSTRU=$(<num.log)
    echo "nstru is $NSTRU"

    cp -f $DIR_FROM/data/$FILE1.dat ./data/$FILE2.dat
    
    for index in $(seq "$NSTRU" )
    do
	if [ "$index" -lt "10" ] 
	then
	    cp -f $DIR_FROM/geoms/${FILE1}_0$index.xyz ./geoms/${FILE2}_0$index.xyz
	    echo "index is 0$index"
	fi

	if [ "$index" -gt "9" ] 
	then
	    cp -f $DIR_FROM/geoms/${FILE1}_$index.xyz ./geoms/${FILE2}_$index.xyz
	    echo "index is 0$index"
	fi
    done 

    cp -f $DIR_FROM/geoms/${FILE1}_l1.log ./geoms/${FILE2}_l1.log
    cp -f $DIR_FROM/geoms/${FILE1}_l1.xyz ./geoms/${FILE2}_l1.xyz
    [[ -e $DIR_FROM/geoms/${FILE1}_l1.molden ]]  &&  cp -f $DIR_FROM/geoms/${FILE1}_l1.molden ./geoms/${FILE2}_l1.molden

#FILE1='ts'
#FILE2='ts'
#echo "file1 is $FILE1"

    if [ "$FILE1" = "ts" ] 
    then
#echo "file1 is $DIR_FROM/geoms/${FILE1}gta_l1.log"
	cp -f $DIR_FROM/geoms/${FILE1}gta_l1.log ./geoms/${FILE2}gta_l1.log
	cp -f $DIR_FROM/geoms/${FILE1}gta_l1.xyz ./geoms/${FILE2}gta_l1.xyz
	cp -f $DIR_FROM/geoms/${FILE1}gta_l1.molden ./geoms/${FILE2}gta_l1.molden
    fi


    cp -f $DIR_FROM/output/hrdata4proj_${FILE1}.dat ./output/hrdata4proj_${FILE2}.dat
    cp -f $DIR_FROM/output/hrproj_freq_${FILE1}.out ./output/hrproj_freq_${FILE2}.out
    cp -f $DIR_FROM/output/${FILE1}_fcmat.log ./output/${FILE2}_fcmat.log
    cp -f $DIR_FROM/output/${FILE1}_hr.out ./output/${FILE2}_hr.out
    cp -f $DIR_FROM/output/${FILE1}_opt.out ./output/${FILE2}_opt.out

    for index in $(seq "$NSTRU" )
    do
	if [ "$index" -lt "10" ] 
	then
	    cp -f $DIR_FROM/output/${FILE1}_opt_0$index.out ./output/${FILE2}_opt_0$index.out
	fi
	
	if [ "$index" -gt "9" ] 
	then
	    cp -f $DIR_FROM/output/${FILE1}_opt_$index.out ./output/${FILE2}_opt_$index.out
	fi
    done 


    if [ "$FILE1" = "ts" ] 
    then
	cp -f $DIR_FROM/output/tauo_${FILE1}.out ./output/tauo_${FILE2}.out
	cp -f $DIR_FROM/output/${FILE1}_zmat.out ./output/${FILE2}_zmat.out
	cp -f $DIR_FROM/irc_files/* ./irc_files
	cp -f $DIR_FROM/md_tunn/* ./md_tunn
    fi

    cp -f $DIR_FROM/hl_logs/${FILE1}_molpro.out ./hl_logs/${FILE2}_molpro.out


    cp -f $DIR_FROM/me_files/${FILE1}_1dge.me ./me_files/${FILE2}_1dge.me
    cp -f $DIR_FROM/me_files/${FILE1}_1dhr.me ./me_files/${FILE2}_1dhr.me
    cp -f $DIR_FROM/me_files/${FILE1}_en.me ./me_files/${FILE2}_en.me
    cp -f $DIR_FROM/me_files/${FILE1}_fr.me ./me_files/${FILE2}_fr.me
    cp -f $DIR_FROM/me_files/${FILE1}_ge.me ./me_files/${FILE2}_ge.me
    cp -f $DIR_FROM/me_files/${FILE1}_hr.me ./me_files/${FILE2}_hr.me
    cp -f $DIR_FROM/me_files/${FILE1}_nofr.me ./me_files/${FILE2}_nofr.me
    cp -f $DIR_FROM/me_files/${FILE1}_unpfr.me ./me_files/${FILE2}_unpfr.me
    cp -f $DIR_FROM/me_files/${FILE1}_zpe.me ./me_files/${FILE2}_zpe.me

    if [ "$FILE1" = "reac1" ] 
    then
	if [ "$FILE2" = "reac2" ] 
	then
#	    sed -ie 's/Bimolecular REACS//g' ./me_files/reac2_ge.me
	    sed -i '1d' ./me_files/reac2_ge.me
	    sed -ie 's/REACT1/REACT2/g' ./me_files/reac2_ge.me
	    echo ' GroundEnergy[kcal/mol] 0.0
 End
 !************************************' > temp.log
	    cp -f ./me_files/reac2_fr.me .
	    sed -i '$ d' reac2_fr.me
	    cat reac2_fr.me temp.log > ./me_files/reac2_fr.me
	fi

	if [ "$FILE2" = "prod1" ] 
	then
	    sed -ie 's/Well REACS/Bimolecular PRODS/g' ./me_files/prod1_ge.me
	    sed -ie 's/Species/Fragment PROD1/g' ./me_files/prod1_ge.me
	    sed -ie 's/REACT1/PROD1/g' ./me_files/prod1_ge.me
	    sed -ie 's/REACS/PRODS/g' ./me_files/prod1_ge.me
	fi
	if [ "$FILE2" = "prod2" ] 
	then
#	    sed -ie 's/Bimolecular REACS//g' ./me_files/prod2_ge.me
	    sed -i '1d' ./me_files/prod2_ge.me
	    sed -ie 's/REACT1/PROD2/g' ./me_files/prod2_ge.me
	    echo ' GroundEnergy[kcal/mol] $proden
 End
 !************************************' > temp.log
	    cp -f ./me_files/prod2_fr.me .
	    sed -i '$ d' prod2_fr.me
	    cat prod2_fr.me temp.log > ./me_files/prod2_fr.me
	fi
	if [ "$FILE2" = "wellp" ] 
	then
	    sed -ie 's/REACT1//g' ./me_files/wellp_ge.me
	    sed -ie 's/REACS/WP/g' ./me_files/wellp_ge.me
	    sed -ie 's/Bimolecular/Well/g' ./me_files/wellp_ge.me
	    sed -ie 's/Fragment/Species/g' ./me_files/wellp_ge.me

	    cp -f ./me_files/wellp_fr.me .
	    echo ' End' > temp.log
	    echo '!*********' >> temp.log
	    sed -i '$ d' wellp_fr.me
	    cat wellp_fr.me temp.log > ./me_files/wellp_fr.me
	fi
    fi
    if [ "$FILE1" = "reac2" ] 
    then
	if [ "$FILE2" = "reac1" ] 
	then
	    echo "Bimolecular REACS" > temp.log
	    cp -f ./me_files/reac1_ge.me .
	    cat temp.log reac1_ge.me > ./me_files/reac1_ge.me
	    sed -ie 's/REACT2/REACT1/g' ./me_files/reac1_ge.me
	    cp -f ./me_files/reac1_fr.me .
	    echo '!*********' > temp.log
	    sed -i '$ d' reac1_fr.me
	    sed -i '$ d' reac1_fr.me
	    sed -i '$ d' reac1_fr.me
	    cat reac1_fr.me temp.log > ./me_files/reac1_fr.me
	fi
	if [ "$FILE2" = "prod1" ] 
	then
	    echo "Bimolecular PRODS" > temp.log
	    cp -f ./me_files/prod1_ge.me .
	    cat temp.log prod1_ge.me > ./me_files/prod1_ge.me
	    sed -ie 's/REACT2/PROD1/g' ./me_files/prod1_ge.me
	    cp -f ./me_files/prod1_fr.me .
	    echo '!*********' > temp.log
	    sed -i '$ d' prod1_fr.me
	    sed -i '$ d' prod1_fr.me
	    sed -i '$ d' prod1_fr.me
	    cat prod1_fr.me temp.log > ./me_files/prod1_fr.me
	fi
	if [ "$FILE2" = "prod2" ] 
	then
	    sed -ie 's/REACT2/PROD2/g' ./me_files/prod2_ge.me
	    sed -ie 's/\] 0.0/\] $proden/g' ./me_files/prod2_fr.me
	fi
    fi

    if [ "$FILE1" = "prod1" ] 
    then
	if [ "$FILE2" = "reac1" ] 
	then
	    sed -ie 's/PROD1/REACT1/g' ./me_files/reac1_ge.me
	    sed -ie 's/PRODS/REACS/g' ./me_files/reac1_ge.me
	fi
	if [ "$FILE2" = "reac2" ] 
	then
	    sed -i '1d' ./me_files/reac2_ge.me
#	    sed -ie 's/Bimolecular PRODS//g' ./me_files/reac2_ge.me
	    sed -ie 's/PROD1/REACT2/g' ./me_files/reac2_ge.me
	fi
	if [ "$FILE2" = "prod2" ] 
	then
	    sed -i '1d' ./me_files/prod2_ge.me
#	    sed -ie 's/Bimolecular REACS//g' ./me_files/prod2_ge.me
	    sed -ie 's/PROD1/PROD2/g' ./me_files/prod2_ge.me
	fi
    fi

    if [ "$FILE1" = "prod2" ] 
    then
	if [ "$FILE2" = "reac1" ] 
	then
	    echo "Bimolecular REACS" > temp.log
	    cp -f ./me_files/reac1_ge.me .
	    cat temp.log reac1_ge.me > ./me_files/reac1_ge.me
	    sed -ie 's/PROD2/REACT1/g' ./me_files/reac1_ge.me
	    cp -f ./me_files/reac1_fr.me .
	    echo '!*********' > temp.log
	    sed -i '$ d' reac1_fr.me
	    sed -i '$ d' reac1_fr.me
	    sed -i '$ d' reac1_fr.me
	    cat reac1_fr.me temp.log > ./me_files/reac1_fr.me
	fi
	if [ "$FILE2" = "prod1" ] 
	then
	    echo "Bimolecular PRODS" > temp.log
	    cp -f ./me_files/prod1_ge.me .
	    cat temp.log prod1_ge.me > ./me_files/prod1_ge.me
	    sed -ie 's/PROD2/PROD1/g' ./me_files/prod1_ge.me
	    cp -f ./me_files/prod1_fr.me .
	    echo '!*********' > temp.log
	    sed -i '$ d' prod1_fr.me
	    sed -i '$ d' prod1_fr.me
	    sed -i '$ d' prod1_fr.me
	    cat prod1_fr.me temp.log > ./me_files/prod1_fr.me
	fi
	if [ "$FILE2" = "reac2" ] 
	then
	    sed -ie 's/PROD2/REACT2/g' ./me_files/reac2_ge.me
	    sed -ie 's/\] $proden/\] 0.0/g' ./me_files/reac2_fr.me
	fi
    fi

    if [ "$FILE1" = "wellr" ] 
    then
	if [ "$FILE2" = "wellp" ] 
	then
	    sed -ie 's/WR/WP/g' ./me_files/wellp_ge.me
	fi
    fi
    if [ "$FILE1" = "wellp" ] 
    then
	if [ "$FILE2" = "wellr" ] 
	then
	    sed -ie 's/WP/WR/g' ./me_files/wellr_ge.me
	fi
    fi

elif [ "$answer" = "no" ]
then
    echo "aborted"
else
    echo "aborted"
fi

#cp -f $DIR_FROM/geoms/$FILE1'_01.dat' ./geoms/
