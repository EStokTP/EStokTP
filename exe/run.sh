#!/bin/bash

numdist=12
dist1=2.0
dist2=2.2
dist3=2.4
dist4=2.6
dist5=2.8
dist6=3.0
dist7=3.2
dist8=3.4
dist9=3.6
dist10=3.8
dist11=4.0
dist12=10.0

numdirs=$(($numdist-1))
irun=-1
if [ $1 == 'all' ]; then
    irun=0
fi
if [ $1 == 'step1' ]; then
    irun=1
fi
if [ $1 == 'step2' ]; then
    irun=2
fi
if [ $1 == 'step3' ]; then
    irun=3
fi
if [ $1 == 'step35' ]; then
    irun=35
fi
if [ $1 == 'step4' ]; then
    irun=4
fi
if [ $1 == 'step5' ]; then
    irun=5
fi
if [ $1 == 'step6' ]; then
    irun=6
fi

echo 'irun' $irun
#exit

if [ $irun -eq 1 ] || [ $irun -eq 0 ]; then
    cp ref/data/estoktp.dat_vrc ref/data/estoktp.dat
    cp -R ref dist_$dist1
    cd  dist_$dist1
    rm -f finished
    rm -f failed
    sed -ie "s/STEPN/Step1/g" data/estoktp.dat
    sed -ie "s/DISTN/$dist1/g" data/ts.dat
    sed -ie "s/nGrid_Opt_ts/Grid_Opt_ts/g" data/estoktp.dat
    sed -ie "ntauo_ts/Tauo_ts/g" data/estoktp.dat
    sed -ie "s/nOpt_ts_1/Opt_ts_1/g" data/estoktp.dat
    sed -ie "s/n1dTau_ts/1dTau_ts/g" data/estoktp.dat
    sed -ie "s/nHL_ts/HL_ts/g" data/estoktp.dat
    sed -ie "s/nSymm_ts/Symm_ts/g" data/estoktp.dat

runes2k2
#estoktp.x
    cd ..

fi

if [ $irun -eq 1 ]; then
    exit 0
fi


#for npoint in $(seq "$numsite" )
#do

if [ $irun -eq 0 ]; then
    while true; do
    # flag represents if the app is ready, 1=ready, 0=not ready
	is_ready=0
	FILE=./dist_$dist1/finished
	if test -f "$FILE"; then
	    echo "$FILE "
	    is_ready=1
	fi
	FILE2=./dist_$dist1/failed
	if test -f "$FILE2"; then
            echo "$FILE2 "
            is_ready=2
	fi

	if [[ $is_ready -eq 2 ]]; then
            echo "Failed at Step 1"
            echo "stopping calculation"
            exit
	fi

	if [[ $is_ready -eq 1 ]]; then
            echo "Proceed to step 2"
            break
	else
        # idle for 60 seconds
            sleep 60
	fi
    done
fi

# Step 2

if [ $irun -eq 2 ] || [ $irun -eq 0 ]; then
    for npoint in $(seq "$numdirs" )
    do
	nprog=$(($npoint+1))
	var="dist$nprog"
	name=dist_"${!var}"
	cp -R dist_$dist1 $name
	cp -f ref/data/estoktp.dat ./$name/data
	cp -f ref/data/ts.dat ./$name/data
	cd ./$name
	sed -ie "s/STEPN/Step2/g" data/estoktp.dat
	sed -ie "s/DISTN/${!var}/g" data/ts.dat
	sed -ie "s/nOpt_ts_1/Opt_ts_1/g" data/estoktp.dat
	sed -ie "s/n1dTau_ts/1dTau_ts/g" data/estoktp.dat
	sed -ie "s/nHL_ts/HL_ts/g" data/estoktp.dat
	sed -ie "s/nSymm_ts/Symm_ts/g" data/estoktp.dat
	rm -f ./finished
	rm -f ./failed
	if [[ $npoint -eq $numdirs ]]; then
    	    sed -ie "s/Step2/Step21/g" data/estoktp.dat
	    sed -ie "s/1dTau_ts/n1dTau_ts/g" data/estoktp.dat
	    sed -ie "s/{frequencies;print,hessian}//g" output/ts_asl1_step2.inp
	fi
runes2k2
#estoktp.x
	cd ..
    done
fi
# now check if step2 calculations finished

if [ $irun -eq 0 ]; then
    while true; do
    # flag represents if the app is ready, 1=ready, 0=not ready
	is_ready=0
	numfinished=0
	for npoint1 in $(seq "$numdirs" )
	do
	    nprog=$(($npoint1+1))
	    var="dist$nprog"
	    name=dist_"${!var}"
	    FILE=./$name/finished
	    if test -f "$FILE"; then
		echo "$FILE "
		is_ready=1
		numfinished=$(($numfinished+1))
		echo "$numfinished"
	    fi
            FILE2=./$name/failed
            if test -f "$FILE2"; then
		echo "$FILE2 "
		is_ready=1
		numfinished=$(($numfinished+1))
		echo "$numfinished"
            fi
	done
	if [[ $numfinished -eq $numdirs ]]; then
            echo "Proceed to step 3"
            break
	else
        # idle for 60 seconds
            sleep 60
	fi
    done
fi

if [ $irun -eq 2 ]; then
    exit 0
fi

# now proceed to step3 calculations 

# Step 3

if [ $irun -eq 3 ] || [ $irun -eq 0 ]; then

    var="dist$numdist"
    name=dist_"${!var}"
    ln -s $name ./100

    for npoint in $(seq "$numdirs" )
    do
	nprog=$(($npoint+1))
	var="dist$npoint"
	name=dist_"${!var}"
	FILE2=./$name/failed
	if [ ! -f "$FILE2" ]; then
#    cp -R dist_$dist1 $name
	    cp -f ref/data/estoktp.dat ./$name/data
	    cp -f ref/data/ts.dat ./$name/data
	    cd ./$name
	    sed -ie "s/STEPN/Step3/g" data/estoktp.dat
	    sed -ie "s/DISTN/${!var}/g" data/ts.dat
#	sed -ie "s/DISTN/$var/g" data/ts.dat
#    sed -ie "s/nOpt_ts_1/Opt_ts_1/g" data/estoktp.dat
#    sed -ie "s/n1dTau_ts/1dTau_ts/g" data/estoktp.dat
            sed -ie "s/nHL_ts/HL_ts/g" data/estoktp.dat
#            sed -ie "s/nSymm_ts/Symm_ts/g" data/estoktp.dat
	    sed -ie "s/nkTP/kTP/g" data/estoktp.dat
	    rm -f ./finished
#estoktp.x	    
runes2k2
	    cd ..
	fi
    done
fi

# now check if step3 calculations finished

if [ $irun -eq 0 ]; then
    while true; do
    # flag represents if the app is ready, 1=ready, 0=not ready
    is_ready=0
    numfinished=0
    for npoint1 in $(seq "$numdirs" )
    do
	nprog=$(($npoint1))
	var="dist$nprog"
	name=dist_"${!var}"
	FILE=./$name/finished
	if test -f "$FILE"; then
	    echo "$FILE "
	    is_ready=1
	    numfinished=$(($numfinished+1))
	    echo "$numfinished"
	fi
        FILE2=./$name/failed
        if test -f "$FILE2"; then
            echo "$FILE2 "
            is_ready=1
            numfinished=$(($numfinished+1))
            echo "$numfinished"
        fi
    done
#  check this conditions, why ge? can't remember but seems needed
    if [[ $numfinished -ge $numdirs ]]; then
        echo "Proceed to step 4"
        break
    else
        # idle for 60 seconds
        sleep 60
    fi
    done
fi

# now proceed to step3.5 calculations 
# all the rates should have been computed

# assemble and do variational analysis

if [ $irun -eq 35 ] || [ $irun -eq 0 ]; then

    cp -R dist_$dist1 variational
    cp -f ref/data/estoktp.dat ./variational/data
    cd variational
    sed -ie "s/nvariational/variational/g" data/estoktp.dat
    sed -ie "s/nkTP/kTP/g" data/estoktp.dat
    sed -ie "s/STEPN/Step3/g" data/estoktp.dat

    name=dist_$dist1
    var=$dist1
    echo $name > blocks.dat
    nfiles=$(($numdirs-1))
    ndata=1

    for npoint in $(seq "$nfiles" )
    do
	nprog=$(($npoint+1))
	var="dist$nprog"
	name=dist_"${!var}"
	FILE2=../$name/failed
	if [ ! -f "$FILE2" ]; then
    	    ndata=$(($ndata+1))
	    echo $name >> blocks.dat
	fi
    done

    echo $ndata >> mep_dist.out
    echo $dist1 >> mep_dist.out 

    for npoint in $(seq "$nfiles" )
    do
	nprog=$(($npoint+1))
	var="dist$nprog"
	name=dist_"${!var}"
	FILE2=../$name/failed
	if [ ! -f "$FILE2" ]; then
            echo "${!var}" >> mep_dist.out
	fi
    done


    cp -f blocks.dat data
    cp -f mep_dist.out output
    /home/chimica2/cavallotti/programs/estoktp/exe/assemble_var.com
runes2k2
#estoktp.x

    cd ..
fi


# now proceed to step 4 to prepare VRC-TST required data for correction potential

if [ $irun -eq 4 ] || [ $irun -eq 0 ]; then

    for npoint in $(seq "$numdist" )
    do
#    nprog=$(($npoint+1))
	var="dist$npoint"
	name=dist_"${!var}"
#    cp -R dist_$dist1 $name
	FILE2=./$name/failed
	if [ ! -f "$FILE2" ]; then
	    cp -f ref/data/estoktp.dat ./$name/data
	    cp -f ref/data/ts.dat ./$name/data
	    cd ./$name
	    sed -ie "s/STEPN/Step4/g" data/estoktp.dat
	    sed -ie "s/DISTN/${!var}/g" data/ts.dat

#	sed -ie "s/DISTN/$var/g" data/ts.dat
#    sed -ie "s/nOpt_ts_1/Opt_ts_1/g" data/estoktp.dat
#    sed -ie "s/n1dTau_ts/1dTau_ts/g" data/estoktp.dat
#    sed -ie "s/nHL_ts/HL_ts/g" data/estoktp.dat
#    sed -ie "s/nSymm_ts/Symm_ts/g" data/estoktp.dat
	    sed -ie "s/npot_corr/pot_corr/g" data/estoktp.dat
	    if [[ $npoint -eq $numdist ]]; then
        	sed -ie "s/pot_corr/pot_corr inf/g" data/estoktp.dat
    	    fi	
	    rm -f ./finished
runes2k2
	    cd ..
	fi
    done
fi

# now check if step 4 calculations finished

if [ $irun -eq 0 ]; then

    while true; do
    # flag represents if the app is ready, 1=ready, 0=not ready
	is_ready=0
	numfinished=0
	for npoint1 in $(seq "$numdist" )
	do
	    nprog=$(($npoint1))
	    var="dist$nprog"
	    name=dist_"${!var}"
	    FILE=./$name/finished
	    if test -f "$FILE"; then
		echo "$FILE "
		is_ready=1
		numfinished=$(($numfinished+1))
		echo "$numfinished"
	    fi
            FILE2=./$name/failed
            if test -f "$FILE2"; then
		echo "$FILE2 "
		is_ready=1
		numfinished=$(($numfinished+1))
		echo "$numfinished"
            fi
	done
#  check this conditions, why ge? can't remember but needed
	if [[ $numfinished -ge $numdist ]]; then
            echo "Proceed to step 5"
            break
	else
        # idle for 60 seconds
            sleep 60
	fi
    done
fi

# now we are ready for step 5 to run the vrc-tsts calculations

if [ $irun -eq 5 ] || [ $irun -eq 0 ]; then

    cp -R variational vrctst
    cp -f ref/data/estoktp.dat ./vrctst/data

    cd vrctst

    sed -ie "s/nvariational/variational rotd_sr rotd_lr/g" data/estoktp.dat
    sed -ie "s/nvrc_tst/vrc_tst/g" data/estoktp.dat
    sed -ie "s/STEPN/Step5/g" data/estoktp.dat
    rm -f ./finished
    estoktp.x
    #runes2k2
    cp -f ./data/machines vrc_tst
    cd vrc_tst
    mkdir lr
    mkdir lr/scratch
    mkdir sr
    mkdir sr/scratch

    cp convert.inp lr
    cp convert.inp sr
    cp -f divsur_lr.inp lr/divsur.inp 
    cp -f divsur_sr.inp sr/divsur.inp 
    cp machines lr
    cp machines sr
    cp mc_flux.inp lr
    cp mc_flux.inp sr
    cp molpro.inp lr
    cp molpro.inp sr
    cp pot.dat lr
    cp pot.dat sr
    cp structure.inp lr
    cp structure.inp sr
    cp -f tst_lr.inp lr/tst.inp 
    cp -f tst_sr.inp sr/tst.inp 
    cp vrc_molpro.tml lr
    cp vrc_molpro.tml sr
 
    cd ..
!    cd ..
!    cd ..

    echo "completed calculation"

fi

# now check if step 5 calculations finished

if [ $irun -eq 0 ]; then

    while true; do
    # flag represents if the app is ready, 1=ready, 0=not ready                                              
                      
	is_ready=0
	FILE=./finished
	if test -f "$FILE"; then
            echo "$FILE "
            is_ready=1
	fi

	if [[ $is_ready -eq 1 ]]; then
            echo "Proceed to step 6"
            break
	else
        # idle for 60 seconds                                                                                
                     
            sleep 60
	fi
    done
fi




