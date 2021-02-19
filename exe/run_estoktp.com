#!/bin/sh

SSH=/usr/bin/ssh

HOST=$1
NODESLOW=$2
NODESHIGH=$3

#while test "x$1" != x ; do
# case $1 in
#  -h | --host           ) shift; HOST="$1"      ;;
#  *                     ) OPTIONS="$OPTIONS $1" ;;
# esac
# shift
#done

# echo $#

DIRECTORY='./data'

if [ ! -d "$DIRECTORY" ] 
then
    echo " "
    echo " it seems the subdirectory data is missing"
    echo " check and run again"
    exit 
fi


if [ "$#" != "1" ] 
then
    echo " usage  ./run_estoktp.com  node "
fi

if [ "$#" != "1" ] 
then
    echo " input file in ./data/estoktp.dat"
fi

if [ "$#" != "1" ] 
then
    exit 
fi


#echo "nodes low level = $NODESLOW"
#echo "nodes high level = $NODESHIGH"
echo "host = $HOST"
echo "$HOST" > host.dat

#cp -f data/estoktp_inp.dat data/estoktp.dat
#sed -ie 's/NODESLOW/'$NODESLOW'/g' data/estoktp.dat
#sed -ie 's/NODESHIGH/'$NODESHIGH'/g' data/estoktp.dat


#if [ "x$HOST" != "x" ]; then

# bash command
# add path to estktp exe subdirectory after 'export PATH=$PATH:~/bin;' if necessary e.g. 'export PATH=$PATH:~/ESTOKTPSUBD/exe;'
#  exec $SSH -n $HOST "cd `pwd`; export PATH=$PATH:~/bin; estoktp.x >& estoktp.log & "

# csh tcsh comman
# add path to estktp exe subdirectory after 'setenv PATH $PATH:~/bin'; if necessary e.g. 'setenv PATH $PATH:~/ESTOKTPSUBD/exe;'
  exec $SSH -n $HOST "cd `pwd`; setenv PATH $PATH:~/bin; estoktp.x >& estoktp.log & "

#else
#  estoktp.x >& estoktp.log &
#fi



