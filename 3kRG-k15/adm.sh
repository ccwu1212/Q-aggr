#!/bin/bash

# This script makes it easier to organize multiple runs of ADMIXTURE for a range of K values

data=$1

minK=$2

maxK=$3

# Number of replications
Nrep=$4

# Additional options for Admixture
adm_options=$5

OUTDIR=$6

### Defaults
DEFAULT_OUTDIR="adm"
DEFAULT_ADM_OPTIONS=" --cv -j4 "

## Usage
if [ $# -lt 3 ];  then
	echo Usage:
	echo "adm.sh <input.bed> <minK> <maxK> [ <Nrep> [ <additional_options> [ <out_dir> ] ] ] "
	echo "Default output directory: $DEFAULT_OUTDIR "
	echo
	exit -1
fi

############  Parse arguments
startSeed=100


base=${data##*/}


#subdir=""

if [ -z "$OUTDIR" ]; then
	OUTDIR=$DEFAULT_OUTDIR
fi

if [ -z "$Nrep" ]; then
	Nrep=1
fi


################################### Start working
mkdir -p $OUTDIR
oldwd=`pwd`
cd $OUTDIR


# Run admixture with 
# 4 threads,
# CV error included
# no boostrapping
for (( rep = 1; rep <= Nrep ; rep ++ )); do
	# Each replication has its own random seed, for reproducibility
	seed=`expr $startSeed + $rep`
	repdir=rep${rep}-s$seed 
	echo "--- Replication $rep ( directory $repdir )"
	mkdir -p $repdir
	cd $repdir
	mkdir -p log
	for (( K=$minK ; K <= $maxK ; K++ )); do
		echo " running K=$K... "
		admixture $data $K -j4 --cv $adm_options \
		   -s $seed   > log/$K.out 2>log/$K.err 
	done
	echo "  -end of run."
	echo
	find log -size 0 -delete
	cd ..
done










