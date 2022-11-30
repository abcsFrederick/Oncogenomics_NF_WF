#!/bin/bash

set -e

SCRIPT_NAME="$0"
SCRIPT_DIRNAME=$(readlink -f $(dirname $0))
SCRIPT_BASENAME=$(basename $0)
WF_HOME=$SCRIPT_DIRNAME

CONFIG_FILE="$WF_HOME/nextflow.config"

# load singularity and nextflow modules
module load singularity nextflow


# set workDir ... by default it goes to `pwd`/work
# this can also be set using "workDir" in nextflow.config
# export NXF_WORK="/data/khanlab2/kopardevn/AWS_MVP_test/work"
export OUTDIR="/data/khanlab3/kopardevn/AWS_MVP_test"
export OUTTAG="5" # workDir will be $OUTDIR/work.$OUTTAG and resultsDir will be $OUTDIR/results.$OUTTAG and singularity cache is set to $OUTDIR/.singularity
export RESULTSDIR="$OUTDIR/results.$OUTTAG"
export WORKDIR="$OUTDIR/work.$OUTTAG"

# set .nextflow dir ... dont want this to go to $HOME/.nextflow
export NXF_HOME="$RESULTSDIR/.nextflow"


printenv|grep NXF

# run
if [ ! -d $RESULTSDIR ]; then mkdir -p $RESULTSDIR;fi
cd $RESULTSDIR
#nextflow run -profile biowulf main.nf -resume
nextflow \
	-C $CONFIG_FILE \
	run \
	-profile test_biowulf \
	$WF_HOME/main.nf -resume
