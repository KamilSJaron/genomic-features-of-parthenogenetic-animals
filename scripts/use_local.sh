#!/bin/bash

# This is a wrapper for a cluster scripts
# It's exected to be called from bsub
# the first argument is a script that should be executed
# the rest of the arguments are arguments of the script
# all all the other arguments will be copied locally if they are valid files
# the last argument got to be the output

set -e

#####################
# PREPARE LOCAL DIR #
#####################

WORKING_DIR=$(pwd)
LONG_INDEX=`printf "%03d" $LSB_JOBINDEX`
LOCAL_DIR="/scratch/local/monthly/$USER/$LSB_JOBID.$LONG_INDEX"
mkdir -p "$LOCAL_DIR/temp"
export TMPDIR="$LOCAL_DIR/temp"

####################
# COPY INPUT FILES #
####################

for arg in "$@"; do
    # if file exists, copy it to
    if [[ -s "$arg" ]]
    then
        RELATIVE_PATH=$(dirname "$arg")
        mkdir -p $LOCAL_DIR/$RELATIVE_PATH
        cp -r $arg $LOCAL_DIR/$RELATIVE_PATH
    fi
done

###########################
# PREPARE FOR OUTPUT FILE #
###########################

# get info about output
OUTPUT=$arg
RELATIVE_PATH=$(dirname "$OUTPUT")
# create directory for the output (both in the work dir and local dir)
mkdir -p $LOCAL_DIR/$RELATIVE_PATH
mkdir -p $WORKING_DIR/$RELATIVE_PATH

# remove script name from arguemnts
SCRIPT=$1
shift 1;

###########
# RUN JOB #
###########

cd $LOCAL_DIR
$SCRIPT $@

###############################
# MOVE RESULTS TO WORKING_DIR #
###############################

if [ -d $OUTPUT ]; then
    mv $OUTPUT $WORKING_DIR/$OUTPUT
else
    mv $OUTPUT* $WORKING_DIR/$RELATIVE_PATH
fi

############
# CLEANING #
############

# remove script and the input files (the output is hopefully already saved)
rm $SCRIPT
for arg in "$@"; do
    # if tghe argument is an existing file, remove it
    if [[ -s "$arg" ]]
    then
        rm -r $arg
    fi
done
# remove all empty direcories
find ./* -depth -type d -exec rmdir {} \;
