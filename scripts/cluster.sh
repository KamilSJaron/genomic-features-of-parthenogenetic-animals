#!/bin/bash

# This is a wrapper for a cluster scripts
# It's exected to be called from bsub
# the first argument is a script that should be executed
# the rest of the arguments are arguments of the script
# all all the other arguments will be copied locally if they are valid files
# the last argument got to be the output

WORKING_DIR=$(pwd)
LOCAL_DIR="/scratch/local/$USER/$JOBID"
mkdir -p "$LOCAL_DIR/temp"
export TMPDIR="$LOCAL_DIR/temp"

cd $LOCAL_DIR

for arg in "$@"; do
    # if file exists, copy it to
    if [[ -s "$arg" ]]
    then
        RELATIVE_PATH=$(dirname $arg)
        mkdir -p $LOCAL_DIR/$RELATIVE_PATH
        cp $arg $LOCAL_DIR/$RELATIVE_PATH
    fi
done

# get info about output
OUTPUT=$arg
RELATIVE_PATH=$(dirname $OUTPUT)
# create directory for the output (both in the work dir and local dir)
mkdir -p $LOCAL_DIR/$RELATIVE_PATH
mkdir -p $WORKING_DIR/$RELATIVE_PATH

# remove script name from arguemnts
SCRIPT=$1
shift 1;

cd $LOCAL_DIR
$SCRIPT $@

mv $OUTPUT $WORKING_DIR/$OUTPUT