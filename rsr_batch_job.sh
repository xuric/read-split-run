#!/bin/bash
# RSR batch job.sh
# Read-Split-Run batch script for running the pipeline comparatively or not...


#Bring in the configuration data
#which had better be in the same directory as this script
#BASEDIR=$(dirname $(realpath $0))
BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )  #moved to this method to avoid platform-specific comamnds
source "${BASEDIR}/rsr_config.sh"

##--------Main-----------

if (( $# < 10 )); then
    yell "Not enough arguments."
    yell "usage: $(basename $0) mode genome readsFile [readsFile2] maxGoodAlignments minSplitSize minSplitdistance maxSplitdistance regionBuffer requiredSupports pathToSaveResults"
    die "params($#): $@"
fi

sanity_check $@

timing_start;

#step 0: Gather metadata
mode=$1
genome=$2
readsGroup1="$3"
readsGroup2="$4"

if [ "$mode" == "comparison" ];then
    #pop this argument off the list
    shift
fi

if [ ! -d "${10}" ]; then 
    mkdir "${10}"
elif [ -f "${10}" ]; then
    die "${10} already exists and is not a directory! Aborting."
fi

#pop mode, genome,readsfile1 off the arguments
shift;shift;shift

#now we set the boundary buffer, for use in the comparative module. it's guaranteed to be #5 now.
boundaryBuffer="$5"

#readsGroup2 will be either a set of reads for a comparative run
#or the minsplitsize. we'll only use readsGroup2 in the event that
#the mode was compare.

log "Job: $mode"
log "reads for batch 1: ${readsGroup1}"
if [ "$mode" == "comparison" ]; then
    log "reads for batch 2: ${readsGroup2}"
fi
log "other params: $@"

#with the variable variables shifted off argv, the remainder is valid input to the batch
#TODO: let this file handle the zip, email and cleanup steps
log "Processing batch for: ${readsGroup1}"
batch1=$(try $PIPELINE "$genome" "${readsGroup1}" "$@")
if (( $? )); then die "Failure."; fi
log "$batch1"

if [ "$mode" == "comparison" ]; then
    log "$mode : Processing batch for: ${readsGroup2}"
    log "$PIPELINE $genome ${readsGroup2} $@"
    batch2=$(${PIPELINE} "$genome" "${readsGroup2}" "$@")
    if (( $? )); then die "Failure."; fi
    log "Comparing $batch1 to $batch2"
    $COMPARE_SCRIPT "$boundaryBuffer" "$batch1" "$batch2"
fi

timing_end;

#Zip here
#Cleanup here
#email here

