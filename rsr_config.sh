#!/bin/bash

#------Only Configure the thing ONCE, please---
if [ -z "$CONFIGURED" ]; then
CONFIGURED=true
BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )  #moved to this method to avoid platform-specific comamnds

#-------USER CONFIGURATION------------
#-------Directories-------------------
BASE_TEMP_DIR="${BASEDIR}/tmp"
BOWTIE_TEMP_DIR="${BASE_TEMP_DIR}/bowtie"
SPLIT_TEMP_DIR="${BASE_TEMP_DIR}/splilt"
RSR_TEMP_DIR="${BASE_TEMP_DIR}/splitpairs"
LOG_DIR="${BASEDIR}/logs"
BASES_TO_TRIM=2
if [ -z "$BOWTIE_INDEXES" ]; then
    BOWTIE_INDEXES=""                            #See INSTALLATION
fi
BOWTIE_INDEX_ROOT=""                          #See INSTALLATION
REFDIR=""                                    #where the refFlats are


#------Programs----------------------
BOWTIE_PROGRAM="/usr/local/bin/bowtie"            # /path/to/bowtie
SPLIT_PROGRAM="${BASEDIR}/srr"                    # Program for splitting reads, compiled from split_read_rsr.c
FORMAT_PROGRAM="${BASEDIR}/sfc"                   # Program for formatting reads, compiled from split_first_column.c
RSR_PROGRAM="${BASEDIR}/sp4"                      # RSR Program ("split pairs"), compiled from splitPairs.cpp
COMPARE_PROGRAM="${BASEDIR}/rsr_compare"          # Program for comparing RSR outputs, compiled from 

#------NOT-USER CONFIGURATION--------
#CHANGE THESE AT YOUR OWN RISK

#Meta defines
DEBUG=false
export -p TERM=vt100                                        # this seesm to be necessary for bowtie...
if [ -z "$LOG_FILE" ] || [ -z "$RUN_ID" ]; then
export -p RUN_ID="$(date +%F.%s)"
export -p LOG_FILE="${LOG_DIR}/RSR_${RUN_ID}.log"
fi

yell() { echo "$(basename $0): $*" >&2; }
die() { yell "$*"; log "$*"; exit 111; kill $$; }
try() { "$@" || die "cannot $*"; }
dprint() { if $DEBUG; then yell "$*"; fi }
log() { echo "$(basename $0): $*" >> "${LOG_FILE}"; dprint "$@"; }
timing_start() { if [ -z "$START_TIME" ]; then START_TIME=$(date +%s); fi }
timing_end() { if [ ! -z "$START_TIME" ]; then log "Duration: $(( $(date +%s) - $START_TIME )) seconds"; fi }



#this sets the location to find the sub-parts of the pipeline
#BASEDIR=$(dirname $(realpath $0))

#rsr_batch_job.sh:  The big guy!
PIPELINE="${BASEDIR}/rsr_pipeline.sh"
COMPARE_SCRIPT="${BASEDIR}/rsr_compare.sh"

#rsr_pipeline.sh: Pipeline Constants
#Define program/script files
ALIGN_SCRIPT="${BASEDIR}/bowtie.sh"
MEASURE_SCRIPT="${BASEDIR}/readlength.sh"
SPLIT_SCRIPT="${BASEDIR}/split.sh"
#FORMAT_SCRIPT="${BASEDIR}/split_1stColumn_RSW.pl"  #there is no format script, only srr
RSR_SCRIPT="${BASEDIR}/splitPairs.sh"   #_SCRIPT'd for consistency

#readlength.sh:  measuring tool
#BASES_TO_TRIM=2   #moved to USER CONFIGURATION section
OFFSET=$(( $BASES_TO_TRIM + 1))      #1 + the number of bases to trim
#the 1 is to account for the newline when reading from the file.

#bowtie.sh:  Alignment constants
#ALIGN_TEMP_DIR="/scr/RSW/bowtie"  #replaced with TEMP_BOWTIE_FILES
ENCODING_GUESSER="${BASEDIR}/guess-encoding.py"
#BOWTIE_INDEXES will be dynamically set in bowtie.sh
QUALITY_TESTS=1000

#split.sh:   reads splitting
OSPLITTER="${BASEDIR}/split_read_RSW.pl"    #Old splitter; not really used
NSPLITTER="${BASEDIR}/srr"                  #New spiltter; actually used

#rsr.sh:   split pair finder
#TMPDIR="/scr/RSW/tmp"  #TEMP_RSW_FILES
#REFDIR=$(realpath "$(dirname $(realpath $0))/../refs")    #where the refFlats are
#REFDIR moved to USER CONFIGURATION section
#RSR="${BASEDIR}/jeff_rsr"  #actual program
RSR="${BASEDIR}/sp4"
RSR_TIMING_FILE="rsr.out"
OUTPUTFILE=""    #dummy value; filled in within the scripts

#rsr_compare.sh :  compare the output of two jobs
COMPARE_PROG="${BASEDIR}/rsr_compare"
LINES_TO_SKIP=22


#cleanup.sh:  gets rid of old files
#none. everything is defined above

#Manage directories
if [ ! -d "$BASE_TEMP_DIR" ]; then
    mkdir "$BASE_TEMP_DIR" || die "Could not make $BASE_TEMP_DIR aborting."
fi
if [ ! -d "$BOWTIE_TEMP_DIR" ]; then
    mkdir "$BOWTIE_TEMP_DIR" || die "Could not make $BOWTIE_TEMP_DIR aborting."
fi
if [ ! -d "$SPLIT_TEMP_DIR" ]; then
    mkdir "$SPLIT_TEMP_DIR" || die "Could not make $SPLIT_TEMP_DIR aborting."
fi
if [ ! -d "$RSR_TEMP_DIR" ]; then
    mkdir "$RSR_TEMP_DIR" || die "Could not make $RSR_TEMP_DIR aborting."
fi
if [ ! -d "$LOG_DIR" ]; then
    mkdir "$LOG_DIR" || die "Could not make $LOG_DIR aborting."
fi


#end if..configured
fi
