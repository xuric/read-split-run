#!/bin/bash

if [ ! $CONFIGURED ] || [[ "$CONFIGURED" == "" ]]; then
    CONFIGURED=true
else
    exit 0
fi

#Meta defines
DEBUG=true

yell() { echo "$(basename $0): $*" >&2; }
die() { yell "$*"; exit 111; kill $$; }
try() { "$@" || die "cannot $*"; }
dprint() { if $DEBUG; then yell "$*"; fi }


#-------USER CONFIGURATION------------
TEMP_BOWTIE_FILES="/scr/RSW/bowtie"
TEMP_SPLIT_FILES="/scr/RSW/tmp"
TEMP_RSW_FILES="/scr/RSW/tmp"
BASES_TO_TRIM=2
BOWTIE="/usr/local/bin/bowtie"
BOWTIE_INDEX_DIR="/net/project/common/data"  #base directory for bowtie indexes
REFDIR="/net/home/bdonham/public_html/dev/public/database/refs"    #where the refFlats are


#------NOT-USER CONFIGURATION--------
#CHANGE THESE AT YOUR OWN RISK

#this sets the location to find the sub-parts of the pipeline
#BASEDIR=$(dirname $(realpath $0))
BASEDIR=$( cd ${0%/*} && pwd -P )  #moved to this method to avoid platform-specific comamnds

#rsw_batch_job.sh:  The big guy!
PIPELINE="${BASEDIR}/rsw_pipeline.sh"
COMPARER="${BASEDIR}/rsw_compare.sh"

#rsw_pipeline.sh: Pipeline Constants
#Define program/script files
ALIGNER="${BASEDIR}/bowtie.sh"
SPLITTER="${BASEDIR}/split.sh"
MEASURER="${BASEDIR}/readlength.sh"
FORMATTER="${BASEDIR}/split_1stColumn_RSW.pl"
RSRER="${BASEDIR}/rsw.sh"   #-ER'd for consistency

#readlength.sh:  measuring tool
#BASES_TO_TRIM=2   #moved to USER CONFIGURATION section
OFFSET=$(( $BASES_TO_TRIM + 1))      #1 + the number of bases to trim
        #the 1 is to account for the newline when reading from the file.

#bowtie.sh:  Alignment constants
#ALIGN_TEMP_DIR="/scr/RSW/bowtie"  #replaced with TEMP_BOWTIE_FILES
ENCODING_GUESSER="${BASEDIR}/guess-encoding.py"
BOWTIE_INDEXES=""   #dummy value
#BOWTIE_INDEXES will be dynamically set in bowtie.sh
QUALITY_TESTS=1000

#split.sh:   reads splitting
OSPLITTER="${BASEDIR}/split_read_RSW.pl"    #Old splitter; not really used
NSPLITTER="${BASEDIR}/srr"                  #New spiltter; actually used

#rsw.sh:   split pair finder
#TMPDIR="/scr/RSW/tmp"  #TEMP_RSW_FILES
#REFDIR=$(realpath "$(dirname $(realpath $0))/../refs")    #where the refFlats are
#REFDIR moved to USER CONFIGURATION section
#RSR="${BASEDIR}/jeff_rsw"  #actual program
RSR="${BASEDIR}/sp4"
RSW_TIMING_FILE="rsw.out"
OUTPUTFILE=""    #dummy value; filled in within the scripts

#rsw_compare.sh :  compare the output of two jobs
COMPARE_PROG="${BASEDIR}/rsw_compare"
LINES_TO_SKIP=22


#cleanup.sh:  gets rid of old files
#none. everything is defined above


