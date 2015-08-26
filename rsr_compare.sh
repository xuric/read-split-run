#!/bin/bash

BASEDIR=$( cd ${0%/*}  >& /dev/null ; pwd -P )
source "${BASEDIR}/rsr_config.sh"

SUFFIX=".results"
OUTNAME="comparison.txt"

function do_compare() {
    try $RSR_COMPARE_PROG "$RSR_LINES_TO_SKIP" $RSR_BOUNDARY $OUTNAME $@ >& $RSR_LOG_FILE
    mv "${OUTNAME}.comparisonSummary.txt" "$RSR_DEST"
    mv "${OUTNAME}.comparedResults.txt" "$RSR_DEST"

}

#---------Main-------

if (( $# < 2 )); then
    yell "usage -- $(basename $0)  <file1> <file2>"
    die  "you had $# : $@"
fi


do_compare $@
