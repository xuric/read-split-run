#!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsw_config.sh"

OPTSFILE="" #to be filled-in
OUTPUTFILE="" #also to be filled-in

# options file needs:
# RNA-seq file: $2
# max split distance: $6
# rna sample length: $4
# known gene reference: $3
# reference boundries: $3.intronBoundries.exonsgaps
# supporting read tolerance: $7
# output basename: $(basename $2)
function make_options_file() {
    if [ -f "$OPTSFILE" ]; then
        rm "$OPTSFILE"
    fi
    if [ -z $REFDIR ];then
        die "Panic! dunno where the ref files are!" 
    fi
    if [ ! -f "$REFDIR/$1/refFlat.txt" ]; then
        die "Error: cannot find knownGene file for $1" 
    fi
    if [ ! -f "$REFDIR/$1/refFlat.txt.intronBoundry.exonsgaps" ]; then
        die "Error: cannot find intron/exon boundry file for $1" 
    fi
	touch "$OPTSFILE"
    echo "$2"  > $OPTSFILE  #reads file
    echo "$6" >> $OPTSFILE  #maxSplitDistance
    echo "$3" >> $OPTSFILE  #sampleLength
    echo "$REFDIR/$1/refFlat.txt" >> $OPTSFILE  #refFlat
    echo "$REFDIR/$1/refFlat.txt.intronBoundry.exonsgaps" >> $OPTSFILE  #intron/exon boundry
    echo "$5" >> $OPTSFILE  #minSplitDistance
    echo "$7" >> $OPTSFILE  #Support tolerance
    if [[ "$OUTPUTFILE" == "" ]]; then
        OUTPUTFILE="$9/$(echo "$(basename $2)" | cut -d. -f1)"
    fi
    echo "$OUTPUTFILE" >> $OPTSFILE  #results base name
    echo "$8" >> $OPTSFILE   #required supports
}

function dry_run() {
    if [ -f "$OPTSFILE" ]; then
        cat "$OPTSFILE"
    fi
}


function run_rsw() {
    if [ -f "$OPTSFILE" ]; then
        dprint "OUTPUTFILE = $OUTPUTFILE"
        $RSR "$OPTSFILE" >& "${OUTPUTFILE}.$RSW_TIMING_FILE"
        if [ ! -f "${OUTPUTFILE}.results" ]; then
            dprint "Panic! rsw failed to generate output file. Check stderr." 
            exit 1
        fi
    fi
}

function cleanup() {
#      OUTPUTFILE=$(tail -2 "$1" | head -1 )
 #     mv ${OUTPUTFILE}.* "$2"
      mv "$OPTSFILE" "$1"
}

if (( $# < 9 )); then
    yell "usage: $0 genome readsFile readLength minSplitSize minSplitdistance maxSplitdistance regionBuffer requiredSupports pathToSaveFesults" 
    die "you had $#" 
    exit 1
fi

date=$(date +%s)

genome=$1
dprint "fyi: \$2 = $2"
OPTSFILE="$TEMP_RSW_FILES/$(basename $2).${date}.options.txt"
OUTPUTFILE="$9/$(echo "$(basename $2)" | cut -d. -f1)"

make_options_file $*

if [ ! -d "${9}" ];then
    mkdir "${9}"
fi


run_rsw
cleanup "${9}"
#dry_run "$OPTSFILE"
echo "${OUTPUTFILE}.results"
