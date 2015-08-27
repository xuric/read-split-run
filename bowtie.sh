#!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsr_config.sh"

function set_indexes() {
    if [ -z "$BOWTIE_INDEXES" ] && [ ! -z "$BOWTIE_INDEX_ROOT" ]; then
        BOWTIE_INDEXES="${BOWTIE_INDEX_ROOT}/$1"
    elif [ ! -z "$BOWTIE_INDEXES" ] && [ -z "$BOWTIE_INDEX_ROOT" ]; then
        : # this is okay. NOOP
    else
        die "Don't know where to look for the bowtie indexes"
    fi
    export -p BOWTIE_INDEXES
}

function run_bowtie() {
    phase=$1
    bowtie_index_file=$2
    params="$3"
    inputs="$4"
    outfiles=()
    log "Using BOWTIE_INDEXES in $(printenv BOWTIE_INDEXES)";
    if [ -z "$(printenv BOWTIE_INDEXES)" ]; then
        die "No BOWTIE_INDEXES"
    fi
    outname="$BOWTIE_TEMP_DIR/$5"
    if [ "$phase" == "phase1" ]; then
        log "bowtie $bowtie_index_file $QUALS  $params -q $inputs --un ${outname}.unmapped.txt ${outname}.bowtie.txt"
        try $BOWTIE_PROGRAM $bowtie_index_file $QUALS $params -q $inputs --un "$outname".unmapped.txt "$outname".bowtie.txt >& "$LOG_FILE"
    elif [ "$phase" == "phase2" ]; then
        log "bowtie $bowtie_index_file $QUALS  $params -q $inputs $outname.bowtie.txt"
        try $BOWTIE_PROGRAM $bowtie_index_file $QUALS $params -q $inputs "$outname".bowtie.txt >& "$LOG_FILE"
    else
        die "Don't know what to do on phase $phase"
    fi

    #yes, it literally passes $4 back out, with some adornment
    echo $outname
}

#------MAIN-------

if (( $# < 4 )); then
yell "usage -- $0 genome phase fileSet maxGoodAlignments"
yell "      fileSet should be a comma-separated list of fastQ files"
yell "              Optionally, a second set of comma-separated files"
yell "              May be added, separated by a pipe | for paired-end"
exit 1
fi


#if [ -z "$BOWTIE_INDEXES" ]; then
#    case $1 in
#        "arabidopsis") result="TAIR10_chrAll";;
#        "mouse") result="mm9sp35";;
#    #    "mouse 79bp") result="mm9sp79";
#        "human") result="hg19sp101";;
#    #    "bacteria") result="GCA_000007425.1_ASM742v1_genomic";
#        *) exit 1
#    esac
#    BOWTIE_INDEXES="$BOWTIE_INDEX_ROOT/$1"
#fi

#try set_indexes $1
genome=$1
if (( $? )) && [ -z "$genome" ]; then
    die "Cannot find genome for $1"
fi

#parameter discovery

if [[ $3 =~ .*\|.* ]]; then # Paired mode
    if $DEBUG; then echo "Paired mode" 1>&2; fi
    OIFS=$IFS
    IFS='|' read pair1 pair2 <<< "$3"
    IFS=$OIFS
    input_params="-1 $pair1 -2 $pair2"
    outname=$(basename "$(echo $pair1 | cut -d, -f1)")
    QUALS=$(awk 'NR % 4 == 0' $(echo "$pair1" | cut -d, -f1) | head -$(( $QUALITY_TESTS * 4)) | python $ENCODING_GUESSER -b -n $QUALITY_TESTS)

else
    input_params=$3
    outname=$(basename "$(echo $3 | cut -d, -f1)")
    QUALS=$(awk 'NR % 4 == 0' $(echo "$3" | cut -d, -f1) | head -$(( $QUALITY_TESTS * 4)) | python $ENCODING_GUESSER -b -n $QUALITY_TESTS)
fi

#number of threads is now 3/4 of total cores on the system.
#watch this value...
bowtie_params="-t -p $(( $(grep -c ^processor /proc/cpuinfo) * 3 / 4))  -k $4 -m $(($4 - 1)) -v "
if [ "$2" == "phase1" ]; then
    bowtie_params+="2"
elif [ "$2" == "phase2" ]; then
    bowtie_params+="0"
fi
log "BOWTIE_INDEX_ROOT: $BOWTIE_INDEX_DIR, BOWTIE_IDNEXES: $BOWTIE_INDEXES"
log "bowtying $2 $genome $bowtie_params $input_params"

echo $(run_bowtie $2 "$genome" "$bowtie_params" "$input_params" "$outname")

