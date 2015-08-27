#!/bin/bash

#rsw_pipeline.sh
#main process for testing data acquired from the website
#for usage, see below

#metadefine
#BASEDIR=$(dirname $(realpath $0))
BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )
#bring in the configuration, which should be in the script's directory
source "${BASEDIR}/rsr_config.sh"

function usage() {
    yell "usage: $0 genome reads maxGoodAlignments minSplitSize minSplitdistance maxSplitdistance regionBuffer requiredSupports pathToSaveFesults" 
    yell "you had $#" 
    die "params: $@" 
}

# align the given reads with the $ALIGNER (bowtie)
# inputs = genome phase readsfiles maxGoodAlignments WhereToPutTheMetadata
#      genome = which bowtie index file to use, see bowtie.sh for list of options
#      phase  = which bowtie alignment to do: phase1 = lax, keep unmapped, phase2 = strict, keep matches
#             Choices are: "phase1" = lax... allows some mismatches; keeps unmapped reads
#                          "phase2" = strict, allows no mismatches ; keeps only matches
#      readsfiles = the files with read-data. see README for information
#      maxGoodAlignments = how many valid alignments are too many? (bowtie parameter -m)
#      whereToPutTheMetadata = the final output directory, bowtie's stdout will be put in a file here
function align() {
    genome=$1;
    phase=$2;
    reads=$3;
    maxGood=$4;
    dest=$5
    log "ALIGN: $genome $phase $reads $maxGood"
    result=$(try $ALIGN_SCRIPT $genome $phase $reads $maxGood)
    if (( $? )); then die "Failed to Align. Aborting."; fi
    #if we made a file for the printed output, move it
    if [ -f "${result}.bowtie.out" ]; then
        cp "${result}".bowtie.out $dest
    fi
    echo $result
}


# split_pairs makes the ... split pairs that will be re-aligned
# inputs = readsFile minSplit maxSplit
#        readsFile = what to split
#        minSplit  = how small do you want the smallest read-segement to be
#        maxSplit  = how big do you want the bigger reads-segement to be
#                    usually this is the whole length of the reads, but
#                     it can be smaller than the reads length to trim off the end.
function split_pairs() {
    #file=$1
    #start=$2
    #len=$3
    log "running $SPLIT_PROGRAM $@ $SPLIT_TEMP_DIR"
    echo $(try $SPLIT_SCRIPT "$@" $SPLIT_TEMP_DIR)
    if (( $? )); then die "Failed to split pairs.Aborting."; fi
}

#Temporarily out of service
#input: bowtie_output_base_name(s)
#function split_columns() {
#    results=()
#    for f in ${@}; do
#        if [ ! -f "${f}.bowtie.txt" ]; then
#            die "No bowtie file for ${f}. Cannot continue."
#        fi
#        $FORMATTER "${f}.bowtie.txt" 
#        if [ ! -f "${f}.bowtie.txt.split1stcolumn" ]; then
#            die "Failed to generate formatted data for ${f}"
#        fi
#        results+=( "${f}.bowtie.txt.split1stcolumn" )
#    done
#    echo "${results[@]}"   
#}

# split_columns
# takes the bowtie output and formats it for use in splitPairs
# inputs = bowtiedFileBaseName  
#         bowtiedFileBaseName = the name of the file bowtie produced, without the .bowtie.txt at the end
# the reason for not passing around the extensions is because it creates a streamline way to keep the original
# filename. this reduces the creep of extensions.
function split_columns() {
    file=$1
    if [ ! -f "${file}.bowtie.txt" ]; then
        die "No bowtie file for ${file}. Cannot continue."
    fi
    try $FORMAT_PROGRAM "${file}.bowtie.txt" >& "${LOG_FILE}"
    if (( $? )); then die "Failed to split columns. Aborting"; fi
    if [ ! -f "${file}.bowtie.txt.split1stcolumn" ]; then
        die "Failed to generate formatted data for ${file}"
    fi
    results="${file}.bowtie.txt.split1stcolumn"
    echo $results
}

# RSR
# this is the thing. it does the heavy lifting
#input: genome readsfile readlength minsplitsize minsplitdist maxsplitdist regionbuffer requiredSuppoerts pathtosaveresults
function rsr() {
    log "rsr.sh $*"
    echo $(try $RSR_SCRIPT "$@")
    if (( $? )); then die "Failed to do candidate matching. Aborting"; fi
}


##--------Main-----------
if (( $# < 9 )); then
    usage
fi

#step 0: Gather metadata
#mode was irrelevant to us, so we removed it
genome=$1
readlength=0
reads=$2
maxGood=$3
destination=$9

log "Preparing read info" 
#pop the genome, reads file, bowtie param (maxgood) off
#rest of params are for rsw
shift
shift
shift

log "Measuring reads... " 
readlength=$($MEASURE_SCRIPT $(echo "$reads" | cut -d, -f1))
log "length of reads $readlength" 

#step 1: Align original reads
log "Aligning reads... " 
results=$(align $genome phase1 $reads $maxGood $destination )
log "bowtie results basename(s)=${results}" 

#step 2: Split the unmapped reads into pieces
log "splitting into pairs..." 
#Note: going to assume that all reads are of the same length.
#TODO: maybe change this to individual splits for individual read lengths?
results=$(split_pairs "${results}" $1 $readlength )
log "split results=${results}" 

#step 3: re-align the split reads
log "Re-aligning reads... " 
results=$(align $genome phase2 ${results} $maxGood $destination)
log "re-align results=${results}" 

#TODO: Possibly integrate split bowtie column here.
# if (( $(du -m ${results}) > 50000 )); then try do_chrom_split ${results}; fi
#step 4: split the column into proper format
log "formatting..." 
results=$(split_columns ${results})
log "done formatting; result=$results" 


#step 5: select candidates
log "running rsr..."
result=$($RSR_SCRIPT $genome "$results" $readlength $@)
echo $result

