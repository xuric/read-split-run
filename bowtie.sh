cleanup.sh                                                                                          0000775 0001753 0001753 00000000555 12563403132 012662  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsr_config.sh"

FIND_MODIFIERS="-L"
FIND_WHAT="-daystart -ignore_readdir_race -mtime +14 "
ACTION="-delete"

find $FIND_MODIFIERS $TMPDIR $FIND_WHAT
find $FIND_MODIFIERS $TEMP_BOWTIE_FILES $FIND_WHAT
find $FIND_MODIFIERS $TEMP_SPLIT_FILES $FIND_WHAT
find $FIND_MODIFIERS $TEMP_RSW_FILES $FIND_WHAT

                                                                                                                                                   guess-encoding.py                                                                                   0000664 0001753 0001753 00000005510 12563403215 014156  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                """
guess-encoding.py: (c) Brent Pedersen

useage: awk 'NR % 4 == 0' your.fastq | python %prog [options]

guess the encoding of a stream of qual lines.
"""
import sys
import optparse

RANGES = {
    'Sanger, Phred+33': (33, 93),
    'Solexa, Solexa+64': (59, 104),
    'Illumina-1.3, Phred+64': (64, 104),
    'Illumina-1.5, Phred+64': (67, 104)
}

BOWTIE_PARAMS = {
    "Sanger, Phred+33": '--phred33-quals',
    "Solexa, Solexa+64": '--solexa-quals',
    "Illumina-1.3, Phred+64": '--phred64-quals',
    "Illumina-1.5, Phred+54": '--phred64-quals'
}

def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89)
    """

    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)

def get_bowtie_option(encoding):
    return BOWTIE_PARAMS[encoding]

def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-n", dest="n", help="number of qual lines to test default:-1"
                 " means test until end of file or until it it possible to "
                 " determine a single file-type",
                 type='int', default=-1)
    p.add_option("-b", dest="b", help="make output the option to use in bowtie"
                 " corresponding to the quality type identified",
                 action="store_true", default=False)

    opts, args = p.parse_args()
    ##print >>sys.stderr, "# reading qualities from stdin"
    gmin, gmax  = 99, 0
    valid = []
    for i, line in enumerate(sys.stdin):
        lmin, lmax = get_qual_range(line.rstrip())
        if lmin < gmin or lmax > gmax:
            gmin, gmax = min(lmin, gmin), max(lmax, gmax)
            valid = get_encodings_in_range(gmin, gmax)
            if len(valid) == 0:
                ##print >>sys.stderr, "no encodings for range: %s" % str((gmin, gmax))
                sys.exit()
            if len(valid) == 1 and opts.n == -1:
                if (opts.b):
                    print get_bowtie_option(valid[0])
                else: 
                    print "\t".join(valid) + "\t" + str((gmin, gmax))
                sys.exit()

        if opts.n > 0 and i > opts.n:
            if (opts.b):
                print get_bowtie_option(valid[0])
            else: 
                print "\t".join(valid) + "\t" + str((gmin, gmax))
            sys.exit()

    if (opts.b):
        print get_bowtie_option(valid[0])
    else: 
        print "\t".join(valid) + "\t" + str((gmin, gmax))


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
                                                                                                                                                                                        gzip_pipe.sh                                                                                        0000775 0001753 0001753 00000001064 12563403132 013215  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash
DEBUG=true

yell() { echo "$(basename $0): $*" >&2; }
die() { yell "$*"; exit 111; kill $$; }
try() { "$@" || die "cannot $*"; }
dprint() { if $DEBUG; then yell "$*"; fi }




#----Main----

if (( $# < 1 )); then
    die "Usage -- $(basename $0) <gzipped file>"
fi

inputfile=$1
outputfile="./$(basename ${inputfile}).pipe"

if [ ! -f "$inputfile" ]; then
    die "No such file $1"
fi

mkfifo $outputfile
if (( $? )) || [ ! -p "$outputfile" ]; then
    die "could not create pipe. aborting"
fi

gunzip -c $inputfile > $outputfile &
echo "$outputfile"

                                                                                                                                                                                                                                                                                                                                                                                                                                                                            INSTALLATION                                                                                        0000664 0001753 0001753 00000007245 12563404513 012527  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                Installing the RSR pipeline

Unzip the contents of the archive (tar -xf rsr.gz).
Change directory to the output directory (cd rsr)
Run 'make' (make); or, compile each .cpp file individually as described below, if you change the output name of the program, you will need to update the rsw_config.sh file to use the name you supplied...  :
    splitPairs.cpp: g++ splitPairs.cpp -std=c++11 -O4 -o sp4
    split_read_rsw.c:  gcc split_read_rsw.c -O4 -o srr
    split_first_column.c: gcc split_first_column.c -O4 -o sfc
    split_bowtie_chrom.c: gcc split_bowtie_chrom.c -O4 -o sbc
    compare.cpp:    g++ compare.cpp -O4 -o rsw_compare

Note: "split bowtie chrom" is sort of a bonus file, see the README for details.

create directories in the location of your choice for holding temporary and permanent output files for the pipeline; you will need these for rsw_config.sh


CONFIGURATION

edit "rsw_config.sh" with the editor of your choice (nano, vim, emacs, etc.).
Navigate to the section labelled "USER CONFIGURATION" where there are two sections to update:
The following fields may need updating, depending on your system. Their default values can work, but 
the TEMP_DIR entries will need to be created.  The entries follow this format:
VARIABLE="VALUE"   . The double-quotes are required, so please keep them.

set the TEMP_DIR values to locations on your system to where the pipeline can store intermediate files. 
It is recommended that these locations have significant allocations of disk space (terabytes, preferrably).

An explanation of each variable follows:

BOWTIE_TEMP_DIR
This is where the unmapped and mapped output files from bowtie will be stored. 

SPLIT_TEMP_DIR
This is where the split unmapped reads will be stored. This directory should have a lot of free disk space,
as these are the largest of the files produced in the process.

RSR_TEMP_DIR
RSR's temporary configuration file and the output files will be stored here before being moved to your 
specified target directory.

LOG_DIR
The RSR pipeline can generate a detailed log of its operation. Specify here where those logs should be
kept.

BASES_TO_TRIM
If your reads have low quality at the right-end (3p), specify here how many to trim off the edge.
This is done at the SPLIT step.

BOWTIE_INDEXES
Where your bowtie index files (ebwt) are located. There is a level of customization here:
If ALL your index files are in ONE directory, set this value to that directory.
if Your index files are in subdirectories of a particular base, set this to an empty string ("").

BOWTIE_INDEX_DIR
If your index files are in subdirectories of a particular base, set THIS variable to that directory.
Otherwise, leave it blank ("").

AVAILABLE_INDEXES
This should be a double-quoted NOT COMMA SEPARATED list of the index prefixes for your available genomes.
These values will be passed to bowtie to specify which index to use. (see bowtie manual:
http://bowtie-bio.sourceforge.net/manual). 

INDEX_NAMES
Like AVAILABLE_INDEXES, this should be a space-delimited, double-quoted list...
This variable does double-duty within the pipeline. It contains the names of the genomes you have available, in
more friendly terms, and IT IS THE LIST OF DIRECTORY NAMES USED TO FIND YOUR INDEXES.

For example: all my index files are stored in directories named after the organism to which they referr.
so I have: "human" "mouse" "arabidopsis" (becasue I can't find a better way to say "mustard weed").

(see also: Index Examples, below)

REFDIR
Where your flat reference files are stored. these are the knownGene.txt files associated with the organisms
for which you have index files.


#------Programs----------------------
BOWTIE_PROGRAM
SPLIT_PROGRAM
FORMAT_PROGRAM
RSR_PROGRAM
COMPARE_PROGRAM

                                                                                                                                                                                                                                                                                                                                                           makefile                                                                                            0000664 0001753 0001753 00000000355 12563406327 012403  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                all: sp sfc srr sbc

sp: 
	g++ -O4 -o sp4 src/splitPairs.cpp -std=c++11

sfc: 
	gcc -O4 -o sfc src/split_columns.c

srr: 
	gcc -O4 -o srr src/split_read_rsw.c

sbc: 
	gcc -O4 -o sbc src/split_on_chrom.c


clean:
	rm -f sp4 sfc srr sbc


                                                                                                                                                                                                                                                                                   readlength.sh                                                                                       0000775 0001753 0001753 00000000210 12563403132 013334  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsr_config.sh"

echo $(( $(head -2 $1 | tail -1 | wc -c) - $OFFSET ))
                                                                                                                                                                                                                                                                                                                                                                                        README                                                                                              0000664 0001753 0001753 00000021156 12564422417 011564  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                README

This is the software package for the Read-Split-Run pipeline. Included are the scripts and software necessary to run the entire process from beginning to end.


INSTALLATION

Installing the pipeline:
    gunzip rsr.gz to the location of your choice; this shall be the INSTALLATION DIRECTORY.

Bowtie 1.0.1 (or later): 
    Ensure that you have bowtie (version 1.0.1 or later) installed in an accessible location (see CONFIGURATION, and DEPENDENCIES, below). It can be downloaded from http://bowtie-bio.sourceforge.net.

Bowtie Indexes, gene reference files:
You will need bowtie indexes and knownGene files for the genome(s) of your choice. 
The index files can be stored in two ways to be compatable with RSR:
    All in one common directory, or each in its own subdirectory of some other common directory.
    ex: (one directory):   /usr/local/bowtie/indexes
        (subdirectories):  /net/projects/genomes/human
                           /net/projects/genomes/mouse

For knownGene reference files:
Create a directory under which you will store your knownGene reference files. In this directory, make subdirectories for your organisms, either by assembly name or common name (of your choice). This is similar to the "common directory" method described for index files, above.
The pipeline requires that the flat reference file be parsed to mark intron/exon boundaries. We have provided a perl script "refflat_parse_RSW.pl" to accomplish this. Once your knownGene files are in place, for each of them run the provided "refflat_parse_RSW.pl" and provide as argument your knownGene file (with full path).
Ex:  perl refflat_parse_RSW.pl /net/projects/refs/human/knownGene.txt



CONFIGURATION:

Edit the configuration file: "rsr_config.sh" in a plain-text editor (vim, nano, etc.)...
change the following variables in the USER CONFIGURATION section to suit your situation:

BOWTIE_TEMP_DIR     location to store intermediate bowtie files.

SPLIT_TEMP_DIR      location to store intermediate split reads files.

RSR_TEMP_DIR        location to store intermediate RSR options and output files.

LOG_DIR             location to store diagnostic and operational logs.

BOWTIE_INDEX_DIR    If your bowtie indexes are separated by directory names under a common
                    directory, specify that place here and leave BOWTIE_INDEXES empty

BOWTIE_INDEXES      If your bowtie indexes are all contained in a single directory, specify
                    that here and leave BOWTIE_INDEX_DIR empty. 

BOWTIE              path to bowtie executable (e.g. /usr/bin/bowtie).

REFDIR              A common directory under which you store knownGene reference files in
                    subdirectories named after their orgnism's name, or common name. 

BASES_TO_TRIM       How many bases should be trimmed off the right-end of your reads? 
                    Set this according ot the quality of your reads. Default: 0
                    This setting is enforced in bowtie, and elsewhere.

AVAILABLE_INDEXES   Here, put the common-names for the organisms of your bowtie indexes and 
                    knownGene files. Inside the parentheses, enter the names with double-quotes
                    round them ( "human" "mouse" ...); do not put commas between them.

INDEX_NAMES         The bowtie indexes available on your system, by name only (the part before
                    the .ebwt portion of the index file). Put the index names in quotes between
                    the parentheses IN THE SAME ORDER AS THE COMMON NAMES.
                    Ex: if your AVAILABLE_INDEXES=( "human "mouse" "arabidopsis" )
                    INDEX_NAMES=( "hg19sp101" "mm9sp35" "TAIR10_ChrAll" )
                    See Also: NAMING CONVENTIONS, below.


All other variables are internal-use and should not be changed. Read the comments in the configuration file
for details as to what function they provide. 


DEPENDENCIES:

All shell scripts (files ending in .sh) rely on "rsw_config.sh," which holds the configuration data.

rsw_batch_job.sh
-   rsw_pipeline.sh
    -   bowtie.sh
        -   bowtie v1.0.1 or newer
    -   split.sh
        -   srr
    -   sfc 
    -   rsw.sh
        -   rsw
-   compare_sh
    -   rsw_comparison

refFlat_parse_RSW.pl and sbc have no dependencies.

RUNNING:

To run the pipeline, execute rsw_batch_job.sh with the following inputs:

mode                Choose "analytic" or "comparison"
                    analytic jobs only produce RSW output files
                    comparative jobs produce RSW output files AND a file which shows the differences
                    between the two

genome              Which genome to which to align reads. See NAMING CONVENTIONS for more info...

readsFile           The file(s) with RNA-Seq data in plain-text FASTQ format.
                    The nature of your run will determine how you should specify your files.
                    if you have ... you should specify your file(s) as...
                        Just one file       the file name, with full path, by itself.
                        Replicates          a double-quoted, comma-separated list
                                            ex: "replicate1.fastq, replicate2.fastq, ..."
                                            ** Yes, include the double-quotes.
                        Paired-end data     place all left-pairs in a comma-separated list, as above, 
                                            then put a vertical-bar | and a second comma-separated 
                                            list with the right-pairs.
                                            ex: "Replicate1_1.fastq,Replicate2_1.fastq|Replicate1_2.fast1,Replicate2_2.fastq"
                                            ** Make sure your pairs are ordered correctly:
                                            "REPLICATE 1 LEFT, REPLICATE 2 LEFT | REPLICATE 1 RIGHT, REPLICATE 2 RIGHT"
                     
           
[readsFile2]        If you are using "comparison" mode, a second set of reads-files goes here,
                    the format is the same as above.

maxGoodAlignments   Maximum number of matches allowed in bowtie (see bowtie -k and -m parameters)

minSplitSize        Smallest length to split your reads into. If you specify more than half the reads' length,
                    the pipeline will exchange it with (readlength - minSplitSize).
                    ** The smaller your split, the more memory, disk space, and time will be needed.

minSplitdistance    Smallest amount of distance allowed between split-reads to be considered a splice-junction candidate.

maxSplitdistance    Largest amount of distance between split-reads to be considered a splice-junction candidate.

regionBuffer        By how much candidate junctions are allowed to vary in their start-position in order to "support" one-another.

requiredSupports    Report only junctions with this many supporting reads, or more.

pathToSaveFesults   Where to put the final output files.


Examples:

Analytic runs:

Normal Run: ./rsw_batch_job.sh analytic mouse "mm9.fastq" 11 11 3 30000 4 2 ~/mm9_results

Run with Replicates:  ./rsw_batch_job analytic human "hg19-1.fastq,hg19-2.fastq,hg19-3.fastq" 11 33 3 50000 5 2 ~/hg19_results

Paired Run: ./rsw_batch_job.sh analytic human "hg19_1.fastq|hg19_2.fastq" 11 33 3 50000 5 2 ~/hg19_paired_results

Paired, replicated: ./rsw_batch_job analytic human "hg19-1_1.fastq,hg19-2_1.fastq|hg19-1_2.fastq,hg19-2_2.fastq" 11 33 3 50000 5 2 ~/hg19_pair_repl_results

Comparative runs:

Normal Run: ./rsw_batch_job.sh comparative mouse "set1.fastq" "set2.fastq" 11 11 3 30000 4 2 ~/mm9_compare
Run with Replicates:  ./rsw_batch_job analytic human "set1replicate1.fastq,set1replicate2.fastq" "set2replicate1.fastq,set2replicate2.fastq" 11 33 3 50000 5 2 ~/hg19_results
Paired Run: ./rsw_batch_job.sh analytic human "set1_1.fastq|set1_2.fastq" "set2_1.fastq|set2_2.fastq" 11 33 3 50000 5 2 ~/hg19_paired_results
...

NAMING CONVENTIONS

Most genome index files are given by their NCBI Assembly names, which we found to be quite inconveninet at times. Therefore we implemented an optional "common name" method of addressing one's genomes (i.e. hg19sp101 = "human", mm9sp35 = "mouse", etc.). These common names are used for the "genome" input parameter to rsw_batch_job.sh, and as the directory name, under your REFDIR, where the knownGene and gene intron/exon boundary files will be located (see INSTALLING).
If you prefer to use the assembly names, simply use them as the directories for your knowngene reference (and possibly for your index subdirectory, if used).

KNOWN ISSUES

The Quality-encoding dection portion of bowtie.sh is known to cause a broken pipe with awk. This is acceptible and does not interfere with the performance of the pipeline.


                                                                                                                                                                                                                                                                                                                                                                                                                  refflat_parse_RSW.pl                                                                                0000755 0001753 0001753 00000002373 12564410473 014611  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/usr/bin/perl

### perl script to parse exons from UCSC refFlat.txt file


if (!defined @ARGV) {
die  "Correct Syntax is: refflat_parse.perl <fileneame> \n\nPlease supply a file to parse.  The output file will be the input file name with \".intronBoundary\".exonsgaps\n\n";
}

$infile = shift @ARGV;

$outfile = $infile . ".intronBoundary.exonsgaps";


open INFILE, "<$infile";
open OUTFILE, ">$outfile";


while (<INFILE>) {
	chomp;
	if ($_ =~ /^\#/) {
		next;
	}
	
	($genename, $name, $chrom, $strand, $txStart, $txEnd, $cdsStart, 
	$cdsEnd, $exonCount, $exonStarts, $exonEnds)=split(/\t/,$_,11);
	
	@exstarts=split(/,/,$exonStarts);
	@exends=split(/,/,$exonEnds);
	
	$count=0;
	while ($count <$exonCount) {
		    $ex_num=$count+1;
			
			if ($ex_num == 1) {
					$introns[$count]="NA";
					$intronsStart[$count]="NA";
					$intronsEnd[$count]="NA";
				}
			if ($ex_num >1) {
					$introns[$count]=$exstarts[$count] - $exends[$count-1];
					$intronsStart[$count]=$exends[$count-1] + 1;
					$intronsEnd[$count]=$exstarts[$count] - 1;
				}
		
			print OUTFILE "$genename\t$name\t$chrom\t$strand\t$txStart\t$txEnd\t$cdsStart\t$cdsEnd\t$ex_num\t$exstarts[$count]\t$exends[$count]\t$introns[$count]\t$intronsStart[$count]--$intronsEnd[$count]\n";
					
			$count++;
	}
}
                                                                                                                                                                                                                                                                     rsr_config.sh                                                                                       0000664 0001753 0001753 00000010562 12563424713 013372  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

#------Only Configure the thing ONCE, please---
if [ -z "$CONFIGURED" ]; then
CONFIGURED=true
BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )  #moved to this method to avoid platform-specific comamnds

#-------USER CONFIGURATION------------
#-------Directories-------------------
BOWTIE_TEMP_DIR="/scr/RSW/bowtie"
SPLIT_TEMP_DIR="/scr/RSW/bowtie"
RSR_TEMP_DIR="/scr/RSW/tmp"
LOG_DIR="${BASEDIR}/logs"
BASES_TO_TRIM=2
BOWTIE_INDEXES=""                            #See INSTALLATION
BOWTIE_INDEX_DIR="/net/project/common/data"  #See INSTALLATION
REFDIR="/net/home/bdonham/public_html/dev/public/database/refs"    #where the refFlats are

AVAILABLE_INDEXES=( "arabidopsis" "mouse" "human" )
INDEX_NAMES=( "TAIR10_chrAll" "mm9sp35" "hg19sp101" )

#------Programs----------------------
BOWTIE_PROGRAM="/usr/local/bin/bowtie" # /path/to/bowtie
SPLIT_PROGRAM="${BASEDIR}/srr"                    # Program for splitting reads, compiled from split_read_rsw.c
FORMAT_PROGRAM="${BASEDIR}/sfc"                   # Program for formatting reads, compiled from split_first_column.c
RSR_PROGRAM="${BASEDIR}/sp4"                      # RSR Program ("split pairs"), compiled from splitPairs.cpp
COMPARE_PROGRAM="${BASEDIR}/rsw_compare"          # Program for comparing RSR outputs, compiled from 

#------NOT-USER CONFIGURATION--------
#CHANGE THESE AT YOUR OWN RISK

#Meta defines
DEBUG=false
if [ -z "$LOG_FILE" ] || [ -z "$RUN_ID" ]; then
export -p RUN_ID="$(date +%F.%s)"
export -p LOG_FILE="RSR_${RUN_ID}.log"
fi

yell() { echo "$(basename $0): $*" >&2; }
die() { yell "$*"; log "$*"; exit 111; kill $$; }
try() { "$@" || die "cannot $*"; }
dprint() { if $DEBUG; then yell "$*"; fi }
log() { echo "$(basename $0): $*" >> "${LOG_DIR}/${LOG_FILE}"; dprint "$@"; }
timing_start() { if [ -z "$START_TIME" ]; then START_TIME=$(date +%s); fi }
timing_end() { if [ ! -z "$START_TIME" ]; then log "Duration: $(( $(date +%s) - $START_TIME )) seconds"; fi }



#this sets the location to find the sub-parts of the pipeline
#BASEDIR=$(dirname $(realpath $0))

#rsw_batch_job.sh:  The big guy!
PIPELINE="${BASEDIR}/rsw_pipeline.sh"
COMPARE_SCRIPT="${BASEDIR}/rsw_compare.sh"

#rsw_pipeline.sh: Pipeline Constants
#Define program/script files
ALIGN_SCRIPT="${BASEDIR}/bowtie.sh"
MEASURE_SCRIPT="${BASEDIR}/readlength.sh"
SPLIT_SCRIPT="${BASEDIR}/split.sh"
#FORMAT_SCRIPT="${BASEDIR}/split_1stColumn_RSW.pl"  #there is no format script, only srr
RSR_SCRIPT="${BASEDIR}/rsr.sh"   #_SCRIPT'd for consistency

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

#Manage directories
if [ ! -d "$BOWTIE_TEMP_DIR" ]; then
    mkdir "$BOWTIE_TEMP_DIR"
    if [ ! -d "$BOWTIE_TEMP_DIR" ]; then
        #couldn't make dir, panic
        die "Could not make BOWTIE_TEMP_DIR. aborting."
    fi
fi
if [ ! -d "$SPLIT_TEMP_DIR" ]; then
    mkdir "$SPLIT_TEMP_DIR"
    if [ ! -d "$SPLIT_TEMP_DIR" ]; then
        #couldn't make dir, panic
        die "Could not make SPLIT_TEMP_DIR. aborting."
    fi
fi
if [ ! -d "$RSR_TEMP_DIR" ]; then
    mkdir "$RSR_TEMP_DIR"
    if [ ! -d "$RSR_TEMP_DIR" ]; then
        #couldn't make dir, panic
        die "Could not make RSR_TEMP_DIR. aborting."
    fi
fi
if [ ! -d "$LOG_DIR" ]; then
    mkdir "$LOG_DIR"
    if [ ! -d "$LOG_DIR" ]; then
        #couldn't make dir, panic
        die "Could not make LOG_DIR. aborting."
    fi
fi


#end if..configured
fi
                                                                                                                                              rsr.sh                                                                                              0000775 0001753 0001753 00000005107 12563403132 012037  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsr_config.sh"

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
    if [ -z "$OUTPUTFILE" ]; then
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
        log "OUTPUTFILE = $OUTPUTFILE"
        if [ -z "$LOG_FILE" ]; then
            logfile="/dev/null"
        else
            logfile="${LOG_DIR}/${LOG_FILE}"
        fi
        try $RSR_PROGRAM "$OPTSFILE" >& $logfile
        if [ ! -f "${OUTPUTFILE}.results" ]; then
            log "Panic! rsw failed to generate output file. Check stderr." 
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
log "fyi: \$2 = $2"
OPTSFILE="$RSR_TEMP_DIR/$(basename $2).${date}.options.txt"
OUTPUTFILE="$9/$(echo "$(basename $2)" | cut -d. -f1)"

make_options_file $*

if [ ! -d "${9}" ];then
    log "${9} didn't exist. something wrong?"
    mkdir "${9}"
fi


run_rsw
cleanup "${9}"
#dry_run "$OPTSFILE"
echo "${OUTPUTFILE}.results"
                                                                                                                                                                                                                                                                                                                                                                                                                                                         rsw_batch_job.sh                                                                                    0000775 0001753 0001753 00000004211 12563424512 014037  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash
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
    yell "usage: $(basename $0) mode genome readsFile [readsFile2] maxGoodAlignments minSplitSize minSplitdistance maxSplitdistance regionBuffer requiredSupports pathToSaveFesults"
    die "params($#): $@"
fi

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

                                                                                                                                                                                                                                                                                                                                                                                       rsw_compare                                                                                         0000775 0001753 0001753 00000057237 12563403140 013153  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                ELF          >    �@     @       �B          @ 8 	 @         @       @ @     @ @     �      �                   8      8@     8@                                          @       @     [5      [5                    �=      �=`     �=`     D      HR                    �=      �=`     �=`                                T      T@     T@     D       D              P�td   �0      �0@     �0@     �       �              Q�td                                                  R�td   �=      �=`     �=`     (      (             /lib64/ld-linux-x86-64.so.2          GNU                        GNU ��?"0���L��e�A[~*�   "              "       (E�LyIk�                                                                                           �                     i                      o                     o                      �                     W                     �                     �                                          '                     �                     3                       �                     �                     O                       �                     �                      �                     �                                          �                                           �                     �                      �                     �                                          !                     �                     �                     �      p@                  P@              libstdc++.so.6 __gmon_start__ _Jv_RegisterClasses _ITM_deregisterTMCloneTable _ITM_registerTMCloneTable _Znam _ZNSs4_Rep10_M_disposeERKSaIcE _ZSt29_Rb_tree_insert_and_rebalancebPSt18_Rb_tree_node_baseS0_RS_ _ZNSt8ios_base4InitD1Ev _ZNSs6appendEPKc _ZNSsC1EPKcRKSaIcE __gxx_personality_v0 _Znwm _ZSt18_Rb_tree_decrementPSt18_Rb_tree_node_base _ZNSt8ios_base4InitC1Ev _ZdlPv libm.so.6 libgcc_s.so.1 _Unwind_Resume libc.so.6 fopen puts strtok strtol fgetc strlen __cxa_atexit strstr fputc memcpy fclose malloc fwrite fprintf memmove strcmp __libc_start_main GCC_3.0 CXXABI_1.3 GLIBCXX_3.4 GLIBC_2.14 GLIBC_2.2.5                                                  �         P&y   ,              0   ӯk   4     t)�   ?        �         ���   K     ui	   V      �?`                   @`                    @`                   (@`                   0@`                   8@`                   @@`                   H@`                   P@`        	           X@`        
           `@`                   h@`                   p@`                   x@`                   �@`        "           �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`                   �@`        #           �@`                    A`                   A`                    A`        !           H��H��3  H��t�3   H���              �5r3  �%t3  @ �%r3  h    ������%j3  h   ������%b3  h   ������%Z3  h   �����%R3  h   �����%J3  h   �����%B3  h   �����%:3  h   �p����%23  h   �`����%*3  h	   �P����%"3  h
   �@����%3  h   �0����%3  h   � ����%
3  h   �����%3  h   � ����%�2  h   ������%�2  h   ������%�2  h   ������%�2  h   ������%�2  h   �����%�2  h   �����%�2  h   �����%�2  h   �����%�2  h   �p����%�2  h   �`����%�2  h   �P����%�2  h   �@����%�2  h   �0����%�2  h   � ����%�2  h   �����%�2  h   � ����%z2  h   �����UH��AWAVAUATSH��   ����p���H��P�����
  ��1��
   ����  H��P���H�x�����Hc�  �׀  H�      PH������H9�wHk�H������H�HH��H�H�S�H��t%H��H�     H�@    H�@    H��H���u�Hc��  I�      �H������H���  H�<�    L9�HG��I���H�r�  HcK�  H�<�    L9�HG��'���H�@�  Hc)�  H�E�    H��   H���H)�H�D$H�E���p�����H��   H��x���H�E�H��P���1��
   L�,@H�\ H�M��L%�  ����H�މ�L����  L��H�  H�XL� L9�H�]�H����   L)�?   L��H��H��H��H��?H�H)�H���  H���   ��   I��$�   L��H���~  H9]�u�    H��H9]�M�>�   L�;L�k�I��D  M�e I�L�m�I�t$�*�����x6u�I�4$I�?������x$u�I�D$I9G|�I�D$ I9G }�f.�     M�&I��L�u��H�u�L����
  �    H�E�H��x���H9E��������~  E1�H�=�~  L�5�/  ����   L�e�O�|m Ic�A��    I��J�?H�
H�BH)�H��H��teL��   1��	H��L��L��L95�/  H�4��c  M��H�tI�H�=D~  H�u/  J�?L�pL�5f/  H�
H�BL�CH)�H��H9�w�H� ~  H��}  B��    B��    A�EI��9�}  �>���H�/  L9��  M��L��H��I)�L��H��H��H��?Hcи?   H)�H� ��
  I���   ��   L���   H��L���	  M9���   L�e�L�u�L����H�E�M�4$H�E�H9E���   L�0H�X�I���    L�+I�~I��I�u�%�����x1u�I�u I�>������xu�I�EI9F| �I�E I9F }�@ M�,$H��M��뫿@A` �y  H�=�|  L�5.  J�?����L��H����  H��P���H�}�H�U�H�p�O���H�E�H�}��i.@ H��H�������H�}�H�E��~.@ H��H�������H�E�H��H�E�H��H�����  H��P���H�}�H��H�s�����H�}���.@ �P���H�}��~.@ �R���H��H��X�����  H�M��   �   ��.@ �:���H��X����   �   ��.@ ������p���H��P�����H��L�d��    H�H�}���.@ 1�����H�H��X�����.@ 1������H��L9�u�H�u��
   � ���H��X����
   �����H�=�,  H�5�,  H��H)�H��H����  Hǅh���   1��H�����` H�H��p���H��z  H�H��u�[ H��H�@H��tH9P s�H�@H��u�H����` t4H9Q w.H��h���H��H)�H��H9�H�S�-  H��H��h����@ H�JH�}���.@ H�1��	�����z  ���   �����1ۉ�d�����x����1fD  H�M��   �   ��.@ �����CH��;;z  �o  H�}�Hc�L�,[�]�L�$�A�$HcȉE�H�Bz  J��H�H�@H)�H��H9�s�H�T+  H��p���L�4�L�<8I�vI���������p���I�6I�?��������]���I�OI;N�O���M�F I�G L)�H�H1�H)�Hc�y  H9��,���A�~0 �g.@ �a.@ M�N(I�V��.@ HD�H�}�H�$1������H��y  Ic$L�=y  J��M��L�,���   I�M �@ I�G�   H��tI��I�W H9�w�I�G1�H��u�@����   H9���   �E�����x����A�$�E��/  ��d����CH��;�x  �����H�u��
   �X�����d����H�5
*  �
  ��x������   H��p���Hc�x���H��x0 ��   H��x  H�=�)  ������A���` fD  L�5)x  M9�t%L������I�M H�P H9��:���M��M���.���I����` ��   I�F I9E A�ǿ(   ����H��t&H�@    H�@    H�@    I�U �     H�P A�����` L��H�������H��w  �����     ��x����]���H��w  ��H�=�(  �����A�   �u���H�U�H�u���/@ 1��/�����.@ �E���H��X����3   �   �80@ �����p0@ � ����:w  1ۅ�~l@ H�9w  L��P�����.@ H��X�����H�,w  I�T� D��1��|���H�w  I�t� ��.@ ��H��v  ��1������CH��9�v  �H�}�����H��X�������H�E�H��H���H�x�H������H�E�H��H�x�����H�e�1�[A\A]A^A_]Ã����d�����x��������0/@ �=���1��v���H��H�E�H��H���H�x��O���H�������H�u���/@ 1�������   �=���H��H�E�H��H���H�x������H�u���/@ ������@ H����` �����8.@ ��` �p@ �N����8.@ ���` �p-@ �}u      H�zu      H��u      H�lu  ��` H�iu  ��` ����H��&      H��&      �8.@ H��&      �@A` ��*@ H������� 1�I��^H��H���PTI�� .@ H���-@ H�Ǡ@ �����f�f.�     f��'A` UH- A` H��H��w]ø    H��t�]� A` ���    � A` UH- A` H��H��H��H��?H�H��u]ú    H��t�]H�ƿ A` ���    �=�%   uUH���~���]��%  ��@ H�=�"   t�    H��tU��=` H����]�{��� �s��� ATUH��SH��H�H�v�Y�����xUA�    t[]D��A\�D  H�3H�} �4�����x0u�H�CH9EA�   |�A�    �H�C H9E A����     [A�   ]D��A\�f�AWAVAUATUSH��H9�H�<$H�t$�   L�oL9�u�  f�H�+I�]H9\$I����   H�$I�m L��L�8L�eI�wL����������   tcM�}�L��@ M�7M��I�v�i�����x5u�I�6H�} �W�����x#u�I�FH9E|�I�F H9E �w����    L�3I��H�}L���I�7H�} ������x u�I�GH9E|�{���I�G H9E �m���I�]L+,$I��M��tJ��    H�4$H��H)�����H9\$H�$I��H�(���� H��[]A\A]A^A_ÐAWH��H)�AVAUATUSH��H��xH=�   H�|$H�T$ H�t$��  H�|$  �;  H�D$H��H�D$(H�T$H�l$ H)�L�rH��H��H��?M�fH�H��H�,�L��H�] L�kL����������  ��  H�D$L��L�x�I�GH��H�$��������X  �  H�4$L�����������  ��  H�T$L�l$H�H�H�E H�H�D$H�\$(H�D$�I@ M�}�M�7L��M��I�v��������   ��   I9ߐ��   H�D$L�3H��I�/H� H�D$L�`�    H�+L��H�$H�}�5�����xAu�H�D$H�} H�0������x*u�H�D$H�HH9M|�k���H�@ H9E �]���D  H���f.�     H�D$I�6H�8�������x,�R���H�D$I�VH9P|�=���I�V H9P �/��� I�������    H�T$ H�t$H�������H��H+D$H=�   ��  H�|$  �;  H�\$� ���H�3I�>�I�����x$�7���H�CI9F|�'���H�C I9F ����H�D$L��L�x�I�GH��H�$��������   ��   H�4$L�����������   �v  H�L$L�t$L�l$H�L�1H�A����I�7I�>������x$�����I�GI9F|�����I�G I9F �����H�T$L�t$L�l$H�L�2H�B�����I�7H�;�c�����x$�b���I�GH9C|�R���I�G H9C �D���H�L$L�l$H�H�H�E H�H�D$�r���I�7H�;������x$�A���I�GH9C|�1���I�G H9C �#���H�T$H�L$H�L�:I��H�A�H�H�D$����M��1�H�D$8�uH��H��H��?H�H��H9���  H����  f�H�l$ H�|$(H�D$I���  H��x[]A\A]A^A_�I�7I�>�\������m����q���I�GI9F�Y����]���I�G I9F �E����J���H�4$H��H�H�H�L$@H��H�|$@��H��H�L$PH�L$@H�D$ H�D$H�\	H�|$ H�L$0I��H��H�\$hH��I��H�T$XH��H�\$`H�T$8H� H�D$HH�D$ I9�|�  �     H;\$ O����   I��I�EH� H��M�H�S�I�*L�T$M�4�H�T$M�&H�}I�t$�X�������   I��L�T$u�I�4$H�} L�T$(H�l$�/�������   L�D$L�T$(�w���I�D$H9E��   �b���I�D$ H9E HL\$ML�ML�H;\$ O���K���H�|$P uH9\$@uL�T$`H�L$XH�\$hI�H�H9\$0|?H�l$8H�|$0 H�D$HI���   H�l$0H�D$8L�l$0�����M��M��H�\$�����H�C�L�t$HI��I��?I�I��K�,�I�vL�m I�}�G�����x0tM���I�6I�} �/�����xu�I�FI9E|�I�F I9E }�I�D$�M�,�H��H��?H�H��L9d$0|I���=���L��I���H�$L�l$H��H�D$ @ H�L$ H�H�D$I�E H�H��L)�H�D$(H��H�D$8H��I��I��?L�H�D$H�|$H�|$ �����E1���    H;\$O�D� �����I��I�GH� H��M�L H�K�I�)L�$M�t� H�L$M�&H�}I�t$�.�����xbI��L�$u�I�4$H�} L�L$0H�,$������x?L�$L�L$0u�I�D$H9E|)�s���I�D$ H9E HL\$ML�ML��W����    M��M��H�\$�@���H�\I�D� I�L�c�L�|$I��K�l� I�wL�u I�~������x3t
M�L� �����I�7I�>�n�����xu�I�GI9F|�I�G I9F }�I�D$�M�t� H��H��?H�H��M��I������L��I���L��L�T$8�=���f�     AUI��ATI��U��S1�H���f.�     9�}HcӃ�A�D ��
tL���t������u߅�~3Hc�I�L��9
t>A�D  ���t%1�9���H��[]A\�D �A]�D  ���A�E  u�H��1�[]A\A]Ð� ��ff.�     ATE1�UH��SH��H�H�v�V�����uH�3H�} �F�����u
H�CH9Et[]D��A\�H�E H+C Hc�h  []H�H1�H)�H9�A��D��A\�f.�     AWAVI��AUATUSH��8�|$H�t$H���@.@ �	���H��H���3  E1�f.�     �8   1��$���H�D$ �D  ��'  HcӃ����h` ��
tOH���'������uڅ�u>��@   E1�A��D;d$ff�     E��u�H������H��8[]A\A]A^A_�@ HcӀ�h` 
�0  Ƃ�h`  ����~  ��'  �>  A��D;d$A�   ~���.@ ��h` ����H��I��t�1�f.�     ��wk���$Ũ0@ f��^.@ L������H����  �  L�L$ L���
   1�H�D$L�L$�i���H�L$L�L$�
   L�|$ 1�H�yI�A �F���I�G(f�1���.@ ���q���H��I���u����������I�FI;F��  H��H�T$ �u  H�I�FH��I�F����@ L������H�PH��H�T$����H�T$L��H������H��H�D$ H��o����    �a.@ L���   �H�D$ ��   �@0�E���D  H�L$ �
   1�L��H�L$�g���H�L$H�H�A�����    H�L$ �
   1�L��H�L$�7���H�L$H�H�A������    L�������H�PH��H�T$�W���H�T$L��H������H��H�D$ H�H����fD  �@0 �����    Ƃh`  �����@ H�D$ H�@(    H�@     �f���H�t$�'  ��.@ 1��`����`���1�����H�t$ L���G   �B���E1��&���H�t$�B.@ 1��)���1�����f�H�?H��t�S��� ��f.�     @ AUATA�   USH��H��H�WH+H��H����   L��H�t$����L�H�KH��H�t$L)�H��H�>H��    H�tvH�8H�;H�KH)�H��H��    L�l H��tH��H������H�;I��H��t����H�+L�L�kH�kH��[]A\A]�fD  H�H9�vI�������Z��� L��E1��H��������H��I������H9�LF��.����    AWAVAUATUSH��H��8H��H�t$(��  H�D$(H�@H��H�D$�\  H�@H��H�D$�'  H�@H��H�D$��   H�@H��H�D$ ��   L�xM����   M�gM��tkI�l$H��tJ@ L�mM��u
�%D  M��I�uH���Q���M�uL������M��u�L�mH���t���M��tL���I�l$L���]���H��tI���I�oL���G���H��tI���s���H�D$ H�hH���)���H��tH�l$ H���C���H�D$H�hH������H��tH�l$H������H�D$H�hH�������H��tH�l$H�������H�D$H�hH�������H��tH�l$H������H�D$(H�hH������H��t
H�l$(�m���H��8[]A\A]A^A_�f.�     D  ATI��USH�_H��u�&fD  H��H�sL���	���H�kH���=���H��u�[]A\� AWA��AVI��AUI��ATL�%  UH�-  SL)�1�H��H������H��t�     L��L��D��A��H��H9�u�H��[]A\A]A^A_�ff.�     ��f�H��H���                 r Error reading from file %s
 -- Novel * .comparedResults.txt w .comparisonSummary.txt Comparison run as:  %s  %s	%s	 %i	%i	%i--%i	%s	 Summary of results... %s	%i	%i
 -	-	-	-	 Error reading results file %s, line exceeded %i characters.
    Usage: ./comp linesToSkip supportPosTolerance outputBasename file1 file2 file3 ...
See splitPairs.cpp file file format of results files.        Unable to open file %s for writing.
    Finished running comparison.  Results written to %s.  Summary writen to %s.
    file 	#known unique to file 	#novel unique to file
     file 	#known unique to file 	#novel unique to file      �(@     �)@     P)@     �)@     0(@     ()@     ;�      �����   ����(  ����  �����   ����  8���P  �����  ����  ����@  (����  ����   ����x  ����`  �����  �����  H���@             zR x�      ���*                  zR x�  $      ����   FJw� ?;*3$"       D   ����           4   \   ����~    B�A�D �]
AEFAGE  L   �   ����O   B�B�B �B(�A0�A8�DP68A0A(B BBB       L   �   �����   B�H�B �B(�A0�A8�G�
8A0A(B BBBA   L   4   ����    B�E�D �C(�F0T
(A ABFFN
(C ABBB   4   �  p���f    B�D�D �o
AEAPAT   <   �  H����    B�B�G �A(�G@�
(A ABBG     L   �  h���~   B�B�E �B(�A0�A8�Dp�
8A0A(B BBBE            zPLR xP@ �  4   $   ����L  45@ A�CP������

A       L   �  `����   B�B�B �B(�A0�A8�Gp�8A0A(B BBB       ,   �  ����=    B�D�A �rAB          $  0����    D�D   <  ����e    B�E�E �E(�H0�H8�M@l8A0A(B BBB    �   ���               ��#w�  �	M� �
�
� �+  �� ��                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              `@     �@     @@                                  v             �             �             h@            $.@            �=`                          �=`                   ���o    �@             @            �@     
       b                                           @`                                        h	@            P	@                   	              ���o    �@     ���o           ���o    �@                                                                                                             �=`                     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@         GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-7) GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-9)  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.gnu.build-id .gnu.hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .gcc_except_table .init_array .fini_array .jcr .dynamic .got .got.plt .data .bss .comment                                                                                8@     8                                    #             T@     T                                     1             t@     t      $                              D   ���o       �@     �      (                             N             �@     �      `                          V              @            b                             ^   ���o       �@     �      H                            k   ���o       �@     �      �                            z             P	@     P	                                  �             h	@     h	                                 �             h@     h                                    �             �@     �                                  �             �@     �      �                             �             $.@     $.      	                              �             0.@     0.      �                             �             �0@     �0      �                              �             h1@     h1      �                             �             45@     45      '                              �             �=`     �=                                    �             �=`     �=                                    �             �=`     �=                                    �             �=`     �=                                  �             �?`     �?                                   �              @`      @                                              A`     A                                                 A`     A       O                                   0               A      X                                                   tA                                                         PJ      (         2                 	                      xU      '	                                                           8@                   T@                   t@                   �@                   �@                    @                   �@                   �@                  	 P	@                  
 h	@                   h@                   �@                   �@                   $.@                   0.@                   �0@                   h1@                   45@                   �=`                   �=`                   �=`                   �=`                   �?`                    @`                   A`                    A`                                       ��                     @     O      �     `@     �          �@     �       7    �`            F   ��                Q    �=`             ^    �@             s     @             �    @@             �     A`            �    �=`             �    `@             �    �=`             F   ��                �    05@                 �=`                  ��                     @`             -     �=`             >     �=`             Q    �=`             Z     A`             e                     y     .@            �    �@             �    �h`     '      �                      �                      �                     �                     �                     �                                          0  "  �+@     �      �    @A`            �    $.@             �                     �     '@     ~                                                (                     G                     �                     �  "  �*@            �     p@             �                          0.@                                 (                     <                      V    A`             c    �`            k                     ~    �@     ~       �    A`             �                     �                     �   8.@             �                     	    �-@     e            �`            &    `A`     '      +    �&@     f       N  "  �*@            u                     �  "  p-@     =       �  "  p-@     =       �    A`             �                     E                     Y    ��`            f                     z    Џ`            �                     �     �`             �                     �  "  �*@     �           ��`                �%@     �       6                     J    A`             Q     P@             r                     �                     �                     �    ��`     0       �                     �  "  �*@     �       	    �@     L      !	    h@              compare.cpp _ZSt16__insertion_sortIN9__gnu_cxx17__normal_iteratorIPP10RSW_resultSt6vectorIS3_SaIS3_EEEEPFbPKS2_SA_EEvT_SD_T0_.constprop.45 _ZSt16__introsort_loopIN9__gnu_cxx17__normal_iteratorIPP10RSW_resultSt6vectorIS3_SaIS3_EEEElPFbPKS2_SA_EEvT_SD_T0_T1_.constprop.40 _GLOBAL__sub_I__Z8get_lineP8_IO_FILEPci _ZStL8__ioinit crtstuff.c __JCR_LIST__ deregister_tm_clones register_tm_clones __do_global_dtors_aux completed.6342 __do_global_dtors_aux_fini_array_entry frame_dummy __frame_dummy_init_array_entry __FRAME_END__ __JCR_END__ _GLOBAL_OFFSET_TABLE_ __init_array_end __init_array_start _DYNAMIC data_start printf@@GLIBC_2.2.5 __libc_csu_fini _start sLine __gmon_start__ _Jv_RegisterClasses puts@@GLIBC_2.2.5 _Znam@@GLIBCXX_3.4 _ZdlPv@@GLIBCXX_3.4 _ZNSs4_Rep10_M_disposeERKSaIcE@@GLIBCXX_3.4 exit@@GLIBC_2.2.5 _ZNSt8_Rb_treeIP10RSW_resultS1_St9_IdentityIS1_ESt4lessIS1_ESaIS1_EE8_M_eraseEPSt13_Rb_tree_nodeIS1_E allResults _fini _ZNSt8ios_base4InitC1Ev@@GLIBCXX_3.4 _Z12read_resultsiPKcRSt6vectorIP10RSW_resultSaIS3_EE malloc@@GLIBC_2.2.5 fopen@@GLIBC_2.2.5 __libc_start_main@@GLIBC_2.2.5 _ZSt18_Rb_tree_decrementPSt18_Rb_tree_node_base@@GLIBCXX_3.4 __cxa_atexit@@GLIBC_2.2.5 _ZNSt6vectorIP10RSW_resultSaIS1_EED1Ev _ZNSt8ios_base4InitD1Ev@@GLIBCXX_3.4 _ITM_deregisterTMCloneTable _IO_stdin_used fputc@@GLIBC_2.2.5 strlen@@GLIBC_2.2.5 _ITM_registerTMCloneTable __data_start results fgetc@@GLIBC_2.2.5 _Z20compare_data_resultsPK10RSW_resultS1_ __TMC_END__ _ZNSsC1EPKcRKSaIcE@@GLIBCXX_3.4 strstr@@GLIBC_2.2.5 __dso_handle strtol@@GLIBC_2.2.5 __libc_csu_init myNovelCount temp _Z13overlapResultPK10RSW_resultS1_ _ZNSt6vectorIP10RSW_resultSaIS1_EED2Ev memmove@@GLIBC_2.2.5 _ZNSt3setIP10RSW_resultSt4lessIS1_ESaIS1_EED2Ev _ZNSt3setIP10RSW_resultSt4lessIS1_ESaIS1_EED1Ev __bss_start _ZSt29_Rb_tree_insert_and_rebalancebPSt18_Rb_tree_node_baseS0_RS_@@GLIBCXX_3.4 strcmp@@GLIBC_2.2.5 myKnownCount strtok@@GLIBC_2.2.5 supportPosTolerance _ZNSs6appendEPKc@@GLIBCXX_3.4 _end fclose@@GLIBC_2.2.5 _ZNSt6vectorIP10RSW_resultSaIS1_EE19_M_emplace_back_auxIIRKS1_EEEvDpOT_ numResultsFiles _Z8get_lineP8_IO_FILEPci fwrite@@GLIBC_2.2.5 _edata __gxx_personality_v0@@CXXABI_1.3 fprintf@@GLIBC_2.2.5 _Znwm@@GLIBCXX_3.4 _Unwind_Resume@@GCC_3.0 alreadyPrinted memcpy@@GLIBC_2.14 _ZNSt6vectorIP10RSW_resultSaIS1_EE19_M_emplace_back_auxIJRKS1_EEEvDpOT_ main _init                                                                                                                                                                                                                                                                                                                                                                  rsw_compare.sh                                                                                      0000775 0001753 0001753 00000001056 12563403132 013551  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

BASEDIR=$( cd ${0%/*}  >& /dev/null ; pwd -P )
source "${BASEDIR}/rsr_config.sh"

SUFFIX=".results"
OUTNAME="comparison.txt"

function do_compare() {
    tolerance=$1
    shift
    whereto=$(dirname $2)
    $COMPARE_PROG "$LINES_TO_SKIP" $tolerance $OUTNAME $@ >& /dev/null
    mv "${OUTNAME}.comparisonSummary.txt" "$whereto"
    mv "${OUTNAME}.comparedResults.txt" "$whereto"

}

#---------Main-------

if (( $# < 3 )); then
    yell "usage -- $(basename $0)  <support tolerance> <file1> <file2>"
    die  "you had $# : $@"
fi


do_compare $@
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  rsw_config.sh                                                                                       0000664 0001753 0001753 00000004657 12563403132 013377  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

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


                                                                                 rsw_pipeline.sh                                                                                     0000775 0001753 0001753 00000012043 12563403132 013726  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

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
    try $FORMAT_PROGRAM "${file}.bowtie.txt" >& "${LOG_DIR}/${LOG_FILE}"
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

#step 4: split the column into proper format
log "formatting..." 
results=$(split_columns ${results})
log "done formatting; result=$results" 


#step 5: select candidates
log "running rsr..."
result=$($RSR_SCRIPT $genome "$results" $readlength $@)
echo $result

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             rsw.sh                                                                                              0000775 0001753 0001753 00000004654 12563403132 012052  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

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
                                                                                    sbc                                                                                                 0000775 0001753 0001753 00000032703 12563406335 011401  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                ELF          >    �	@     @       "          @ 8 	 @         @       @ @     @ @     �      �                   8      8@     8@                                          @       @     �      �                          `     `     �      �                    (      (`     (`     �      �                   T      T@     T@     D       D              P�td   �      �@     �@     L       L              Q�td                                                  R�td         `     `     �      �             /lib64/ld-linux-x86-64.so.2          GNU                        GNU  Cc��-"rӎYM���(tX                         9�                            �                      )                      B                      o                      I                                            ~                      n                      =                      �                       6                      O                                            "                      /                      v                                            ]                      d                      V     � `             libc.so.6 sprintf fopen __strdup perror clock strtok strtol feof strlen fputs malloc stderr fwrite fcloseall fprintf getline __libc_start_main free __gmon_start__ GLIBC_2.2.5                                          ui	   �       �`        
           � `                    `                     `                   ( `                   0 `                   8 `                   @ `                   H `                   P `                   X `        	           ` `        
           h `                   p `                   x `                   � `                   � `                   � `                   � `                   � `                   � `                   H��H�=  H��t�   H���      �52  �%4  @ �%2  h    ������%*  h   ������%"  h   ������%  h   �����%  h   �����%
  h   �����%  h   �����%�  h   �p����%�  h   �`����%�  h	   �P����%�  h
   �@����%�  h   �0����%�  h   � ����%�  h   �����%�  h   � ����%�  h   ������%�  h   ������%�  h   ������%�  h   �����USH��H����~y������� ` H�Ź4  1�H���H�H�{�����H�{�9@ H�?!  �J���H��tYH���-  1���������H)�q@ �H*��   �Y?  ����1�H��[]�H�H�=  �;@ 1������   ��H�SH�=  �S@ �����6@ ������   �f�1�I��^H��H���PTI���@ H���@ H��	@ �;����f��     �� ` UH-� ` H��H��w]ø    H��t�]�� ` ���    �� ` UH-� ` H��H��H��H��?H�H��u]ú    H��t�]H�ƿ� ` ���    �=A   uUH���~���]�.  ��@ H�=x   t�    H��tU� ` H����]�{��� �s��� AUATUH��SH��H��tH�1��
   ������Lc�u'H��  H����   �    H��H��[]A\A]�f�Hc�H��� ` H��u�=5  ��   H������L�-F  H��L������H�|�1���A�@ H��H��H��L��!@ 1��c���H�߾)@ ����H��H��J��� ` �r���H�=  ��@ 1������6@ �����H��  �	   �   �+@ �����5���D  H�-�  A�5  H������H�x����A�@ H��H�ù@ H��!@ 1�������Z���ffff.�     U�5@ �7@ S�   H���v���H��u�/�    1��5@ �\���H��tH��H��H���u�H��[]� H��1�[]��    AUI��ATE1�U�   SH��(H�$    L����������   H�t$L��H������H�����   M��tL�������H�<$�7@ �>����5@ H��I������H��I��u�d@ 1��5@ ����H��I��tLL��H��H���u�L������H���k���H�<$H�������L���
������_���H��([]A\A]��    E1��f.�     �AWA��AVI��AUI��ATL�%x  UH�-x  SL)�1�H��H�������H��t�     L��L��D��A��H��H9�u�H��[]A\A]A^A_�ff.�     ��f�H��H���                 .txt Chr_unknown %s.%s%s w Ignoring. 	 Chr Usage: %s <reads file>
 Could not open %s for reading Elapsed time: %lf seconds
     Panic! Cannot open chrom file %s for writing!
  @     �����ư>;L       ����   @���x  ���h    ����   `���   ����8  �����   ����             zR x�      ����*                  zR x�  $      `���@   FJw� ?;*3$"    <   D   8���S   B�B�A �D(�D0w
(D ABBC      4   �   X���Y    A�K�I u
AADDCA     <   �   �����    B�E�D �F(�DP�
(A ABBH     ,   �   �����    A�A�G {
AAA     D   ,   ���e    B�E�E �E(�H0�H8�M@l8A0A(B BBB    t  (���                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   �
@     �
@                                  �@            �@            `                          `                   ���o    �@            �@            �@     
       �                                             `            �                           �@            �@            0       	              ���o    �@     ���o           ���o    h@                                                                                                             (`                     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     	@         GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-7) GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-9)  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.gnu.build-id .gnu.hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .init_array .fini_array .jcr .dynamic .got .got.plt .data .bss .comment                                                                                  8@     8                                    #             T@     T                                     1             t@     t      $                              D   ���o       �@     �      $                             N             �@     �      �                          V             �@     �      �                              ^   ���o       h@     h      *                            k   ���o       �@     �                                   z             �@     �      0                            �             �@     �      �                          �             �@     �                                    �             �@     �      @                            �             	@     	      �                             �             �@     �      	                              �              @            �                              �             �@     �      L                              �              @            �                             �             `                                         �             `                                         �              `                                          �             (`     (      �                           �             �`     �                                   �               `             �                             �             � `     �                                     �             � `     �       �	                              �      0               �       X                                                   !                                                         �)      X         -                 	                      �1      �                                                           8@                   T@                   t@                   �@                   �@                   �@                   h@                   �@                  	 �@                  
 �@                   �@                   �@                   	@                   �@                    @                   �@                    @                   `                   `                    `                   (`                   �`                     `                   � `                   � `                                       ��                    ��                      `             *     
@             ?     @
@             R     �
@             h     � `            w     `             �     �
@             �     `                 ��                �     �@             �      `                  ��                �      `             �     (`             �      `                   `             &    �@            6                     H                      d     � `             o    � `             v    �@            ~                     �    �@             �                     �                     �                     �                     �                         � `                                  )                     ;                      J                     ^   @             k     @            z    � `     �	      �    �@     e       �                     �    �*`             �    �	@             �    � `            �    � `             �    	@     �       �    �*`            �                     �                     �                                           $                     9    �@     �       ?                     T                     h   � `             t    0@     Y       ~                      �    �@             �                     �    �
@     S      �    � `             split_on_chrom.c crtstuff.c __JCR_LIST__ deregister_tm_clones register_tm_clones __do_global_dtors_aux completed.6342 __do_global_dtors_aux_fini_array_entry frame_dummy __frame_dummy_init_array_entry __FRAME_END__ __JCR_END__ __init_array_end _DYNAMIC __init_array_start _GLOBAL_OFFSET_TABLE_ __libc_csu_fini free@@GLIBC_2.2.5 _ITM_deregisterTMCloneTable data_start _edata file310 clock@@GLIBC_2.2.5 _fini strlen@@GLIBC_2.2.5 printf@@GLIBC_2.2.5 fputs@@GLIBC_2.2.5 __strdup@@GLIBC_2.2.5 __libc_start_main@@GLIBC_2.2.5 __data_start fprintf@@GLIBC_2.2.5 feof@@GLIBC_2.2.5 __gmon_start__ strtol@@GLIBC_2.2.5 __dso_handle _IO_stdin_used chrFiles __libc_csu_init malloc@@GLIBC_2.2.5 _end _start openFiles __bss_start main prefix fopen@@GLIBC_2.2.5 perror@@GLIBC_2.2.5 strtok@@GLIBC_2.2.5 _Jv_RegisterClasses getline@@GLIBC_2.2.5 parse sprintf@@GLIBC_2.2.5 fwrite@@GLIBC_2.2.5 __TMC_END__ get_chrom _ITM_registerTMCloneTable _init fcloseall@@GLIBC_2.2.5 pick_file stderr@@GLIBC_2.2.5                                                              sfc                                                                                                 0000775 0001753 0001753 00000032520 12563406334 011401  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                ELF          >    �	@     @        "          @ 8 	 @         @       @ @     @ @     �      �                   8      8@     8@                                          @       @     �      �                          `     `     �      �                    (      (`     (`     �      �                   T      T@     T@     D       D              P�td           @      @     D       D              Q�td                                                  R�td         `     `     �      �             /lib64/ld-linux-x86-64.so.2          GNU                        GNU ��Fb����7Oe�� �~�                         9�                            �                      G                      3                      c                      A                      :                      H                                            �                      j                      .                      �                       N                      z                                            '                      r                                                                  \                      U     � `             libc.so.6 exit sprintf fopen __strdup perror feof strlen memset fputs fclose malloc stderr fwrite strchr fprintf getline memmove __libc_start_main free __gmon_start__ GLIBC_2.2.5                                     ui	   �       �`                   � `                    `                     `                   ( `                   0 `                   8 `                   @ `                   H `                   P `                   X `        	           ` `        
           h `                   p `                   x `                   � `                   � `                   � `                   � `                   � `                   � `                   � `                   H��H�  H��t��   H���      �5  �%  @ �%  h    ������%�  h   ������%�  h   ������%�  h   �����%�  h   �����%�  h   �����%�  h   �����%�  h   �p����%�  h   �`����%�  h	   �P����%�  h
   �@����%�  h   �0����%�  h   � ����%�  h   �����%�  h   � ����%�  h   ������%�  h   ������%z  h   ������%r  h   ������%j  h   ����USH����~/��H�^H�l�f�     H�;H����  H9�u�1�H��[]�H�H�=,  �p@ 1�� ����   �ܐ1�I��^H��H���PTI�� @ H���@ H��P	@ �����f�f.�     f��� ` UH-� ` H��H��w]ø    H��t�]�� ` ���    �� ` UH-� ` H��H��H��H��?H�H��u]ú    H��t�]H�ƿ� ` ���    �=q   uUH���~���]�^  ��@ H�=�   t�    H��tU� ` H����]�{��� �s��� ATH��USH��tg���}���Hc�H9�vX�|�/t��Hc��D+	�+	[]A\�f.�     D�d�A�D$�<w�H�t�H�|��   ������E�H��/D�d�H��  �@@ �+   �   �����   ������     AU��@ ATUSH��H��(H�$    ����H��I����   H������L�hL���Q���L��H��1�H���������@ H�ھ�@ H��1��z���H��������@ �D H���0���H��I����   1��XH�t$L��H���0���H���tNH��tH������H�<$����H�,$�	   H��H���0���H)�H�߉�����L��H���(���L��������t�L�������L�������H��([]A\A]ÐH�=i  H�ھ�@ 1��:�����@ ����H��([]A\A]�H�=>  H���@ 1�������@ �e���L��1�������f�     AWA��AVI��AUI��ATL�%H  UH�-H  SL)�1�H��H�������H��t�     L��L��D��A��H��H9�u�H��[]A\A]A^A_�ff.�     ��f�H��H���                 No field supplied to split_field. aborting
     usage -- %s <file to split> [additonal files...]
       r Could not open %s for reading .split1stcolumn %s%s w Could not open %s for writing    ;D       ����   P���8  ����`   �����   @����   ����h   ����             zR x�      @���*                  zR x�  $      h���P   FJw� ?;*3$"    ,   D   �����    B�D�A �k
ABK    L   t   P���g   B�G�A �A(�GP�
(A ABBBd
(A ABBA   ,   �   ���W    A�A�D q
AAA      D   �   @���e    B�E�E �E(�H0�H8�M@l8A0A(B BBB    <  h���                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   p
@     P
@                                  �@            $@            `                          `                   ���o    �@            �@            �@     
       �                                             `            �                            @            �@            0       	              ���o    �@     ���o           ���o    �@                                                                                                             (`                     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     	@     	@     &	@     6	@     F	@         GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-7) GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-9)  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.gnu.build-id .gnu.hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .init_array .fini_array .jcr .dynamic .got .got.plt .data .bss .comment                                                                                  8@     8                                    #             T@     T                                     1             t@     t      $                              D   ���o       �@     �      $                             N             �@     �                                V             �@     �      �                              ^   ���o       �@     �      ,                            k   ���o       �@     �                                   z             �@     �      0                            �              @            �                          �             �@     �                                    �              @            P                            �             P	@     P	      �                             �             $@     $      	                              �             0@     0      �                              �              @            D                              �             H@     H      �                             �             `                                         �             `                                         �              `                                          �             (`     (      �                           �             �`     �                                   �               `             �                             �             � `     �                                     �             � `     �                                     �      0               �       X                                                   !                                                         �)      �         -                 	                      �1      �                                                           8@                   T@                   t@                   �@                   �@                   �@                   �@                   �@                  	 �@                  
  @                   �@                    @                   P	@                   $@                   0@                    @                   H@                   `                   `                    `                   (`                   �`                     `                   � `                   � `                                       ��                    ��                      `             )     �	@             >     
@             Q     P
@             g     � `            v     `             �     p
@             �     `                 ��                �     �@             �      `                  ��                �      `             �     (`             �      `                   `             %     @            5                     G                      c     � `             n    � `             u                     �    $@             �                     �    �
@     �       �                     �                     �                     �                     �                                          2    � `             ?                     T                     f                      u   8@             �    0@            �    �@     e       �                     �    � `             �    �	@             �    � `             �    P	@     W       �                     �                     �                                           "                     7    @@     g      =                     R                     d                     x   � `             �                      �    �@             �    � `             split_columns.c crtstuff.c __JCR_LIST__ deregister_tm_clones register_tm_clones __do_global_dtors_aux completed.6342 __do_global_dtors_aux_fini_array_entry frame_dummy __frame_dummy_init_array_entry __FRAME_END__ __JCR_END__ __init_array_end _DYNAMIC __init_array_start _GLOBAL_OFFSET_TABLE_ __libc_csu_fini free@@GLIBC_2.2.5 _ITM_deregisterTMCloneTable data_start _edata fclose@@GLIBC_2.2.5 _fini strlen@@GLIBC_2.2.5 split_field strchr@@GLIBC_2.2.5 fputs@@GLIBC_2.2.5 memset@@GLIBC_2.2.5 close@@GLIBC_2.2.5 __strdup@@GLIBC_2.2.5 __libc_start_main@@GLIBC_2.2.5 __data_start fprintf@@GLIBC_2.2.5 feof@@GLIBC_2.2.5 __gmon_start__ __dso_handle _IO_stdin_used __libc_csu_init malloc@@GLIBC_2.2.5 _end _start __bss_start main memmove@@GLIBC_2.2.5 fopen@@GLIBC_2.2.5 perror@@GLIBC_2.2.5 _Jv_RegisterClasses getline@@GLIBC_2.2.5 parse sprintf@@GLIBC_2.2.5 exit@@GLIBC_2.2.5 fwrite@@GLIBC_2.2.5 __TMC_END__ _ITM_registerTMCloneTable _init stderr@@GLIBC_2.2.5                                                                                                                                                                                 sp4                                                                                                 0000775 0001753 0001753 00000153043 12563406334 011340  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                ELF          >    l>@     @        �          @ 8 	 @         @       @ @     @ @     �      �                   8      8@     8@                                          @       @     ?�      ?�                    ؝      ؝`     ؝`     �      ��                    ��      ��`     ��`                                T      T@     T@     D       D              P�td   $�      $�@     $�@     �      �             Q�td                                                  R�td   ؝      ؝`     ؝`     (      (             /lib64/ld-linux-x86-64.so.2          GNU                        GNU xܩC�hO�K�XO�#A>�   1            `1   2   3   )��Y@�(E�LyIk�                                                 �                                                                   �                     >                                          �                      �                     �                      &                     �                     �                                          	                     �                     3                       �                     /                     �                     O                       ~                      }                     �                     #                                          �                     �                     	                     v                     i                       �                     �                      �                     �                     �                     �                     �                      �                     _                     �                                          �                     ]                     �                     �                     �                     %                     �    ��`            q    ��`             �       @             �     p@              libstdc++.so.6 __gmon_start__ _Jv_RegisterClasses _ITM_deregisterTMCloneTable _ITM_registerTMCloneTable __pthread_key_create _ZNSs4_Rep10_M_destroyERKSaIcE _ZSt11_Hash_bytesPKvmm _ZSt17__throw_bad_allocv _ZNSs4_Rep10_M_disposeERKSaIcE __cxa_rethrow _ZNSt8ios_base4InitD1Ev _ZNSsC1EPKcRKSaIcE _ZNKSt8__detail20_Prime_rehash_policy14_M_need_rehashEmmm __cxa_begin_catch _ZNSs4_Rep20_S_empty_rep_storageE __gxx_personality_v0 _ZNSs12_M_leak_hardEv _Znwm _ZNKSt8__detail20_Prime_rehash_policy11_M_next_bktEm __cxa_end_catch _ZNSsC1ERKSs _ZNSs7reserveEm _ZNSt8ios_base4InitC1Ev _ZdlPv libm.so.6 libgcc_s.so.1 _Unwind_Resume libc.so.6 strcpy sprintf fopen puts time popen strtok strtol fgetc strlen __cxa_atexit memset strstr memcmp stdout fputc fclose malloc fscanf fwrite fprintf memmove __libc_start_main setpriority free GCC_3.0 GLIBCXX_3.4.18 CXXABI_1.3.5 CXXABI_1.3 GLIBCXX_3.4 GLIBC_2.2.5                                                                 O         P&y   4              P   h��   <     uѯ   K     ӯk   X     t)�   c        l         ui	   o      ��`                   ��`        2           ��`        1           �`                    �`                   (�`                   0�`                   8�`                   @�`                   H�`                   P�`        	           X�`        
           `�`                   h�`                   p�`                   x�`                   ��`                   ��`                   ��`        3           ��`                   ��`                   ��`                   ��`                   ��`                   ��`                   Ƞ`                   Р`                   ؠ`                   �`                   �`                   �`                   ��`                     �`        !           �`        "           �`        #           �`        $            �`        %           (�`        &           0�`        '           8�`        (           @�`        )           H�`        4           P�`        *           X�`        +           `�`        ,           h�`        -           p�`        .           x�`        /           ��`        0           H��H��  H��t�C   H���              �5�  �%�  @ �%�  h    ������%��  h   ������%�  h   ������%�  h   �����%�  h   �����%ڎ  h   �����%Ҏ  h   �����%ʎ  h   �p����%  h   �`����%��  h	   �P����%��  h
   �@����%��  h   �0����%��  h   � ����%��  h   �����%��  h   � ����%��  h   ������%��  h   ������%z�  h   ������%r�  h   ������%j�  h   �����%b�  h   �����%Z�  h   �����%R�  h   �����%J�  h   �p����%B�  h   �`����%:�  h   �P����%2�  h   �@����%*�  h   �0����%"�  h   � ����%�  h   �����%�  h   � ����%
�  h   ������%�  h    ������%��  h!   ������%�  h"   ������%�  h#   �����%�  h$   �����%ڍ  h%   �����%ҍ  h&   �����%ʍ  h'   �p����%  h(   �`����%��  h)   �P����%��  h*   �@����%��  h+   �0����%��  h,   � ����%��  h-   ����AW�   AVAUATUH��1�S��1�H���   �\���1�������H�s�  �#$  H�}��-  �0  H�=��  �0  H�=ԍ  ��/  H�=��  ��/  H�=��  ��/  H�=��  ��/  H�=d�  ��/  H�=x+! �:  H�D�  L�=5�  L9�H�D$8�`  H��I��H���������L)�A�?   H��H��H��H��H��?H�I)�M�H��/  ��%  M����  I�G0H�t$8H�D$XI��M)�I�O@I��I��L��L��H��?I�I��K�dH��L�H�P9���  ��  H�~�9���   ��  9��"  f���  L�M�M�OM�GA�I�O M�L�XI�W(M�_L�XM�_L�XM�_L�X M�_ L�X(M�_(L�L�HL�@@�xH�H H�P(I��M���   f�I�@�;PI����   ��   L9���   H�M�$I��0M�T$�E�L$�I�|$�I�T$�I�\$�H�XI�\$�H�XI�\$�H�XI�\$�H�X I�\$�H�X(I�\$�L�L�PH�HD�HH�x H�P(I�W@ I�L$L��9�x�^���I�G I9D$ �O���I��0��f�     H�x I9 �H���H��0�-���f�     �`?@ L��L���U  L��L)�H=/  �f  M����  L���!���H�x I9P�P���H�~�9��i�����   9�xyf�     tdM�W0M�M�GI�A�OI�W M�M�W8I�G(M�O0M�G8M�WM�W@I�@M�WM�WHA�OHM�WM�WPI�WPM�W M�WXI�GXM�W(�M���H�F�I9GP}�L�V�M�M�GI�A�OI�W M�L�V�I�G(M�WL�V�M�WL�V�M�WL�V�M�W L�V�M�W(L�N�L�F�H�~��N�H�V�H�F������H�V�H9P �u�������H�N�I9OP�����M���1������H+��  �@�@ H��1������H�=��  �.  H�=�'! �:  1�����H+L�  ���@ H��1�����H�=^�  �i.  H�=�'! �1  1��V���H+�  �Ђ@ H��1�����H�=)�  �4.  H�U�  L�=F�  L9�H�D$H�h  H��H��I���������L)��?   H��H��I��H��H��?H�H)�H�H��/  ��  H���4  I�G0H�t$HL�t$0H�D$@L)�I�O0H��H��H�\$0H��H��?H�H��H�[H��L�H�9��(  ��  H�~�9��8  f��D  9���  f���  L� I�E�_E�WE�OM�GM�'L�`I� A�O(A�W,M�gL�`M�gL�`M�gL�` M�g L�`(M�g(H�D�XD�PD�HL�@H�x �H(�P,H��H�\$@�   �H�G�;H���  ��   H9��  L�0D�kH��0D�c�D�[�L�S�L�K�L�s�L�pD�C��K�L�s�L�pL�s�L�pL�s�L�p L�s�L�p(L�s�H�D�hD�`D�XL�PL�H D�@(�H,I��H�I��9�x6�`���A�G8C|&�P���A�G9Cx�A���I�G H9C �3��� H��0�f.�     D�PE8W|%�$���D�XE9_x����L�P M9W �����H��0������    ��?@ H��H��� W  H��L)�H=/  ��  H���  H��������xA8@|'������xA9Hx@ �����H�x I9P�����H�~�9�������b  9���   f���  I�_0M�E�WE�OE�GI�I�I�_8I�O A�W(A�G,M�_0I�_I�_@E�W8E�O<E�G@I�_I�_HI�HI�_I�_PI�OPI�_ I�_XA�WXA�G\I�_(������~�@8x|"�[����~�9xx�M���H�~�H9x �?���H�^�M�E�WE�OE�GI�I�H�^�I�O A�W(A�G,I�_H�^�I�_H�^�I�_H�^�I�_ H�^�I�_(L�^�D�V�D�N�D�F�H�~�H�N��V��F��8����N�A8O@�����������N�A9OH����������H�N�I9OP���������� �V�8P�����������V�9P�o����y���H�V�H9P �[����f���1�����H+j�  ��@ H��1������H�=|�  �)  H���  H+��  H�����������      H��H���,
  �x���H�D$00   1ɉD$XH�T$hH�|$PfD  ��tH�E�  H�|$0H�T8�H9T8�t���  H;L$P��	  H�|$0L��  H�o�L�I�4)H�H;�  D�fD:g�   H�VH;W��  D�D$XH�\$0A)�H��H�L$HL�L�<II���L�    L9���  L���  H��0I�4)I�<H�H;��  D�fD:g��  H�VH;W��  L�G H�N D��!! M��I)�Mc�M9��r  D�^D:_t��vw;5�!! u�A��L�_  A��RH��H���g  I�|9 E�d1 Ed1��H�|$8Lc�D)�A9��B���;54!! �6���H�=��  L���  H�|$@L)T$@H�|$@H�|$@H��tKL�_Mc�   E1��D  H�wL9�t,I��H��K�4�L9n(u�L9v0u�H;u�H9Vu�L9L$@�����H�=T�  L�U�  I)�I��L�L$hM����  I���   E1��f�     H�FI��L9���   H��K�Dm H��H�H9Pu�L�P L9�|�M9�|�H�@(H9��I9���D$`   �ZD  H�AH�D$H�    H�L$HH�D$00����A��+tH��H�������    A��-�����H��H�������D$`    ��   �1���H�xhH�@P    H�@X    H�@`    �@h  �?�
   H�@p    I���y���H��������I�FPH9��  H��H��H�D$@�����H�T$@1�H��H�D$@�M����D$`H�L$@L��$�   ��I�NH��   �   ���@ I�FH��H��  A�VMc�I�~HH��$�   H�H�HI�I�N�@M�f8M�f(A�F HcD$8I�F@I�F0H��$�   L��$�   ��5  H�6�  H;7�  �1  H��H��$�   �   H�H��  H��H��  ����K�Dm 1�H��H$�  H� �F���L�D$XH�L$P��uH��H��H��?H�H��H9���  H����  f.�     I��0H�|$0_�|$8L� D�hD�`D�XH�hH�X �x(D�p,��  M�W`M��   H�l$@M���Af�     D�L$0H�H��0D�hD�X@�pH�xL�@ D�P(D�H,L��I��0I9���  H�] A;H����   ��   D�M,D�mH�U��uD�]H�}L�E D�U(D�L$0@ ;H��x1�{���@8r%D  �j���;zxD  �Z���L9B �P���L�
H��0L�L�J8L�HL�J@L�HL�JHL�HL�JPL�H L�JXL�H(H���A�@8}|#�V���A�9}x�G���I� H9} �9����E,H���������D�]D�UD�ML�EH�M �D$0H��D�m(L)�H��H��H��tGH��L��L��H��H)�H�L$`L�D$XD�L$PD�T$@D�\$8�����H�L$`L�D$XD�L$PD�T$@D�\$8�D$0H��0I�E�_E�WE�OM�GA�G,L��I��0I9�I�O E�o(�g���L��H9L$HL�t$Hu2����� H��0L� �XI9ΉhD�HL�PL�X D�h(D�`,������YD�IH�Q�L��iH��L�QL�Y D�i(D�a,�     D;H��x u�D8J|�D;Rxu�L;Z f�}�fD  H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���L��1�H�L$@��uH��H��H��?H�H��H9���  H����  I��0H�|$0_L�(L�`L�@@�hL�X H�X(��  M��   M�g`H�l$X�1fD  L�L�PH��0H�XD�HL�@ H�p(L��I��0I9��  H�]A;_H����  tiD�ML�] H�U�L�UL�E H�u(f�     ;ZH��xu�L;B }�H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���L�u M;w }�H��H���������L�U L)�L�MD�EH��H�M(H��H��t=H��L��L��H��H)�H�L$PD�D$HL�L$@L�T$0�X���H�L$PD�D$HL�L$@L�T$0L��H��0I��0I9�M�M�OI�_E�GM�w I�O(�����L9l$8H�l$8u1�a����     I��0L�H�ZL9�H�zD�ZL�B L�J(�5���E�]M�U I�E�I�]I�}L��M�E M�M(@ ;xH��xu�L9@ ~�H�0H��0H�2H�p8H�rH�p@H�rH�pHH�rH�pPH�r H�pXH�r(H���L�u ����H�~�H9x �*��������F�A8G@������'����F�A9GH�����@ ����H�F�I9GP�����������    1������1��    ����H+��  �@�@ H��1��E���H�=�y  ��  H�5��  H�=��  ��?@ ��P  1������H+��  ���@ H��1�����H�=�y  �  H�=y�  H�j�  H��H)�H��H����  H�D$0   L�T$01��    H��L��H)�H��H��H9��d  N�<�    H��L�$�    J�l9�H�UH9S�=  M���   �    H;W8}H�W8H�P0H;W@~H�W@H�H��$�   H��HA��H��$�   H��$�   �`.  A����   H���  H�=��  H��H)�H��L9���   J�,9J�!I��I��H�UH9S��   L�K(H�C0L�E(L)�H��H��?H1�H��H�E0H)�L)�H��H��?H1�H)�H9�ucH�EH9CuYL��Hc5�! L)�H�H1�H)�H9��o���H�E H9�b���E1�E���  H��H��H�` H�P(�����H�W8H�P0�����H��H)�H��H�T$0H9�H�rsH�t$0I���g���1��   �����H+��  ���@ H��1�����H�=�w  �  H�5z�  H�=k�  �0@@ �N  1�����H+[�  � �@ H��1������H�=mw  �x  1��q���H+2�  �H�@ H��1�����H�=Dw  �O  1��H���1�H���  �:���H+��  ���@ H��1��d���H�=w  �  H�=w  �  H�=Uw  �   H�=)w  ��  H�=-w  ��  H�=w  ��  H�=�v  ��  H��v  �R   �   �؄@ �E���H�n�  H�o�  1�H)�H��H����   D  H��H�=�v  �]{@ H�P@D�H L�@H�HH�T$ H�P8H�T$H�P`H�T$H�P0H�T$H�P(H�$H�1������H���  H��H�XXH��t+f�     H�CH�=%v  �y{@ H�1�����H�H��u�H�5v  �
   �-���H���  H���  H�MH)�H��H9�vH��H���4���H�v  �`   �   �0�@ A�   �3���H��u  �`   �   �0�@ ����H��u  �`   �   �0�@ �����H��u  �`   �   �0�@ �����L��  L��  1�L��L)�H��H����  �    M�4�L�u�  1�L�,�    I�H��I��H�S�  H�4�H��H���l  L�I�zL��f�     H9��\  H�	H��tH�y1�H��I��H9�t�Hc�! I;F`�2  I�BH9�t/I�
H���"  H�q1�H��I��H9��  H��L��I��H9�u�H�> ��  L�D�  H�5E�  H��  L)�N�(H��H��m۶m۶mH��A�@yH����   M�H0Hct! L�^M+H(�   1�� H�AL9�t_H��H��H�<�    H��H��H)�L�L9H u�I�x8H+x(H��H��?H1�H)�H9�|�I�x@H+x0H��H��?H1�H)�H9�|�H9�vA�@y A�x �=
  H�-�s  H��s  L��H���$  H��
   ����H�0�  H��J�4(�  H�޿	   ����H��  J�(L�pXM����   I�F�~{@ H��H�1������M�~L�v�  1�I�/H��I��H�\�  H��I��H��t.H�L�Af�     L9�tIH�	H��tL�A1�L��I��I9�t�   ����H��tH�     I�H�PH��H��L�� �` �8  M�6M���W���H�޿
   �����L�B�  L�C�  L��I�T$L)�H��L9�vL��I���4���H�=�r  ����H�=yr  �t���H�=}r  �h���H�=Qr  �\���H�=5r  �P���1�����H+��  ���@ H��1��C���H�=�q  ��  H���  H���  H��H)�H��H����   fD  H�h�H��teH�}XH��u	�@ H��H��-���H��u�H�EPH�}H1�H��    �����H�}HH�E`    H�EX    �����H�������H�)�  H�*�  H��H��H��  H)�H��H���t���H���   1�[]A\A]A^A_�H��$�   ���` �~-  ����I�G0H�\$8H�D$XL��L)�H��_�/  H��H���������H�\$`H��H�x�H�|$0H�|$0H��L�t$0H����H�D$@H�|$@H�|$PIk�0L�I��H�D$hK�D6H�D$xHk�0L�H�D$pL�l$@L��M� I�hM�HA�XM�P M9�M�X(�  L��L�t$H�WuL�p L;w HL�HL�H�IH�0H��L�L9�H�1H�pH�qH�pH�qH�pH�qH�p H�q H�p(H�q(�|  H��H�AH� H�r�H��H�<vH��L�H��L�D�wD9p�y���H��H���H�\$`H��L)�H��_�$���L�S��     I�M�*H���������M�bM�BA�jM�Z I�I�GI�Z(I�BI�GI�BI�GI�BI�G I�B I�G(I�B(L��L)�H�D$0H��H��H�G�H�|$@I��I��?I�I��M���=���1��af.�     uL�p L;w HL�HL�H�IH�0H��L�L9�H�1H�pH�qH�pH�qH�pH�qH�p H�q H�p(H�q(�����H��L�qK�6H�r�H�vH��I�<J��H��D�wL�D9p�x���H��H���L�t$HH�|$P uGH9T$0u@H�|$pH�L$hH�T$xH�H�H�GH�AH�GH�AH�GH�AH�G H�A H�G(H�A(H��I9���   I��0M��L� H�hL�H�XL�P L�X(�?���I���t���H�TH�RH��L�H�1H�0H�qH�pH�qH�pH�qH�pH�q H�p H�I(H�H(H�J�H��H�IH��L�D9@xtH�RH��L������L9X �`  H�RH�0H��L�H�2H�pH�rH�pH�rH�pH�rH�p H�r H�p(H�r(H�Q�H��H��?H�H��H�������H��H���u���H�B�H��H��?H�H��H�IH��L�D9HxtH�RH��L������L9P ��  H�RH�0H��L�H�2H�pH�rH�pH�rH�pH�rH�p H�r H�p(H�r(H�Q�H��H��?H�H��I9��s���H��H���u���I�G0L�d$HH�D$@L)�H��_��   H��H���������I��L��I�m�H��Hk�0L�� H��L�H��0L��$�   H�{8H��$�   H�s@H��$�   H�KHH��$�   H�SPH��$�   H�CXH�|$H�t$L��H�L$H�T$ ��?@ H��$�   L�$L��H�D$(H���i9  H���x���L��L)�H��_����M�T$� A�B(M�H���������E�jE�bE�ZI�j�D$8I�I�Z E�r,I�I�GI�BI�GI�BI�GI�BI�G I�B I�G(I�B(L��L)�H�D$0H��H��H�G�H�|$PI��I��?I�I��M����  1�L�D$X�EH�IH�0H��L�I9�H�1H�pH�qH�pH�qH�pH�qH�p H�q H�p(H�q(�����H��H�AH� H�r�H��H�<vH��L�H��L�D�D9 x9u�D�GD8@|,�D�GD9@x �r���L�@ L;G HL�HL��]���D  H��H���M���H�TH�RH��L�H�1H�0H�qH�pH�qH�pH�qH�pH�q H�p H�I(H�H(H�J�H��H�IH��L�D9 x.tH�RH��L��&���D8X|�9hxu�H9X  �D  H�RH�0H��L�H�2H�pH�rH�pH�rH�pH�rH�p H�r H�p(H�r(H�Q�H��H��?H�H��H�������H��H���d���L��1��k����m  �����E1��D$`    �����L�������H�-�i  H��i  ����I�o0H9l$H�����M�g`L�t$0�8 H�D�pD�PD�@H�pL�H D�X(D�h,L��I��0H��0H9D$H�K���H�] A;H����   ��   D�uD�EH�U�D�UH�uL�M D�](D�m, ;H��x)u�D8B!��z���;rxD  �j���L9J �`���H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���A�@8}|#�c���A�9}x�T���I� H9} �F���H��H�T$0D�]L)�D�UD�MH��L�EH�M H��D�u(D�m,H��tGH��L��L��H��H)�H�L$`L�D$XD�L$PD�T$@D�\$8����H�L$`L�D$XD�L$PD�T$@D�\$8I�E�_L��E�WE�OM�GI�O E�w(E�o,�y���Hc�! I;F`������   �,���H��tH�     I�H�PH��H��H�� �` �V-  �����Hk�0L�����Hk�0L��j���H��L������H��������q���H��H�������M�g0L9d$8�0���M�o`H�l$0�1@ L�H�pH�XD�HL�@ L�X(L��I��0I��0H9D$8�����I�\$A;_L����   tnE�L$M�$I�T$�I�t$M�D$ M�\$(f�     ;ZH��xu�L9B ~�H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���M�t$ M;w }�L��H�T$0M�$L)�M�L$E�D$H��I�l$(H��H��t3H��L��L��H��H)�D�D$PL�L$HL�T$@����D�D$PL�L$HL�T$@M�M�OL��I�_E�GM�w I�o(�����M�t$ �x���Hk�0L�����f.�     AT�T@a USH���������H����(z@ �T@a � @ �����
   ���` H��      H��      H��      �ߴ    �?H�ܴ      �_���H9�H���  ��  L�$�    L������1�L��H��H���?����(z@ �`�` �`Y@ H�-i�  ������
   �@�` H��      H��      H��      ��    �?H��      �����H9�H�ݳ  �  H�,�    H���/���H��1�H��H�������(z@ � �` � Y@ H���  �d����(z@ ���` ��X@ H�J�      H�G�      H�D�      �/����(z@ ���` ��X@ H���      H��      H��      ������(z@ ���` ��X@ H���      H���      H���      �����[]A\�(z@ ���` ��X@ H�G�      H�D�      H�A�      ��������� 1�I��^H��H���PTI��z@ H���y@ H���@ �;����f��     ���` UH-��` H��H��w]ø    H��t�]���` ���    ���` UH-��` H��H��H��H��?H�H��u]ú    H��t�]H�ƿ��` ���    �=�b   uUH���~���]��b  ��@ H�=�^   t�    H��tU��` H����]�{��� �s��� �F9Gx �    t	���    H�F H9G ���@ �   �f��9x:�    t�Ð�V�   8W|�    �N9Oxu�H�F H9G ���D  �   �f.�     �F9Gx8�    t	���    H�N(�W0L�G(�D)�;V0xu�I9����f.�     �   �f.�     �F9Gx0�    t	���    H�W0V(+W(;V0xu�H�G`H9F`��� �   �f�AUI��ATI��U��S1�H���f.�     9�}HcӃ�A�D ��
tL���������u߅�~3Hc�I�L��9
t>A�D  ���t%1�9���H��[]A\�D �A]�D  ���A�E  u�H��1�[]A\A]Ð� ��ff.�     H�7�����   H���   f�     ��H��Hc�Hi� ���H�� ���)ȉ���)����V�)Ȅ�u�H�w���t=H��D  ��H��Hc�Hi� ���H�� ���)ȉ���)����V�)Ȅ�u��W��Hc�Hi� ���H�� ���)ȉ���)���)�ø   넋9x2�    t�Ð�V�   8W|�    �H�GH+F����    �   �f�H�G`H9F`���@ AU�0z@ I��ATUSH��X�y���H��I���T  H�$��` ���` 1��@ �C�H��L����������t4H���a ��   ��
uڃ��C� ����   Hc�H����f�     L���H�������   H�$H�|$�
   1�H�8�  �����H�|$�
   1��a�  �����F�  H�D$�
   H�|$(1�H���  H�D$ H���  �}���H�|$0�
   1����  �f���H�|$@���  �
   H�D$81�H���  �C����m�  H��X[]A\A]ÿ�{@ �'  L��1�����1������L���o����+�����{@ ����1�������2z@ L��1��Z���1������ H���|@ H�<�  Iz@ �r�  @�  �X�  !   H��  Vz@ H���  �{@ �(�     ��     H���  bz@ ���     �����H|@ �
�����|@ H�������fff.�     U��|@ SH��H��8�l�  L���  D���  ���  H���  �D$ H�V�  H�D$���  �D$���  �D$H�F�  H�$1������H�-T�  H���l���H=�&  wH��8[]�H��H�꾘~@ 1�����1�����@ H��H���  �jz@ ��a 1��D����uz@ ��a ����H��H��\  �-  ��a ��~@ 1������H���  �wz@ ��a 1�������uz@ ��a �]���H��H��\  ��   ��a �0@ 1�����H�U�  ��z@ ��a 1������uz@ ��a ����H��H�[\  ��   ��a �x@ 1��D���H��  ��z@ ��a 1��l����uz@ ��a �����H��H��[  tY��a ��@ 1�� ���H���  ��z@ ��a 1��(����uz@ ��a ����H��H��[  t��a ��@ 1�H��������~@ ��a 1�����1�����H��(�~y ��z@ ��z@ L�N0L�F`H�NHD�L+N(H�VH�D$H�F@H�D$H�F8��z@ H�$1������H��(�ffffff.�     AU�0z@ I����z@ ATA�   USH��8  �����H�\$@H��H�$    H�D$    �D$    H�D$     H�D$(    �D$0    �fD  ��z@ H��L����u  1�H�ھ{@ H���������t�L�d$ H��H��H�������������1��M���L��4   �   �8�@ H���  ����H�P�  H+A�  L��p�@ 1�H��H���������  L�ﾠ�@ 1�����H���  H+�  L��Ѐ@ 1�H��H������H���  H+��  H��m۶m۶mL�� �@ H��H��H��1��N���H�_�  H+P�  L��0�@ 1�H���-���H��  L��`�@ 1�����L��H��L�ﾐ�@ 1�����H��  H+�  L��ȁ@ 1������L��
   �h���H��8  []A\A]�f.�     L�d$ H��{@ H��1�H��L���c�������fffff.�     AVI���i�AUATUSH�?H�w������H�-U�  1�I��H��H�>�  H��I��H��t(H�H�K I9�tkH�H��tH�K1�H��H��I9�t�   �0���H��H��tH�xH�     L���e���H�C    L��H��L��`�` �z  []A\A]�   A^�@ I�>H�sH�W�H;V�u��2������w���H��1�[]A\A]A^�H���t���H���l�������H���/���H�������    AWAVAUI��ATUSH��8������~ZHcп{@ �   I�t���H  ��t;I�t��${@ �   �u(L��){@ �@�` 1�������0z@ �@�` �@���H����0z@ L���>���H��H���  A�    @ �0   E1������H���%fD  A��'  Ic�A����`�` ��
��  H���������u�E����  �t~   E1�X{@ �`�` �����H���l  E1�fD  H�T$H�|$ H���^���H�D$ H�p�H9p��k  A��w%D���$���@  H�|$ �
   1��G����C(@ 1��X{@ A���`���H���w  H�T$ H�z�H����` t�M��H�O��  ����������g���H�t$H�D$����H�D$�N���fD  H�D$ �P���xH�|$ �����H�D$ � �C�r���f�H�|$ �
   1�����H�C �V���fD  H�|$ �����H�@H�C�9���f�     H�D$ �H���xH�|$ �z���H�D$ � �C�
���f.�     H�|$ ����H�@H������f.�     H�|$ �
   1������C������    H�|$ 1����������    IcԀ�_�` 
�|  Ƃ`�`  ����  A��'  �z  A�   �����H�D$ H�x�H����` ��  A���c  H�3�   �H����t5D  ��H��HcșHi� �H�� ���)щ���)�)Љ��F���u�H�s���t=H��f�     ��H��HcșHi� �H�� ���)щ���)�)Љ��F���u��C��Hc���Hi� ���H�� ���)Љ���)�H�)�  )�H;(�  �K,tsH��t.H�H�H�SH�PH�SH�PH�SH�PH�S H�P H�S(H�P(H��0H�ݢ  H���-���A���k���H�������H��8[]A\A]A^A_�Ƃ_�`  ����H�޿��` ��  뼺'  L�� �@ 1�����H��������L��{@ ������J��q��r��������<{@ L��1��g���1������M��H�Wt(����������<���H�t$�l����-���E1������P��J��H�����H���Z����%����5���H��H�D$ H�t$H�x��Z���H������f�AW�0z@ AVI��AUATUSH��8����H��I���   A�     �0   1��4���H���$�    ��'  HcӃ���`�` ��
��  L����������uօ���  ��y   E1��X{@ �`�` ����H����  1�f�H�T$H�|$ H������H�D$ H�H�H9H��;  ��w&���$�(�@ D  H�|$ �
   1�����H�E  1��X{@ ������H���H  H�T$ H�z�H����` t�M��H�O�  ����������h���H�t$H�D$�����H�D$�O����    H�|$ �n���H�@H�E�@ H�D$ �P���xH�|$ ����H�D$ � �E�Z���f.�     H�|$ �&���H�@H�E �9���f�     H�|$ �
   1�����H�E(����fD  H�|$ �����H�@H�E�����f�     H�|$ 1��4��������    HcӀ�_�` 
��   Ƃ`�`  ����N  ��'  ��   A�   �5����H�D$ H�x�H����` ��   ��
��   H�E H�U(H9�~
H�H�U H�E(H�+�  H;,�  ttH��t/H�U H�H�UH�PH�UH�PH�UH�PH�U H�P H�U(H�P(H��0H��  H���S���A���I���L�������H��8[]A\A]A^A_�Ƃ_�`  �%���H����` �Y	  뼺'  L��� �@ 1��þ��H�������밋J��q��r���������<{@ L��1�蚾��1�����M��H�Wt(��������������H�t$蟿�������E1�������P��J��H�����H�������X�������H��H�D$ H�t$H�x�荾��H������D  AW�0z@ AVAUATUSH��8H�<$贾��H��I����  �     �8   1��d���H���$�    ��'  HcӃ���`�` ��
�+  L���������uօ��  ��u   E1��X{@ �`�` �M���H���  1�f�H�T$H�|$ H���޾��H�D$ H�H�H9H���  ��w~���$�X�@ D  L�|$ I�G�H�x�ֽ��L��I��H���ؾ���Z{@ L��蛾��H��I���g  �
   1�L���  荾��I�H�E(�
   1��y���H�E0L�������D  1��X{@ ��艾��H���p  H�T$ H�z�H����` �*����    H�OH���!  �������������H�t$H�D$訽��H�D$�����fD  H�D$ �P���xH�|$ ����H�D$ � �E�j���f�H�|$ �
   1�追��H�H�E �L���@ H�|$ �����H�@H�E�1���f�     H�|$ �����H�@H�E ����f�     H�|$ ����H�@H�E�����f�     H�|$ 1������T����    H�E0    H�E(    ���� HcӀ�_�` 
��   Ƃ`�`  ����N  ��'  ��   A�   ����f�     H�D$ H�x�H����` ��   ����   H�Қ  H;Ӛ  t|H��t7H�U H�H�UH�PH�UH�PH�UH�PH�U H�P H�U(H�P(H�U0H�P0H��8H���  H������A�������L��谼��H��8[]A\A]A^A_�Ƃ_�`  �,���H��` �h  �H�4$�'  � �@ 1�聺��H��蹺��믋J��q��r��������H�4$�<{@ 1��W���1�������    H�WH��t(��������������H�t$�W��������E1��~����P��J��H�����H���E�����������H��H�D$ H�t$H�x��E���H���m����fff.�     AT1�I��USH��H�L�OH��I��H�H��H��H��t%H�L�AL9�tkH�	H��tL�A1�L��I��H9�t�   H�t$����H��tH�t$H�     L�L�VL�HL�PH��H��L��H���R
  H���   []A\�@ H��H��1�[]A\�f�ATI��UH��S1�H����tfL��H��H�x` H�J(tKH;H8}H�H8H�J0H;H@~H�H@H�
H�xHH���H�T$H�$�������u�H��[]A\�f.�     H�H8H�J0�H��L���f.�     @ H�?H��t�ø�� ��f.�     @ H�?H��t飸�� ��f.�     @ H�?H��t郸�� ��f.�     @ H�?H��t�c��� ��f.�     @ UH��SH��H�H��u	�@ H��H��-���H��u�H�EH�} 1�H��    ����H�E    H�E    H�} H��[]����f�ATI��USH��H�_H��t0�    H��u�zH��H�CH�+H�x�H����` uIH��豷��H��u�I�D$I�<$1�H��    �d���I�<$I�D$    I�D$    �y���H��[]A\ù������H����H�t$�H����H�CH�+H�x�H����` uH���:���H��t�H���ڋP��J��҉H��H�t$�	������    AUH���������ATU�0   SH��H��H�WH+H��H��H����   H��H�t$����H�;H�SI��L��H���������H)�H��I��H��L��H���   H�t$H�H�H�VH�QH�VH�QH�VH�QH�V H�Q H�V(H�Q(H�;H�SH)�H��H��H��M�,H��tH��L��蠷��H�;I��0H��t�/���L�L�#L�kH�kH��[]A\A]�fD  H�H9�vH���������� L��E1��Hk�0H�UUUUUUUH������H9�HF������f.�     D  AUH���������ATU�0   SH��H��H�WH+H��H��H����   H��H�t$迷��H�;H�SI��L��H���������H)�H��I��H��L��H���   H�t$H�H�H�VH�QH�VH�QH�VH�QH�V H�Q H�V(H�Q(H�;H�SH)�H��H��H��M�,H��tH��L���P���H�;I��0H��t�ߴ��L�L�#L�kH�kH��[]A\A]�fD  H�H9�vH���������� L��E1��Hk�0H�UUUUUUUH������H9�HF������f.�     D  AUH��m۶m۶mATU�8   SH��H��H�WH+H��H��H����   H��H�t$�o���H�;H�SI��L��H��m۶m۶mH)�H��I��H��L��H���   H�t$H�H�H�VH�QH�VH�QH�VH�QH�V H�Q H�V(H�Q(H�V0H�Q0H�;H�SH)�H��H��H��M�,H��tH��L�������H�;I��8H��t至��L�L�#L�kH�kH��[]A\A]�fD  H�H9�vH���������� L��E1��Hk�8H��$I�$I�H������H9�HF�������    AUATA�   USH��H��H�WH+H��H����   L��H�t$�,���L�H�KH��H�t$L)�H��H�>H��    H�tvH�8H�;H�KH)�H��H��    L�l H��tH��H�������H�;I��H��t臲��H�+L�L�kH�kH��[]A\A]�fD  H�H9�vI�������Z��� L��E1��H��������H��I������H9�LF��.����    AWAVAUATUSH��H9�H�<$H�t$��   L�oI9���   H��L�w�Ff.�     L��H+$L��H��H��tH�4$H��L��H)�����I��I��H9\$H�$L� tRH�$I�} L��H�0�Մ�M�e u�M�}��fD  I�H�L��I��I�7L���Մ�u�L�#I��L��I��H9\$u�H��[]A\A]A^A_� AWI��AVAUI��ATI�̹   UH��H� SH��H�W�H�w�L�w蜳������   H��������H��H9��.  L�<�    L������1�H��L��I��螰��H�u1�H�E    L�UH��u�UD  I�H�I�H�2H��t?H��H�F1�H�H��M��M�M��u�H�EH�H�uM�H�> ��   H��I�4�H��u�H�} �W���L��1�H�]H��L�u I��I��M�l$K�>H�H��t/H�I�$H� L� H�EH��L��[]A\A]A^A_� L�u �f�H�UI�$L�eI�$H��tH�B1�H�uL��M�$�HE H�UH��D  H������蛱��H��賱��L�u(����H���r���H��蚱��I�D$H�x�H����` uL��耯��軯���    H�WH��t'��������H�t$�H�����H������H��薱���P��J��H�����f�     AWI��AVAUI��ATI�̹   UH��H� SH��H�W�H�w�L�w茱������   H��������H��H9��6  L�<�    L������1�H��L��I��莮��H�u1�H�E    L�UH��u�UD  I�H�I�H�2H��t?H��H�F1�H�H��M��M�M��u�H�EH�H�uM�H�> ��   H��I�4�H��u�H�} �G���1�L��L�u H��H�]I��I��K�>H�H��t8H�I�$H� L� H�EH��L��[]A\A]A^A_�L�u I��K�>H�H��u�H�UI�$L�eI�$H��tH�B1�H�uL��M�$�HE H�UH��f�     H������胯��H��蛯��L�u(�ҭ��H���Z���H��肯��L���z���赭��H���=���H��赯��D  AVM��AUI��ATI��H��UH��SH��L��A�Є�L����   H��A�ք��  L��L��A�ք�L�H�{H�s�KH�S H�C(��   M�M L�M�ML�KM�ML�KM�ML�KM�M L�K M�M(L�K(I�}I�uA�MI�U I�E([]A\M�E A]A^��    L��A�ք�u~L��H��A�ք�L�H�{H�s�KH�S H�C(�w���L�M L�L�ML�KL�ML�KL�ML�KL�M L�K L�M(L�K(L�E H�}H�u�MH�U H�E([]A\A]A^��    �KL�H�{H�sH�S H�C(M�$L�M�L$L�KM�L$L�KM�L$L�KM�L$ L�K M�L$(L�K(I�|$I�t$A�L$I�T$ I�D$([]M�$A\A]A^��    L�H�{H�s�KH�S H�C(����f.�     f�AVM��AUI��ATI��H��UH��SH��L��A�Є�L����   H��A�ք���  L��L��A�ք�L�D�KD�C�{H�sH�K �S(�C,�  M�] L�M�]L�[M�]L�[M�]L�[M�] L�[ M�](L�[(E�ME�EA�}I�uI�M A�U(A�E,[]A\M�U A]A^�L��A�ք���   L��H��A�ք�L�D�KD�C�{H�sH�K �S(�C,�k���L�] L�L�]L�[L�]L�[L�]L�[L�] L�[ L�](L�[(L�U D�MD�E@�}H�uH�M �U(�E,[]A\A]A^�D  D�K�{L�D�CH�sH�K �S(�C,M�$L�M�\$L�[M�\$L�[M�\$L�[M�\$ L�[ M�\$(L�[(E�L$E�D$A�|$I�t$I�L$ A�T$(A�D$,[]M�$A\A]A^�fD  L�D�KD�C�{H�sH�K �S(�C,�����f.�      AWI��AVAUI��ATI�̹   UH��H� SH��H�W�H�w�L�w�̫������   H��������H��H9��6  L�<�    L���N���1�H��L��I���Ψ��H�u1�H�E    L�UH��u�UD  I�H�I�H�2H��t?H��H�F1�H�H��M��M�M��u�H�EH�H�uM�H�> ��   H��I�4�H��u�H�} 臨��1�L��L�u H��H�]I��I��K�>H�H��t8H�I�$H� L� H�EH��L��[]A\A]A^A_�L�u I��K�>H�H��u�H�UI�$L�eI�$H��tH�B1�H�uL��M�$�HE H�UH��f�     H�������é��H���۩��L�u(����H��蚩��H���©��L��躧�������H���}���H�������D  AWH��H��AVAUI��ATI��I��?UL�SH��XH�D$H�|$H�D$H�t$H�T$H�$H9��c  I��� I��I�FL�< I�_�I�<�H�$H�4[H��M�d= H��I�l5 L��H����L���K�vHE�ID�H�H��L�H;\$H�H�VH�PH�VH�PH�VH�PH�V H�P H�V(H�P(�z����D$uL�t$I��L��H��?I�I��I9���   H�[H��L�H��$�   H;\$H�T$ H��$�   H�T$(H��$�   H�T$0H��$�   H�T$8H��$�   H�T$@H��$�   H�T$H��   H�T$ H�H�T$(H�PH�T$0H�PH�T$8H�PH�T$@H�P H�T$HH�P(H��X[]A\A]A^A_�@ H�\K�vH�[H��L�H��L�H�H�
H�HH�JH�HH�JH�HH�JH�H H�J H�H(H�J(����f.�     H�C�I��I��?I�I��O�$vH�t$ H�$I��M�L���Є�H�[uH��L��"���f�I�$H��L�H�I�T$H�PI�T$H�PI�T$H�PI�T$ H�P I�T$(H�P(I�F�H��H��?H�H��L9t$|L�������D  L��I���f���H���#���f.�      AWAVI��AUATI��USH��xH�T$0H��H)�H��/  ��  H�|$0 I���B  H�O0H��H�L$8H��L)�H���������M��H��H��L��H��I�N�H�l$0L��H��H��?H�H��H�@H��I������H��f.�     H��H��L��A�Մ�H�S0u�L�}��     L��L��L��I��0A�Մ�u�H9�vfL�M L�H�{H�s�KH�S L�L�MH�C(L�KL�ML�KL�ML�KL�M L�K L�M(L�K(H�U H�S0L�E H�}H�u�MH�E(�a���H�T$0L��L��H������H��H��L)�H��/  �s  H�|$0 tI��H�\$8�����H��H��H���������I��L��M�w�I��K�,vH��L��I��H�}H�uH��0H�MHH�UPL�E0H�EXH�|$HH�t$PH�L$XH�T$`H�|$H�t$L��H�L$H�T$ L��L�D$@H�D$hL��L�$H�D$(L���+���M��u�H��0I���������I�$L�H��L�CH�{L)�H�sH�K H��0H�S0I�T$H�CXH�|$PH�t$XH�L$`H�S8I�T$L�L$@L�D$HH�D$hH�S@I�T$H�SHI�T$ H�SPI�T$(H�SXH��H�|$H��H�t$H�L$ I��1�L�$L�D$H�D$(L��L���o���H��_�M���H��x[]A\A]A^A_�fD  AWH��H��AVAUI��ATI��I��?UL�SH��XH�D$H�|$H�D$H�t$H�T$H�$H9��c  I��� I��I�FL�< I�_�I�<�H�$H�4[H��M�d= H��I�l5 L��H����L���K�vHE�ID�H�H��L�H;\$H�H�VH�PH�VH�PH�VH�PH�V H�P H�V(H�P(�z����D$uL�t$I��L��H��?I�I��I9���   H�[H��L�H��$�   H;\$H�T$ H��$�   H�T$(H��$�   H�T$0H��$�   H�T$8H��$�   H�T$@H��$�   H�T$H��   H�T$ H�H�T$(H�PH�T$0H�PH�T$8H�PH�T$@H�P H�T$HH�P(H��X[]A\A]A^A_�@ H�\K�vH�[H��L�H��L�H�H�
H�HH�JH�HH�JH�HH�JH�H H�J H�H(H�J(����f.�     H�C�I��I��?I�I��O�$vH�t$ H�$I��M�L���Є�H�[uH��L��"���f�I�$H��L�H�I�T$H�PI�T$H�PI�T$H�PI�T$ H�P I�T$(H�P(I�F�H��H��?H�H��L9t$|L�������D  L��I���f���H���#���f.�      AWAVAUI��ATUH��SH��xH�T$0H��H)�H��/  ��  H�|$0 I���Q  H�O0H��H�L$8I��H)�H���������M��H��L��H��H��I�M�H�l$0L��H��H��?H�H��H�@H��H�T ����L��f�     I��H��H��A�Ԅ�I�W0u�L�s��     L��L��H��I��0A�Ԅ�u�L9�vuL�M�E�OE�GA�I�wM�L�[A�W(I�O A�G,M�_L�[M�_L�[M�_L�[ M�_ L�[(M�_(�S(I�W0L�D�KD�C@�{H�sH�K �C,�R���H�T$0L��L��L������L��L��H)�H��/  �o  H�|$0 tI��L�t$8����I��H��H���������I��L��M�n�I��K�\m H��H��I��H�{H�sH��0H�KHH�SPL�C0H�CXH�|$HH�t$PH�L$XH�T$`H�|$H�t$H��H�L$H�T$ L��L�D$@H�D$hL��L�$H�D$(L������M��u�I�_�I���������H�U L�I��L�CH�{I)�H�sH�K H��0H�S0H�UH�CXH�|$PH�t$XH�L$`H�S8H�UL�L$@L�D$HH�D$hH�S@H�UH�SHH�U H�SPH�U(H�SXL��H�|$H��H�t$H�L$ I��1�L�$L�D$H�D$(L��H���d���I��_�R���H��x[]A\A]A^A_�f.�     �AWH��H��AVI��AUATI��I��?UL�SH��8H�D$H�|$H�D$H�t$ H�T$(H�L$L�D$H9��*  I���D  I��I�OH�D$L�$	H��M�,I�\$�I�} I�,�H�u ��L���HE�ID�H;\$H�K��|��D$(uL�|$(I��L��H��?I�I��I9�t0H;\$ I��<D  H�T$H�H��8[]A\A]A^A_�f�     H�\H;\$ I��H�K��~�H�C�I��I��?I�I��O�$�H�t$H�D$I�<$�Є�uI���f.�     I�$I��I�G�H��H��?H�H��L9|$ |L���h����     L��I���H������AWAVI��AUATUH��SH��H�$H��H)�H���   ��  H�<$ I���e  H�OH��H�L$H)�H�,$H�}H��H��H��?H�H��H�\� H�3A�Ԅ�I�v���   H�;A�Ԅ���   H�E H�H�U H�H�u M��H�T$�	D  H�u H��H�:A�Ԅ�H�Su�M�}��    I�7M��H�} I��A�Ԅ�u�I9�vI�U H�H�I�E H�SH�u �H�$L��L��H�������H��H��H)�H���   ��   H�<$ txI������H�}A�Ԅ�t=H�E H�uH�EH�u �H���I�v�H�}A�Ԅ�H�E t�I�V�H�U I�F�H�u � ���I�v�H�;A�Ԅ�H�E u�H�H�U H�H�u �����H��H��L�j�I��I���I��J�L� M��L��L��H������M��u�H��H�E I��H�I)�1�M��L��H��H��H�CH������I���H��[]A\A]A^A_�f.�     @ AWI��AVAUATUSH��H��H9���   H��I�պ?   H)�L��H��H��H��H��?H�H)�H�����H���   ~dL���   L��H��L������M9�t; L��M�&I�^��@ H�H�E H��H��H�3L��A�Մ�u�I��L�e M9�u�H��[]A\A]A^A_�L��L��H������H��[]A\A]A^A_�f.�     f�AWA��AVI��AUI��ATL�% $  UH�-($  SL)�1�H��H������H��t�     L��L��D��A��H��H9�u�H��[]A\A]A^A_�ff.�     ��f�H��H���                 r Error opening file %s
 RSW_test.txt refFlat.txt RSW_tst %s.results w %s.results.full %s.results.unknown %s.results.unknown.full %s.results.splitPairs Novel * %s	%s	%li	%li	%li--%li	%s /proc/self/status VmRSS: %19s %19s %999s .gz gunzip -c %s .lrz cat %s | ./lrunzip Error reading from file %s
 	 -- %s	%s	%s	%c	%i-%i	%i	%i-%i
 %s,  , %s      Options file %s is more than the max of %i bytes.
      Not enough lines in options file.       refFlat.txt.intronBoundary.exonsgaps    Not enough arguments given, using default values.       Usage is to load options from file: ./splitPairs optionsFile.txt        And make sure options file has options in order, each on their own line with no extra lines.    Running with options...
  file with read data          %s
  max distance between matches %i
  length of samples            %i
  refFlat file                 %s
  refFlat intron boundary file %s
  minimum splice length        %i
  tolerance of difference in position for supporting reads  %i
  base of file name for writing results                     %s
  minimum number of supporting reads                        %i

      Error, filename %s is too long.
        Error opening file %s for writing.
     Will write summary results that match in known genes to file
   %s
     Will write full results that match in known genes to file
   %s
        Will write summary results that do NOT match in known genes to file
   %s
      Will write full results that do NOT match in known genes to file
   %s
 Will write split pairs to file
   %s
   Finished processing data, results written to files.
    Number of entries in data file:             %i
 Number of different reads:                  %i
 Number of entries in refFlat file:          %i
 Number of entries in refFlat boundary file: %i
 Number of matches:                          %i
 String table size:                          %i
 VmRSS, memory resident set size:            %s %s
      Total time to process:                      %i seconds
 Error reading data file %s, line exceeded %i characters.
       Done reading/sorting refFlat, total time elapsed %i seconds
    Done reading refFlat intron/exon boundaries, total time elapsed %i seconds
     Done reading read data, total time elapsed %i seconds
  Done sorting read data, total time elapsed %i seconds
  Done finding matched pairs, total time elapsed %i seconds
      Done sorting matched pairs, total time elapsed %i seconds
      Done computing supporting reads, total time elapsed %i seconds
 Done sorting matched pairs again, total time elapsed %i seconds
        Done filtering matched pairs, total time elapsed %i seconds
    Done calculating/filtering supporting reads, total time elapsed %i seconds
     Line No	Id	Gene	Chr	Strand	Splice region	Supporting reads	Supporting splice range
      GeneName	Chromosome	# supporting reads	splice length	range of supporting reads	Novel or not (*)
        Done saving results, total time elapsed %i seconds
                     `L@     0L@     �L@     �K@     L@     �K@     `K@     `K@     HK@     �P@      Q@     xP@     �P@     �O@     �P@     U@     �T@     0U@     �T@     8T@     8T@     8T@     8T@     8T@     8T@     8T@     �T@     �S@     UNFOUND_                                                                                            ;�  1   ܉���  ̌���  ���<  H����  <���  l����  �����  ����  L���$  ���t  �����  �����  ����  �����  ,���  ̽��D  L���d  ����|  ����t  �����  ����  ����l  ���  ����L  \���  |���4  ����L  ����d  ����|  <����  ����  l����  ����<  ���|  ����  ����$  �����  �����  �����  ����T	  �����	  L����	  <���D
  �����
  �����
  |���4  �����  |���l  �����         zR x�      ����*                  zR x�  $      �����   FJw� ?;*3$"       D   0���.              \   8���              t   @���              �   H���              �   P���           $   �   X���^    A�D�D NAA   �   ����F              �   ����F                0���>           4   ,  H����    B�D�A �D0p
 AABA     L   d   ����    B�E�D �C(�F0T
(A ABFFN
(C ABBB      �  p����              �  (���>              �  P���           <   �  H����   B�J�A �A(�D�%
(A ABBA      <  �����    D{,   T  ����    A�F�GPr
AAA        �  �����   D_
E       �  ���Q    D0L<   �  (���   B�O�G �A(�G��
(A ABBK   <   �  X���A   B�L�A �F(�G@�
(A ABBG     <   <  h���A   B�L�A �F(�G@�
(A ABBG     <   |  x���I   B�L�A �F(�G@�
(A ABBG     <   �  �����    B�B�G �A(�G@�
(A ABBG     D   �  H����    B�B�B �B(�A0�A8�DP�8A0A(B BBB       zPLR xp@ �  L   $   ����  ��@ B�E�B �E(�I0�H8�DP
8D0A(B BBBDT   t   P���	  �@ B�J�B �A(�A0��
(A BBGEd
(A BBBA     L   �   ���  D�@ B�B�B �E(�A0�A8�Dp(
8A0A(B BBBAL     �����  t�@ B�G�E �B(�A0�A8�Dp�
8A0A(B BBBAL   l  X���C  ��@ B�G�B �B(�A0�A8�DpQ
8A0A(B BBBAL   �  8����  Ԕ@ B�E�B �E(�I0�H8�D@�
8D0A(B BBBA <   L  ����    B�F�A �D0�
 FABED FAB 4   �  x����    B�D�D �F0Q
 AABK     d   �  @����   B�E�E �G(�D0��
(A BFBHy
(A BBBHe
(A FBBH     d   ,  ����   B�E�E �G(�D0��
(A BFBA�
(A BBBFv
(A FBBG     L   T  `����  ��@ B�E�B �E(�I0�H8�D@�
8D0A(B BBBA L   �  �����   B�I�B �E(�H0�D8�D�w
8A0A(B BBBE   L   4  P����   B�B�E �B(�D0�A8�D��8A0A(B BBB      L   �  �����   B�I�B �E(�H0�D8�D�w
8A0A(B BBBE   L   �  `����   B�B�B �E(�A0�D8�D��8A0A(B BBB      L   $	  ����   B�I�E �B(�H0�D8�Dp�
8A0A(B BBBJ     L   t	  @���   B�B�E �B(�A0�D8�DP�8A0A(B BBB       \   �	  ����    B�E�B �B(�A0�A8�G@�
8A0A(B BBBAR8A0A(B BBBT   �  ����F(  $�@ B�G�B �B(�A0�F8�K��
8C0A(B BBBA       ,   |
  ب��)   B�F�A ��
ABu  D   �
  ���e    B�E�E �E(�H0�H8�M@l8A0A(B BBB    �
  0���               �%/  ]������ �    }    �%k  ���  �� �          �-"k�  ���
 ��	�y  �	�
 �
          �-"�  ��� ���l  �� �          �-"�  ��� ���n  �� �          �%/  ]������ �    }    �%/  ]������ �    }    ��?�  �-�M ��2  �M�M                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0?@     @<@     ?@                                  E             O             l             �@            z@            ؝`                          �`                   ���o    �@            �@            �@     
       {                                           �`            P                           �@            @@            H       	              ���o    �@     ���o           ���o    D@                                                                                                             ��`                     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@         GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-7) GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-9)  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.gnu.build-id .gnu.hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .gcc_except_table .init_array .fini_array .jcr .dynamic .got .got.plt .data .bss .comment                                                                                8@     8                                    #             T@     T                                     1             t@     t      $                              D   ���o       �@     �      4                             N             �@     �      �                          V             �@     �      {                             ^   ���o       D@     D      j                            k   ���o       �@     �      �                            z             @@     @      H                            �             �@     �      P                          �             �@     �                                    �              @            �                            �             �@     �      $f                             �             z@     z      	                              �              z@      z                                    �             $�@     $�      �                             �             ��@     ��      <                             �             ��@     ��      K                             �             ؝`     ؝                                    �             �`     �                                    �             �`     �                                    �             ��`     ��                                  �             ��`     ��                                   �              �`      �      �                                        ��`     ��                                                ��`     ��      ��                                   0               ��      X                                                   �                                                         ��      H         3                 	                      �                                                                 8@                   T@                   t@                   �@                   �@                   �@                   D@                   �@                  	 @@                  
 �@                   �@                    @                   �@                   z@                    z@                   $�@                   ��@                   ��@                   ؝`                   �`                   �`                   ��`                   ��`                    �`                   ��`                   ��`                                       ��                     �H@     	      =    @W@     �       �    ��@     d       �    @<@     )      �    T@a            �   ��                �    �`             �    �>@             �    �>@                 ?@             "    ȡ`            1    �`             X    0?@             d    ؝`             �   ��                �    �@             �    �`                  ��                �     �`             �     �`             �     ؝`             �    ��`             �     ��`             �    A@     �           0@@     >       .                     B  "  �X@            `    @@a            m  "  0^@     �       �                     �    z@            �    l>@             �  "   r@     �      M    `�`     '      S  "   Y@     ^       �                      �                      �                     �     �`            �    P@a            �                     �                                          9                     K  "  �X@            i    z@             o                     �                     �  "   u@     �      %                     9                     L                     k  "  �X@            �  "  0^@     �       �     �`            �    �A@     >       �                       "  �i@     �      �                     �    ��`            �      @             �    ��`            �                      �     z@            	                     	    �F@           7	    �?a            K	    ��`            P	                     b	                     v	     �`     0       �	                      �	    ��`             �	                     �	  "  0_@     �       I
    ��`     '      Q
    �R@     C      h
                     }
                     �
    �`            �
     @a            �
    `?@     .       �
    �a     '      �
                     �
   ��`             �
                         0B@     �      4                     H  "  �\@     I      �    pF@     Q       �  "  `@           �   (z@             �    �W@     �       �  "  �[@     A      ,                     @    �y@     e       P  "  `Y@     �       �    �?@     F       �  "  �c@     �      5    �D@     �      J    @�`     '      O  "  �X@            y  "   b@     �      �                     �    ��`            �  "  �X@            �    �?@     F       �    �?a                 B@            0  "  �x@     �       �    ��`             �  "  �\@     I      �    �?a            �    `�`     0       	  "  pl@     �      |     O@     �      �                     �                      �                     �  "   Y@     ^                            5     @a            I                     ]    X@a             b                     v    0@a            �  "  �X@            �                     �                                          7    ��`             f                     z    @a            �    p@@     �       �                     �  "  �X@            �    ��`            �    �`            �                         �C@     �       (  "  �g@     �      &  "  �[@     A      k    ��`            w    ��`             ~     J@           �     p@             �                     �  "  �v@           8                     K  "  �e@           �                     �  "  `o@     �      H                     [    ��`            m  "  @Z@     A      �  "  @Z@     A      �    �`            �                         PD@     �       <  "  `Y@     �       y    ��`            �  "  �X@            �                     �    �@     F(      �    �@             �                      splitPairs.cpp _ZNSt10_HashtableISsSsSaISsENSt8__detail9_IdentityESt8equal_toISsESt4hashISsENS1_18_Mod_range_hashingENS1_20_Default_ranged_hashENS1_20_Prime_rehash_policyENS1_17_Hashtable_traitsILb1ELb1ELb1EEEE9_M_insertIRKSsEESt4pairINS1_14_Node_iteratorISsLb1ELb1EEEbEOT_St17integral_constantIbLb1EE.constprop.175 _ZNSt10_HashtableIPKcSt4pairIKS1_P10RSW_spliceESaIS6_ENSt8__detail10_Select1stESt8equal_toIS1_ESt4hashIS1_ENS8_18_Mod_range_hashingENS8_20_Default_ranged_hashENS8_20_Prime_rehash_policyENS8_17_Hashtable_traitsILb0ELb0ELb1EEEE9_M_insertIRKS6_EES2_INS8_14_Node_iteratorIS6_Lb0ELb0EEEbEOT_St17integral_constantIbLb1EE.isra.157 _ZL14unfound_string _GLOBAL__sub_I__Z8get_lineP8_IO_FILEPci _ZStL8__ioinit crtstuff.c __JCR_LIST__ deregister_tm_clones register_tm_clones __do_global_dtors_aux completed.6342 __do_global_dtors_aux_fini_array_entry frame_dummy __frame_dummy_init_array_entry __FRAME_END__ __JCR_END__ _GLOBAL_OFFSET_TABLE_ __init_array_end __init_array_start _DYNAMIC data_start _Z12compute_hashP3RSW _Z24compare_spliceByChromLenP10RSW_spliceS0_ printf@@GLIBC_2.2.5 _ZNSt6vectorI3RSWSaIS0_EED1Ev sampleLength _ZNSt6vectorIP10RSW_spliceSaIS1_EE19_M_emplace_back_auxIJRKS1_EEEvDpOT_ memset@@GLIBC_2.2.5 __libc_csu_fini _start _ZSt16__introsort_loopIN9__gnu_cxx17__normal_iteratorIP3RSWSt6vectorIS2_SaIS2_EEEElPFbRKS2_S9_EEvT_SC_T0_T1_ sLine _ZNSt13unordered_setIPKcSt4hashIS1_ESt8equal_toIS1_ESaIS1_EED2Ev __gmon_start__ _Jv_RegisterClasses puts@@GLIBC_2.2.5 fKnown maxDistance _ZdlPv@@GLIBCXX_3.4 _ZNSs7reserveEm@@GLIBCXX_3.4 _ZNSs4_Rep10_M_disposeERKSaIcE@@GLIBCXX_3.4 exit@@GLIBC_2.2.5 _ZNSt6vectorI3RSWSaIS0_EED2Ev _fini __cxa_rethrow@@CXXABI_1.3 _ZNSt8ios_base4InitC1Ev@@GLIBCXX_3.4 _ZSt13__adjust_heapIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEElS3_PFbS3_S3_EEvT_T0_SC_T1_T2_ malloc@@GLIBC_2.2.5 fopen@@GLIBC_2.2.5 __libc_start_main@@GLIBC_2.2.5 _ZNSt6vectorI9RSW_KnownSaIS0_EED1Ev _ZNSt6vectorIP10RSW_spliceSaIS1_EE19_M_emplace_back_auxIIRKS1_EEEvDpOT_ fKnownFull _Z12compare_dataRK3RSWS1_ _ZNSsC1ERKSs@@GLIBCXX_3.4 _ZSt13__adjust_heapIN9__gnu_cxx17__normal_iteratorIP9RSW_KnownSt6vectorIS2_SaIS2_EEEElS2_PFbRKS2_S9_EEvT_T0_SD_T1_T2_ __cxa_atexit@@GLIBC_2.2.5 data_known _ZNSt8ios_base4InitD1Ev@@GLIBCXX_3.4 beginTime _ITM_deregisterTMCloneTable _IO_stdin_used fputc@@GLIBC_2.2.5 _Z10printStatsP8_IO_FILE refFlatBoundaryFile data free@@GLIBC_2.2.5 strlen@@GLIBC_2.2.5 readIdsReported _ITM_registerTMCloneTable __data_start _ZNSs4_Rep10_M_destroyERKSaIcE@@GLIBCXX_3.4 _ZSt16__insertion_sortIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEEPFbS3_S3_EEvT_SB_T0_ options _Z15read_boundariesPKc sprintf@@GLIBC_2.2.5 fgetc@@GLIBC_2.2.5 fUnknownFull refFlatFile _Z18compare_data_knownRK9RSW_KnownS1_ buff setpriority@@GLIBC_2.2.5 __TMC_END__ _ZNSsC1EPKcRKSaIcE@@GLIBCXX_3.4 _Z19readOptionsFromFilePKc strstr@@GLIBC_2.2.5 _ZNSt6vectorI14RSW_BoundariesSaIS0_EE19_M_emplace_back_auxIIRKS0_EEEvDpOT_ _Z11printSpliceP8_IO_FILEP10RSW_splice _ZNSt10_HashtableISsSsSaISsENSt8__detail9_IdentityESt8equal_toISsESt4hashISsENS1_18_Mod_range_hashingENS1_20_Default_ranged_hashENS1_20_Prime_rehash_policyENS1_17_Hashtable_traitsILb1ELb1ELb1EEEE21_M_insert_unique_nodeEmmPNS1_10_Hash_nodeISsLb1EEE __dso_handle _Z19updateSpliceSupportP10RSW_spliceS0_ _ZNSt6vectorI9RSW_KnownSaIS0_EE19_M_emplace_back_auxIIRKS0_EEEvDpOT_ strtol@@GLIBC_2.2.5 __libc_csu_init _ZNSt13unordered_setISsSt4hashISsESt8equal_toISsESaISsEED1Ev _Z24compare_spliceByChromPosP10RSW_spliceS0_ _ZSt22__move_median_to_firstIN9__gnu_cxx17__normal_iteratorIP9RSW_KnownSt6vectorIS2_SaIS2_EEEEPFbRKS2_S9_EEvT_SC_SC_SC_T0_ _Z15openOutputFilesv temp _ZNSt6vectorI14RSW_BoundariesSaIS0_EED2Ev _ZNSt10_HashtableIPKcSt4pairIKS1_P10RSW_spliceESaIS6_ENSt8__detail10_Select1stESt8equal_toIS1_ESt4hashIS1_ENS8_18_Mod_range_hashingENS8_20_Default_ranged_hashENS8_20_Prime_rehash_policyENS8_17_Hashtable_traitsILb0ELb0ELb1EEEE21_M_insert_unique_nodeEmmPNS8_10_Hash_nodeIS6_Lb0EEE memmove@@GLIBC_2.2.5 endTime _ZNSt6vectorIP10RSW_spliceSaIS1_EED1Ev _Z18compare_dataToSortRK3RSWS1_ resultsBaseName _Z23compare_spliceBySupportP10RSW_spliceS0_ _ZSt4sortIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEEPFbS3_S3_EEvT_SB_T0_ __bss_start _ZNSt6vectorI14RSW_BoundariesSaIS0_EE19_M_emplace_back_auxIJRKS0_EEEvDpOT_ minSupportingReads stringTable _ZSt16__introsort_loopIN9__gnu_cxx17__normal_iteratorIP9RSW_KnownSt6vectorIS2_SaIS2_EEEElPFbRKS2_S9_EEvT_SC_T0_T1_ _Z14read_knownGenePKc strcpy@@GLIBC_2.2.5 __pthread_key_create strtok@@GLIBC_2.2.5 _ZNSt13unordered_setIPKcSt4hashIS1_ESt8equal_toIS1_ESaIS1_EED1Ev _ZSt11_Hash_bytesPKvmm@@CXXABI_1.3.5 supportPosTolerance memcmp@@GLIBC_2.2.5 _end fclose@@GLIBC_2.2.5 minSpliceLength _ZNSt6vectorI14RSW_BoundariesSaIS0_EED1Ev _ZNKSt8__detail20_Prime_rehash_policy11_M_next_bktEm@@GLIBCXX_3.4.18 __cxa_end_catch@@CXXABI_1.3 _ZSt17__throw_bad_allocv@@GLIBCXX_3.4 _ZNSs4_Rep20_S_empty_rep_storageE@@GLIBCXX_3.4 fscanf@@GLIBC_2.2.5 sampleDataFile _Z8get_lineP8_IO_FILEPci __cxa_begin_catch@@CXXABI_1.3 _ZNSt6vectorI9RSW_KnownSaIS0_EED2Ev data_boundaries fUnknown fwrite@@GLIBC_2.2.5 _Z17setDefaultOptionsv _ZNSt10_HashtableIPKcS1_SaIS1_ENSt8__detail9_IdentityESt8equal_toIS1_ESt4hashIS1_ENS3_18_Mod_range_hashingENS3_20_Default_ranged_hashENS3_20_Prime_rehash_policyENS3_17_Hashtable_traitsILb0ELb1ELb1EEEE21_M_insert_unique_nodeEmmPNS3_10_Hash_nodeIS1_Lb0EEE _ZNSt6vectorI9RSW_KnownSaIS0_EE19_M_emplace_back_auxIJRKS0_EEEvDpOT_ data_splice _edata _Z9read_dataPKc __gxx_personality_v0@@CXXABI_1.3 fprintf@@GLIBC_2.2.5 _ZSt16__introsort_loopIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEElPFbS3_S3_EEvT_SB_T0_T1_ _Znwm@@GLIBCXX_3.4 _ZSt22__move_median_to_firstIN9__gnu_cxx17__normal_iteratorIP3RSWSt6vectorIS2_SaIS2_EEEEPFbRKS2_S9_EEvT_SC_SC_SC_T0_ _Unwind_Resume@@GCC_3.0 _ZSt13__adjust_heapIN9__gnu_cxx17__normal_iteratorIP3RSWSt6vectorIS2_SaIS2_EEEElS2_PFbRKS2_S9_EEvT_T0_SD_T1_T2_ popen@@GLIBC_2.2.5 numDifferentReads _ZNSt6vectorI3RSWSaIS0_EE19_M_emplace_back_auxIJRKS0_EEEvDpOT_ _ZNSt6vectorI3RSWSaIS0_EE19_M_emplace_back_auxIIRKS0_EEEvDpOT_ fSplitPairs _ZNSs12_M_leak_hardEv@@GLIBCXX_3.4 _Z19printCurrentOptionsP8_IO_FILE _ZNSt13unordered_setISsSt4hashISsESt8equal_toISsESaISsEED2Ev stdout@@GLIBC_2.2.5 _ZNSt6vectorIP10RSW_spliceSaIS1_EED2Ev time@@GLIBC_2.2.5 main _init _ZNKSt8__detail20_Prime_rehash_policy14_M_need_rehashEmmm@@GLIBCXX_3.4.18                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              split.sh                                                                                            0000775 0001753 0001753 00000005066 12563403132 012370  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #!/bin/bash

BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )
#bring in the configuration, which should be in the script's directory
source "${BASEDIR}/rsr_config.sh"


#inputs: reads_file left_cut reads_length
function n_split() {
    if [ -z "$LOG_FILE" ]; then
        output="/dev/null"
    else
        output=$LOG_FILE
    fi
    $SPLIT_PROGRAM "$1" "$2" "$3" "$SPLIT_TEMP_DIR" >& $output
    log "n_split:: ${SPLIT_TEMP_DIR}/$(basename $1).split"
    echo "${SPLIT_TEMP_DIR}/$(basename $1).split"
}

#this function does not obery $SPLIT_TEMP_DIR, so you get
#your results where your data started.
function o_split() {
    file=$1
    start=$2
    stop=$3
    log "o_split::file=$1,start=$2,stop=$3"
    for ((i = $start;i <= (($stop - $start));i++));do
        log "./split_read_RSW.pl $file $i $(($stop - $i))"
        $SPLITTER "$file" $i $(($stop - $i)) >& /dev/null
    done
    log $(ls ${file}.*)
    if [ -f "${file}.split" ];then
        rm -f "${file}.split"
    fi
    cat $file.split.left* >> "${file}.split"
    rm -f $file.split.left*
    result="${file}.split"
    echo "$result"
}

#input: resultsFileBase minimum-cut-size readsLength 
function split_pairs() {
    base=$1
    start=$2
    len=$3
    if (( $start >= $len / 2 )); then
        start=$(( $len - $start ))
    fi
    log "split_pairs::file=$1,start=$2,len=$3"
    if [ -f "${base}.unmapped.txt" ]; then
        result=$(n_split "${base}.unmapped.txt" $start $len)
    elif [ -f "${base}.unmapped_1.txt" ] && [ -f "${base}.unmapped_2.txt" ]; then
        tmp=$(n_split "${base}.unmapped_1.txt" $start $(( $len - $start ))  )
        mv "$tmp" "${tmp}-1.txt"
        result="${tmp}-1.txt"
        tmp=$(n_split "${base}.unmapped_2.txt" $start $(( $len - $start ))  )
        mv "$tmp" "${tmp}-2.txt"
        result+="|${tmp}-2.txt"
    else
        die "split_pairs::Don't know what do do with $base"
    fi
    log "split_pairs:: $result"
    echo "$result"
}


#-----Main---

if (( $# < 3 )); then
    yell "Usage: $(basename $0) <reads filename base> <min split> <reads length>"
    yell "The reads filename should be the part before the \"unmapped.txt\" or \"unmapped_1.txt\""
    yell "The script will detect whether there is a pair and split both."
    yell "examples: "
    yell "if your file to split is \"Dtt_KO.txt.unmapped.txt\" and you wanted a minimum split of 11bp with a reads length of 35bp, then..."
    yell "         $(basename $0) Dtt_KO.txt 11 35"
    yell "if your files to aplit are \"Hg.txt.unmapped_1.txt\" and \"Hg.unmapped_2.txt\" with min split 33 and read length 101."
else
    echo $(split_pairs $@)
fi
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          sr4                                                                                                 0000775 0001753 0001753 00000153043 12563406232 011337  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                ELF          >    l>@     @        �          @ 8 	 @         @       @ @     @ @     �      �                   8      8@     8@                                          @       @     ?�      ?�                    ؝      ؝`     ؝`     �      ��                    ��      ��`     ��`                                T      T@     T@     D       D              P�td   $�      $�@     $�@     �      �             Q�td                                                  R�td   ؝      ؝`     ؝`     (      (             /lib64/ld-linux-x86-64.so.2          GNU                        GNU xܩC�hO�K�XO�#A>�   1            `1   2   3   )��Y@�(E�LyIk�                                                 �                                                                   �                     >                                          �                      �                     �                      &                     �                     �                                          	                     �                     3                       �                     /                     �                     O                       ~                      }                     �                     #                                          �                     �                     	                     v                     i                       �                     �                      �                     �                     �                     �                     �                      �                     _                     �                                          �                     ]                     �                     �                     �                     %                     �    ��`            q    ��`             �       @             �     p@              libstdc++.so.6 __gmon_start__ _Jv_RegisterClasses _ITM_deregisterTMCloneTable _ITM_registerTMCloneTable __pthread_key_create _ZNSs4_Rep10_M_destroyERKSaIcE _ZSt11_Hash_bytesPKvmm _ZSt17__throw_bad_allocv _ZNSs4_Rep10_M_disposeERKSaIcE __cxa_rethrow _ZNSt8ios_base4InitD1Ev _ZNSsC1EPKcRKSaIcE _ZNKSt8__detail20_Prime_rehash_policy14_M_need_rehashEmmm __cxa_begin_catch _ZNSs4_Rep20_S_empty_rep_storageE __gxx_personality_v0 _ZNSs12_M_leak_hardEv _Znwm _ZNKSt8__detail20_Prime_rehash_policy11_M_next_bktEm __cxa_end_catch _ZNSsC1ERKSs _ZNSs7reserveEm _ZNSt8ios_base4InitC1Ev _ZdlPv libm.so.6 libgcc_s.so.1 _Unwind_Resume libc.so.6 strcpy sprintf fopen puts time popen strtok strtol fgetc strlen __cxa_atexit memset strstr memcmp stdout fputc fclose malloc fscanf fwrite fprintf memmove __libc_start_main setpriority free GCC_3.0 GLIBCXX_3.4.18 CXXABI_1.3.5 CXXABI_1.3 GLIBCXX_3.4 GLIBC_2.2.5                                                                 O         P&y   4              P   h��   <     uѯ   K     ӯk   X     t)�   c        l         ui	   o      ��`                   ��`        2           ��`        1           �`                    �`                   (�`                   0�`                   8�`                   @�`                   H�`                   P�`        	           X�`        
           `�`                   h�`                   p�`                   x�`                   ��`                   ��`                   ��`        3           ��`                   ��`                   ��`                   ��`                   ��`                   ��`                   Ƞ`                   Р`                   ؠ`                   �`                   �`                   �`                   ��`                     �`        !           �`        "           �`        #           �`        $            �`        %           (�`        &           0�`        '           8�`        (           @�`        )           H�`        4           P�`        *           X�`        +           `�`        ,           h�`        -           p�`        .           x�`        /           ��`        0           H��H��  H��t�C   H���              �5�  �%�  @ �%�  h    ������%��  h   ������%�  h   ������%�  h   �����%�  h   �����%ڎ  h   �����%Ҏ  h   �����%ʎ  h   �p����%  h   �`����%��  h	   �P����%��  h
   �@����%��  h   �0����%��  h   � ����%��  h   �����%��  h   � ����%��  h   ������%��  h   ������%z�  h   ������%r�  h   ������%j�  h   �����%b�  h   �����%Z�  h   �����%R�  h   �����%J�  h   �p����%B�  h   �`����%:�  h   �P����%2�  h   �@����%*�  h   �0����%"�  h   � ����%�  h   �����%�  h   � ����%
�  h   ������%�  h    ������%��  h!   ������%�  h"   ������%�  h#   �����%�  h$   �����%ڍ  h%   �����%ҍ  h&   �����%ʍ  h'   �p����%  h(   �`����%��  h)   �P����%��  h*   �@����%��  h+   �0����%��  h,   � ����%��  h-   ����AW�   AVAUATUH��1�S��1�H���   �\���1�������H�s�  �#$  H�}��-  �0  H�=��  �0  H�=ԍ  ��/  H�=��  ��/  H�=��  ��/  H�=��  ��/  H�=d�  ��/  H�=x+! �:  H�D�  L�=5�  L9�H�D$8�`  H��I��H���������L)�A�?   H��H��H��H��H��?H�I)�M�H��/  ��%  M����  I�G0H�t$8H�D$XI��M)�I�O@I��I��L��L��H��?I�I��K�dH��L�H�P9���  ��  H�~�9���   ��  9��"  f���  L�M�M�OM�GA�I�O M�L�XI�W(M�_L�XM�_L�XM�_L�X M�_ L�X(M�_(L�L�HL�@@�xH�H H�P(I��M���   f�I�@�;PI����   ��   L9���   H�M�$I��0M�T$�E�L$�I�|$�I�T$�I�\$�H�XI�\$�H�XI�\$�H�XI�\$�H�X I�\$�H�X(I�\$�L�L�PH�HD�HH�x H�P(I�W@ I�L$L��9�x�^���I�G I9D$ �O���I��0��f�     H�x I9 �H���H��0�-���f�     �`?@ L��L���U  L��L)�H=/  �f  M����  L���!���H�x I9P�P���H�~�9��i�����   9�xyf�     tdM�W0M�M�GI�A�OI�W M�M�W8I�G(M�O0M�G8M�WM�W@I�@M�WM�WHA�OHM�WM�WPI�WPM�W M�WXI�GXM�W(�M���H�F�I9GP}�L�V�M�M�GI�A�OI�W M�L�V�I�G(M�WL�V�M�WL�V�M�WL�V�M�W L�V�M�W(L�N�L�F�H�~��N�H�V�H�F������H�V�H9P �u�������H�N�I9OP�����M���1������H+��  �@�@ H��1������H�=��  �.  H�=�'! �:  1�����H+L�  ���@ H��1�����H�=^�  �i.  H�=�'! �1  1��V���H+�  �Ђ@ H��1�����H�=)�  �4.  H�U�  L�=F�  L9�H�D$H�h  H��H��I���������L)��?   H��H��I��H��H��?H�H)�H�H��/  ��  H���4  I�G0H�t$HL�t$0H�D$@L)�I�O0H��H��H�\$0H��H��?H�H��H�[H��L�H�9��(  ��  H�~�9��8  f��D  9���  f���  L� I�E�_E�WE�OM�GM�'L�`I� A�O(A�W,M�gL�`M�gL�`M�gL�` M�g L�`(M�g(H�D�XD�PD�HL�@H�x �H(�P,H��H�\$@�   �H�G�;H���  ��   H9��  L�0D�kH��0D�c�D�[�L�S�L�K�L�s�L�pD�C��K�L�s�L�pL�s�L�pL�s�L�p L�s�L�p(L�s�H�D�hD�`D�XL�PL�H D�@(�H,I��H�I��9�x6�`���A�G8C|&�P���A�G9Cx�A���I�G H9C �3��� H��0�f.�     D�PE8W|%�$���D�XE9_x����L�P M9W �����H��0������    ��?@ H��H��� W  H��L)�H=/  ��  H���  H��������xA8@|'������xA9Hx@ �����H�x I9P�����H�~�9�������b  9���   f���  I�_0M�E�WE�OE�GI�I�I�_8I�O A�W(A�G,M�_0I�_I�_@E�W8E�O<E�G@I�_I�_HI�HI�_I�_PI�OPI�_ I�_XA�WXA�G\I�_(������~�@8x|"�[����~�9xx�M���H�~�H9x �?���H�^�M�E�WE�OE�GI�I�H�^�I�O A�W(A�G,I�_H�^�I�_H�^�I�_H�^�I�_ H�^�I�_(L�^�D�V�D�N�D�F�H�~�H�N��V��F��8����N�A8O@�����������N�A9OH����������H�N�I9OP���������� �V�8P�����������V�9P�o����y���H�V�H9P �[����f���1�����H+j�  ��@ H��1������H�=|�  �)  H���  H+��  H�����������      H��H���,
  �x���H�D$00   1ɉD$XH�T$hH�|$PfD  ��tH�E�  H�|$0H�T8�H9T8�t���  H;L$P��	  H�|$0L��  H�o�L�I�4)H�H;�  D�fD:g�   H�VH;W��  D�D$XH�\$0A)�H��H�L$HL�L�<II���L�    L9���  L���  H��0I�4)I�<H�H;��  D�fD:g��  H�VH;W��  L�G H�N D��!! M��I)�Mc�M9��r  D�^D:_t��vw;5�!! u�A��L�_  A��RH��H���g  I�|9 E�d1 Ed1��H�|$8Lc�D)�A9��B���;54!! �6���H�=��  L���  H�|$@L)T$@H�|$@H�|$@H��tKL�_Mc�   E1��D  H�wL9�t,I��H��K�4�L9n(u�L9v0u�H;u�H9Vu�L9L$@�����H�=T�  L�U�  I)�I��L�L$hM����  I���   E1��f�     H�FI��L9���   H��K�Dm H��H�H9Pu�L�P L9�|�M9�|�H�@(H9��I9���D$`   �ZD  H�AH�D$H�    H�L$HH�D$00����A��+tH��H�������    A��-�����H��H�������D$`    ��   �1���H�xhH�@P    H�@X    H�@`    �@h  �?�
   H�@p    I���y���H��������I�FPH9��  H��H��H�D$@�����H�T$@1�H��H�D$@�M����D$`H�L$@L��$�   ��I�NH��   �   ���@ I�FH��H��  A�VMc�I�~HH��$�   H�H�HI�I�N�@M�f8M�f(A�F HcD$8I�F@I�F0H��$�   L��$�   ��5  H�6�  H;7�  �1  H��H��$�   �   H�H��  H��H��  ����K�Dm 1�H��H$�  H� �F���L�D$XH�L$P��uH��H��H��?H�H��H9���  H����  f.�     I��0H�|$0_�|$8L� D�hD�`D�XH�hH�X �x(D�p,��  M�W`M��   H�l$@M���Af�     D�L$0H�H��0D�hD�X@�pH�xL�@ D�P(D�H,L��I��0I9���  H�] A;H����   ��   D�M,D�mH�U��uD�]H�}L�E D�U(D�L$0@ ;H��x1�{���@8r%D  �j���;zxD  �Z���L9B �P���L�
H��0L�L�J8L�HL�J@L�HL�JHL�HL�JPL�H L�JXL�H(H���A�@8}|#�V���A�9}x�G���I� H9} �9����E,H���������D�]D�UD�ML�EH�M �D$0H��D�m(L)�H��H��H��tGH��L��L��H��H)�H�L$`L�D$XD�L$PD�T$@D�\$8�����H�L$`L�D$XD�L$PD�T$@D�\$8�D$0H��0I�E�_E�WE�OM�GA�G,L��I��0I9�I�O E�o(�g���L��H9L$HL�t$Hu2����� H��0L� �XI9ΉhD�HL�PL�X D�h(D�`,������YD�IH�Q�L��iH��L�QL�Y D�i(D�a,�     D;H��x u�D8J|�D;Rxu�L;Z f�}�fD  H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���L��1�H�L$@��uH��H��H��?H�H��H9���  H����  I��0H�|$0_L�(L�`L�@@�hL�X H�X(��  M��   M�g`H�l$X�1fD  L�L�PH��0H�XD�HL�@ H�p(L��I��0I9��  H�]A;_H����  tiD�ML�] H�U�L�UL�E H�u(f�     ;ZH��xu�L;B }�H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���L�u M;w }�H��H���������L�U L)�L�MD�EH��H�M(H��H��t=H��L��L��H��H)�H�L$PD�D$HL�L$@L�T$0�X���H�L$PD�D$HL�L$@L�T$0L��H��0I��0I9�M�M�OI�_E�GM�w I�O(�����L9l$8H�l$8u1�a����     I��0L�H�ZL9�H�zD�ZL�B L�J(�5���E�]M�U I�E�I�]I�}L��M�E M�M(@ ;xH��xu�L9@ ~�H�0H��0H�2H�p8H�rH�p@H�rH�pHH�rH�pPH�r H�pXH�r(H���L�u ����H�~�H9x �*��������F�A8G@������'����F�A9GH�����@ ����H�F�I9GP�����������    1������1��    ����H+��  �@�@ H��1��E���H�=�y  ��  H�5��  H�=��  ��?@ ��P  1������H+��  ���@ H��1�����H�=�y  �  H�=y�  H�j�  H��H)�H��H����  H�D$0   L�T$01��    H��L��H)�H��H��H9��d  N�<�    H��L�$�    J�l9�H�UH9S�=  M���   �    H;W8}H�W8H�P0H;W@~H�W@H�H��$�   H��HA��H��$�   H��$�   �`.  A����   H���  H�=��  H��H)�H��L9���   J�,9J�!I��I��H�UH9S��   L�K(H�C0L�E(L)�H��H��?H1�H��H�E0H)�L)�H��H��?H1�H)�H9�ucH�EH9CuYL��Hc5�! L)�H�H1�H)�H9��o���H�E H9�b���E1�E���  H��H��H�` H�P(�����H�W8H�P0�����H��H)�H��H�T$0H9�H�rsH�t$0I���g���1��   �����H+��  ���@ H��1�����H�=�w  �  H�5z�  H�=k�  �0@@ �N  1�����H+[�  � �@ H��1������H�=mw  �x  1��q���H+2�  �H�@ H��1�����H�=Dw  �O  1��H���1�H���  �:���H+��  ���@ H��1��d���H�=w  �  H�=w  �  H�=Uw  �   H�=)w  ��  H�=-w  ��  H�=w  ��  H�=�v  ��  H��v  �R   �   �؄@ �E���H�n�  H�o�  1�H)�H��H����   D  H��H�=�v  �]{@ H�P@D�H L�@H�HH�T$ H�P8H�T$H�P`H�T$H�P0H�T$H�P(H�$H�1������H���  H��H�XXH��t+f�     H�CH�=%v  �y{@ H�1�����H�H��u�H�5v  �
   �-���H���  H���  H�MH)�H��H9�vH��H���4���H�v  �`   �   �0�@ A�   �3���H��u  �`   �   �0�@ ����H��u  �`   �   �0�@ �����H��u  �`   �   �0�@ �����L��  L��  1�L��L)�H��H����  �    M�4�L�u�  1�L�,�    I�H��I��H�S�  H�4�H��H���l  L�I�zL��f�     H9��\  H�	H��tH�y1�H��I��H9�t�Hc�! I;F`�2  I�BH9�t/I�
H���"  H�q1�H��I��H9��  H��L��I��H9�u�H�> ��  L�D�  H�5E�  H��  L)�N�(H��H��m۶m۶mH��A�@yH����   M�H0Hct! L�^M+H(�   1�� H�AL9�t_H��H��H�<�    H��H��H)�L�L9H u�I�x8H+x(H��H��?H1�H)�H9�|�I�x@H+x0H��H��?H1�H)�H9�|�H9�vA�@y A�x �=
  H�-�s  H��s  L��H���$  H��
   ����H�0�  H��J�4(�  H�޿	   ����H��  J�(L�pXM����   I�F�~{@ H��H�1������M�~L�v�  1�I�/H��I��H�\�  H��I��H��t.H�L�Af�     L9�tIH�	H��tL�A1�L��I��I9�t�   ����H��tH�     I�H�PH��H��L�� �` �8  M�6M���W���H�޿
   �����L�B�  L�C�  L��I�T$L)�H��L9�vL��I���4���H�=�r  ����H�=yr  �t���H�=}r  �h���H�=Qr  �\���H�=5r  �P���1�����H+��  ���@ H��1��C���H�=�q  ��  H���  H���  H��H)�H��H����   fD  H�h�H��teH�}XH��u	�@ H��H��-���H��u�H�EPH�}H1�H��    �����H�}HH�E`    H�EX    �����H�������H�)�  H�*�  H��H��H��  H)�H��H���t���H���   1�[]A\A]A^A_�H��$�   ���` �~-  ����I�G0H�\$8H�D$XL��L)�H��_�/  H��H���������H�\$`H��H�x�H�|$0H�|$0H��L�t$0H����H�D$@H�|$@H�|$PIk�0L�I��H�D$hK�D6H�D$xHk�0L�H�D$pL�l$@L��M� I�hM�HA�XM�P M9�M�X(�  L��L�t$H�WuL�p L;w HL�HL�H�IH�0H��L�L9�H�1H�pH�qH�pH�qH�pH�qH�p H�q H�p(H�q(�|  H��H�AH� H�r�H��H�<vH��L�H��L�D�wD9p�y���H��H���H�\$`H��L)�H��_�$���L�S��     I�M�*H���������M�bM�BA�jM�Z I�I�GI�Z(I�BI�GI�BI�GI�BI�G I�B I�G(I�B(L��L)�H�D$0H��H��H�G�H�|$@I��I��?I�I��M���=���1��af.�     uL�p L;w HL�HL�H�IH�0H��L�L9�H�1H�pH�qH�pH�qH�pH�qH�p H�q H�p(H�q(�����H��L�qK�6H�r�H�vH��I�<J��H��D�wL�D9p�x���H��H���L�t$HH�|$P uGH9T$0u@H�|$pH�L$hH�T$xH�H�H�GH�AH�GH�AH�GH�AH�G H�A H�G(H�A(H��I9���   I��0M��L� H�hL�H�XL�P L�X(�?���I���t���H�TH�RH��L�H�1H�0H�qH�pH�qH�pH�qH�pH�q H�p H�I(H�H(H�J�H��H�IH��L�D9@xtH�RH��L������L9X �`  H�RH�0H��L�H�2H�pH�rH�pH�rH�pH�rH�p H�r H�p(H�r(H�Q�H��H��?H�H��H�������H��H���u���H�B�H��H��?H�H��H�IH��L�D9HxtH�RH��L������L9P ��  H�RH�0H��L�H�2H�pH�rH�pH�rH�pH�rH�p H�r H�p(H�r(H�Q�H��H��?H�H��I9��s���H��H���u���I�G0L�d$HH�D$@L)�H��_��   H��H���������I��L��I�m�H��Hk�0L�� H��L�H��0L��$�   H�{8H��$�   H�s@H��$�   H�KHH��$�   H�SPH��$�   H�CXH�|$H�t$L��H�L$H�T$ ��?@ H��$�   L�$L��H�D$(H���i9  H���x���L��L)�H��_����M�T$� A�B(M�H���������E�jE�bE�ZI�j�D$8I�I�Z E�r,I�I�GI�BI�GI�BI�GI�BI�G I�B I�G(I�B(L��L)�H�D$0H��H��H�G�H�|$PI��I��?I�I��M����  1�L�D$X�EH�IH�0H��L�I9�H�1H�pH�qH�pH�qH�pH�qH�p H�q H�p(H�q(�����H��H�AH� H�r�H��H�<vH��L�H��L�D�D9 x9u�D�GD8@|,�D�GD9@x �r���L�@ L;G HL�HL��]���D  H��H���M���H�TH�RH��L�H�1H�0H�qH�pH�qH�pH�qH�pH�q H�p H�I(H�H(H�J�H��H�IH��L�D9 x.tH�RH��L��&���D8X|�9hxu�H9X  �D  H�RH�0H��L�H�2H�pH�rH�pH�rH�pH�rH�p H�r H�p(H�r(H�Q�H��H��?H�H��H�������H��H���d���L��1��k����m  �����E1��D$`    �����L�������H�-�i  H��i  ����I�o0H9l$H�����M�g`L�t$0�8 H�D�pD�PD�@H�pL�H D�X(D�h,L��I��0H��0H9D$H�K���H�] A;H����   ��   D�uD�EH�U�D�UH�uL�M D�](D�m, ;H��x)u�D8B!��z���;rxD  �j���L9J �`���H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���A�@8}|#�c���A�9}x�T���I� H9} �F���H��H�T$0D�]L)�D�UD�MH��L�EH�M H��D�u(D�m,H��tGH��L��L��H��H)�H�L$`L�D$XD�L$PD�T$@D�\$8����H�L$`L�D$XD�L$PD�T$@D�\$8I�E�_L��E�WE�OM�GI�O E�w(E�o,�y���Hc�! I;F`������   �,���H��tH�     I�H�PH��H��H�� �` �V-  �����Hk�0L�����Hk�0L��j���H��L������H��������q���H��H�������M�g0L9d$8�0���M�o`H�l$0�1@ L�H�pH�XD�HL�@ L�X(L��I��0I��0H9D$8�����I�\$A;_L����   tnE�L$M�$I�T$�I�t$M�D$ M�\$(f�     ;ZH��xu�L9B ~�H�:H��0H�8H�z8H�xH�z@H�xH�zHH�xH�zPH�x H�zXH�x(H���M�t$ M;w }�L��H�T$0M�$L)�M�L$E�D$H��I�l$(H��H��t3H��L��L��H��H)�D�D$PL�L$HL�T$@����D�D$PL�L$HL�T$@M�M�OL��I�_E�GM�w I�o(�����M�t$ �x���Hk�0L�����f.�     AT�T@a USH���������H����(z@ �T@a � @ �����
   ���` H��      H��      H��      �ߴ    �?H�ܴ      �_���H9�H���  ��  L�$�    L������1�L��H��H���?����(z@ �`�` �`Y@ H�-i�  ������
   �@�` H��      H��      H��      ��    �?H��      �����H9�H�ݳ  �  H�,�    H���/���H��1�H��H�������(z@ � �` � Y@ H���  �d����(z@ ���` ��X@ H�J�      H�G�      H�D�      �/����(z@ ���` ��X@ H���      H��      H��      ������(z@ ���` ��X@ H���      H���      H���      �����[]A\�(z@ ���` ��X@ H�G�      H�D�      H�A�      ��������� 1�I��^H��H���PTI��z@ H���y@ H���@ �;����f��     ���` UH-��` H��H��w]ø    H��t�]���` ���    ���` UH-��` H��H��H��H��?H�H��u]ú    H��t�]H�ƿ��` ���    �=�b   uUH���~���]��b  ��@ H�=�^   t�    H��tU��` H����]�{��� �s��� �F9Gx �    t	���    H�F H9G ���@ �   �f��9x:�    t�Ð�V�   8W|�    �N9Oxu�H�F H9G ���D  �   �f.�     �F9Gx8�    t	���    H�N(�W0L�G(�D)�;V0xu�I9����f.�     �   �f.�     �F9Gx0�    t	���    H�W0V(+W(;V0xu�H�G`H9F`��� �   �f�AUI��ATI��U��S1�H���f.�     9�}HcӃ�A�D ��
tL���������u߅�~3Hc�I�L��9
t>A�D  ���t%1�9���H��[]A\�D �A]�D  ���A�E  u�H��1�[]A\A]Ð� ��ff.�     H�7�����   H���   f�     ��H��Hc�Hi� ���H�� ���)ȉ���)����V�)Ȅ�u�H�w���t=H��D  ��H��Hc�Hi� ���H�� ���)ȉ���)����V�)Ȅ�u��W��Hc�Hi� ���H�� ���)ȉ���)���)�ø   넋9x2�    t�Ð�V�   8W|�    �H�GH+F����    �   �f�H�G`H9F`���@ AU�0z@ I��ATUSH��X�y���H��I���T  H�$��` ���` 1��@ �C�H��L����������t4H���a ��   ��
uڃ��C� ����   Hc�H����f�     L���H�������   H�$H�|$�
   1�H�8�  �����H�|$�
   1��a�  �����F�  H�D$�
   H�|$(1�H���  H�D$ H���  �}���H�|$0�
   1����  �f���H�|$@���  �
   H�D$81�H���  �C����m�  H��X[]A\A]ÿ�{@ �'  L��1�����1������L���o����+�����{@ ����1�������2z@ L��1��Z���1������ H���|@ H�<�  Iz@ �r�  @�  �X�  !   H��  Vz@ H���  �{@ �(�     ��     H���  bz@ ���     �����H|@ �
�����|@ H�������fff.�     U��|@ SH��H��8�l�  L���  D���  ���  H���  �D$ H�V�  H�D$���  �D$���  �D$H�F�  H�$1������H�-T�  H���l���H=�&  wH��8[]�H��H�꾘~@ 1�����1�����@ H��H���  �jz@ ��a 1��D����uz@ ��a ����H��H��\  �-  ��a ��~@ 1������H���  �wz@ ��a 1�������uz@ ��a �]���H��H��\  ��   ��a �0@ 1�����H�U�  ��z@ ��a 1������uz@ ��a ����H��H�[\  ��   ��a �x@ 1��D���H��  ��z@ ��a 1��l����uz@ ��a �����H��H��[  tY��a ��@ 1�� ���H���  ��z@ ��a 1��(����uz@ ��a ����H��H��[  t��a ��@ 1�H��������~@ ��a 1�����1�����H��(�~y ��z@ ��z@ L�N0L�F`H�NHD�L+N(H�VH�D$H�F@H�D$H�F8��z@ H�$1������H��(�ffffff.�     AU�0z@ I����z@ ATA�   USH��8  �����H�\$@H��H�$    H�D$    �D$    H�D$     H�D$(    �D$0    �fD  ��z@ H��L����u  1�H�ھ{@ H���������t�L�d$ H��H��H�������������1��M���L��4   �   �8�@ H���  ����H�P�  H+A�  L��p�@ 1�H��H���������  L�ﾠ�@ 1�����H���  H+�  L��Ѐ@ 1�H��H������H���  H+��  H��m۶m۶mL�� �@ H��H��H��1��N���H�_�  H+P�  L��0�@ 1�H���-���H��  L��`�@ 1�����L��H��L�ﾐ�@ 1�����H��  H+�  L��ȁ@ 1������L��
   �h���H��8  []A\A]�f.�     L�d$ H��{@ H��1�H��L���c�������fffff.�     AVI���i�AUATUSH�?H�w������H�-U�  1�I��H��H�>�  H��I��H��t(H�H�K I9�tkH�H��tH�K1�H��H��I9�t�   �0���H��H��tH�xH�     L���e���H�C    L��H��L��`�` �z  []A\A]�   A^�@ I�>H�sH�W�H;V�u��2������w���H��1�[]A\A]A^�H���t���H���l�������H���/���H�������    AWAVAUI��ATUSH��8������~ZHcп{@ �   I�t���H  ��t;I�t��${@ �   �u(L��){@ �@�` 1�������0z@ �@�` �@���H����0z@ L���>���H��H���  A�    @ �0   E1������H���%fD  A��'  Ic�A����`�` ��
��  H���������u�E����  �t~   E1�X{@ �`�` �����H���l  E1�fD  H�T$H�|$ H���^���H�D$ H�p�H9p��k  A��w%D���$���@  H�|$ �
   1��G����C(@ 1��X{@ A���`���H���w  H�T$ H�z�H����` t�M��H�O��  ����������g���H�t$H�D$����H�D$�N���fD  H�D$ �P���xH�|$ �����H�D$ � �C�r���f�H�|$ �
   1�����H�C �V���fD  H�|$ �����H�@H�C�9���f�     H�D$ �H���xH�|$ �z���H�D$ � �C�
���f.�     H�|$ ����H�@H������f.�     H�|$ �
   1������C������    H�|$ 1����������    IcԀ�_�` 
�|  Ƃ`�`  ����  A��'  �z  A�   �����H�D$ H�x�H����` ��  A���c  H�3�   �H����t5D  ��H��HcșHi� �H�� ���)щ���)�)Љ��F���u�H�s���t=H��f�     ��H��HcșHi� �H�� ���)щ���)�)Љ��F���u��C��Hc���Hi� ���H�� ���)Љ���)�H�)�  )�H;(�  �K,tsH��t.H�H�H�SH�PH�SH�PH�SH�PH�S H�P H�S(H�P(H��0H�ݢ  H���-���A���k���H�������H��8[]A\A]A^A_�Ƃ_�`  ����H�޿��` ��  뼺'  L�� �@ 1�����H��������L��{@ ������J��q��r��������<{@ L��1��g���1������M��H�Wt(����������<���H�t$�l����-���E1������P��J��H�����H���Z����%����5���H��H�D$ H�t$H�x��Z���H������f�AW�0z@ AVI��AUATUSH��8����H��I���   A�     �0   1��4���H���$�    ��'  HcӃ���`�` ��
��  L����������uօ���  ��y   E1��X{@ �`�` ����H����  1�f�H�T$H�|$ H������H�D$ H�H�H9H��;  ��w&���$�(�@ D  H�|$ �
   1�����H�E  1��X{@ ������H���H  H�T$ H�z�H����` t�M��H�O�  ����������h���H�t$H�D$�����H�D$�O����    H�|$ �n���H�@H�E�@ H�D$ �P���xH�|$ ����H�D$ � �E�Z���f.�     H�|$ �&���H�@H�E �9���f�     H�|$ �
   1�����H�E(����fD  H�|$ �����H�@H�E�����f�     H�|$ 1��4��������    HcӀ�_�` 
��   Ƃ`�`  ����N  ��'  ��   A�   �5����H�D$ H�x�H����` ��   ��
��   H�E H�U(H9�~
H�H�U H�E(H�+�  H;,�  ttH��t/H�U H�H�UH�PH�UH�PH�UH�PH�U H�P H�U(H�P(H��0H��  H���S���A���I���L�������H��8[]A\A]A^A_�Ƃ_�`  �%���H����` �Y	  뼺'  L��� �@ 1��þ��H�������밋J��q��r���������<{@ L��1�蚾��1�����M��H�Wt(��������������H�t$蟿�������E1�������P��J��H�����H�������X�������H��H�D$ H�t$H�x�荾��H������D  AW�0z@ AVAUATUSH��8H�<$贾��H��I����  �     �8   1��d���H���$�    ��'  HcӃ���`�` ��
�+  L���������uօ��  ��u   E1��X{@ �`�` �M���H���  1�f�H�T$H�|$ H���޾��H�D$ H�H�H9H���  ��w~���$�X�@ D  L�|$ I�G�H�x�ֽ��L��I��H���ؾ���Z{@ L��蛾��H��I���g  �
   1�L���  荾��I�H�E(�
   1��y���H�E0L�������D  1��X{@ ��艾��H���p  H�T$ H�z�H����` �*����    H�OH���!  �������������H�t$H�D$訽��H�D$�����fD  H�D$ �P���xH�|$ ����H�D$ � �E�j���f�H�|$ �
   1�追��H�H�E �L���@ H�|$ �����H�@H�E�1���f�     H�|$ �����H�@H�E ����f�     H�|$ ����H�@H�E�����f�     H�|$ 1������T����    H�E0    H�E(    ���� HcӀ�_�` 
��   Ƃ`�`  ����N  ��'  ��   A�   ����f�     H�D$ H�x�H����` ��   ����   H�Қ  H;Ӛ  t|H��t7H�U H�H�UH�PH�UH�PH�UH�PH�U H�P H�U(H�P(H�U0H�P0H��8H���  H������A�������L��谼��H��8[]A\A]A^A_�Ƃ_�`  �,���H��` �h  �H�4$�'  � �@ 1�聺��H��蹺��믋J��q��r��������H�4$�<{@ 1��W���1�������    H�WH��t(��������������H�t$�W��������E1��~����P��J��H�����H���E�����������H��H�D$ H�t$H�x��E���H���m����fff.�     AT1�I��USH��H�L�OH��I��H�H��H��H��t%H�L�AL9�tkH�	H��tL�A1�L��I��H9�t�   H�t$����H��tH�t$H�     L�L�VL�HL�PH��H��L��H���R
  H���   []A\�@ H��H��1�[]A\�f�ATI��UH��S1�H����tfL��H��H�x` H�J(tKH;H8}H�H8H�J0H;H@~H�H@H�
H�xHH���H�T$H�$�������u�H��[]A\�f.�     H�H8H�J0�H��L���f.�     @ H�?H��t�ø�� ��f.�     @ H�?H��t飸�� ��f.�     @ H�?H��t郸�� ��f.�     @ H�?H��t�c��� ��f.�     @ UH��SH��H�H��u	�@ H��H��-���H��u�H�EH�} 1�H��    ����H�E    H�E    H�} H��[]����f�ATI��USH��H�_H��t0�    H��u�zH��H�CH�+H�x�H����` uIH��豷��H��u�I�D$I�<$1�H��    �d���I�<$I�D$    I�D$    �y���H��[]A\ù������H����H�t$�H����H�CH�+H�x�H����` uH���:���H��t�H���ڋP��J��҉H��H�t$�	������    AUH���������ATU�0   SH��H��H�WH+H��H��H����   H��H�t$����H�;H�SI��L��H���������H)�H��I��H��L��H���   H�t$H�H�H�VH�QH�VH�QH�VH�QH�V H�Q H�V(H�Q(H�;H�SH)�H��H��H��M�,H��tH��L��蠷��H�;I��0H��t�/���L�L�#L�kH�kH��[]A\A]�fD  H�H9�vH���������� L��E1��Hk�0H�UUUUUUUH������H9�HF������f.�     D  AUH���������ATU�0   SH��H��H�WH+H��H��H����   H��H�t$迷��H�;H�SI��L��H���������H)�H��I��H��L��H���   H�t$H�H�H�VH�QH�VH�QH�VH�QH�V H�Q H�V(H�Q(H�;H�SH)�H��H��H��M�,H��tH��L���P���H�;I��0H��t�ߴ��L�L�#L�kH�kH��[]A\A]�fD  H�H9�vH���������� L��E1��Hk�0H�UUUUUUUH������H9�HF������f.�     D  AUH��m۶m۶mATU�8   SH��H��H�WH+H��H��H����   H��H�t$�o���H�;H�SI��L��H��m۶m۶mH)�H��I��H��L��H���   H�t$H�H�H�VH�QH�VH�QH�VH�QH�V H�Q H�V(H�Q(H�V0H�Q0H�;H�SH)�H��H��H��M�,H��tH��L�������H�;I��8H��t至��L�L�#L�kH�kH��[]A\A]�fD  H�H9�vH���������� L��E1��Hk�8H��$I�$I�H������H9�HF�������    AUATA�   USH��H��H�WH+H��H����   L��H�t$�,���L�H�KH��H�t$L)�H��H�>H��    H�tvH�8H�;H�KH)�H��H��    L�l H��tH��H�������H�;I��H��t臲��H�+L�L�kH�kH��[]A\A]�fD  H�H9�vI�������Z��� L��E1��H��������H��I������H9�LF��.����    AWAVAUATUSH��H9�H�<$H�t$��   L�oI9���   H��L�w�Ff.�     L��H+$L��H��H��tH�4$H��L��H)�����I��I��H9\$H�$L� tRH�$I�} L��H�0�Մ�M�e u�M�}��fD  I�H�L��I��I�7L���Մ�u�L�#I��L��I��H9\$u�H��[]A\A]A^A_� AWI��AVAUI��ATI�̹   UH��H� SH��H�W�H�w�L�w蜳������   H��������H��H9��.  L�<�    L������1�H��L��I��螰��H�u1�H�E    L�UH��u�UD  I�H�I�H�2H��t?H��H�F1�H�H��M��M�M��u�H�EH�H�uM�H�> ��   H��I�4�H��u�H�} �W���L��1�H�]H��L�u I��I��M�l$K�>H�H��t/H�I�$H� L� H�EH��L��[]A\A]A^A_� L�u �f�H�UI�$L�eI�$H��tH�B1�H�uL��M�$�HE H�UH��D  H������蛱��H��賱��L�u(����H���r���H��蚱��I�D$H�x�H����` uL��耯��軯���    H�WH��t'��������H�t$�H�����H������H��薱���P��J��H�����f�     AWI��AVAUI��ATI�̹   UH��H� SH��H�W�H�w�L�w茱������   H��������H��H9��6  L�<�    L������1�H��L��I��莮��H�u1�H�E    L�UH��u�UD  I�H�I�H�2H��t?H��H�F1�H�H��M��M�M��u�H�EH�H�uM�H�> ��   H��I�4�H��u�H�} �G���1�L��L�u H��H�]I��I��K�>H�H��t8H�I�$H� L� H�EH��L��[]A\A]A^A_�L�u I��K�>H�H��u�H�UI�$L�eI�$H��tH�B1�H�uL��M�$�HE H�UH��f�     H������胯��H��蛯��L�u(�ҭ��H���Z���H��肯��L���z���赭��H���=���H��赯��D  AVM��AUI��ATI��H��UH��SH��L��A�Є�L����   H��A�ք��  L��L��A�ք�L�H�{H�s�KH�S H�C(��   M�M L�M�ML�KM�ML�KM�ML�KM�M L�K M�M(L�K(I�}I�uA�MI�U I�E([]A\M�E A]A^��    L��A�ք�u~L��H��A�ք�L�H�{H�s�KH�S H�C(�w���L�M L�L�ML�KL�ML�KL�ML�KL�M L�K L�M(L�K(L�E H�}H�u�MH�U H�E([]A\A]A^��    �KL�H�{H�sH�S H�C(M�$L�M�L$L�KM�L$L�KM�L$L�KM�L$ L�K M�L$(L�K(I�|$I�t$A�L$I�T$ I�D$([]M�$A\A]A^��    L�H�{H�s�KH�S H�C(����f.�     f�AVM��AUI��ATI��H��UH��SH��L��A�Є�L����   H��A�ք���  L��L��A�ք�L�D�KD�C�{H�sH�K �S(�C,�  M�] L�M�]L�[M�]L�[M�]L�[M�] L�[ M�](L�[(E�ME�EA�}I�uI�M A�U(A�E,[]A\M�U A]A^�L��A�ք���   L��H��A�ք�L�D�KD�C�{H�sH�K �S(�C,�k���L�] L�L�]L�[L�]L�[L�]L�[L�] L�[ L�](L�[(L�U D�MD�E@�}H�uH�M �U(�E,[]A\A]A^�D  D�K�{L�D�CH�sH�K �S(�C,M�$L�M�\$L�[M�\$L�[M�\$L�[M�\$ L�[ M�\$(L�[(E�L$E�D$A�|$I�t$I�L$ A�T$(A�D$,[]M�$A\A]A^�fD  L�D�KD�C�{H�sH�K �S(�C,�����f.�      AWI��AVAUI��ATI�̹   UH��H� SH��H�W�H�w�L�w�̫������   H��������H��H9��6  L�<�    L���N���1�H��L��I���Ψ��H�u1�H�E    L�UH��u�UD  I�H�I�H�2H��t?H��H�F1�H�H��M��M�M��u�H�EH�H�uM�H�> ��   H��I�4�H��u�H�} 臨��1�L��L�u H��H�]I��I��K�>H�H��t8H�I�$H� L� H�EH��L��[]A\A]A^A_�L�u I��K�>H�H��u�H�UI�$L�eI�$H��tH�B1�H�uL��M�$�HE H�UH��f�     H�������é��H���۩��L�u(����H��蚩��H���©��L��躧�������H���}���H�������D  AWH��H��AVAUI��ATI��I��?UL�SH��XH�D$H�|$H�D$H�t$H�T$H�$H9��c  I��� I��I�FL�< I�_�I�<�H�$H�4[H��M�d= H��I�l5 L��H����L���K�vHE�ID�H�H��L�H;\$H�H�VH�PH�VH�PH�VH�PH�V H�P H�V(H�P(�z����D$uL�t$I��L��H��?I�I��I9���   H�[H��L�H��$�   H;\$H�T$ H��$�   H�T$(H��$�   H�T$0H��$�   H�T$8H��$�   H�T$@H��$�   H�T$H��   H�T$ H�H�T$(H�PH�T$0H�PH�T$8H�PH�T$@H�P H�T$HH�P(H��X[]A\A]A^A_�@ H�\K�vH�[H��L�H��L�H�H�
H�HH�JH�HH�JH�HH�JH�H H�J H�H(H�J(����f.�     H�C�I��I��?I�I��O�$vH�t$ H�$I��M�L���Є�H�[uH��L��"���f�I�$H��L�H�I�T$H�PI�T$H�PI�T$H�PI�T$ H�P I�T$(H�P(I�F�H��H��?H�H��L9t$|L�������D  L��I���f���H���#���f.�      AWAVI��AUATI��USH��xH�T$0H��H)�H��/  ��  H�|$0 I���B  H�O0H��H�L$8H��L)�H���������M��H��H��L��H��I�N�H�l$0L��H��H��?H�H��H�@H��I������H��f.�     H��H��L��A�Մ�H�S0u�L�}��     L��L��L��I��0A�Մ�u�H9�vfL�M L�H�{H�s�KH�S L�L�MH�C(L�KL�ML�KL�ML�KL�M L�K L�M(L�K(H�U H�S0L�E H�}H�u�MH�E(�a���H�T$0L��L��H������H��H��L)�H��/  �s  H�|$0 tI��H�\$8�����H��H��H���������I��L��M�w�I��K�,vH��L��I��H�}H�uH��0H�MHH�UPL�E0H�EXH�|$HH�t$PH�L$XH�T$`H�|$H�t$L��H�L$H�T$ L��L�D$@H�D$hL��L�$H�D$(L���+���M��u�H��0I���������I�$L�H��L�CH�{L)�H�sH�K H��0H�S0I�T$H�CXH�|$PH�t$XH�L$`H�S8I�T$L�L$@L�D$HH�D$hH�S@I�T$H�SHI�T$ H�SPI�T$(H�SXH��H�|$H��H�t$H�L$ I��1�L�$L�D$H�D$(L��L���o���H��_�M���H��x[]A\A]A^A_�fD  AWH��H��AVAUI��ATI��I��?UL�SH��XH�D$H�|$H�D$H�t$H�T$H�$H9��c  I��� I��I�FL�< I�_�I�<�H�$H�4[H��M�d= H��I�l5 L��H����L���K�vHE�ID�H�H��L�H;\$H�H�VH�PH�VH�PH�VH�PH�V H�P H�V(H�P(�z����D$uL�t$I��L��H��?I�I��I9���   H�[H��L�H��$�   H;\$H�T$ H��$�   H�T$(H��$�   H�T$0H��$�   H�T$8H��$�   H�T$@H��$�   H�T$H��   H�T$ H�H�T$(H�PH�T$0H�PH�T$8H�PH�T$@H�P H�T$HH�P(H��X[]A\A]A^A_�@ H�\K�vH�[H��L�H��L�H�H�
H�HH�JH�HH�JH�HH�JH�H H�J H�H(H�J(����f.�     H�C�I��I��?I�I��O�$vH�t$ H�$I��M�L���Є�H�[uH��L��"���f�I�$H��L�H�I�T$H�PI�T$H�PI�T$H�PI�T$ H�P I�T$(H�P(I�F�H��H��?H�H��L9t$|L�������D  L��I���f���H���#���f.�      AWAVAUI��ATUH��SH��xH�T$0H��H)�H��/  ��  H�|$0 I���Q  H�O0H��H�L$8I��H)�H���������M��H��L��H��H��I�M�H�l$0L��H��H��?H�H��H�@H��H�T ����L��f�     I��H��H��A�Ԅ�I�W0u�L�s��     L��L��H��I��0A�Ԅ�u�L9�vuL�M�E�OE�GA�I�wM�L�[A�W(I�O A�G,M�_L�[M�_L�[M�_L�[ M�_ L�[(M�_(�S(I�W0L�D�KD�C@�{H�sH�K �C,�R���H�T$0L��L��L������L��L��H)�H��/  �o  H�|$0 tI��L�t$8����I��H��H���������I��L��M�n�I��K�\m H��H��I��H�{H�sH��0H�KHH�SPL�C0H�CXH�|$HH�t$PH�L$XH�T$`H�|$H�t$H��H�L$H�T$ L��L�D$@H�D$hL��L�$H�D$(L������M��u�I�_�I���������H�U L�I��L�CH�{I)�H�sH�K H��0H�S0H�UH�CXH�|$PH�t$XH�L$`H�S8H�UL�L$@L�D$HH�D$hH�S@H�UH�SHH�U H�SPH�U(H�SXL��H�|$H��H�t$H�L$ I��1�L�$L�D$H�D$(L��H���d���I��_�R���H��x[]A\A]A^A_�f.�     �AWH��H��AVI��AUATI��I��?UL�SH��8H�D$H�|$H�D$H�t$ H�T$(H�L$L�D$H9��*  I���D  I��I�OH�D$L�$	H��M�,I�\$�I�} I�,�H�u ��L���HE�ID�H;\$H�K��|��D$(uL�|$(I��L��H��?I�I��I9�t0H;\$ I��<D  H�T$H�H��8[]A\A]A^A_�f�     H�\H;\$ I��H�K��~�H�C�I��I��?I�I��O�$�H�t$H�D$I�<$�Є�uI���f.�     I�$I��I�G�H��H��?H�H��L9|$ |L���h����     L��I���H������AWAVI��AUATUH��SH��H�$H��H)�H���   ��  H�<$ I���e  H�OH��H�L$H)�H�,$H�}H��H��H��?H�H��H�\� H�3A�Ԅ�I�v���   H�;A�Ԅ���   H�E H�H�U H�H�u M��H�T$�	D  H�u H��H�:A�Ԅ�H�Su�M�}��    I�7M��H�} I��A�Ԅ�u�I9�vI�U H�H�I�E H�SH�u �H�$L��L��H�������H��H��H)�H���   ��   H�<$ txI������H�}A�Ԅ�t=H�E H�uH�EH�u �H���I�v�H�}A�Ԅ�H�E t�I�V�H�U I�F�H�u � ���I�v�H�;A�Ԅ�H�E u�H�H�U H�H�u �����H��H��L�j�I��I���I��J�L� M��L��L��H������M��u�H��H�E I��H�I)�1�M��L��H��H��H�CH������I���H��[]A\A]A^A_�f.�     @ AWI��AVAUATUSH��H��H9���   H��I�պ?   H)�L��H��H��H��H��?H�H)�H�����H���   ~dL���   L��H��L������M9�t; L��M�&I�^��@ H�H�E H��H��H�3L��A�Մ�u�I��L�e M9�u�H��[]A\A]A^A_�L��L��H������H��[]A\A]A^A_�f.�     f�AWA��AVI��AUI��ATL�% $  UH�-($  SL)�1�H��H������H��t�     L��L��D��A��H��H9�u�H��[]A\A]A^A_�ff.�     ��f�H��H���                 r Error opening file %s
 RSW_test.txt refFlat.txt RSW_tst %s.results w %s.results.full %s.results.unknown %s.results.unknown.full %s.results.splitPairs Novel * %s	%s	%li	%li	%li--%li	%s /proc/self/status VmRSS: %19s %19s %999s .gz gunzip -c %s .lrz cat %s | ./lrunzip Error reading from file %s
 	 -- %s	%s	%s	%c	%i-%i	%i	%i-%i
 %s,  , %s      Options file %s is more than the max of %i bytes.
      Not enough lines in options file.       refFlat.txt.intronBoundary.exonsgaps    Not enough arguments given, using default values.       Usage is to load options from file: ./splitPairs optionsFile.txt        And make sure options file has options in order, each on their own line with no extra lines.    Running with options...
  file with read data          %s
  max distance between matches %i
  length of samples            %i
  refFlat file                 %s
  refFlat intron boundary file %s
  minimum splice length        %i
  tolerance of difference in position for supporting reads  %i
  base of file name for writing results                     %s
  minimum number of supporting reads                        %i

      Error, filename %s is too long.
        Error opening file %s for writing.
     Will write summary results that match in known genes to file
   %s
     Will write full results that match in known genes to file
   %s
        Will write summary results that do NOT match in known genes to file
   %s
      Will write full results that do NOT match in known genes to file
   %s
 Will write split pairs to file
   %s
   Finished processing data, results written to files.
    Number of entries in data file:             %i
 Number of different reads:                  %i
 Number of entries in refFlat file:          %i
 Number of entries in refFlat boundary file: %i
 Number of matches:                          %i
 String table size:                          %i
 VmRSS, memory resident set size:            %s %s
      Total time to process:                      %i seconds
 Error reading data file %s, line exceeded %i characters.
       Done reading/sorting refFlat, total time elapsed %i seconds
    Done reading refFlat intron/exon boundaries, total time elapsed %i seconds
     Done reading read data, total time elapsed %i seconds
  Done sorting read data, total time elapsed %i seconds
  Done finding matched pairs, total time elapsed %i seconds
      Done sorting matched pairs, total time elapsed %i seconds
      Done computing supporting reads, total time elapsed %i seconds
 Done sorting matched pairs again, total time elapsed %i seconds
        Done filtering matched pairs, total time elapsed %i seconds
    Done calculating/filtering supporting reads, total time elapsed %i seconds
     Line No	Id	Gene	Chr	Strand	Splice region	Supporting reads	Supporting splice range
      GeneName	Chromosome	# supporting reads	splice length	range of supporting reads	Novel or not (*)
        Done saving results, total time elapsed %i seconds
                     `L@     0L@     �L@     �K@     L@     �K@     `K@     `K@     HK@     �P@      Q@     xP@     �P@     �O@     �P@     U@     �T@     0U@     �T@     8T@     8T@     8T@     8T@     8T@     8T@     8T@     �T@     �S@     UNFOUND_                                                                                            ;�  1   ܉���  ̌���  ���<  H����  <���  l����  �����  ����  L���$  ���t  �����  �����  ����  �����  ,���  ̽��D  L���d  ����|  ����t  �����  ����  ����l  ���  ����L  \���  |���4  ����L  ����d  ����|  <����  ����  l����  ����<  ���|  ����  ����$  �����  �����  �����  ����T	  �����	  L����	  <���D
  �����
  �����
  |���4  �����  |���l  �����         zR x�      ����*                  zR x�  $      �����   FJw� ?;*3$"       D   0���.              \   8���              t   @���              �   H���              �   P���           $   �   X���^    A�D�D NAA   �   ����F              �   ����F                0���>           4   ,  H����    B�D�A �D0p
 AABA     L   d   ����    B�E�D �C(�F0T
(A ABFFN
(C ABBB      �  p����              �  (���>              �  P���           <   �  H����   B�J�A �A(�D�%
(A ABBA      <  �����    D{,   T  ����    A�F�GPr
AAA        �  �����   D_
E       �  ���Q    D0L<   �  (���   B�O�G �A(�G��
(A ABBK   <   �  X���A   B�L�A �F(�G@�
(A ABBG     <   <  h���A   B�L�A �F(�G@�
(A ABBG     <   |  x���I   B�L�A �F(�G@�
(A ABBG     <   �  �����    B�B�G �A(�G@�
(A ABBG     D   �  H����    B�B�B �B(�A0�A8�DP�8A0A(B BBB       zPLR xp@ �  L   $   ����  ��@ B�E�B �E(�I0�H8�DP
8D0A(B BBBDT   t   P���	  �@ B�J�B �A(�A0��
(A BBGEd
(A BBBA     L   �   ���  D�@ B�B�B �E(�A0�A8�Dp(
8A0A(B BBBAL     �����  t�@ B�G�E �B(�A0�A8�Dp�
8A0A(B BBBAL   l  X���C  ��@ B�G�B �B(�A0�A8�DpQ
8A0A(B BBBAL   �  8����  Ԕ@ B�E�B �E(�I0�H8�D@�
8D0A(B BBBA <   L  ����    B�F�A �D0�
 FABED FAB 4   �  x����    B�D�D �F0Q
 AABK     d   �  @����   B�E�E �G(�D0��
(A BFBHy
(A BBBHe
(A FBBH     d   ,  ����   B�E�E �G(�D0��
(A BFBA�
(A BBBFv
(A FBBG     L   T  `����  ��@ B�E�B �E(�I0�H8�D@�
8D0A(B BBBA L   �  �����   B�I�B �E(�H0�D8�D�w
8A0A(B BBBE   L   4  P����   B�B�E �B(�D0�A8�D��8A0A(B BBB      L   �  �����   B�I�B �E(�H0�D8�D�w
8A0A(B BBBE   L   �  `����   B�B�B �E(�A0�D8�D��8A0A(B BBB      L   $	  ����   B�I�E �B(�H0�D8�Dp�
8A0A(B BBBJ     L   t	  @���   B�B�E �B(�A0�D8�DP�8A0A(B BBB       \   �	  ����    B�E�B �B(�A0�A8�G@�
8A0A(B BBBAR8A0A(B BBBT   �  ����F(  $�@ B�G�B �B(�A0�F8�K��
8C0A(B BBBA       ,   |
  ب��)   B�F�A ��
ABu  D   �
  ���e    B�E�E �E(�H0�H8�M@l8A0A(B BBB    �
  0���               �%/  ]������ �    }    �%k  ���  �� �          �-"k�  ���
 ��	�y  �	�
 �
          �-"�  ��� ���l  �� �          �-"�  ��� ���n  �� �          �%/  ]������ �    }    �%/  ]������ �    }    ��?�  �-�M ��2  �M�M                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0?@     @<@     ?@                                  E             O             l             �@            z@            ؝`                          �`                   ���o    �@            �@            �@     
       {                                           �`            P                           �@            @@            H       	              ���o    �@     ���o           ���o    D@                                                                                                             ��`                     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@         GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-7) GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-9)  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.gnu.build-id .gnu.hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .gcc_except_table .init_array .fini_array .jcr .dynamic .got .got.plt .data .bss .comment                                                                                8@     8                                    #             T@     T                                     1             t@     t      $                              D   ���o       �@     �      4                             N             �@     �      �                          V             �@     �      {                             ^   ���o       D@     D      j                            k   ���o       �@     �      �                            z             @@     @      H                            �             �@     �      P                          �             �@     �                                    �              @            �                            �             �@     �      $f                             �             z@     z      	                              �              z@      z                                    �             $�@     $�      �                             �             ��@     ��      <                             �             ��@     ��      K                             �             ؝`     ؝                                    �             �`     �                                    �             �`     �                                    �             ��`     ��                                  �             ��`     ��                                   �              �`      �      �                                        ��`     ��                                                ��`     ��      ��                                   0               ��      X                                                   �                                                         ��      H         3                 	                      �                                                                 8@                   T@                   t@                   �@                   �@                   �@                   D@                   �@                  	 @@                  
 �@                   �@                    @                   �@                   z@                    z@                   $�@                   ��@                   ��@                   ؝`                   �`                   �`                   ��`                   ��`                    �`                   ��`                   ��`                                       ��                     �H@     	      =    @W@     �       �    ��@     d       �    @<@     )      �    T@a            �   ��                �    �`             �    �>@             �    �>@                 ?@             "    ȡ`            1    �`             X    0?@             d    ؝`             �   ��                �    �@             �    �`                  ��                �     �`             �     �`             �     ؝`             �    ��`             �     ��`             �    A@     �           0@@     >       .                     B  "  �X@            `    @@a            m  "  0^@     �       �                     �    z@            �    l>@             �  "   r@     �      M    `�`     '      S  "   Y@     ^       �                      �                      �                     �     �`            �    P@a            �                     �                                          9                     K  "  �X@            i    z@             o                     �                     �  "   u@     �      %                     9                     L                     k  "  �X@            �  "  0^@     �       �     �`            �    �A@     >       �                       "  �i@     �      �                     �    ��`            �      @             �    ��`            �                      �     z@            	                     	    �F@           7	    �?a            K	    ��`            P	                     b	                     v	     �`     0       �	                      �	    ��`             �	                     �	  "  0_@     �       I
    ��`     '      Q
    �R@     C      h
                     }
                     �
    �`            �
     @a            �
    `?@     .       �
    �a     '      �
                     �
   ��`             �
                         0B@     �      4                     H  "  �\@     I      �    pF@     Q       �  "  `@           �   (z@             �    �W@     �       �  "  �[@     A      ,                     @    �y@     e       P  "  `Y@     �       �    �?@     F       �  "  �c@     �      5    �D@     �      J    @�`     '      O  "  �X@            y  "   b@     �      �                     �    ��`            �  "  �X@            �    �?@     F       �    �?a                 B@            0  "  �x@     �       �    ��`             �  "  �\@     I      �    �?a            �    `�`     0       	  "  pl@     �      |     O@     �      �                     �                      �                     �  "   Y@     ^                            5     @a            I                     ]    X@a             b                     v    0@a            �  "  �X@            �                     �                                          7    ��`             f                     z    @a            �    p@@     �       �                     �  "  �X@            �    ��`            �    �`            �                         �C@     �       (  "  �g@     �      &  "  �[@     A      k    ��`            w    ��`             ~     J@           �     p@             �                     �  "  �v@           8                     K  "  �e@           �                     �  "  `o@     �      H                     [    ��`            m  "  @Z@     A      �  "  @Z@     A      �    �`            �                         PD@     �       <  "  `Y@     �       y    ��`            �  "  �X@            �                     �    �@     F(      �    �@             �                      splitPairs.cpp _ZNSt10_HashtableISsSsSaISsENSt8__detail9_IdentityESt8equal_toISsESt4hashISsENS1_18_Mod_range_hashingENS1_20_Default_ranged_hashENS1_20_Prime_rehash_policyENS1_17_Hashtable_traitsILb1ELb1ELb1EEEE9_M_insertIRKSsEESt4pairINS1_14_Node_iteratorISsLb1ELb1EEEbEOT_St17integral_constantIbLb1EE.constprop.175 _ZNSt10_HashtableIPKcSt4pairIKS1_P10RSW_spliceESaIS6_ENSt8__detail10_Select1stESt8equal_toIS1_ESt4hashIS1_ENS8_18_Mod_range_hashingENS8_20_Default_ranged_hashENS8_20_Prime_rehash_policyENS8_17_Hashtable_traitsILb0ELb0ELb1EEEE9_M_insertIRKS6_EES2_INS8_14_Node_iteratorIS6_Lb0ELb0EEEbEOT_St17integral_constantIbLb1EE.isra.157 _ZL14unfound_string _GLOBAL__sub_I__Z8get_lineP8_IO_FILEPci _ZStL8__ioinit crtstuff.c __JCR_LIST__ deregister_tm_clones register_tm_clones __do_global_dtors_aux completed.6342 __do_global_dtors_aux_fini_array_entry frame_dummy __frame_dummy_init_array_entry __FRAME_END__ __JCR_END__ _GLOBAL_OFFSET_TABLE_ __init_array_end __init_array_start _DYNAMIC data_start _Z12compute_hashP3RSW _Z24compare_spliceByChromLenP10RSW_spliceS0_ printf@@GLIBC_2.2.5 _ZNSt6vectorI3RSWSaIS0_EED1Ev sampleLength _ZNSt6vectorIP10RSW_spliceSaIS1_EE19_M_emplace_back_auxIJRKS1_EEEvDpOT_ memset@@GLIBC_2.2.5 __libc_csu_fini _start _ZSt16__introsort_loopIN9__gnu_cxx17__normal_iteratorIP3RSWSt6vectorIS2_SaIS2_EEEElPFbRKS2_S9_EEvT_SC_T0_T1_ sLine _ZNSt13unordered_setIPKcSt4hashIS1_ESt8equal_toIS1_ESaIS1_EED2Ev __gmon_start__ _Jv_RegisterClasses puts@@GLIBC_2.2.5 fKnown maxDistance _ZdlPv@@GLIBCXX_3.4 _ZNSs7reserveEm@@GLIBCXX_3.4 _ZNSs4_Rep10_M_disposeERKSaIcE@@GLIBCXX_3.4 exit@@GLIBC_2.2.5 _ZNSt6vectorI3RSWSaIS0_EED2Ev _fini __cxa_rethrow@@CXXABI_1.3 _ZNSt8ios_base4InitC1Ev@@GLIBCXX_3.4 _ZSt13__adjust_heapIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEElS3_PFbS3_S3_EEvT_T0_SC_T1_T2_ malloc@@GLIBC_2.2.5 fopen@@GLIBC_2.2.5 __libc_start_main@@GLIBC_2.2.5 _ZNSt6vectorI9RSW_KnownSaIS0_EED1Ev _ZNSt6vectorIP10RSW_spliceSaIS1_EE19_M_emplace_back_auxIIRKS1_EEEvDpOT_ fKnownFull _Z12compare_dataRK3RSWS1_ _ZNSsC1ERKSs@@GLIBCXX_3.4 _ZSt13__adjust_heapIN9__gnu_cxx17__normal_iteratorIP9RSW_KnownSt6vectorIS2_SaIS2_EEEElS2_PFbRKS2_S9_EEvT_T0_SD_T1_T2_ __cxa_atexit@@GLIBC_2.2.5 data_known _ZNSt8ios_base4InitD1Ev@@GLIBCXX_3.4 beginTime _ITM_deregisterTMCloneTable _IO_stdin_used fputc@@GLIBC_2.2.5 _Z10printStatsP8_IO_FILE refFlatBoundaryFile data free@@GLIBC_2.2.5 strlen@@GLIBC_2.2.5 readIdsReported _ITM_registerTMCloneTable __data_start _ZNSs4_Rep10_M_destroyERKSaIcE@@GLIBCXX_3.4 _ZSt16__insertion_sortIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEEPFbS3_S3_EEvT_SB_T0_ options _Z15read_boundariesPKc sprintf@@GLIBC_2.2.5 fgetc@@GLIBC_2.2.5 fUnknownFull refFlatFile _Z18compare_data_knownRK9RSW_KnownS1_ buff setpriority@@GLIBC_2.2.5 __TMC_END__ _ZNSsC1EPKcRKSaIcE@@GLIBCXX_3.4 _Z19readOptionsFromFilePKc strstr@@GLIBC_2.2.5 _ZNSt6vectorI14RSW_BoundariesSaIS0_EE19_M_emplace_back_auxIIRKS0_EEEvDpOT_ _Z11printSpliceP8_IO_FILEP10RSW_splice _ZNSt10_HashtableISsSsSaISsENSt8__detail9_IdentityESt8equal_toISsESt4hashISsENS1_18_Mod_range_hashingENS1_20_Default_ranged_hashENS1_20_Prime_rehash_policyENS1_17_Hashtable_traitsILb1ELb1ELb1EEEE21_M_insert_unique_nodeEmmPNS1_10_Hash_nodeISsLb1EEE __dso_handle _Z19updateSpliceSupportP10RSW_spliceS0_ _ZNSt6vectorI9RSW_KnownSaIS0_EE19_M_emplace_back_auxIIRKS0_EEEvDpOT_ strtol@@GLIBC_2.2.5 __libc_csu_init _ZNSt13unordered_setISsSt4hashISsESt8equal_toISsESaISsEED1Ev _Z24compare_spliceByChromPosP10RSW_spliceS0_ _ZSt22__move_median_to_firstIN9__gnu_cxx17__normal_iteratorIP9RSW_KnownSt6vectorIS2_SaIS2_EEEEPFbRKS2_S9_EEvT_SC_SC_SC_T0_ _Z15openOutputFilesv temp _ZNSt6vectorI14RSW_BoundariesSaIS0_EED2Ev _ZNSt10_HashtableIPKcSt4pairIKS1_P10RSW_spliceESaIS6_ENSt8__detail10_Select1stESt8equal_toIS1_ESt4hashIS1_ENS8_18_Mod_range_hashingENS8_20_Default_ranged_hashENS8_20_Prime_rehash_policyENS8_17_Hashtable_traitsILb0ELb0ELb1EEEE21_M_insert_unique_nodeEmmPNS8_10_Hash_nodeIS6_Lb0EEE memmove@@GLIBC_2.2.5 endTime _ZNSt6vectorIP10RSW_spliceSaIS1_EED1Ev _Z18compare_dataToSortRK3RSWS1_ resultsBaseName _Z23compare_spliceBySupportP10RSW_spliceS0_ _ZSt4sortIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEEPFbS3_S3_EEvT_SB_T0_ __bss_start _ZNSt6vectorI14RSW_BoundariesSaIS0_EE19_M_emplace_back_auxIJRKS0_EEEvDpOT_ minSupportingReads stringTable _ZSt16__introsort_loopIN9__gnu_cxx17__normal_iteratorIP9RSW_KnownSt6vectorIS2_SaIS2_EEEElPFbRKS2_S9_EEvT_SC_T0_T1_ _Z14read_knownGenePKc strcpy@@GLIBC_2.2.5 __pthread_key_create strtok@@GLIBC_2.2.5 _ZNSt13unordered_setIPKcSt4hashIS1_ESt8equal_toIS1_ESaIS1_EED1Ev _ZSt11_Hash_bytesPKvmm@@CXXABI_1.3.5 supportPosTolerance memcmp@@GLIBC_2.2.5 _end fclose@@GLIBC_2.2.5 minSpliceLength _ZNSt6vectorI14RSW_BoundariesSaIS0_EED1Ev _ZNKSt8__detail20_Prime_rehash_policy11_M_next_bktEm@@GLIBCXX_3.4.18 __cxa_end_catch@@CXXABI_1.3 _ZSt17__throw_bad_allocv@@GLIBCXX_3.4 _ZNSs4_Rep20_S_empty_rep_storageE@@GLIBCXX_3.4 fscanf@@GLIBC_2.2.5 sampleDataFile _Z8get_lineP8_IO_FILEPci __cxa_begin_catch@@CXXABI_1.3 _ZNSt6vectorI9RSW_KnownSaIS0_EED2Ev data_boundaries fUnknown fwrite@@GLIBC_2.2.5 _Z17setDefaultOptionsv _ZNSt10_HashtableIPKcS1_SaIS1_ENSt8__detail9_IdentityESt8equal_toIS1_ESt4hashIS1_ENS3_18_Mod_range_hashingENS3_20_Default_ranged_hashENS3_20_Prime_rehash_policyENS3_17_Hashtable_traitsILb0ELb1ELb1EEEE21_M_insert_unique_nodeEmmPNS3_10_Hash_nodeIS1_Lb0EEE _ZNSt6vectorI9RSW_KnownSaIS0_EE19_M_emplace_back_auxIJRKS0_EEEvDpOT_ data_splice _edata _Z9read_dataPKc __gxx_personality_v0@@CXXABI_1.3 fprintf@@GLIBC_2.2.5 _ZSt16__introsort_loopIN9__gnu_cxx17__normal_iteratorIPP10RSW_spliceSt6vectorIS3_SaIS3_EEEElPFbS3_S3_EEvT_SB_T0_T1_ _Znwm@@GLIBCXX_3.4 _ZSt22__move_median_to_firstIN9__gnu_cxx17__normal_iteratorIP3RSWSt6vectorIS2_SaIS2_EEEEPFbRKS2_S9_EEvT_SC_SC_SC_T0_ _Unwind_Resume@@GCC_3.0 _ZSt13__adjust_heapIN9__gnu_cxx17__normal_iteratorIP3RSWSt6vectorIS2_SaIS2_EEEElS2_PFbRKS2_S9_EEvT_T0_SD_T1_T2_ popen@@GLIBC_2.2.5 numDifferentReads _ZNSt6vectorI3RSWSaIS0_EE19_M_emplace_back_auxIJRKS0_EEEvDpOT_ _ZNSt6vectorI3RSWSaIS0_EE19_M_emplace_back_auxIIRKS0_EEEvDpOT_ fSplitPairs _ZNSs12_M_leak_hardEv@@GLIBCXX_3.4 _Z19printCurrentOptionsP8_IO_FILE _ZNSt13unordered_setISsSt4hashISsESt8equal_toISsESaISsEED2Ev stdout@@GLIBC_2.2.5 _ZNSt6vectorIP10RSW_spliceSaIS1_EED2Ev time@@GLIBC_2.2.5 main _init _ZNKSt8__detail20_Prime_rehash_policy14_M_need_rehashEmmm@@GLIBCXX_3.4.18                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              src/                                                                                                0000775 0001753 0001753 00000000000 12563406227 011466  5                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                src/RSW.h                                                                                           0000654 0001753 0001753 00000003625 12563405654 012322  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #ifndef RSW_H_
#define RSW_H_

class RSW {
 public:
  const char *id;
  char side;
  int length;
  char direction;
  const char *chromosome;
  long int position;
  int count;
  int hash;
};



class RSW_Known {
 public:
  const char *id1;
  const char *id2;
  const char *chromosome;
  char direction;
  long int position1;
  long int position2;
};

class RSW_Boundaries {
 public:
  const char *id1;
  const char *id2;
  const char *chromosome;
  char direction;
  long int length;
  long int position1;
  long int position2;
};



class RSW_splice {
 public:
  const char * id;
  const char * geneName;
  int geneUnknown;
  const char * chromosome;
  char direction;
  long int positionSmaller;
  long int positionLarger;
  long int minSmallSupport;
  long int maxLargeSupport;
  unordered_map <const char *, RSW_splice *> supported_reads;
  //vector <RSW_splice *> supported_reads;
  bool reported; // has this result already been reported as supporting another read
  bool novel; // does this show up in refFlatBoundary file as already known
};

class RSW_result {
 public:
  const char * geneName;
  const char * chromosome;
  long int supportCount;
  long int spliceLength;
  long int minSmallSupport;
  long int maxLargeSupport;
  bool novel; // does this show up in refFlatBoundary file as already known
};



// return 1 if line received, return 0 if EOF was encountered, -1 if error
// remove newline character from string.
int get_line(FILE *f, char s[], int maxChars) {
  int i = 0;
  int ch;
  do {
    ch = fgetc(f);
    if (ch == EOF) 
      break;
    if (i < maxChars) {
      s[i] = ch; i++;
    }
  } while (ch != '\n');
  if (i > 0) {
    if (s[i-1] == '\n')
      s[i-1] = '\0';
    else 
      s[i] = '\0';
  }
  else 
    s[0] = '\0';

  if (ch == EOF) return 0;
  else if (i >= maxChars) return -1; // error, not enough room in s for line
  else return 1;
}

#define MAX_LINE 10000
#define MAX_STR_LEN 10000


#endif
                                                                                                           src/split_read_rsw.c                                                                                0000664 0001753 0001753 00000007553 12563406227 014665  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <libgen.h>

const unsigned int K = 1024;
int LEFT, RIGHT, READLEN;
const char *SUFFIX = ".split";

void parse(FILE *in, FILE *out) {
    //char id[K],seq[K],str[K],score[K], tmp[K];
    char *id,*seq,*str,*score,tmp[K];
    size_t n;
    unsigned long i = 0u, k = 0u, l = 0u;
    unsigned int j;
    unsigned int left_cut, right_cut; 
    id = seq = str = score = 0;
    while(!feof(in)) {
        left_cut = LEFT, right_cut = RIGHT;
        if (getline(&id, &n, in) == -1) break;
            if (id[strlen(id)-1]== '\n') id[strlen(id)-1] = 0;
        if (getline(&seq, &n, in) == -1) break;
            if (seq[strlen(seq)-1] == '\n') seq[strlen(seq)-1] = 0;
        if (getline(&str, &n, in) == -1) break;
            if (str[strlen(str)-1] == '\n') str[strlen(str)-1] = 0;
        if (getline(&score, &n, in) == -1) break;
            if (score[strlen(score)-1] == '\n') score[strlen(score)-1] = 0;
        k+=4;
        l = 0u;
        //must have all 4 to contnue...
        while(left_cut <= right_cut) {
            ++l;
            //left
            fprintf(out, "%s-L-%d\n",id,left_cut);
            strncpy(tmp,seq,left_cut);
            tmp[left_cut] = 0;
            fprintf(out, "%s\n",tmp);
            fprintf(out, "%s-L-%d\n",str,left_cut);
            strncpy(tmp,score,left_cut);
            tmp[left_cut] = 0;
            fprintf(out, "%s\n",tmp);
            //right
            fprintf(out, "%s-R-%d\n",id,right_cut);
            fprintf(out, "%.*s\n",right_cut,seq+left_cut);
            fprintf(out, "%s-R-%d\n",str,right_cut);
            fprintf(out, "%.*s\n",right_cut,score+left_cut);
            //REVERSED!
            //left
            fprintf(out, "%s-L-%d\n",id,right_cut);
            strncpy(tmp,seq,right_cut);
            tmp[right_cut] = 0;
            fprintf(out, "%s\n",tmp);
            fprintf(out, "%s-L-%d\n",str,right_cut);
            strncpy(tmp,score,right_cut);
            tmp[right_cut] = 0;
            fprintf(out, "%s\n",tmp);
            //right
            fprintf(out, "%s-R-%d\n",id,left_cut);
            fprintf(out, "%.*s\n",left_cut,seq+right_cut);
            fprintf(out, "%s-R-%d\n",str,left_cut);
            fprintf(out, "%.*s\n",left_cut,score+right_cut);
            //next
            ++left_cut; --right_cut;
            i+=16;
        }
        free(id); id = 0;
        free(seq); seq = 0;
        free(str); str = 0;
        free(score); score = 0;
    }
    fprintf(stderr,"in: %lu, out: %lu, passes each line: %lu, total passes: %lu\n", k, i, l, l*k);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr,"Usage; %s <file> <min size> <total size>\n",argv[0]);
        return 1;
    }
    clock();
    clock_t after, before;
    before = clock();
    FILE *fr, *fw;
    char *newfn, *name;
    int sp = 0, sn, ss = strlen(SUFFIX);
    char flag = 0;
    if (argc >= 5) {
        sp = strlen(argv[4]);
        if ((flag = ('/' == argv[4][sp-1]))) ++sp;
        name = basename(argv[1]);
    } else {
        name = argv[1];
    }
    sn = strlen(name);
    newfn = (char *)malloc(sizeof(char) * (sp+sn+ss));
    memset(newfn,0,sizeof(char)*(sp+sn+ss));
    if (sp) strcpy(newfn,argv[4]);
    if (sp && !flag) strcat(newfn,"/");
    newfn = strcat(strcat(newfn, name),SUFFIX);

    LEFT = atoi(argv[2]);
    READLEN = atoi(argv[3]);
    RIGHT = READLEN - LEFT;
    if (!(fr = fopen(argv[1], "r"))) {
        fprintf(stderr,"Unable to open %s for reading", argv[1]);
        perror("");
        return 1;
    }

    if (!(fw = fopen(newfn, "w"))) {
        fprintf(stderr,"Unable to open %s for reading", argv[1]);
        perror("");
        fclose(fr);
        return 1;
    }
    
    parse(fr,fw);
    fclose(fr);
    fclose(fw);
    after = clock();
    fprintf(stderr,"time elapsed: %lf seconds.\n",(after-before)*(1.0/CLOCKS_PER_SEC));
    return 0;
}
                                                                                                                                                     src/split_bowtie_chrom.c                                                                            0000664 0001753 0001753 00000004721 12563404716 015533  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHR "chr"
#define SUFFIX ".txt"
#define MAX_CHROMOSOMES (308 + 1)
#define MAX_FILES  ( MAX_CHROMOSOMES + 5 )

#define chrX "chrX"
#define chrY "chrY"
#define chrM "chrM"
#define UC "chr_uc"
#define UNKNOWN "chr_unknown"

FILE *chrFiles[MAX_FILES];
int openFiles = 0;
char *prefix;

FILE *pick_file(char *chr) {
    int id;
    char *fn_chr = chr;
    if (!chr) id = MAX_FILES -1;
    else {
        id = atoi((chr+3));
        if (0 == id) {
            if (strcmp(chr, chrX) == 0)  id = MAX_CHROMOSOMES +1;
            else if (strcmp(chr, chrY) == 0) id = MAX_CHROMOSOMES +2;
            else if (strcmp(chr, chrM) == 0) id = MAX_CHROMOSOMES +3;
            else if (strncmp(chr, "uc_", 3) == 0) {
                id = MAX_CHROMOSOMES +4;
                fn_chr = UC;
            } else {
                //id = 0
                fn_chr = UNKNOWN;
            }
       }
    }
    //fprintf(stderr,"pick_file:: guessing file %i for %s\n",id,chr+3);
    if (!chrFiles[id]) {
        char *fn;
        int sz = sizeof(char) * (strlen(fn_chr)+strlen(prefix)+strlen(SUFFIX)+1);
        fn = (char *)malloc(sizeof(char) * sz);
        memset(fn, 0, sz);
        sprintf(fn,"%s.%s%s",prefix,fn_chr,SUFFIX);
        if (!(chrFiles[id] = fopen(fn,"w"))) {
            fprintf(stderr,"Panic! Cannot open chrom file %s for writing!\n",fn);
            perror("");
            fprintf(stderr,"Ignoring.");
            return 0;
        }
    }
    return chrFiles[id];
}

char *get_chrom(char *line) {
    char *tok;
    tok = strtok(line,"\t");
    tok = strtok(NULL,"\t");
    tok = strtok(NULL,"\t");
    return tok;
}

void parse(FILE *in) {
    char *line = 0;
    char *tmp, *chrom;
    size_t n;
    FILE *outfile;

    tmp = 0;
    while(!(feof(in)) && (getline(&line,&n,in) != -1)) {
        if (tmp) { free(tmp); tmp=0; }
        tmp = strdup(line);
        chrom = get_chrom(tmp);
        outfile = pick_file(chrom);
        if (!outfile) continue;
        fprintf(outfile,"%s",line);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <reads file>\n", argv[0]);
        return 1;
    }
    FILE *fr;

    memset(chrFiles, 0, MAX_CHROMOSOMES * sizeof(FILE *));

    prefix = strdup(argv[1]);
    if (!(fr = fopen(argv[1], "r"))) {
        fprintf(stderr, "Could not open %s for reading", argv[1]);
        perror("");
        return 1;
    }
    
    parse(fr);
    fcloseall();
    return 0;
}
                                               src/splitPairs.cpp                                                                                  0000640 0001753 0001753 00000071673 12563404616 014334  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                /*
  To compile: g++ splitPairs.cpp -o sp -O4 -std=c++11
  To run:     ./sp options.txt
  where options.txt is an options file.  If the program is run
  with no command-line arguments it by default processes 
  RSW_test.txt with the same parameters as the scripts downloaded
  from Yongsheng Bai's website.  When you run the program it
  prints which files output is written to.

Modification history...
* sp4 - 7/30/2015
* 7/30/2015 - fix bug in computing supporting reads that appeared in sp3.
* 7/30/2015 - fixed bug where compare_data was being used in main loop that computes
              matched pairs; compare_data was written to be used in sorting, but was
	      being used looking for exact match between lines.  Just stopped using
	      it in the main loop and put in the if tests there (to avoid repeating
	      this mistake in the future).
* 7/30/2015 - had been working on using openmp to parallelize the main loops.  was
              getting seg fault for unknown reason, so commented out all #pragma
	      lines that were part of that effort.
* sp3 - 7/29/2015
* 7/29/2015 - Jeff Kinne - fixed bug that was reporting the number of supporting reads
              as one two small for each junction.
* 7/9/2015 - Jeff Kinne - store strings in stringTable to avoid duplicating the 
             strings.  Should save memory and speed up (because compare string pointers
	     rather than strcmp).
* 7/8/2015 - started tracking mod history...
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <time.h>
#include <sys/resource.h>
//#include <omp.h>
using namespace std;

#include "RSW.h"


int maxDistance;
int sampleLength;
int minSpliceLength;
int supportPosTolerance;
char const *sampleDataFile;
char const *refFlatFile;
char const *refFlatBoundaryFile;
char const *resultsBaseName;
int minSupportingReads;

char buff[MAX_STR_LEN];
char options[MAX_STR_LEN];

time_t beginTime, endTime;

unordered_set<string> stringTable;
unordered_set<const char *> readIdsReported;

vector<struct RSW> data;
vector<struct RSW_Known> data_known;
vector<struct RSW_Boundaries> data_boundaries;
const char unfound_string[100] = "UNFOUND_";
vector<RSW_splice *> data_splice;
int numDifferentReads;

int compute_hash(RSW *d) {
  int h = 1, i;
  for(i=0; d->id[i] != '\0'; i++) h = h * d->id[i] % SHRT_MAX;
  for(i=0; d->chromosome[i] != '\0'; i++) h = h * d->chromosome[i] % SHRT_MAX;
  h = h * d->direction % SHRT_MAX;
  return h;
}


char sLine[MAX_LINE+1];
char temp[MAX_LINE+1];

void read_data(const char *filename) {
  FILE *f;
  int len = strlen(filename);
  if (len > 3 && strcmp(filename+len-3,".gz")==0) {
    sprintf(temp,"gunzip -c %s", filename);
    f = popen(temp, "r");
  }
  else if (len > 4 && strcmp(filename+len-4,".lrz")==0) {
    sprintf(temp,"cat %s | ./lrunzip", filename);
    f = popen(temp, "r");
  }
  else f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW *r = new RSW;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete r;
      break;
    }
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0:
	result = stringTable.insert(temp);
	r->id = (result.first)->c_str();

	//r->id = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(r->id, temp);
	break;
      case 1:
	r->side = temp[0];
	break;
      case 2:
	r->length = atoi(temp.c_str());
	break;
      case 3:
	r->direction = temp[0];
	break;
      case 4:
	result = stringTable.insert(temp);
	r->chromosome = (result.first)->c_str();
	//r->chromosome = (char *) malloc(sizeof(char) * (strlen(temp)+1));
	//strcpy(r->chromosome, temp);
	break;
      case 5:
	r->position = atol(temp.c_str());
	break;
      case 8:
	r->count = atoi(temp.c_str());
	break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 9) {delete r; break;}

    r->hash = compute_hash(r);
    data.push_back(*r);
    delete r;

    //if (data.size() >= 15000000) break; // cut off early, for debugging.
  }

  fclose(f);
}

void read_knownGene(const char *filename) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW_Known * rk = new RSW_Known;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete rk;
      break;
    }
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0:
	result = stringTable.insert(temp);
	rk->id1 = (result.first)->c_str();

	//rk->id1 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id1, temp);
	break;
      case 1:
	result = stringTable.insert(temp);
	rk->id2 = (result.first)->c_str();

	//rk->id2 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id2, temp);
	break;
      case 2:
	result = stringTable.insert(temp);
	rk->chromosome = (result.first)->c_str();

	//rk->chromosome = (char *) malloc(sizeof(char) * (strlen(temp)+1));
	//strcpy(rk->chromosome, temp);
	break;
      case 3:
	rk->direction = temp[0];
	break;
      case 4:
	rk->position1 = atol(temp.c_str());
	break;
      case 5:
	rk->position2 = atol(temp.c_str());
	break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 11) {delete rk; break;}

    if (rk->position1 > rk->position2) {
      int t = rk->position1;
      rk->position1 = rk->position2;
      rk->position2 = t;
    }
    data_known.push_back(*rk);
    delete rk;
  }

  fclose(f);
}

void read_boundaries(const char *filename) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW_Boundaries * rk = new RSW_Boundaries;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete rk;
      break;
    }
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0:
	result = stringTable.insert(temp);
	rk->id1 = (result.first)->c_str();

	//rk->id1 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id1, temp);
	break;
      case 1:
	result = stringTable.insert(temp);
	rk->id2 = (result.first)->c_str();

	//rk->id2 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id2, temp);
	break;
      case 2:
	result = stringTable.insert(temp);
	rk->chromosome = (result.first)->c_str();

	//rk->chromosome = (char *) malloc(sizeof(char) * (strlen(temp)+1));
	//strcpy(rk->chromosome, temp);
	break;
      case 3:
	rk->direction = temp[0];
	break;
      case 11:
	rk->length = atoi(temp.c_str());
	break;
      case 12:
	char *tempB = (char *) malloc(sizeof(char)* (temp.size()+1));
	strcpy(tempB,temp.c_str());
	char *temp1 = strstr(tempB, "--");
	if (temp1 == NULL) {
	  rk->position1 = rk->position2 = 0;
	}
	else {
	  temp1[0] = '\0';
	  rk->position1 = atol(tempB);
	  rk->position2 = atol(temp1+2);
	}
	free(tempB);
	break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 13) {delete rk; break;}

    data_boundaries.push_back(*rk);
    delete rk;
  }

  fclose(f);
}

bool compare_data(RSW const &aa, RSW const &bb) {
  //int temp = strcmp(aa.id, bb.id);
  int temp = aa.id-bb.id;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.direction < bb.direction) 
    return true;
  else if (bb.direction < aa.direction)
    return false;

  //temp = strcmp(aa.chromosome, bb.chromosome);
  temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  return false;
}

bool compare_dataToSort(RSW const &aa, RSW const &bb) {
  //int temp = strcmp(aa.id, bb.id);
  int temp = aa.id-bb.id;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.direction < bb.direction) 
    return true;
  else if (bb.direction < aa.direction)
    return false;

  //temp = strcmp(aa.chromosome, bb.chromosome);
  temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.position < bb.position) return true;
  else return false;
}


bool compare_data_known(RSW_Known const &aa, RSW_Known const &bb) {
  //int temp = strcmp(aa.chromosome, bb.chromosome);
  int temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  return aa.position1 < bb.position1;
}

bool compare_spliceBySupport(RSW_splice *aa, RSW_splice *bb) {
  if (bb->supported_reads.size() < aa->supported_reads.size()) return true;
  else return false;
}

// sort by chromosome, length fo splice, # supporting reads - used in deciding which splices to print
bool compare_spliceByChromLen(RSW_splice *aa, RSW_splice *bb) {
  //int temp = strcmp(aa->chromosome, bb->chromosome);
  int temp = aa->chromosome-bb->chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  temp = (aa->positionLarger - aa->positionSmaller) - 
    (bb->positionLarger - bb->positionSmaller);
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (bb->supported_reads.size() < aa->supported_reads.size()) return true;
  //if (aa->positionSmaller < bb->positionSmaller) return true;
  return false;
}

// sort by chromosome, length of splice, position - used for calculating supporting reads
bool compare_spliceByChromPos(RSW_splice *aa, RSW_splice *bb) {
  //int temp = strcmp(aa->chromosome, bb->chromosome);
  int temp = aa->chromosome-bb->chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  temp = (aa->positionLarger - aa->positionSmaller) - 
    (bb->positionLarger - bb->positionSmaller);
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa->positionSmaller < bb->positionSmaller) return true;
  return false;
}



void readOptionsFromFile(const char *filename) {
  FILE * fOptions = fopen(filename, "r");
  if (fOptions == NULL) { printf("Error opening file %s\n", filename); exit(0); }

  int pos = 0; int count = 0;
  char *fields[9]; fields[0] = &options[0];
  int ch;
  while ((ch = fgetc(fOptions)) != EOF) {
    if (pos >= MAX_STR_LEN) { printf("Options file %s is more than the max of %i bytes.\n", filename, MAX_STR_LEN); exit(0); }
    if (ch == '\n') {
      options[pos++] = '\0'; count++;
      if (count < 9)
	fields[count] = &options[pos];
      else break;
    }
    else options[pos++] = ch;
  }

  fclose(fOptions);

  if (count < 9) { printf("Not enough lines in options file.\n"); exit(0); }

  sampleDataFile = fields[0];
  maxDistance = atoi(fields[1]);
  sampleLength = atoi(fields[2]);
  refFlatFile = fields[3];
  refFlatBoundaryFile = fields[4];
  minSpliceLength = atoi(fields[5]);
  supportPosTolerance = atoi(fields[6]);
  resultsBaseName = fields[7];
  minSupportingReads = atoi(fields[8]);
}

void setDefaultOptions() {
  sampleDataFile = "RSW_test.txt";
  maxDistance = 40000;
  sampleLength = 33;
  refFlatFile = "refFlat.txt";
  refFlatBoundaryFile = "refFlat.txt.intronBoundary.exonsgaps";
  minSpliceLength = 2;
  supportPosTolerance = 5;
  resultsBaseName = "RSW_tst";
  minSupportingReads = 2;

  printf("Not enough arguments given, using default values.\n");
  printf("Usage is to load options from file: ./splitPairs optionsFile.txt \n");
  printf("And make sure options file has options in order, each on their own line with no extra lines.\n");
}

void printCurrentOptions(FILE *f) {
  fprintf(f, "Running with options...\n"
	  "  file with read data          %s\n"
	  "  max distance between matches %i\n"
	  "  length of samples            %i\n"
	  "  refFlat file                 %s\n"
	  "  refFlat intron boundary file %s\n"
	  "  minimum splice length        %i\n"
	  "  tolerance of difference in position for supporting reads  %i\n"
	  "  base of file name for writing results                     %s\n"
	  "  minimum number of supporting reads                        %i\n\n",
	  sampleDataFile, maxDistance, sampleLength, refFlatFile, refFlatBoundaryFile, minSpliceLength, supportPosTolerance, resultsBaseName, minSupportingReads);

  if (strlen(sampleDataFile) > MAX_STR_LEN - 100) {
    fprintf(f,"Error, filename %s is too long.\n", sampleDataFile); exit(0);
  }
}

FILE * fKnown, *fUnknown, *fKnownFull, *fUnknownFull, *fSplitPairs;

void openOutputFiles() {
  sprintf(buff,"%s.results", resultsBaseName);
  fKnown = fopen(buff, "w");
  if (fKnown == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write summary results that match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.full", resultsBaseName);
  fKnownFull = fopen(buff, "w");
  if (fKnownFull == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write full results that match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.unknown", resultsBaseName);
  fUnknown = fopen(buff, "w");
  if (fUnknown == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write summary results that do NOT match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.unknown.full", resultsBaseName);
  fUnknownFull = fopen(buff, "w");
  if (fUnknownFull == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write full results that do NOT match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.splitPairs", resultsBaseName);
  fSplitPairs = fopen(buff, "w");
  if (fSplitPairs == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write split pairs to file\n"
	 "   %s\n", buff);
}

void printSplice(FILE *f, RSW_splice *sp) {
  fprintf(f,
	  "%s\t%s\t%li\t%li\t%li--%li\t%s",
	  sp->geneName, sp->chromosome, 
	  sp->supported_reads.size(),
	  sp->positionLarger-sp->positionSmaller, 
	  sp->minSmallSupport,sp->maxLargeSupport,
	  //sp->positionSmaller, sp->positionLarger, 
	  sp->novel ? "Novel" : "*"
	  );
}


void updateSpliceSupport(RSW_splice *splice1, RSW_splice *splice2) {

  for(int i=0; i < 2; i++) {
    RSW_splice * sp1, *sp2;
    if (i == 0) {sp1 = splice1; sp2 = splice2;}
    else {sp1 = splice2; sp2 = splice1;}

    //#pragma omp critical
    {
      if (sp1->supported_reads.size() <= 0 ||
	  sp2->positionSmaller < sp1->minSmallSupport)
	sp1->minSmallSupport = sp2->positionSmaller;
      if (sp1->supported_reads.size() <= 0 ||
	  sp2->positionLarger > sp1->maxLargeSupport)
	sp1->maxLargeSupport = sp2->positionLarger;
      sp1->supported_reads.insert({sp2->id, sp2});
    }
  }
}

void printStats(FILE *f) {
  FILE *fStatus = fopen("/proc/self/status","r");
  char s[1000], mem[20]="", units[20]="";
  while (fscanf(fStatus,"%999s",s) == 1) {
    if (strcmp(s,"VmRSS:") == 0) {
      fscanf(fStatus,"%19s %19s",mem, units);
      break;
    }
  }
  fclose(fStatus);

  endTime = time(NULL);
  fprintf(f, "Finished processing data, results written to files.\n");
  fprintf(f, "Number of entries in data file:             %i\n", data.size());
  fprintf(f, "Number of different reads:                  %i\n", numDifferentReads);
  fprintf(f, "Number of entries in refFlat file:          %i\n", data_known.size());
  fprintf(f, "Number of entries in refFlat boundary file: %i\n", data_boundaries.size());
  fprintf(f, "Number of matches:                          %i\n", data_splice.size());
  fprintf(f, "String table size:                          %i\n", stringTable.size());
  fprintf(f, "VmRSS, memory resident set size:            %s %s\n", mem, units);
  fprintf(f, "Total time to process:                      %i seconds\n", endTime-beginTime);
  fprintf(f, "\n");
}

int main(int argc, char *argv[]) {
  setpriority(0, 0, 20); // so other processes get priority over this one

  beginTime = time(NULL);

  if (argc > 1) 
    readOptionsFromFile(argv[1]);
  else 
    setDefaultOptions();

  openOutputFiles();

  printCurrentOptions(stdout);
  printCurrentOptions(fKnown);
  printCurrentOptions(fKnownFull);
  printCurrentOptions(fUnknown);
  printCurrentOptions(fUnknownFull);
  printCurrentOptions(fSplitPairs);


  // read from refFlat file into data_known array, 
  read_knownGene(refFlatFile);
  sort(data_known.begin(), data_known.end(), compare_data_known);
  printf("Done reading/sorting refFlat, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read from refFlat boundary file into data_boundaries array, 
  read_boundaries(refFlatBoundaryFile);
  printf("Done reading refFlat intron/exon boundaries, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read the read data
  read_data(sampleDataFile);
  printf("Done reading read data, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort the read data
  sort(data.begin(), data.end(), compare_dataToSort);
  printf("Done sorting read data, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // look at all pairs of read segments, looking for matches
  numDifferentReads = 0;
  const int inputSize = data.size();
  //#pragma omp parallel for
  for(int left=0; left < inputSize; left++) {
    int right;
    if (left == 0 || //strcmp(data[left].id,data[left-1].id) != 0) 
	(data[left].id  != data[left-1].id) ) {
      //#pragma omp atomic
      numDifferentReads++;
    }

    for(right=left+1; right < inputSize; right++) {
      // read segments are ordered by id/chromosome/strand, so if 
      // there isn't a match we can skip the rest of the read segments
      // for "right", and go to the next iteration of the "left" loop
      if (data[left].id != data[right].id ||
	  data[left].direction != data[right].direction ||
	  data[left].chromosome != data[right].chromosome) break;

      if (data[right].position - data[left].position > maxDistance) break;

      // want it to be from two sides of the same segment
      if (data[left].side == data[right].side) continue;

      // total length should be correct
      if (data[left].length + data[right].length != sampleLength) continue;

      // calculate the end of the segments, since what is given
      // in the data is the beginning of the segments
      int endSmaller, endLarger; // splice is between endSmaller and endLarger
      int first, second;
      if (data[left].side == 'L' && data[left].direction == '+' ||
	  data[left].side == 'R' && data[left].direction == '-') { first = left; second = right; }
      else { first = right; second = left; }
      endSmaller = data[first].position + data[first].length;
      endLarger = data[second].position;

      // splice length, and check that it is within specified bounds
      int spliceLength = endLarger - endSmaller;
      if (spliceLength > maxDistance) continue;
      if (spliceLength < minSpliceLength) continue;

      // check if we already have this splice from this read...
      int i;
      for(i=0; i < data_splice.size(); i++) {
	if (data_splice[i]->positionSmaller != endSmaller ||
	    data_splice[i]->positionLarger != endLarger ||
	    data_splice[i]->id != data[left].id ||
	    //strcmp(data_splice[i]->id, data[left].id) != 0 ||
	    data_splice[i]->chromosome != data[left].chromosome) 
	  //strcmp(data_splice[i]->chromosome, data[left].chromosome) != 0) 
	  continue;
	break;
      }
      // if already have this exact splice for this chromosome from this read, don't include it again.
      if (i < data_splice.size()) continue;

      // note: could print this match here, step 5 done.

      // look for this in the known gene...
      int k; int foundInGene = 0;
      for(k=0; k < data_known.size(); k++) {
	if ((((data[left].chromosome == data_known[k].chromosome) &&
	      //if ((((strcmp(data[left].chromosome, data_known[k].chromosome) == 0) &&
	      (data[left].position >= data_known[k].position1 &&
	       data[right].position >= data_known[k].position1) &&
	      (data[left].position <= data_known[k].position2 &&
	       data[right].position <= data_known[k].position2)))) {
	  if (foundInGene) ;//printf("DUPLICATE_"); // duplicate
	  foundInGene = 1;
	  break; // just cut off search, don't look for duplicates
	}
      }
      int geneIndex = k;
      //if (!found) printf("UNFOUND_\t"); // not found in any known gene.
      //else printf("%s\t", data_known[k].id1);

      // add this splice match to the data_splice vector.
      RSW_splice *sp = new RSW_splice;
      sp;
      if (foundInGene) {
	sp->geneName = data_known[geneIndex].id1;
	sp->geneUnknown = 0;
      }	 
      else {
	sp->geneName = unfound_string;
	sp->geneUnknown = 1;
      }
      sp->id = data[left].id;
      sp->chromosome = data[left].chromosome;
      sp->direction = data[left].direction;
      sp->positionSmaller = sp->minSmallSupport = endSmaller;
      sp->positionLarger = sp->maxLargeSupport = endLarger;
      //sp->reported = false;
      sp->supported_reads.insert({sp->id, sp});
      //#pragma omp critical
      {
	data_splice.push_back(sp);
      }
    }
  }

  printf("Done finding matched pairs, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  //fclose(fKnown); fclose(fKnownFull); fclose(fUnknown); fclose(fUnknownFull); fclose(fSplitPairs);
  //exit(0);


  // sort splices by chromosome and length
  sort(data_splice.begin(), data_splice.end(), compare_spliceByChromPos);

  printf("Done sorting matched pairs, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);


  // compute supporting reads

  //#pragma omp parallel for
  for(int sp1=0; sp1 < data_splice.size(); sp1++) {
    int sp2;
    for(sp2=sp1; sp2 < data_splice.size(); sp2++) {
      // see if these reads support each other

      // because splices are sorted by chromosome and length, if these don't match then 
      // can break sp2 loop.
      if (data_splice[sp1]->chromosome != data_splice[sp2]->chromosome) break;
      //if (strcmp(data_splice[sp1]->chromosome, data_splice[sp2]->chromosome) != 0) break;
      if (abs(data_splice[sp1]->positionLarger - data_splice[sp1]->positionSmaller) != 
	  abs(data_splice[sp2]->positionLarger - data_splice[sp2]->positionSmaller)) break; // splice length must be the same
      if (data_splice[sp1]->geneName != data_splice[sp2]->geneName) break; // ???
      //if (strcmp(data_splice[sp1]->geneName, data_splice[sp2]->geneName) != 0) break;

      // difference in where the splice is, 0 for exact match.  if doesn't match, 
      // continue to next sp2.
      if (abs(data_splice[sp1]->positionSmaller-data_splice[sp2]->positionSmaller) > supportPosTolerance) 
	continue;
	//break;
      if (data_splice[sp1]->id == data_splice[sp2]->id) continue; // read can't support itself

      // they are matches for each other, so add each to each other's list, as 
      // long as not a match already there for the same read id.
      updateSpliceSupport(data_splice[sp1], data_splice[sp2]);
    }
  }

  printf("Done computing supporting reads, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort splices by chromosome and length, and # supporting reads
  sort(data_splice.begin(), data_splice.end(), compare_spliceByChromLen);

  printf("Done sorting matched pairs again, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // set anything that is included within another splice region, or the same splice region,
  // as already reported.  if saying something is reported, take one with larger supporting reads as
  // one to print.
  /*  for(sp1=0; sp1 < data_splice.size(); sp1++) {
    if (data_splice[sp1]->reported) continue; // already reported sp1, so go to next.

    for(sp2=sp1+1; sp2 < data_splice.size(); sp2++) {
      // sorted by chromosome and length, if those don't match then not supporting each other, so 
      // go to next sp1.
      if (data_splice[sp1]->chromosome != data_splice[sp2]->chromosome) break;
      //if (strcmp(data_splice[sp1]->chromosome, data_splice[sp2]->chromosome) != 0) break;
      if (abs(data_splice[sp1]->positionLarger - data_splice[sp1]->positionSmaller) != 
	  abs(data_splice[sp2]->positionLarger - data_splice[sp2]->positionSmaller)) break; // splice length must be the sam

      if (data_splice[sp2]->reported) continue; // already not reporting sp2, so go to next.

      int m1 = data_splice[sp1]->minSmallSupport - data_splice[sp2]->minSmallSupport;
      int m2 = data_splice[sp1]->maxLargeSupport - data_splice[sp2]->maxLargeSupport;
      if (m1 <= 0 && m2 >= 0) {// splice sp2 region included in splice sp1 region
	// and # reads of sp2 should be <= # reads of sp1 since sorted by that ordering
	// so choose not to report sp2.
	data_splice[sp2]->reported = true;
      }
    }
    }*/

  printf("Done filtering matched pairs, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  endTime = time(NULL);


  printf("Done calculating/filtering supporting reads, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // print the splice results.  start with the statistics
  printStats(stdout);
  printStats(fKnown);
  printStats(fKnownFull);
  printStats(fUnknown);
  printStats(fUnknownFull);
  printStats(fSplitPairs);

  int k;


  // save all the splices, and reads that support each one.
  fprintf(fSplitPairs, "Line No\tId\tGene\tChr\tStrand\tSplice region\tSupporting reads\tSupporting splice range\n");
  for(k=0; k < data_splice.size(); k++) {
    fprintf(fSplitPairs, "%s\t%s\t%s\t%c\t%i-%i\t%i\t%i-%i\n",
	    data_splice[k]->id, data_splice[k]->geneName, 
	    data_splice[k]->chromosome, data_splice[k]->direction,
	    data_splice[k]->positionSmaller,data_splice[k]->positionLarger,
	    data_splice[k]->supported_reads.size(),
	    data_splice[k]->minSmallSupport,data_splice[k]->maxLargeSupport);

    for(auto &x: data_splice[k]->supported_reads)
      fprintf(fSplitPairs, "%s, ", x.second->id);
    fprintf(fSplitPairs, "\n");
  }

  // save the tabulated results.
  fprintf(fKnown, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  fprintf(fUnknown, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  fprintf(fKnownFull, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  fprintf(fUnknownFull, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  FILE * f, *fFull;
  for(k=0; k < data_splice.size(); k++) {
    if (readIdsReported.find(data_splice[k]->id) != readIdsReported.end()) {continue;}

    /*    if (data_splice[k]->reported) {continue;}  // if already reported as supporting another read, don't print it.
	  data_splice[k]->reported = true;*/

    // if doesn't have enough supporting reads, don't print it
    if (data_splice[k]->supported_reads.size() < minSupportingReads) continue;

    readIdsReported.insert(data_splice[k]->id);

    // check if novel or not.
    int j;
    data_splice[k]->novel = true;
    for(j=0; j < data_boundaries.size(); j++) {
      if (data_splice[k]->positionLarger-data_splice[k]->positionSmaller != 
	  data_boundaries[j].length) continue;
      if (abs(data_splice[k]->minSmallSupport-data_boundaries[j].position1) <= supportPosTolerance &&
	  abs(data_splice[k]->maxLargeSupport-data_boundaries[j].position2) <= supportPosTolerance)
	break;
    }
    if (j < data_boundaries.size()) 
      data_splice[k]->novel = false;

    // going into known file or unknown
    if (data_splice[k]->geneUnknown) {
      f = fUnknown; fFull = fUnknownFull;
    }
    else {
      f = fKnown; fFull = fKnownFull;
    }

    // print out results
    printSplice(f, data_splice[k]);
    fprintf(f, "\n");

    // print out full results, that includes id's of supporting reads
    printSplice(fFull, data_splice[k]);

    fprintf(fFull, "\t");
    for(auto& x: data_splice[k]->supported_reads) {
      fprintf(fFull,", %s", x.second->id);
      //x.second->reported = true; // already reported this splice as supporting another read that has higher support, so don't report it again
      readIdsReported.insert(x.second->id);
    }
    fprintf(fFull, "\n");
  }

  fclose(fKnown); fclose(fKnownFull); fclose(fUnknown); fclose(fUnknownFull); fclose(fSplitPairs);

  printf("Done saving results, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // free memory.  good to do so we can run a memory checker and verify
  // we don't have any memory leaks.
  while (data_splice.size() > 0) {
    delete data_splice.back(); data_splice.pop_back();
  }
  
  return 0;
}
                                                                     src/split_columns.c                                                                                 0000664 0001753 0001753 00000004573 12563404716 014537  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>


#define SUFFIX ".split1stcolumn"

//D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:-R-68/2	

void split_field(char *readName, int len) {
    int i, offset;

    if (!readName || strlen(readName) <= len) {
        fprintf(stderr,"No field supplied to split_field. aborting\n");
        exit(1);
    }
//actual output
    if (readName[len-2] == '/' && (readName[len-1] == '1' || readName[len-1] == '2')) { //paired
        //before:
        //D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:-R-68/2	
        //let's do some horse-trading...
        //put the /1 or /2 at the beginning of the substring...
        offset = 7;
        char t1,t2;
        t1 = readName[len - 2];
        t2 = readName[len - 1];
        //shift everything right
        for (i = (len - 3);i > (len - offset);i--) {
             readName[i+2] = readName[i];
        }
        readName[i] = t1;
        readName[i+1] = t2;
    }
    //after, OR, the /2 is irrelevant
    //D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:/2-R-68
    i = len - 5;
    readName[i] = '\t';
    readName[i+2] = '\t';
}

void parse(char *filename) {
    char *line = 0;
    char *readName,*tmp;
    size_t n;
    int pos;
    FILE *in, *out;
    char *outname;

    
    if (!(in = fopen(filename,"r"))) {
        fprintf(stderr, "Could not open %s for reading", filename);
        perror("");
        return;
    }
    
    outname = (char *)malloc(sizeof(char) * (strlen(filename)+strlen(SUFFIX)));
    memset(outname, 0, sizeof(char) * (strlen(filename)+strlen(SUFFIX)));
    sprintf(outname, "%s%s",filename,SUFFIX);
    outname[strlen(filename)+strlen(SUFFIX)] = 0;
    if (!(out = fopen(outname, "w"))) {
        fprintf(stderr,"Could not open %s for writing", outname);
        perror("");
        close(in);
        return;
    }

    tmp = 0;
    while((!feof(in)) && (getline(&line,&n,in) != -1)) {
        if (tmp) { free(tmp); tmp=0; }
        tmp = strdup(line);
        char *found = strchr(line, '\t');
        pos = (found - line);
        split_field(tmp, pos);
        fputs(tmp,out);
    }

    fclose(in);
    fclose(out);
}


int main(int argc, char *argv[]) {
    int i;

    if (argc < 2) {
        fprintf(stderr, "usage -- %s <file to split> [additonal files...]\n",argv[0]);
        return 1;
    }

    for (i = 1;i < argc;i++)
        parse(argv[i]);

    return 0;
}
                                                                                                                                     src/split_on_chrom.c                                                                                0000664 0001753 0001753 00000004342 12563405032 014645  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define CHR "Chr"
#define SUFFIX ".txt"
#define MAX_CHROMOSOMES 308


FILE *chrFiles[MAX_CHROMOSOMES + 2];
int openFiles = 0;
char *prefix;
const char *const file310 = "Chr_unknown";

FILE *pick_file(char *chr) {
    int id;
    if (!chr) id = MAX_CHROMOSOMES+1;
    else id = atoi((chr+3));
    if (0 == id) id = MAX_CHROMOSOMES+1;
    //fprintf(stderr,"pick_file:: guessing file %i for %s\n",id,chr+3);
    if (!chrFiles[id]) {
        char *fn;
        if (MAX_CHROMOSOMES+1 == id) {
            fn = (char *)malloc(strlen(file310)+strlen(prefix)+strlen(SUFFIX)+3);
            sprintf(fn,"%s.%s%s",prefix,file310,SUFFIX);
        }
        else {
            fn = (char *)malloc(strlen(chr)+strlen(prefix)+strlen(SUFFIX)+3);
            sprintf(fn,"%s.%s%s",prefix,chr,SUFFIX);
        }
        if (!(chrFiles[id] = fopen(fn,"w"))) {
            fprintf(stderr,"Panic! Cannot open chrom file %s for writing!\n");
            perror("");
            fprintf(stderr,"Ignoring.");
            return 0;
        }
    }
    return chrFiles[id];
}

char *get_chrom(char *line) {
    char *tok;
    tok = strtok(line,"\t");
    while(tok && strncmp(tok,CHR,strlen(CHR)) != 0)
        tok = strtok(NULL,"\t");
    return tok;
}

void parse(FILE *in) {
    char *line = 0;
    char *tmp, *chrom;
    size_t n;
    FILE *outfile;

    tmp = 0;
    while(!(feof(in)) && (getline(&line,&n,in) != -1)) {
        if (tmp) { free(tmp); tmp=0; }
        tmp = strdup(line);
        chrom = get_chrom(tmp);
        outfile = pick_file(chrom);
        if (!outfile) continue;
        fprintf(outfile,"%s",line);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <reads file>\n", argv[0]);
        return 1;
    }
    FILE *fr;
    clock_t after, before = clock();

    memset(chrFiles, 0, MAX_CHROMOSOMES * sizeof(FILE *));

    prefix = strdup(argv[1]);
    if (!(fr = fopen(argv[1], "r"))) {
        fprintf(stderr, "Could not open %s for reading", argv[1]);
        perror("");
        return 1;
    }
    
    parse(fr);
    fcloseall();
    after = clock();
    printf("Elapsed time: %lf seconds\n", (after-before)*(1.0/CLOCKS_PER_SEC));
    return 0;
}
                                                                                                                                                                                                                                                                                              srr                                                                                                 0000775 0001753 0001753 00000032571 12563406335 011443  0                                                                                                    ustar   bdonham                         bdonham                                                                                                                                                                                                                ELF          >    L@     @       ("          @ 8 	 @         @       @ @     @ @     �      �                   8      8@     8@                                          @       @     d      d                          `     `     �      �                    (      (`     (`     �      �                   T      T@     T@     D       D              P�td   �      �@     �@     <       <              Q�td                                                  R�td         `     `     �      �             /lib64/ld-linux-x86-64.so.2          GNU                        GNU aQ����F��-��!�����                         9�                            �                      '                      /                      O                      �                      A                      H                      {                      k                      <                      �                       5                                            V                                                                   s                      ]                      d     � `             libc.so.6 __xpg_basename fopen perror strncpy clock strtol feof strlen memset fclose malloc strcat stderr fprintf getline __libc_start_main stpcpy free __gmon_start__ GLIBC_2.2.5                                       ui	   �       �`                   � `                    `                     `                   ( `                   0 `                   8 `                   @ `                   H `                   P `                   X `        	           ` `        
           h `                   p `                   x `                   � `                   � `                   � `                   � `                   � `                   H��H�m  H��t�   H���      �5b  �%d  @ �%b  h    ������%Z  h   ������%R  h   ������%J  h   �����%B  h   �����%:  h   �����%2  h   �����%*  h   �p����%"  h   �`����%  h	   �P����%  h
   �@����%
  h   �0����%  h   � ����%�  h   �����%�  h   � ����%�  h   ������%�  h   ������%�  h   �����AWAVAUATU��SH��H������  ����������H�=�  I���E1�I��L��D����H��J�1�  H�kD��L��H���H��D�t�Mc�L���?���L��1�H��I������L�5h  H��L���]���H��L���R���H�{1��
   I�������H�{1��
   �J  �������+=  H�{�|@ �&  �$  �����H��H���)  L�羜@ ����H��I���;  H��H���W  H�������L������������L)�H�=�  ��@ �H*��   �Y�  ����1�H��[]A\A]A^A_�H�s H�$L��H���H�{H��L�Hc�D�y�|�/�����D$DD�����L��H��H��D��H�$�D��)�D�t�Mc�L�������1�L��H��I���l���E�������H�s L���7����|$ �����f� / ����H�H�=�  �0@ 1��L����   �8���H�SH�=�  �~@ 1��+�����@ �����   ����H�SH�=�  �~@ 1�� �����@ �f���H�������   �����1�I��^H��H���PTI���@ H��`@ H���@ �����f��     �� ` UH-� ` H��H��w]ø    H��t�]�� ` ���    �� ` UH-� ` H��H��H��H��?H�H��u]ú    H��t�]H�ƿ� ` ���    �=�   uUH���~���]��  ��@ H�=   t�    H��tU� ` H����]�{��� �s��� UH��AWAVAUATSH��H��hH��x���H�E�    H��   H�E�    H�E�    H�E�    H�E�    I��Hǅp���    H�E�    �������;   H��x���H�u�H�}�D�-  D�5  �����H����  L�}�L������I�D��8
�*  H��x���H�u�H�}�����H�����  L�}�L�������I�D��8
�  H��x���H�u�H�}��n���H�����  L�}�L������I�D��8
��  H��x���H�u�H�}��6���H����d  L�}�L���p���I�D��8
��  H��p���E9�H�E�    ��  f�     H�U�D��`@ H��1�H�E�E���R���H�u�L��L�������L��i@ H��1�C�< �,���H�U�D��`@ H��1�����H�u�L��L������L��i@ H��1�C�< �����H�U�D��m@ H��1������L��HM�D��v@ H��1������H�U�D��m@ H��1�����L��HM�D��v@ H��1�E������H�U�D��`@ H��1��y���H�u�L��L�������L��i@ H��1�C�< �S���H�U�D��`@ H��1��=���H�u�L��L������L��i@ H��1�C�< ����H�U�D��m@ H��1�����L��HM�D��v@ H��1�A�������H�U�D��m@ H��1������L��HM�D��1��v@ H��A������H�E�E9��3���H�}�����H�}�H�E�    �	���H�}�H�E�    �����H�}�H�E�    �����H��x���H�E�    �c���������� L�E�H��p�����@ H�M�H�=�  1�M��L������H�e�[A\A]A^A_]��  ������     �  ������     �  �����     �  �F����     AWA��AVI��AUI��ATL�%�  UH�-�  SL)�1�H��H�������H��t�     L��L��D��A��H��H9�u�H��[]A\A]A^A_�ff.�     ��f�H��H���                 in: %lu, out: %lu, passes each line: %lu, total passes: %lu
    Usage; %s <file> <min size> <total size>
       %s-L-%d
 %s
 %s-R-%d
 %.*s
 r Unable to open %s for reading w time elapsed: %lf seconds.
 .split       �����ư>;<      �����    ����   |���X   p����   ����0   ���x             zR x�      ���*                  zR x�  $      @���0   FJw� ?;*3$"    ,   D   ����   A�CI������
A   L   t   ���|   B�B�B �B(�A0�C8�GPG
8A0A(B BBBA    D   �   X���e    B�E�E �E(�H0�H8�M@l8A0A(B BBB      ����                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           @     �@                                  �@            �@            `                          `                   ���o    �@            �@            �@     
       �                                             `            �                           �@            �@            0       	              ���o    �@     ���o           ���o    T@                                                                                                             (`                     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@                             �@     GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-7) GCC: (GNU) 4.8.3 20140911 (Red Hat 4.8.3-9)  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.gnu.build-id .gnu.hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .init_array .fini_array .jcr .dynamic .got .got.plt .data .bss .comment                                                                              8@     8                                    #             T@     T                                     1             t@     t      $                              D   ���o       �@     �      $                             N             �@     �      �                          V             �@     �      �                              ^   ���o       T@     T      (                            k   ���o       �@     �                                   z             �@     �      0                            �             �@     �      �                          �             �@     �                                    �             �@     �      0                            �             �@     �                                   �             �@     �      	                              �             �@     �      �                              �             �@     �      <                              �             @           T                             �             `                                         �             `                                         �              `                                          �             (`     (      �                           �             �`     �                                   �               `             �                             �             � `     �                                     �             � `     �                                     �      0               �       X                                                    !                                                         �)      (         -                 	                      �1      �                                                           8@                   T@                   t@                   �@                   �@                   �@                   T@                   �@                  	 �@                  
 �@                   �@                   �@                   �@                   �@                   �@                   �@                   @                   `                   `                    `                   (`                   �`                     `                   � `                   � `                                       ��                    ��                      `             *     �@             ?     �@             R     �@             h     � `            w     `             �     @             �     `                 ��                �     `@             �      `                  ��                �      `             �     (`             �      `                   `             &    �@            6                     H    � `            P                     e                      �     � `             �    � `             �                     �                     �                     �    � `            �    �@             �                     �    � `            �    � `            �                                          -    � `             :                     O                     a                      p                     �   �@             �    �@            �                     �    `@     e       �                     �    � `             �    �@            �    L@             �    � `             �    �@     |      �                                          &                      :                     O                     c    @@           i   � `             u                      �    �@             �    � `             split_read_rsw.c crtstuff.c __JCR_LIST__ deregister_tm_clones register_tm_clones __do_global_dtors_aux completed.6342 __do_global_dtors_aux_fini_array_entry frame_dummy __frame_dummy_init_array_entry __FRAME_END__ __JCR_END__ __init_array_end _DYNAMIC __init_array_start _GLOBAL_OFFSET_TABLE_ __libc_csu_fini free@@GLIBC_2.2.5 READLEN strncpy@@GLIBC_2.2.5 _ITM_deregisterTMCloneTable data_start _edata clock@@GLIBC_2.2.5 fclose@@GLIBC_2.2.5 stpcpy@@GLIBC_2.2.5 RIGHT _fini strlen@@GLIBC_2.2.5 SUFFIX LEFT memset@@GLIBC_2.2.5 __libc_start_main@@GLIBC_2.2.5 __data_start fprintf@@GLIBC_2.2.5 feof@@GLIBC_2.2.5 __gmon_start__ strtol@@GLIBC_2.2.5 __dso_handle _IO_stdin_used __xpg_basename@@GLIBC_2.2.5 __libc_csu_init malloc@@GLIBC_2.2.5 _end K _start __bss_start main fopen@@GLIBC_2.2.5 perror@@GLIBC_2.2.5 _Jv_RegisterClasses getline@@GLIBC_2.2.5 strcat@@GLIBC_2.2.5 parse __TMC_END__ _ITM_registerTMCloneTable _init stderr@@GLIBC_2.2.5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        