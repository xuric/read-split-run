#!/bin/bash
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

