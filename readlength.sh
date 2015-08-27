#!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsr_config.sh"

if (( $# < 1 )); then
    die "usage  -- $0 <reads file>"
else
    echo $(( $(head -2 $1 | tail -1 | wc -c) - $OFFSET ))
fi
