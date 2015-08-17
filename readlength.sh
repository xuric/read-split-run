#!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsr_config.sh"

echo $(( $(head -2 $1 | tail -1 | wc -c) - $OFFSET ))
