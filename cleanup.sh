#!/bin/bash

BASEDIR=$(dirname $(realpath $0))
source "${BASEDIR}/rsr_config.sh"

FIND_MODIFIERS="-L"
FIND_WHAT="-daystart -ignore_readdir_race -mtime +14 "
ACTION="-delete"

find $FIND_MODIFIERS $TMPDIR $FIND_WHAT
find $FIND_MODIFIERS $TEMP_BOWTIE_FILES $FIND_WHAT
find $FIND_MODIFIERS $TEMP_SPLIT_FILES $FIND_WHAT
find $FIND_MODIFIERS $TEMP_RSR_FILES $FIND_WHAT

