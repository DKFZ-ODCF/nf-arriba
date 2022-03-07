#!/bin/bash
#
# Copyright (c) 2021 DKFZ.
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/nf-arriba/blob/master/LICENSE).
#

# This script is to be used as an argument to the STAR parameter `--readFilesCommand`.
# It converts paired-end/single-end BAM to SAM.

# abort on error
set -e -u -o pipefail

if [ $# -ne 2 ]; then
	echo "Usage: $(basename $0) SE|PE BAM_FILE" >&2
	exit 1
fi

# read command-line arguments
LAYOUT="$1"
BAM="$2"

if [ "$LAYOUT" = "PE" ]; then
	samtools collate -u -f -r 1000000 -O "$BAM" | # STAR expects BAM file to be collated
	samtools view - |
	awk -F '\t' '$1 == name && $2 != flag { print prev; print $0 } { name = $1; flag = $2; prev = $0 }' # skip singletons and invalid flags
elif [ "$LAYOUT" = "SE" ]; then
	samtools view -F 2304 "$BAM" # discard supplementary/secondary alignments
fi |
cut -f 1-11 # remove SAM tags, because STAR would keep them

