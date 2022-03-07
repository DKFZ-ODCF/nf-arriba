#!/bin/bash
#
# Copyright (c) 2021 DKFZ.
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/nf-arriba/blob/master/LICENSE).
#

# enable logging and abort on error
set -x -e -u -o pipefail

# convert SOPHIA SVs into a format that Arriba understands
if [ -n "$SOPHIA_SVS" ]; then
	convert_sophia_svs.awk "$SOPHIA_SVS"
fi > "structural_variants.tsv"

# determine how to extract input files for STAR (BAM/gzipped FastQ)
if [ -n "$FASTQ1" ]; then
	READ_FILES_TYPE="Fastx"
	READ_FILES_COMMAND="zcat"
elif [ -n "$BAM" ]; then
	# auto-detect layout (single-end vs. paired-end)
	LAYOUT=$(samtools view "$BAM" | head -n1 | awk '{print ($2 % 2) ? "PE" : "SE"}' || true)
	READ_FILES_TYPE="SAM $LAYOUT"
	READ_FILES_COMMAND="bam_file_reader.sh $LAYOUT"
fi

# run Arriba in parallel to STAR using named pipe
mkfifo arriba_pipe
arriba \
	-x arriba_pipe \
	-a "$ASSEMBLY" -g "$ANNOTATION" \
	-b "$BLACKLIST" -k "$KNOWN_FUSIONS" -t "$KNOWN_FUSIONS" -p "$PROTEIN_DOMAINS" \
	-d "structural_variants.tsv" -f "no_genomic_support,genomic_support" \
	-o "fusions.tsv" -O "fusions.discarded.tsv" & ARRIBA_PID=$!

# alignment
STAR \
	--runThreadN "$THREADS" \
	--genomeDir "$STAR_INDEX" --genomeLoad NoSharedMemory \
	--readFilesIn "$FASTQ1" "$FASTQ2" "$BAM" --readFilesType $READ_FILES_TYPE --readFilesCommand "$READ_FILES_COMMAND" \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 \
	--alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0.5 \
	--chimSegmentMin 10 --chimScoreDropMax 30 --chimScoreSeparation 1 --chimScoreJunctionNonGTAG 0 \
	--chimJunctionOverhangMin 10 --chimOutType WithinBAM HardClip --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
tee arriba_pipe |

# create a reduced BAM file lacking sequence and quality information,
# because we only need it to draw the coverage track in the fusion plots
samtools view -h - |
awk -F '\t' '{
	print # pass on full SAM record to next command
	if ($1 !~ /^@/) # this is not a SAM header line
		$0 = "0\t0\t" $3 "\t" $4 "\t0\t" $6 "\t*\t0\t0\t*\t*" # reduce SAM record to essential info
	print > "reduced.sam"
}' |

# quantify expression of viral genomes
quantify_virus_expression.sh /dev/stdin /dev/stdout > "virus_expression.tsv"

# draw fusion plots
samtools sort -m "$((SORT_MEMORY/THREADS))M" -l 1 -@ "$THREADS" "reduced.sam" |
tee "coverage.bam" | samtools index -@ "$THREADS" - "coverage.bam.bai"
if ! wait $ARRIBA_PID; then # wait for fusions.tsv
	echo "Arriba process exited with error" >&2
	exit 1
fi
draw_fusions.R \
	--fusions="fusions.tsv" --output="fusions.pdf" --alignments="coverage.bam" \
	--annotation="$ANNOTATION" --cytobands="$CYTOBANDS" --proteinDomains="$PROTEIN_DOMAINS"

# save some disk space
gzip -9 "fusions.discarded.tsv"

