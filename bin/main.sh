#!/bin/bash
#
# Copyright (c) 2022 DKFZ.
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/nf-arriba/blob/master/LICENSE).
#

# enable logging and abort on error
set -x -e -u -o pipefail

# convert SOPHIA SVs into a format that Arriba understands
if [ -n "$SOPHIA_SVS" ]; then
	convert_sophia_svs.awk "$SOPHIA_SVS"
fi > structural_variants.tsv

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

# alignment
STAR \
	--runThreadN "$THREADS" \
	--genomeDir "$STAR_INDEX" --genomeLoad NoSharedMemory \
	--readFilesIn "$FASTQ1" "$FASTQ2" "$BAM" --readFilesType $READ_FILES_TYPE --readFilesCommand "$READ_FILES_COMMAND" \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outBAMcompression 1 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 \
	--alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0.5 \
	--chimSegmentMin 10 --chimScoreDropMax 30 --chimScoreSeparation 1 --chimScoreJunctionNonGTAG 0 \
	--chimJunctionOverhangMin 10 --chimOutType WithinBAM HardClip --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
tee alignments.bam |

# fusion calling
arriba \
	-x /dev/stdin \
	-a "$ASSEMBLY" -g "$ANNOTATION" \
	-b "$BLACKLIST" -k "$KNOWN_FUSIONS" -t "$KNOWN_FUSIONS" -p "$PROTEIN_DOMAINS" \
	-d structural_variants.tsv -f "no_genomic_support,genomic_support" \
	-o fusions.tsv -O fusions_discarded.tsv

# sort and index alignments
samtools sort -@ "$THREADS" -m "$((SORT_MEMORY/THREADS))M" -l 1 --write-index -o sorted.bam##idx##sorted.bam.bai alignments.bam
rm -f alignments.bam

# quantify expression of viral genomes
quantify_virus_expression.sh sorted.bam virus_expression.tsv

# find all viruses with mapped reads and extract their alignments
samtools idxstats sorted.bam |
awk -v contigs="$VIRAL_CONTIGS" '$3 > 0 && match($1, contigs) {print $1 "\t0\t" $2}' |
samtools view -@ "$THREADS" -b -M -L /dev/stdin -O "BAM,level=9" --write-index -o virus_alignments.bam##idx##virus_alignments.bam.bai sorted.bam

# extract fusion-supporting alignments
extract_fusion-supporting_alignments.sh fusions.tsv sorted.bam extracted
samtools merge -@ "$THREADS" -O SAM - extracted*.bam |
awk '!duplicate[$0]++' |
samtools view -@ "$THREADS" -O "BAM,level=9" --write-index -o fusions_alignments.bam##idx##fusions_alignments.bam.bai
rm -f extracted*.bam*

# draw fusion plots
draw_fusions.R \
	--fusions=fusions.tsv --output=fusions.pdf --alignments=sorted.bam \
	--annotation="$ANNOTATION" --cytobands="$CYTOBANDS" --proteinDomains="$PROTEIN_DOMAINS"

# final cleanup
gzip -9 fusions_discarded.tsv
rm -f sorted.bam*

