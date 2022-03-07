/**
 *  Copyright (c) 2022 DKFZ.
 *
 *  Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/nf-arriba/blob/master/LICENSE).
 *
 *  Author: Sebastian Uhrig
 */


/** Check whether parameters are correct (names and values)
 *
 * @param parameters
 * @param allowedParameters
 */
void checkParameters(parameters, List<String> allowedParameters) {
	Set<String> unknownParameters = parameters.
			keySet().
			grep {
			  !it.contains('-') // Nextflow creates hyphenated versions of camel-cased parameters.
			}.
			minus(allowedParameters)
	if (!unknownParameters.empty) {
		log.error("There are unrecognized parameters: ${unknownParameters}")
		exit(1)
	}
}
allowedParameters = ['arribaVersion','threads','memory','fastq1','fastq2','bam','outputDir','sophiaSVs','preset','dataDir','starIndex','assembly','annotation','blacklist','knownFusions','proteinDomains','cytobands','viralContigs','ignoreViruses']
checkParameters(params, allowedParameters)


process arriba {
	cpus params.threads
	memory { 1.MB * params.memory }
	time { 48.hours * (2**task.attempt - 1) }
	maxRetries 2

	publishDir params.outputDir, mode: "move"

	input:
		file fastq1 from Channel.fromPath(params.fastq1.split(',').toList()).collect()
		file fastq2 from Channel.fromPath(params.fastq2.split(',').toList()).collect()
		file bam from Channel.fromPath(params.bam)
		file sophiaSVs from Channel.fromPath(params.sophiaSVs)
		file starIndex from Channel.fromPath(params.starIndex)
		file assembly from Channel.fromPath(params.assembly)
		file annotation from Channel.fromPath(params.annotation)
		file blacklist from Channel.fromPath(params.blacklist)
		file knownFusions from Channel.fromPath(params.knownFusions)
		file proteinDomains from Channel.fromPath(params.proteinDomains)
		file cytobands from Channel.fromPath(params.cytobands)

	output:
		file 'fusions.tsv' into fusions
		file 'fusions_discarded.tsv.gz' into fusions_discarded
		file 'fusions.pdf' into fusions_pdf
		file 'fusions_alignments.bam' into fusions_alignments
		file 'fusions_alignments.bam.bai' into fusions_alignments_index
		file 'virus_expression.tsv' into virus_expression
		file 'virus_alignments.bam' into virus_alignments
		file 'virus_alignments.bam.bai' into virus_alignments_index

	shell:
		"""
		FASTQ1="${params.fastq1.contains('dummy_files') ? '' : fastq1.join(',')}" \
		FASTQ2="${params.fastq2.contains('dummy_files') ? '' : fastq2.join(',')}" \
		BAM="${params.bam.contains('dummy_files') ? '' : bam}" \
		SOPHIA_SVS="${params.sophiaSVs.contains('dummy_files') ? '' : sophiaSVs}" \
		STAR_INDEX="$starIndex" \
		ASSEMBLY="$assembly" \
		ANNOTATION="$annotation" \
		BLACKLIST="$blacklist" \
		KNOWN_FUSIONS="$knownFusions" \
		PROTEIN_DOMAINS="$proteinDomains" \
		CYTOBANDS="$cytobands" \
		THREADS="$params.threads" \
		MEMORY="$params.memory" \
		OUTPUT_DIR="$params.outputDir" \
		VIRAL_CONTIGS="$params.viralContigs" \
		IGNORE_VIRUSES="$params.ignoreViruses" \
		main.sh 2>&1
		"""

}

