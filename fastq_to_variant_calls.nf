#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
==================================================================================================
    Nextflow DSL2 implementation of data processing based on GATK best practices: raw FASTQ to VCFs
==================================================================================================
 Author: Netsanet Gebremedhin and Michael DeRan
 Diamond Age Data Science <opensource@diamondage.com>
------------------------------------------------------------------------------------------------------------
*/


include { formatConversion } from './paired-fastq-to-unmapped-bam.nf'
include { preprocessing_mapping; quality_recalibraiton; variant_discovery  } from './processing-for-variant-discovery-gatk4.nf'

fastq_files_list = file(params.input_fofn)


// Convert input manifest to a channel.

fastq_params = channel.fromPath(fastq_files_list)
                   .splitText(keepHeader:false)
                    .map { line ->
                           cols = line.tokenize('\t')
                                [
                                    file(cols[2]),
                                    file(cols[3]),
                                    cols[6],
                                    cols[1],
                                    cols[4],
                                    cols[7],
                                    cols[5],
                                    cols[0],
                                    cols[8]

                    ]
}


workflow {
    formatConversion({fastq_params})
    preprocessing_mapping(formatConversion.out)
    quality_recalibraiton(preprocessing_mapping.out)
    variant_discovery(quality_recalibraiton.out)
    variant_discovery.out
	        .map {
                sampleId, vcfFile, VcfIndex ->
                "${sampleId}\t${vcfFile}"
            }
            .collectFile(
                name: "${params.projectId}_merged_vcfs.tsv", newLine: true, storeDir: "${params.outdir}/${params.projectId}/vcfs"
            )

}
