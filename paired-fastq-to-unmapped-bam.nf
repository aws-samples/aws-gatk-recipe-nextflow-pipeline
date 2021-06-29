#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
================================================================================================
    Nextflow DSL2 implemenation of fastq to unmapped bam conversion based on GATK best practices
================================================================================================
 Author: Netsanet Gebremedhin and Michael DeRan
 Diamond Age Data Science <opensource@diamondage.com>
------------------------------------------------------------------------------------------------------------
*/


fastq_files_list = file(params.input_fofn)
gatk_docker= params.gatk_docker

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



process pairedFastqToUnmappedBam {
    tag "${sample_name}"


    memory '8 GB'
    cpus 2

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:

    tuple   file(fastq_1),
            file(fastq_2),
            val(run_date),
            val(sample_name),
            val(library_name),
            val(platform_name),
            val(platform_unit),
            val(readgroup_name),
            val(sequencing_center)

    output:
    tuple val(sample_name), path("${readgroup_name}.unmapped.bam")

    script:

    gatk_path = "/gatk/gatk"

    """
    ${gatk_path}  \
    FastqToSam \
    --FASTQ ${fastq_1} \
    --FASTQ2 ${fastq_2} \
    --OUTPUT ${readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ${readgroup_name} \
    --SAMPLE_NAME ${sample_name} \
    --LIBRARY_NAME ${library_name} \
    --PLATFORM_UNIT ${platform_unit} \
    --RUN_DATE ${run_date} \
    --PLATFORM ${platform_name} \
    --SEQUENCING_CENTER ${sequencing_center}
    """
}


workflow formatConversion{
    take: data

    main:
        pairedFastqToUnmappedBam(data)
    emit:
        pairedFastqToUnmappedBam.out
}


workflow {
    formatConversion({fastq_params})
    formatConversion.out
        .map { sample_name, unmapped_bam ->
            "${sample_name}\ts3:/${unmapped_bam}"
        }
        .collectFile(
            name: 'unmapped_bams.tsv', newLine: true, storeDir: "${params.outdir}"
        )
}
