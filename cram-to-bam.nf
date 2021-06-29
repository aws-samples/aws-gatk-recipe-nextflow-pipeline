#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
============================================================================================================
Nextflow DSL2 pipeline for conversion of cram files to unmapped bam files
============================================================================================================
 Author: Netsanet Gebremedhin and Michael DeRan
 Diamond Age Data Science <opensource@diamondage.com>
------------------------------------------------------------------------------------------------------------
*/


params.cram_list = file("cramToBam.tsv")
params.java_opt_validateBAM = "-Xmx3G"

cram_list = file(params.input_crams)

cram_params = channel.fromPath(cram_list)
                   .splitText(keepHeader:false)
                    .map { line ->
                           cols = line.tokenize('\t')
                                [
                                    file(cols[1].trim()),
                                    cols[0]
                                ]
                    }

ref_fasta = file(params.fasta)
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta",".dict") )



process cramToBam {
    tag "${sample_name}"

    memory '32 GB'
    cpus 32

    errorStrategy 'retry'
    maxRetries 6

    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"


    publishDir "${params.outdir}/BAMs/", mode: 'copy', overwrite: true

    input:
    tuple   path(input_cram),
            val(sample_name) 
    path(ref_fasta)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bai")

    script:
    """
    set -eo pipefail

    samtools view -h -T ${ref_fasta} ${input_cram} | samtools view -b -@ 32 -o ${sample_name}.bam -
    samtools index -b ${sample_name}.bam
    mv ${sample_name}.bam.bai ${sample_name}.bai

    """

}



process revertSam {

    tag "${sample_name}"

    memory '4 GB'
    cpus 4

    errorStrategy 'retry'
    maxRetries 6

    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"

    publishDir "${params.outdir}/uBAMs/", mode: 'copy', overwrite: true

    input:
    tuple   val(sample_name),
            path(bam),
            path(bam_index)

    output:
    tuple   val(sample_name), 
            path("${sample_name}.u.bam")

    script:
    output_name = "${sample_name}.u.bam"

    """
    java ${params.java_opt_validateBAM} -jar /usr/gitc/picard.jar \
      RevertSam \
      INPUT=${bam} \
      OUTPUT=${output_name}
    """

}




process validateBam {
    tag "${sample_name}"

    memory '4 GB'
    cpus 4

    errorStrategy 'retry'
    maxRetries 6

    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"

    publishDir "${params.outdir}/reports/", mode: 'copy', overwrite: true

    input:
    tuple   val(sample_name),
            path(unmapped_bam)

    output:
    tuple   val(sample_name), 
            path("${unmapped_bam}.validation_report")

    script:
    output_name = "${unmapped_bam}.validation_report"

    """
    java ${params.java_opt_validateBAM} -jar /usr/gitc/picard.jar \
      ValidateSamFile \
      INPUT=${unmapped_bam} \
      OUTPUT=${output_name} \
      MODE=SUMMARY \
      IS_BISULFITE_SEQUENCED=false 
    """


}




workflow formatConversion{
    take: data

    main:
        cramToBam(data, ref_fasta)
        
        revertSam(
            cramToBam.out
        )
        
        validateBam(
            revertSam.out
        )

    emit:
        cramToBam.out
        revertSam.out
        validateBam.out

}




workflow {
    formatConversion({cram_params})

    Channel.fromPath("${params.outdir}/BAMs/*.bam")
	        .map { bam ->
		        sampleId = bam.baseName.tokenize('.')[0]
		        bamFile = "s3:/${bam}"
		        [ sampleId, bamFile]
	        }

	        .map {
                sampleId, bamFile ->
                "${sampleId}\t${bamFile}"
            }

            .collectFile(
                name: "bams.tsv", newLine: true, storeDir: "${params.outdir}"
            )

    Channel.fromPath("${params.outdir}/uBAMs/*.u.bam")
	        .map { bam ->
		        sampleId = bam.baseName.tokenize('.')[0]
		        bamFile = "s3:/${bam}"
		        [ sampleId, bamFile]
	        }

	        .map {
                sampleId, bamFile ->
                "${sampleId}\t${bamFile}"
            }

            .collectFile(
                name: "unmapped_bams.tsv", newLine: true, storeDir: "${params.outdir}"
            )


    Channel.fromPath("${params.outdir}/reports/*.validation_report")
	        .map { report ->
		        sampleId = report.baseName.tokenize('.')[0]
		        validationReport = "s3:/${report}"
		        [ sampleId, validationReport]
	        }

	        .map {
                sampleId, validationReport ->
                "${sampleId}\t${validationReport}"
            }

            .collectFile(
                name: "validation_reports.tsv", newLine: true, storeDir: "${params.outdir}"
            )

}