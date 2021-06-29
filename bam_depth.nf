#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
============================================================================================================
    Nextflow DSL2 implementation of data processing based on GATK best practices: Germline Variant Discovery
============================================================================================================
 Author: Netsanet Gebremedhin and Michael DeRan
 Diamond Age Data Science <opensource@diamondage.com>
------------------------------------------------------------------------------------------------------------
*/


/*
*
################################################ WORKFLOW INPUTS #######################################################
*
*/

unmapped_bams = file(params.unmapped_bams_list)

ref_fasta = file(params.fasta)
ref_alt =  file("${params.fasta}.64.alt")
ref_amb = file("${params.fasta}.64.amb")
ref_ann = file("${params.fasta}.64.ann")
ref_bwt = file("${params.fasta}.64.bwt")
ref_pac = file("${params.fasta}.64.pac")
ref_sa = file("${params.fasta}.64.sa")
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta",".dict") )
dbSNP_vcf = file(params.dbSNP_vcf)
dbSNP_vcf_index = file("${params.dbSNP_vcf}.idx")
known_indels_mills = file(params.known_indels_mills)
known_indels_mills_index = file("${params.known_indels_mills}.tbi")
known_indels_dbSNP = file(params.known_indels_dbSNP)
known_indels_dbSNP_index = file("${params.known_indels_dbSNP}.tbi")
sequence_grouping = file(params.sequence_grouping)
sequence_grouping_unmapped = file(params.sequence_grouping_unmapped)
scattered_calling_interval = file(params.scattered_calling_interval)

gatk_docker= params.gatk_docker
gotc_docker= params.gotc_docker
gatk4110_docker= params.gatk4110_docker

/*
*
################################################ PROCESS DECLARATION ###################################################
*
*/

process getBwaVersion {
    tag "BWA version"

    container gotc_docker

    output:
    stdout

    script:

    bwa_path = "/usr/gitc"
    """

    ${bwa_path}/bwa 2>&1  \
    | grep -e '^Version'  \
    | sed 's/Version: //' \
    | tr -d '\n'
    """
}

process samToFastqBwaMem {
    tag "${sampleId}"

    memory '64 GB'
    cpus 32

    errorStrategy 'retry'
    maxRetries 3

    container gotc_docker

    input:
    tuple val(sampleId), path(input_unmapped_bam)

    path(ref_alt)
    path(ref_amb)
    path(ref_ann)
    path(ref_bwt)
    path(ref_pac)
    path(ref_sa)
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    output:
    val(sampleId)
    path "*.mapped.bam"
    path(input_unmapped_bam)

    script:
    
    gotc_path = "/usr/gitc"
    bwa_path = gotc_path
    """
    set -o pipefail
    set -e

    java -Dsamjdk.compression_level=${params.compression_level} ${params.java_opt_samtofastq} \
        -jar ${gotc_path}/picard.jar \
        SamToFastq \
        INPUT=${input_unmapped_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true \
    | \
        ${bwa_path}/bwa mem \
         -K 100000000 -p -v 3 -t 32 -Y ${ref_fasta} /dev/stdin -  2> >(tee ${sampleId}.bwa.stderr.log >&2) \
    | \
        samtools view -1 - > ${sampleId}.mapped.bam
    """
}


process mergeBamAlignment {
    tag "${sampleId}"

    memory '64 GB'
    cpus 32

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:

    val(sampleId)
    path(input_mapped_bam)
    path(input_unmapped_bam)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    val(bwa_version)

    output:
    val(sampleId)
    path "*.mapped.merged.bam"

    script:
    // DRY: Move this to config or params file
    gatk_path = "/gatk/gatk"
    bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 32 -Y ${ref_fasta}"

    """
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opt_mergebams}" \
    MergeBamAlignment  \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ALIGNED_BAM ${input_mapped_bam} \
    --UNMAPPED_BAM ${input_unmapped_bam} \
    --OUTPUT ${sampleId}.mapped.merged.bam \
    --REFERENCE_SEQUENCE ${ref_fasta} \
    --PAIRED_RUN true \
    --SORT_ORDER "unsorted" \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --PROGRAM_RECORD_ID "bwamem" \
    --PROGRAM_GROUP_VERSION "${bwa_version}" \
    --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
    --PROGRAM_GROUP_NAME "bwamem" \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true
    """
}


process markDuplicates {
    tag "${sampleId}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:
    val(sampleId)
    path(input_mapped_merged_bam)

    output:
    val(sampleId)
    path("*_merged.deduped.bam")
    path("*_merged.deduped.metrics.txt")

    script:
    gatk_path = "/gatk/gatk"

    """
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opt_markdups}" \
    MarkDuplicates \
    --INPUT ${input_mapped_merged_bam} \
    --OUTPUT ${sampleId}_merged.deduped.bam \
    --METRICS_FILE ${sampleId}_merged.deduped.metrics.txt \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER "queryname" \
    --CREATE_MD5_FILE true
    """

}

process sortAndFixTags {
    tag "${sampleId}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    publishDir "${params.outdir}/sorted.bams/", mode: 'copy', overwrite: true

    input:
    val(sampleId)
    path(input_mapped_merged_marked_bam)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    output:
    val(sampleId)
    path("${sampleId}.mapped.merged.duplicate_marked.sorted.bam")
    path("${sampleId}.mapped.merged.duplicate_marked.sorted.bai")
    path("${sampleId}.mapped.merged.duplicate_marked.sorted.bam.md5")


    script:
    gatk_path = "/gatk/gatk"

    """
    set -o pipefail

    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opt_sort}" \
      SortSam \
    --INPUT ${input_mapped_merged_marked_bam} \
    --OUTPUT /dev/stdout \
    --SORT_ORDER "coordinate" \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    | \
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opt_fix}" \
    SetNmMdAndUqTags \
    --INPUT /dev/stdin \
    --OUTPUT ${sampleId}.mapped.merged.duplicate_marked.sorted.bam \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true \
    --REFERENCE_SEQUENCE ${ref_fasta}
    """
}


process mosdepth {
    tag "mosdepth_${sampleId}"

    memory '16 GB'
    cpus 8

    errorStrategy 'retry'
    maxRetries 3

    container "quay.io/biocontainers/mosdepth:0.3.1--ha7ba039_0"

    publishDir "${params.outdir}/mosdepth/", mode: 'copy', overwrite: true

    input:
    val(sampleId)
    path(input_mapped_bam)
    path(input_mapped_bam_bai)

    output:
    val(sampleId)
    path("*.mosdepth.global.dist.txt")
    path("*.mosdepth.region.dist.txt")
    path("*.mosdepth.summary.txt")
    path("*.regions.bed.gz")
    path("*.regions.bed.gz.csi")

    script:
    """
    mosdepth -n -t 4 --fast-mode --by 500 ${sampleId} ${input_mapped_bam}
    """
}

/*
*
################################################ STAGING ###############################################################
*
*/


unmapped_bams_channel = channel.fromPath(unmapped_bams)
                                .splitText()
                                .map { line ->[ line.tokenize("\t")[0],file(line.tokenize("\t")[1].trim())] }

/*
*
################################################ EXECUTION #############################################################
*
*/

workflow bam_depth {
    take: data

    main:
        samToFastqBwaMem(data,ref_alt,ref_amb,ref_ann,ref_bwt,ref_pac,ref_sa,ref_dict,ref_fasta,ref_fasta_fai)

        bwa_version = getBwaVersion()

        mergeBamAlignment(
            samToFastqBwaMem.out[0],
            samToFastqBwaMem.out[1],
            samToFastqBwaMem.out[2],
            ref_dict,
            ref_fasta,
            ref_fasta_fai,
            bwa_version
        )

        markDuplicates(mergeBamAlignment.out[0],mergeBamAlignment.out[1])

        sortAndFixTags(markDuplicates.out[0],markDuplicates.out[1],ref_dict,ref_fasta,ref_fasta_fai)

        mosdepth(
            sortAndFixTags.out[0],
            sortAndFixTags.out[1],
            sortAndFixTags.out[2],
        )

}


workflow {

        bam_depth(unmapped_bams_channel)
        
}
