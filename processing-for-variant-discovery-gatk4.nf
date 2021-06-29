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

if( !params["unmapped_bams_list"] ){
    Channel
        .from(1)
        .map{ chr -> 
            "sample_${chr}\tfile_${chr}"
        }

        .collectFile(
            name: "unmapped_bams.tsv", newLine: true, storeDir: "${params.outdir}/${params.projectId}/temp"
        )

    params["unmapped_bams_list"] = "${params.outdir}/${params.projectId}/temp/unmapped_bams.tsv"
}

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

    input:
    val(sampleId)
    path(input_mapped_merged_marked_bam)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    output:
        tuple val(sampleId),
        path("${sampleId}.mapped.merged.duplicate_marked.sorted.bam"),
        path("${sampleId}.mapped.merged.duplicate_marked.sorted.bai"),
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

process createSequenceGrouping {
    tag { "Create sequence grouping" }

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gotc_docker

    input:
    path(ref_dict)

    output:
    path("*txt")
    path("*unmapped.txt")

    script:
    """
    sequenceGrouping.py ${ref_dict} "sequence_grouping.txt" "sequence_grouping_with_unmapped.txt"
    """
}



process baseRecalibrator {
    tag "${sampleId}_${subgroup_name}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:
    tuple val(sampleId),
           path(input_mapped_merged_marked_sorted_bam),
           path(input_mapped_merged_marked_sorted_bai),
           path(input_mapped_merged_marked_sorted_md5),
           val(subgroup_name),
           val(subgroup)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    path(dbSNP_vcf)
    path(dbSNP_vcf_index)
    path(known_indels_mills)
    path(known_indels_mills_index)
    path(known_indels_dbSNP)
    path(known_indels_dbSNP_index)

    output:
        tuple val(sampleId),
        path("${sampleId}_recalibration_report_${subgroup_name}.recal_data.csv")

    script:
    gatk_path = "/gatk/gatk"
    subgroup_trimmed = subgroup.trim().split("\t").join(" -L ")

    """
    ${gatk_path} --java-options ${params.java_opt_baserecal} \
    BaseRecalibrator \
    -R ${ref_fasta} \
    -I ${input_mapped_merged_marked_sorted_bam} \
    --use-original-qualities \
    -O "${sampleId}_recalibration_report_${subgroup_name}.recal_data.csv" \
    --known-sites ${dbSNP_vcf} \
    --known-sites ${known_indels_mills} \
    --known-sites ${known_indels_dbSNP} \
    -L ${subgroup_trimmed}
    """
}


process gatherBqsrReports {
    tag "${sampleId}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:
    tuple val(sampleId), path(input_bqsr_reports)

    output:
        tuple val(sampleId),
        path("${sampleId}.recal_data.csv")

    script:
    gatk_path = "/gatk/gatk"
    input_bqsr_params = input_bqsr_reports.join(" -I ")

    """
    ${gatk_path} --java-options ${params.java_opt_bqsrreport} \
    GatherBQSRReports \
    -I ${input_bqsr_params} \
    -O "${sampleId}.recal_data.csv"
    """
}


process applyBQSR {
    tag "${sampleId}_${subgroup_unmapped_name}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:
    tuple val(sampleId),
           path(input_mapped_merged_marked_sorted_bam),
           path(input_mapped_merged_marked_sorted_bai),
           path(input_mapped_merged_marked_sorted_md5),
           path(input_merged_bqsr_report),
           val(scatter_id),
           val(subgroup_unmapped_name),
           val(subgroup_unmapped)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)


    output:
    tuple val(sampleId),
        path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${subgroup_unmapped_name}.recal.bam")

    script:

    gatk_path = "/gatk/gatk"
    subgroup_unmapped_trimmed = subgroup_unmapped.trim().split("\t").join(" -L ")

    """
    ${gatk_path} --java-options "${params.java_opt_applybqsr}" \
    ApplyBQSR \
    -R ${ref_fasta}  \
    -I ${input_mapped_merged_marked_sorted_bam}  \
    -O "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${subgroup_unmapped_name}.recal.bam"  \
    -L ${subgroup_unmapped_trimmed} \
    -bqsr ${input_merged_bqsr_report} \
    --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities
    """
}


process gatherBamFiles {
    tag "${sampleId}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:
    tuple val(sampleId), path(input_recalibrated_bams)

    output:
    tuple   val(sampleId),
            path("${sampleId}.recal.merged.bam"),
            path("${sampleId}.recal.merged.bai"),
            path("${sampleId}.recal.merged.bam.md5")

    script:

    gatk_path = "/gatk/gatk"
    inputs_bams_to_merge = input_recalibrated_bams
                            .sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --INPUT ")

    """
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level}  \
    ${params.java_opt_gatherbams} " \
    GatherBamFiles \
    --INPUT ${inputs_bams_to_merge} \
    --OUTPUT "${sampleId}.recal.merged.bam" \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true
    """
}


process haplotypeCaller {
    tag "${sampleId}_${interval_chunk_name}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    input:

    tuple val(sampleId),
           path(input_recal_merged_bam),
           path(input_recal_merged_bai),
           path(input_recal_merged_md5),
           val(scatter_id),
           val(interval_chunk_name),
           path(interval_list_file)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)


    output:

    tuple   val(sampleId),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf"),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf.idx")


    script:

    gatk_path = "/gatk/gatk"

    """
    set -e
    ${gatk_path} --java-options "${params.java_opt_haplotype}" \
    HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${input_recal_merged_bam} \
    --output "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf" \
    -contamination 0 \
    -ERC GVCF \
    -L ${interval_list_file}
    """
}


process mergeVCFs {
    tag "${sampleId}"

    memory '32 GB'
    cpus 16

    errorStrategy 'retry'
    maxRetries 3

    container gatk_docker

    publishDir "${params.outdir}/${params.projectId}/vcfs/", mode: 'copy', overwrite: true

    input:

    tuple val(sampleId),
           path(input_vcfs_to_merge),
           path(inputs_vcf_indices)


    output:
    tuple   val(sampleId),
            path("${sampleId}.merged.vcf.gz"),
            path("${sampleId}.merged.vcf.gz.tbi")

    script:

    gatk_path = "/gatk/gatk"
    input_vcfs_params = input_vcfs_to_merge
                            .sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --INPUT ")

    """
    set -e

    ${gatk_path} --java-options "${params.java_opt_mergevcfs}"  \
    MergeVcfs \
    --INPUT ${input_vcfs_params} \
    --OUTPUT ${sampleId}.merged.vcf

    bgzip -c ${sampleId}.merged.vcf>${sampleId}.merged.vcf.gz
    tabix -fp vcf ${sampleId}.merged.vcf.gz
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


// Prepare channel inputs for Base Recalibration

subgrouping = channel.fromPath(sequence_grouping)
                    .splitText()
                    .map { line ->[ line.tokenize(':')[0].trim(), line.trim() ] }

// Prepare channel inputs for  applying BQSR

recal_scatter_counter = 0
subgrouping_unmapped = channel.fromPath(sequence_grouping_unmapped)
                            .splitText()
                            .map { line ->   recal_scatter_counter = recal_scatter_counter + 1
                                [ recal_scatter_counter, line.tokenize(':')[0].trim(), line.trim() ]
                            }


// Prepare channel for scattered calling intervals: input to haplotypecaller

calling_scatter_counter = 0

calling_intervals = channel.fromPath(scattered_calling_interval)
                        .splitText()
                        .map { line -> calling_scatter_counter = calling_scatter_counter + 1
                            [calling_scatter_counter, line.split("/")[6], line.trim() ]
                        }



/*
*
################################################ EXECUTION #############################################################
*
*/

workflow preprocessing_mapping {
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

    emit:
        sortAndFixTags.out
}


workflow quality_recalibraiton {
    take: data

    main:
        baseRecalibrator(
                data.combine(subgrouping), \
                ref_dict,   \
                ref_fasta,  \
                ref_fasta_fai,  \
                dbSNP_vcf,  \
                dbSNP_vcf_index,    \
                known_indels_mills, \
                known_indels_mills_index,   \
                known_indels_dbSNP,     \
                known_indels_dbSNP_index   \
        )

        gatherBqsrReports(
                baseRecalibrator.out.groupTuple()
        )


        applyBQSR(
                data.join(gatherBqsrReports.out).combine(subgrouping_unmapped), \
                ref_dict,   \
                ref_fasta,  \
                ref_fasta_fai  \
        )

        gatherBamFiles(
                applyBQSR.out.groupTuple()
        )

    emit:
        gatherBamFiles.out

}


workflow variant_discovery {
    take: data

    main:

        haplotypeCaller(
                data.combine(calling_intervals), \
                ref_dict,   \
                ref_fasta,  \
                ref_fasta_fai \
        )

        mergeVCFs(
            haplotypeCaller.out.groupTuple()
        )
    emit:
        mergeVCFs.out


}


workflow {

        preprocessing_mapping(unmapped_bams_channel)

        quality_recalibraiton(preprocessing_mapping.out)

        variant_discovery(quality_recalibraiton.out)

        Channel.fromPath("${params.outdir}/${params.projectId}/vcfs/*.vcf.gz")
	        .map { vcf ->
		        sampleId = vcf.baseName.tokenize('.')[0]
		        vcfFile = "s3:/${vcf}"
		        [ sampleId, vcfFile]
	        }
	        .map {
                sampleId, vcfFile ->
                "${sampleId}\t${vcfFile}"
            }
            .collectFile(
                name: "${params.projectId}_merged_vcfs.tsv", newLine: true, storeDir: "${params.outdir}/${params.projectId}/vcfs"
            )
}
