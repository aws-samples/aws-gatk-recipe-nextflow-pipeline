#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
==========================================================================================
    Nextflow DSL2 implementation of joint discovery and VQSR based on GATK best practices
==========================================================================================
 Author: Netsanet Gebremedhin and Michael DeRan
 Diamond Age Data Science <opensource@diamondage.com>
------------------------------------------------------------------------------------------------------------
*/


callset_name = params.callset_name
input_vcfs_path = params.input_vcf_path
scatter_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"

ref_fasta = file(params.fasta)
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta",".dict") )
dbSNP_vcf = file(params.dbSNP_vcf)
dbSNP_vcf_index = file("${params.dbSNP_vcf}.idx")
known_indels_mills = file(params.known_indels_mills)
known_indels_mills_index = file("${params.known_indels_mills}.tbi")
axiompoly_vcf = file(params.axiompoly_vcf)
axiompoly_vcf_index = file("${params.axiompoly_vcf}.tbi")
hapmap_vcf = file(params.hapmap_vcf)
hapmap_vcf_index = file("${params.hapmap_vcf}.tbi")
omni_vcf = file(params.omni_vcf)
omni_vcf_index = file("${params.omni_vcf}.tbi")
onekg_vcf = file(params.onekg_vcf)
onekg_vcf_index = file("${params.onekg_vcf}.tbi")

interval_list = file(params.unpadded_intervals_file)
eval_interval_list  = file(params.eval_interval_list)
scatter_count = 23
workspace_dir_name = "genomicsdb"
batch_size = 50
excess_het_threshold = 54.69
max_gaussians_indels_recal = 4
max_gaussians_snps_recal = 6
snp_filter_level = 99.7
indel_filter_level = 99.0
snp_recalibration_tranche_values = params.snp_recalibration_tranche_values
snp_recalibration_annotation_values = params.snp_recalibration_annotation_values
indel_recalibration_tranche_values = params.indel_recalibration_tranche_values
indel_recalibration_annotation_values = params.indel_recalibration_annotation_values

gatk_joint_docker= params.gatk_joint_docker
gotc_docker= params.gotc_docker
gatk4110_docker= params.gatk4110_docker

sample_vcf_ch = Channel.fromPath("${input_vcfs_path}/*.vcf.gz", followLinks: false)
                                .map { vcf ->
                                    sampleId = vcf.baseName.tokenize('.')[0]
                                    vcfFile = "s3:/${vcf}"
                                    vcfIndex = "s3:/${vcf}.tbi"
                                    [ sampleId, vcfFile, vcfIndex ]
                                }


process splitIntervalList {
    tag { "Split interval list" }
	memory '4 GB'
    cpus 1

    container gatk4110_docker


    input:

    path(interval_list)
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)
    val(scatter_count)

    output:
    path("scatterDir/*-scattered.interval_list")


    script:

    """
    gatk --java-options -Xms3g SplitIntervals \
    -L ${interval_list} -O  scatterDir -scatter ${scatter_count} -R ${ref_fasta} \
    -mode ${scatter_mode} --interval-merging-rule OVERLAPPING_ONLY
    """

}

process normalizeVcfs {
    tag { "Normalize_VCF_${sampleId}" }

    errorStrategy 'retry'
    maxRetries 3


	memory '16 GB'
    cpus 4

    container gatk_joint_docker

    beforeScript 'conda install -c bioconda bcftools'

    //publishDir "${params.outdir}/${callset_name}/normalized_vcfs/", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId), path(vcf), path(vcfIndex)

    output:
    tuple path("${sampleId}.normalized.vcf.gz"), path("${sampleId}.normalized.vcf.gz.tbi")

    script:
    """
    /opt/miniconda/bin/bcftools norm -m +any -O 'z' -o ${sampleId}.normalized.vcf.gz ${vcf}
    tabix -fp vcf ${sampleId}.normalized.vcf.gz
    """
}


process GenomicsDBImport {
    tag { "GenomicsDBImport_${scattered_interval_id}" }

    errorStrategy 'retry'
    maxRetries 3


	memory '16 GB'
    cpus 4

    container gatk_joint_docker

    input:

	each path(interval_file)
	path(input_vcfs_and_indices)
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)
	path(sample_name_map)

	output:
    tuple val(scattered_interval_id), path(interval_file), path("${workspace_dir_name}.${scattered_interval_id}.tar")

    script:

    scattered_interval_id = interval_file.getBaseName().split('-')[0]

	"""
    set -euo pipefail

    rm -rf "${workspace_dir_name}.${scattered_interval_id}"

	gatk --java-options -Xms8g \
    GenomicsDBImport \
    --genomicsdb-workspace-path "${workspace_dir_name}.${scattered_interval_id}" \
    --batch-size ${batch_size} \
    -L ${interval_file} \
    --sample-name-map ${sample_name_map} \
    --reader-threads 5 \
    --merge-input-intervals \
    --consolidate

    tar -cf "${workspace_dir_name}.${scattered_interval_id}.tar" "${workspace_dir_name}.${scattered_interval_id}"

	"""
}


process genotypeGVCFs {
    tag { "Genotype_${scattered_interval_id}" }

    errorStrategy 'retry'
    maxRetries 3


	memory '16 GB'
    cpus 2

    container gatk_joint_docker

    input:
	tuple val(scattered_interval_id), path(interval_file), path(genomicsDB)
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)
    path(dbSNP_vcf)
    path(dbSNP_vcf_index)

	output:
    tuple val(scattered_interval_id), path("${callset_name}.${scattered_interval_id}.vcf.gz"), path("${callset_name}.${scattered_interval_id}.vcf.gz.tbi")

    script:
    dbName = genomicsDB.getBaseName()

	"""
    set -euo pipefail

    tar -xf $genomicsDB
    WORKSPACE=${dbName}

    gatk --java-options -Xms8g  \
    GenotypeGVCFs  \
    -R ${ref_fasta}  \
    -O "${callset_name}.${scattered_interval_id}.vcf.gz"  \
    -D ${dbSNP_vcf}  \
    -G StandardAnnotation -G AS_StandardAnnotation  \
    --only-output-calls-starting-in-intervals  \
    --use-new-qual-calculator  \
    -V gendb://\$WORKSPACE  \
    -L ${interval_file}  \
    --merge-input-intervals
	"""

}


process hardFilterAndMakeSitesOnlyVcf {
    tag { "Genotype_${scattered_interval_id}" }

    errorStrategy 'retry'
    maxRetries 3

	memory '4 GB'
    cpus 1

    container gatk_joint_docker

    input:
	tuple val(scattered_interval_id), path(genotyped_vcf), path(genotyped_vcf_idx)


	output:
    path("${callset_name}.${scattered_interval_id}.variant_filtered.vcf.gz")
    path("${callset_name}.${scattered_interval_id}.variant_filtered.vcf.gz.tbi")
    path("${callset_name}.${scattered_interval_id}.sites_only.variant_filtered.vcf.gz")
    path("${callset_name}.${scattered_interval_id}.sites_only.variant_filtered.vcf.gz.tbi")


    script:
	"""
    set -euo pipefail

    gatk --java-options -Xms3g \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O "${callset_name}.${scattered_interval_id}.variant_filtered.vcf.gz" \
      -V ${genotyped_vcf}

    gatk --java-options -Xms3g \
      MakeSitesOnlyVcf \
      -I "${callset_name}.${scattered_interval_id}.variant_filtered.vcf.gz" \
      -O "${callset_name}.${scattered_interval_id}.sites_only.variant_filtered.vcf.gz"
	"""
}


process gatherVcfs {
    tag { "GatherVCFs" }

    errorStrategy 'retry'
    maxRetries 3

	memory '8 GB'
    cpus 2

    container gatk4110_docker

    input:
    path(sites_only_filtered_vcfs)
    path(sites_only_filtered_vcf_indices)



    output:
    path("${callset_name}.sites_only.vcf.gz")
    path("${callset_name}.sites_only.vcf.gz.tbi")

    script:
    input_vcfs_params = sites_only_filtered_vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --input ")
    """
    set -euo pipefail

    gatk --java-options -Xms6g \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ${input_vcfs_params} \
      --output "${callset_name}.sites_only.vcf.gz"

    tabix "${callset_name}.sites_only.vcf.gz"
    """

}

process indelsVariantRecalibrator {
    tag { "IndelsVariantRecalibrator" }

    errorStrategy 'retry'
    maxRetries 3

	memory '16 GB'
    cpus 2

    container gatk4110_docker

    input:

    path(sites_only_variant_filtered_vcf)
    path(sites_only_variant_filtered_vcf_index)

    path(dbSNP_vcf)
    path(dbSNP_vcf_index)
    path(known_indels_mills)
    path(known_indels_mills_index)
    path(axiompoly_vcf)
    path(axiompoly_vcf_index)

    output:
    tuple  path("${callset_name}.indels.recal"),
           path("${callset_name}.indels.recal.idx"),
           path("${callset_name}.indels.tranches")



    script:

    recalibration_annotation_inputs_params = indel_recalibration_annotation_values.join(" -an ")
    recalibration_tranche_input_params = indel_recalibration_tranche_values.join(" -tranche ")

    """
    set -euo pipefail

    gatk --java-options -Xms14g \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O "${callset_name}.indels.recal" \
      --tranches-file "${callset_name}.indels.tranches" \
      --trust-all-polymorphic \
      -tranche ${recalibration_tranche_input_params} \
      -an ${recalibration_annotation_inputs_params} \
      --use-allele-specific-annotations \
      -mode INDEL \
      --max-gaussians ${max_gaussians_indels_recal} \
      -resource:mills,known=false,training=true,truth=true,prior=12 ${known_indels_mills} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiompoly_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbSNP_vcf}
    """
}


process snpsVariantRecalibrator {
    tag { "SNPsVariantRecalibrator" }

    errorStrategy 'retry'
    maxRetries 3

	memory '16 GB'
    cpus 2

    container gatk4110_docker

    input:

    path(sites_only_variant_filtered_vcf)
    path(sites_only_variant_filtered_vcf_index)

    path(dbSNP_vcf)
    path(dbSNP_vcf_index)
    path(omni_vcf)
    path(omni_vcf_index)
    path(hapmap_vcf)
    path(hapmap_vcf_index)
    path(onekg_vcf)
    path(onekg_vcf_index)

    output:
    tuple   path("${callset_name}.snps.recal"),
            path("${callset_name}.snps.recal.idx"),
            path("${callset_name}.snps.tranches")



    script:

    recalibration_annotation_inputs_params = snp_recalibration_annotation_values.join(" -an ")
    recalibration_tranche_input_params = snp_recalibration_tranche_values.join(" -tranche ")

    """
    set -euo pipefail

    gatk --java-options -Xms14g \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O "${callset_name}.snps.recal" \
      --tranches-file "${callset_name}.snps.tranches" \
      --trust-all-polymorphic \
      -tranche ${recalibration_tranche_input_params} \
      -an ${recalibration_annotation_inputs_params} \
      --use-allele-specific-annotations \
      -mode SNP \
      --max-gaussians ${max_gaussians_snps_recal} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ${omni_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ${onekg_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbSNP_vcf}
    """
}


// Scatter and apply recalibration to both snps and indels vcfs

process applyRecalibration {
    tag { "ApplyRecalibration" }

    errorStrategy 'retry'
    maxRetries 3

	memory '8 GB'
    cpus 2

    container gatk4110_docker

    input:
    path(input_vcf)
    path(input_vcf_indices)

    tuple   path(input_indels_recal),
            path(input_indels_recal_index),
            path(input_indels_tranches)

    tuple   path(input_snps_recal),
            path(input_snps_recal_index),
            path(input_snps_tranches)

    output:
    path("${callset_name}.${scattered_interval_id}.filtered.vcf.gz")
    path("${callset_name}.${scattered_interval_id}.filtered.vcf.gz.tbi")


    script:
    scattered_interval_id = input_vcf.getBaseName().tokenize('.')[1]
    """
    set -euo pipefail

    gatk --java-options -Xms5g \
      ApplyVQSR \
      -O 'tmp.indel.recalibrated.vcf' \
      -V ${input_vcf} \
      --recal-file ${input_indels_recal} \
      --use-allele-specific-annotations  \
      --tranches-file ${input_indels_tranches} \
      --truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    gatk --java-options -Xms5g \
      ApplyVQSR \
      -O "${callset_name}.${scattered_interval_id}.filtered.vcf.gz" \
      -V 'tmp.indel.recalibrated.vcf' \
      --recal-file ${input_snps_recal} \
      --use-allele-specific-annotations \
      --tranches-file ${input_snps_tranches} \
      --truth-sensitivity-filter-level ${snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
    """
}



process gatherFinalVcf {
    tag { "GatherFinalRecalibratedVCFs" }

    errorStrategy 'retry'
    maxRetries 3

	memory '8 GB'
    cpus 2

    container gatk4110_docker

    publishDir "${params.outdir}/${callset_name}/", mode: 'copy', overwrite: true

    input:
    path(input_recal_vcfs)
    path(input_recal_vcf_indices)



    output:
    tuple path("${callset_name}.vcf.gz"), path("${callset_name}.vcf.gz.tbi")

    script:
    input_recal_vcfs_params = input_recal_vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --input ")
    """
    set -euo pipefail

    gatk --java-options -Xms6g \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ${input_recal_vcfs_params} \
      --output "${callset_name}.vcf.gz"

    tabix "${callset_name}.vcf.gz"
    """

}


process collectVariantCallingMetrics {
    tag { "CollectVariantCallingMetrics" }

    errorStrategy 'retry'
    maxRetries 3

	memory '8 GB'
    cpus 2

    container gatk4110_docker

    publishDir "${params.outdir}/${callset_name}/", mode: 'copy', overwrite: true

    input:
    tuple path(final_recal_full_vcf), path(final_recal_full_vcf_index)
    path(ref_dict)
    path(dbSNP_vcf)
    path(dbSNP_vcf_index)
    path(eval_interval_list)



    output:
    tuple   path("${callset_name}.variant_calling_detail_metrics"),
            path("${callset_name}.variant_calling_summary_metrics")

    script:
    """
    set -euo pipefail

    gatk --java-options -Xms6g \
      CollectVariantCallingMetrics \
      --INPUT ${final_recal_full_vcf} \
      --DBSNP ${dbSNP_vcf} \
      --SEQUENCE_DICTIONARY ${ref_dict} \
      --OUTPUT ${callset_name} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ${eval_interval_list}
    """
}



workflow Genotyping {
    take:   data

    main:
        splitIntervalList(interval_list,
                        ref_dict,
                        ref_fasta,
                        ref_fasta_fai,
                        scatter_count
                    )

        normalizeVcfs(data)

        normalizeVcfs.out
            .map {
                normalizedVcf, normalizedVcfIndex ->
                sampleId = normalizedVcf.baseName.tokenize('.')[0]
                normalizedVcfName = normalizedVcf.baseName
                [sampleId, normalizedVcfName]
            }
            .map {
                sampleId, vcfFile ->
                "${sampleId}\t${vcfFile}.gz"
            }
            .collectFile(
                name: "${callset_name}_normalized_vcf_inputs.tsv", newLine: true, storeDir: "${params.outdir}/${callset_name}/"
            )

        sample_normalized_vcf_map = "${params.outdir}/${callset_name}/${callset_name}_normalized_vcf_inputs.tsv"

        GenomicsDBImport(splitIntervalList.out,
                        normalizeVcfs.out.collect(),
                        ref_dict,
                        ref_fasta,
                        ref_fasta_fai,
	                    sample_normalized_vcf_map
                    )

        genotypeGVCFs(GenomicsDBImport.out,
                    ref_dict,
                    ref_fasta,
                    ref_fasta_fai,
                    dbSNP_vcf,
                    dbSNP_vcf_index
                )


    emit:
        genotypeGVCFs.out

}

workflow Recalibration {
    take:
        data

    main:

        hardFilterAndMakeSitesOnlyVcf(data)

        gatherVcfs(
                    hardFilterAndMakeSitesOnlyVcf.out[2].collect(),
                    hardFilterAndMakeSitesOnlyVcf.out[3].collect()
                 )

        indelsVariantRecalibrator(gatherVcfs.out[0],
                                 gatherVcfs.out[1],
                                dbSNP_vcf,
                                dbSNP_vcf_index,
                                known_indels_mills,
                                known_indels_mills_index,
                                axiompoly_vcf,
                                axiompoly_vcf_index
                            )

        snpsVariantRecalibrator(gatherVcfs.out[0],
                                gatherVcfs.out[1],
                                dbSNP_vcf,
                                dbSNP_vcf_index,
                                omni_vcf,
                                omni_vcf_index,
                                hapmap_vcf,
                                hapmap_vcf_index,
                                onekg_vcf,
                                onekg_vcf_index
                            )

        applyRecalibration(
                            hardFilterAndMakeSitesOnlyVcf.out[0],
                            hardFilterAndMakeSitesOnlyVcf.out[1],
                            indelsVariantRecalibrator.out,
                            snpsVariantRecalibrator.out
                        )

        gatherFinalVcf(
                    applyRecalibration.out[0].collect(),
                    applyRecalibration.out[1].collect()
                 )

    emit:
        gatherFinalVcf.out

}


workflow CollectMetrics {
    take:   data

    main:
        collectVariantCallingMetrics(
                    data,
                    ref_dict,
                    dbSNP_vcf,
                    dbSNP_vcf_index,
                    eval_interval_list
            )
}



workflow {

        Genotyping(sample_vcf_ch)
        Recalibration(Genotyping.out)
        CollectMetrics(Recalibration.out)
}
