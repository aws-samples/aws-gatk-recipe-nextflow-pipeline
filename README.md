# Nextflow Pipelines for GATK4 #

## Paired-fastqs to uBAM

```
nextflow run aws-nextflow-gatk/paired-fastq-to-unmapped-bam.nf \
  -work-dir s3://<your_bucket>/work \
  --outdir s3://<your_bucket>/results/fastq-to-ubam/ \
  --input_fofn s3://<your_bucket>/input_files/fastq_manifest.txt
```

### Inputs
- work-dir: nextflow work directory. Must be a location in s3
- outdir: location in s3 where outputs are saved
- input_fofn: an input manifest file of tab-delimited value of the following parameters:
```
readgroup_name
sample_name
fastq_1
fastq_2
library_name
platform_unit
run_date
platform_name
sequencing_center
```

#### Example Input:

```
NA19725_A	NA19725	s3://1000genomes/phase3/data/NA19725/sequence_read/SRR032764_1.filt.fastq.gz	s3://1000genomes/phase3/data/NA19725/sequence_read/SRR032764_2.filt.fastq.gz	Solexa-16044	BI.PE.091118_SL-XBD_0005_FC43265AAXX.1	2010-01-05T00:00:00Z	ILLUMINA	BI
NA19723_A	NA19723	s3://1000genomes/phase3/data/NA19723/sequence_read/SRR032770_1.filt.fastq.gz	s3://1000genomes/phase3/data/NA19723/sequence_read/SRR032770_2.filt.fastq.gz	Solexa-16043	BI.PE.091118_SL-XBD_0005_FC43265AAXX.7	2010-01-05T00:00:00Z	ILLUMINA	BI
NA19722_A	NA19722	s3://1000genomes/phase3/data/NA19722/sequence_read/SRR032772_1.filt.fastq.gz	s3://1000genomes/phase3/data/NA19722/sequence_read/SRR032772_2.filt.fastq.gz	Solexa-16042	BI.PE.091203_SL-XAM_0006_FC4328HAAXX.1	2010-01-05T00:00:00Z	ILLUMINA	BI
```


## Per-Sample Variant Calling

```
nextflow run aws-nextflow-gatk/processing-for-variant-discovery-gatk4.nf \
  -work-dir s3://nf-work-bucket-us-east-2/work \
  --outdir s3://<your_bucket>/results/variant-discovery/ \
  --unmapped_bams_list s3://<your_bucket>/results/fastq-to-ubam/unmapped_bams.tsv \
  --projectId my-project
```

### Inputs
- work-dir: nextflow work directory. Must be a location in s3
- outdir: location in s3 where outputs are saved
- projectId: used as directory name for output vcf files
- unmapped_bams_list: a tab-delimited file with sample names and ubam locations

#### Example Input: 

```
NA19771	s3://<your_bucket>/unmapped_bams/NA19771_A.unmapped.bam
NA19654	s3://<your_bucket>/unmapped_bams/NA19654_A.unmapped.bam
NA19731	s3://<your_bucket>/unmapped_bams/NA19731_A.unmapped.bam
```

## Joint Genotyping

```
nextflow run aws-nextflow-gatk/JointGenotyping.nf \
  -work-dir s3://<your_bucket>/work \
  --outdir s3://<your_bucket>/results/jointGenotyping/ \
  --input_vcf_path s3://<your_bucket>/results/variant-discovery/my-project/vcfs \
  --callset_name my-callset \
  --manifest s3://<your_bucket>/results/variant-discovery/merged_vcf_out.txt
```

### Inputs:
- work-dir: nextflow work directory. Must be a location in s3
- outdir: location in s3 where outputs are saved
- input_vcf_path: directory containing vcfs from per-sample variant calling
- callset_name: name for directory and output vcf
- manifest: tab-delimited file with sample names and vcf file locations from per-sample variant calling step


## Format Conversion + Variant calling in one-step (BAM-to-VCF)
A single workflow that combines `Paired-fastqs to uBAM` and `Per-Sample Variant Calling`.

Takes fastq files as inputs and outputs vcf files. 

```
nextflow run aws-nextflow-gatk/fastq_to_variant_calls.nf \
  -work-dir s3://<your_bucket>/work \
  --outdir s3://<your_bucket>/results/fastq-to-ubam/ \
  --input_fofn s3://<your_bucket>/input_files/fastq_manifest.txt \
  --projectId my-project
```

### Inputs:
- work-dir: nextflow work directory. Must be a location in s3
- outdir: location in s3 where outputs are saved
- input_fofn: an input manifest file of tab-delimited values as in `Paired-fastqs to uBAM`
- projectId: used as directory name for output vcf files
