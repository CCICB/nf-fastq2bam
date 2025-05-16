## Introduction

Many bioinformatic tools require a BAM file as input. **nf_fastq2bam** is a nextflow workflow that takes FASTQ files as input and outputs a single mapped BAM for both Whole genome sequencing(WGS) as well as RNA sequencing experiments.  
**nf_fastq2bam** runs the following steps from FASTQ files:

- TUMOR/NORMAL WGS FASTQ -> BWAMEM2_ALIGN -> SAMBAMBA_MERGE -> SAMTOOLS_SORT -> GATK4_MARKDUPLICATES 

- RNA FASTQ -> STAR_ALIGN -> SAMBAMBA_MERGE -> SAMTOOLS_SORT -> GATK4_MARKDUPLICATES 

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.


Create a samplesheet with your inputs (WGS/WTS BAMs in this example):

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
PATIENT1_WGTS,PATIENT1,PATIENT1_tumoursample,tumor,dna,fastq,library_id:HH5FYCCXY_library;lane:6,/path/2/data/patient1_tumoursample/HH5FYCCXY_6_180321_FD01114733_R1.fastq.gz;/path/2/data/patient1_tumoursample/HH5FYCCXY_6_180321_FD01114733_R2.fastq.gz
PATIENT1_WGTS,PATIENT1,PATIENT1_tumoursample,tumor,dna,fastq,library_id:HH5FYCCXY_library;lane:5,/path/2/data/patient1_tumoursample/HH5FYCCXY_5_180321_FD01114733_R1.fastq.gz;/path/2/data/patient1_tumoursample/HH5FYCCXY_5_180321_FD01114733_R2.fastq.gz
PATIENT1_WGTS,PATIENT1,PATIENT1_normalsample,normal,dna,fastq,library_id:HH5FYCCXZ_library;lane:7,/path/2/data/patient1_normalsample/HH5FYCCXZ_7_180321_FD01114735_R1.fastq.gz;/path/2/data/patient1_normalsample/HH5FYCCXZ_7_180321_FD01114735_R2.fastq.gz
PATIENT1_WGTS,PATIENT1,PATIENT1_rnasample,tumor,rna,fastq,library_id:HH5FYCCXZ_library;lane:9,/path/2/data/patient1_rnasample/HH5FYCCXZ_9_180321_FD01114735_R1.fastq.gz;/path/2/data/patient1_rnasample/HH5FYCCXZ_9_180321_FD01114735_R2.fastq.gz
```

Launch `nf_fastq2bam`:

```bash
git clone git@bitbucket.org:cciacb/nf_fastq2bam.git
```

```bash
nextflow run nf_fastq2bam \
  -profile <docker|singularity|...> \
  -c <your nf-fastq2bam config files> \
  -params-file <your nf-fastq2bam params file> \
  --mode <wgts|targeted> \
  --genome_type <alt|no_alt> \
  --genome <GRCh37_hmf|GRCh38_hmf> \
  --input samplesheet.csv \
  --outdir output/
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).


## Pipeline output

The output of this pipeline contains processed BAM files for different sample types. This includes aligned, merged, sorted, and duplicate-marked BAM files from tumour, normal and RNA-Seq samples.



## Citations


You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia,
> Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).