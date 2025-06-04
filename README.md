## Introduction

Many bioinformatic tools require a BAM file as input. **nf-fastq2bam** is a nextflow workflow that takes FASTQ files as input and outputs a single mapped BAM for both Whole genome sequencing(WGS) as well as RNA sequencing experiments.  
**nf-fastq2bam** runs the following steps from FASTQ files:

- TUMOR/NORMAL WGS FASTQ -> BWAMEM2_ALIGN -> SAMBAMBA_MERGE -> SAMTOOLS_SORT -> GATK4_MARKDUPLICATES 

- RNA FASTQ -> STAR_ALIGN -> SAMBAMBA_MERGE -> SAMTOOLS_SORT -> GATK4_MARKDUPLICATES 

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.


Create a `samplesheet.csv` with your inputs (WGS/WTS BAMs in this example):

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
PATIENT1_WGTS,PATIENT1,PATIENT1_tumoursample,tumor,dna,fastq,library_id:HH5FYCCXY_library;lane:6,/path/2/data/patient1_tumoursample/HH5FYCCXY_6_180321_FD01114733_R1.fastq.gz;/path/2/data/patient1_tumoursample/HH5FYCCXY_6_180321_FD01114733_R2.fastq.gz
PATIENT1_WGTS,PATIENT1,PATIENT1_tumoursample,tumor,dna,fastq,library_id:HH5FYCCXY_library;lane:5,/path/2/data/patient1_tumoursample/HH5FYCCXY_5_180321_FD01114733_R1.fastq.gz;/path/2/data/patient1_tumoursample/HH5FYCCXY_5_180321_FD01114733_R2.fastq.gz
PATIENT1_WGTS,PATIENT1,PATIENT1_normalsample,normal,dna,fastq,library_id:HH5FYCCXZ_library;lane:7,/path/2/data/patient1_normalsample/HH5FYCCXZ_7_180321_FD01114735_R1.fastq.gz;/path/2/data/patient1_normalsample/HH5FYCCXZ_7_180321_FD01114735_R2.fastq.gz
PATIENT1_WGTS,PATIENT1,PATIENT1_rnasample,tumor,rna,fastq,library_id:HH5FYCCXZ_library;lane:9,/path/2/data/patient1_rnasample/HH5FYCCXZ_9_180321_FD01114735_R1.fastq.gz;/path/2/data/patient1_rnasample/HH5FYCCXZ_9_180321_FD01114735_R2.fastq.gz
```

Your `params.yaml` file:

```
ref_data_genome_fasta           : /scratch/df13/resources/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta
ref_data_genome_fai             : /scratch/df13/resources/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.fai
ref_data_genome_dict            : /scratch/df13/resources/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.dict
ref_data_genome_bwa_index_image : /scratch/df13/resources/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.img
ref_data_genome_gridss_index    : /scratch/df13/resources/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.gridsscache
```

Your `run.config` file for HPC environments, e.g. NCI Gadi:

```
params {
  // Tip: use conf/hmf_genomes.config as a guide
  config_profile_description = 'NCI Gadi HPC profile'
  config_profile_url = 'https://opus.nci.org.au/display/Help/Gadi+User+Guide'
  genomes {
		'GRCh38_hmf' {
			fasta         = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta"
			fai           = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.fai"
			dict          = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.dict"
			img           = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.img"
			gridss_index  = "/path/2/hmf_resources/38/gridss_index-2.13.2.tar.gz"
            bwamem2_index = "/path/2/hmf_resources/38/bwa-mem2_index-2.2.1.tar.gz"
            star_index    = "/path/2/hmf_resources/38/star_index-gencode_38-2.7.3a.tar.gz"
		}
	}

	ref_data_hmf_data_path        = "/path/2/hmf_resources/38/hmf_pipeline_resources.38_v2.0--3"
	ref_data_panel_data_path      = "/path/2/hmf_resources/38/panel/tso500/hmf_panel_resources.tso500.38_v2.1.0--3"
}

singularity {
    enabled = true
    cacheDir = '/path/2/cache'
    autoMounts = true
}

executor {
  queueSize = 200
  pollInterval = '5 min'
  queueStatInterval = '5 min'
  submitRateLimit = '20 min'
}

process {
  resourceLimits = [ cpus: 32, memory: 1020.GB, time: 48.h ]
  errorStrategy = 'retry'
  maxRetries = 3
  executor = 'pbspro'
  project = '<YOUR PROJECT ID>'
  module = 'singularity'
  cache = 'lenient'
  stageInMode = 'symlink'
}

```

Your `run.config` file if running in a local environment:

```
params {
	genomes {
		'GRCh38_hmf' {
			fasta         = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta"
			fai           = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.fai"
			dict          = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.dict"
			img           = "/path/2/hmf_resources/38/hg38_alts_decoys_phiX_masked_GRC.fasta.img"
			gridss_index  = "/path/2/hmf_resources/38/gridss_index-2.13.2.tar.gz"
            bwamem2_index = "/path/2/hmf_resources/38/bwa-mem2_index-2.2.1.tar.gz"
            star_index    = "/path/2/hmf_resources/38/star_index-gencode_38-2.7.3a.tar.gz"
		}
	}
	ref_data_hmf_data_path        = "/path/2/hmf_resources/38/hmf_pipeline_resources.38_v2.0--3"
	ref_data_panel_data_path      = "/path/2/hmf_resources/38/panel/tso500/hmf_panel_resources.tso500.38_v2.1.0--3"
}

process {
    //  defaults for all processes.
    resourceLimits = [ cpus: 32, memory: 1020.GB, time: 48.h ]
    cpus   = { check_max( 1    * task.attempt, 'cpus'   ) }
    memory = { check_max( 6.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h  * task.attempt, 'time'   ) }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    // The 'local' executor is used by default. Please see here for other executor options https://www.nextflow.io/docs/latest/executor.html#executors
    // such as AWS batch, AZURE batch, PBS, Slurm, etc.
    executor = 'local' 
}

```

Recommended `resources.config`:
```
process {
  withName: BAM_TOOLS {
    memory =  { 24.GB * task.attempt }  
    cpus   =  { 16    * task.attempt }
  }
  withName: 'BWAMEM2_ALIGN.*' {
    cpus = 4
    memory =  { 64.GB * task.attempt }  

  }
  withName: 'FASTP.*' {
    memory = { 64.GB * task.attempt }  
    cpus = 16
  }
}
```

Launch `nf-fastq2bam`:

```bash
git clone git@github.com:CCICB/nf-fastq2bam.git
```

```bash
nextflow run nf-fastq2bam \
  -profile <docker|singularity|...> \
  -c <your run.config and resources.config files here, separated by commas. You may also list other config files.> \
  -params-file <your params.yaml> \
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
