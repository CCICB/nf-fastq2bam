//
// Align DNA reads
//

import Constants
import Utils

include { BWAMEM2_ALIGN  } from '../../../modules/local/bwa-mem2/mem/main'
include { FASTP          } from '../../../modules/local/fastp/main'
include { GATK4_MARKDUPLICATES as GATK4_MARKDUPLICATES_TUMOR } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_MARKDUPLICATES as GATK4_MARKDUPLICATES_NORMAL } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_MARKDUPLICATES as GATK4_MARKDUPLICATES_DONOR } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TUMOR } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NORMAL } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DONOR } from '../../../modules/nf-core/samtools/sort/main'
include { SAMBAMBA_MERGE as SAMBAMBA_MERGE_TUMOR } from '../../../modules/local/sambamba/merge/main'
include { SAMBAMBA_MERGE as SAMBAMBA_MERGE_NORMAL } from '../../../modules/local/sambamba/merge/main'
include { SAMBAMBA_MERGE as SAMBAMBA_MERGE_DONOR } from '../../../modules/local/sambamba/merge/main'

workflow READ_ALIGNMENT_DNA {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/

    // Params
    max_fastq_records    // numeric: [optional]  max number of FASTQ records per split
    umi_enable           // boolean: [mandatory] enable UMI processing
    umi_location         //  string: [optional]  fastp UMI location argument (--umi_loc)
    umi_length           // numeric: [optional]  fastp UMI length argument (--umi_len)
    umi_skip             // numeric: [optional]  fastp UMI skip argument (--umi_skip)

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs, separate by tumor and normal
    // channel: [ meta ]
    ch_inputs_tumor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_TUMOR)
            runnable: Utils.hasTumorDnaFastq(meta) && !has_existing
            skip: true
        }

    ch_inputs_normal_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_NORMAL)
            runnable: Utils.hasNormalDnaFastq(meta) && !has_existing
            skip: true
        }

    ch_inputs_donor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_DONOR)
            runnable: Utils.hasDonorDnaFastq(meta) && !has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta -> [meta, Utils.getTumorDnaSample(meta), 'tumor'] },
            ch_inputs_normal_sorted.runnable.map { meta -> [meta, Utils.getNormalDnaSample(meta), 'normal'] },
            ch_inputs_donor_sorted.runnable.map { meta -> [meta, Utils.getDonorDnaSample(meta), 'donor'] },
        )
        .flatMap { meta, meta_sample, sample_type ->
            meta_sample
                .getAt(Constants.FileType.FASTQ)
                .collect { key, fps ->
                    def (library_id, lane) = key

                    def meta_fastq = [
                        key: meta.group_id,
                        id: "${meta.group_id}_${meta_sample.sample_id}",
                        sample_id: meta_sample.sample_id,
                        library_id: library_id,
                        lane: lane,
                        sample_type: sample_type,
                    ]

                    return [meta_fastq, fps['fwd'], fps['rev']]
                }
        }

    //
    // MODULE: fastp
    //
    // Split FASTQ into chunks if requested for distributed processing
    // channel: [ meta_fastq_ready, fastq_fwd, fastq_fwd ]
    ch_fastqs_ready = Channel.empty()
    if (max_fastq_records > 0 || umi_enable) {

        // Run process
        FASTP(
            ch_fastq_inputs,
            max_fastq_records,
            umi_location,
            umi_length,
            umi_skip,
        )

        ch_versions = ch_versions.mix(FASTP.out.versions)

    }

    // Now prepare according to FASTQs splitting
    if (max_fastq_records > 0) {

        ch_fastqs_ready = FASTP.out.fastq
            .flatMap { meta_fastq, reads_fwd, reads_rev ->

                def data = [reads_fwd, reads_rev]
                    .transpose()
                    .collect { fwd, rev ->

                        def split_fwd = fwd.name.replaceAll('\\..+$', '')
                        def split_rev = rev.name.replaceAll('\\..+$', '')

                        assert split_fwd == split_rev

                        // NOTE(SW): split allows meta_fastq_ready to be unique, which is required during reunite below
                        def meta_fastq_ready = [
                            *:meta_fastq,
                            id: "${meta_fastq.id}_${split_fwd}",
                            split: split_fwd,
                        ]

                        return [meta_fastq_ready, fwd, rev]
                    }

                return data
            }

    } else {

        // Select appropriate source
        ch_fastq_source = umi_enable ? FASTP.out.fastq : ch_fastq_inputs

        ch_fastqs_ready = ch_fastq_source
            .map { meta_fastq, fastq_fwd, fastq_rev ->

                def meta_fastq_ready = [
                    *:meta_fastq,
                    split: null,
                ]

                return [meta_fastq_ready, fastq_fwd, fastq_rev]
            }

    }

    //
    // MODULE: BWA-MEM2
    //
    // Create process input channel
    // channel: [ meta_bwamem2, fastq_fwd, fastq_rev ]
    ch_bwamem2_inputs = ch_fastqs_ready
        .map { meta_fastq_ready, fastq_fwd, fastq_rev ->

            def meta_bwamem2 = [
                *:meta_fastq_ready,
                read_group: "${meta_fastq_ready.sample_id}.${meta_fastq_ready.library_id}.${meta_fastq_ready.lane}",
            ]

            return [meta_bwamem2, fastq_fwd, fastq_rev]
        }

    // Run process
    BWAMEM2_ALIGN(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    ch_versions = ch_versions.mix(BWAMEM2_ALIGN.out.versions)

    // Branch the output by sample type for further processing
    // channel: [ meta_bwamem2, bam ]
    ch_bwamem2_output_by_type = BWAMEM2_ALIGN.out.bam
        .branch { meta_bwamem2, bam ->
            tumor:  meta_bwamem2.sample_type == 'tumor'
            normal: meta_bwamem2.sample_type == 'normal'
            donor:  meta_bwamem2.sample_type == 'donor'
        }

    //
    // MODULE: Sambamba merge
    //
    // First, count expected BAMs per sample for non-blocking groupTuple op
    // Count for tumor samples
    // channel: [ meta_count, group_size ]
    ch_tumor_fastq_counts = ch_bwamem2_inputs
        .filter { meta_bwamem2, reads_fwd, reads_rev -> meta_bwamem2.sample_type == 'tumor' }
        .map { meta_bwamem2, reads_fwd, reads_rev ->
            def meta_count = [key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type]
            return [meta_count, meta_bwamem2]
        }
        .groupTuple()
        .map { meta_count, meta_bwamem2s -> return [meta_count, meta_bwamem2s.size()] }

    // Count for normal samples
    ch_normal_fastq_counts = ch_bwamem2_inputs
        .filter { meta_bwamem2, reads_fwd, reads_rev -> meta_bwamem2.sample_type == 'normal' }
        .map { meta_bwamem2, reads_fwd, reads_rev ->
            def meta_count = [key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type]
            return [meta_count, meta_bwamem2]
        }
        .groupTuple()
        .map { meta_count, meta_bwamem2s -> return [meta_count, meta_bwamem2s.size()] }

    // Count for donor samples
    ch_donor_fastq_counts = ch_bwamem2_inputs
        .filter { meta_bwamem2, reads_fwd, reads_rev -> meta_bwamem2.sample_type == 'donor' }
        .map { meta_bwamem2, reads_fwd, reads_rev ->
            def meta_count = [key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type]
            return [meta_count, meta_bwamem2]
        }
        .groupTuple()
        .map { meta_count, meta_bwamem2s -> return [meta_count, meta_bwamem2s.size()] }

    // Now, group with expected size for tumor
    // channel: [ meta_group, [bam, ...] ]
    ch_tumor_bams_united = ch_tumor_fastq_counts
        .cross(
            // First element to match meta_count above for `cross`
            ch_bwamem2_output_by_type.tumor.map { meta_bwamem2, bam -> 
                [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], bam] 
            }
        )
        .map { count_tuple, bam_tuple ->
            def group_size = count_tuple[1]
            def (meta_bam, bam) = bam_tuple

            def meta_group = [
                *:meta_bam,
            ]

            return tuple(groupKey(meta_group, group_size), bam)
        }
        .groupTuple()

    // For normal samples
    ch_normal_bams_united = ch_normal_fastq_counts
        .cross(
            ch_bwamem2_output_by_type.normal.map { meta_bwamem2, bam -> 
                [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], bam] 
            }
        )
        .map { count_tuple, bam_tuple ->
            def group_size = count_tuple[1]
            def (meta_bam, bam) = bam_tuple

            def meta_group = [
                *:meta_bam,
            ]

            return tuple(groupKey(meta_group, group_size), bam)
        }
        .groupTuple()

    // For donor samples
    ch_donor_bams_united = ch_donor_fastq_counts
        .cross(
            ch_bwamem2_output_by_type.donor.map { meta_bwamem2, bam -> 
                [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], bam] 
            }
        )
        .map { count_tuple, bam_tuple ->
            def group_size = count_tuple[1]
            def (meta_bam, bam) = bam_tuple

            def meta_group = [
                *:meta_bam,
            ]

            return tuple(groupKey(meta_group, group_size), bam)
        }
        .groupTuple()

    // Sort into merge-eligible BAMs (at least two BAMs required)
    // channel: runnable: [ meta_group, [bam, ...] ]
    // channel: skip: [ meta_group, bam ]
    ch_tumor_bams_united_sorted = ch_tumor_bams_united
        .branch { meta_group, bams ->
            runnable: bams.size() > 1
            skip: true
                return [meta_group, bams[0]]
        }

    ch_normal_bams_united_sorted = ch_normal_bams_united
        .branch { meta_group, bams ->
            runnable: bams.size() > 1
            skip: true
                return [meta_group, bams[0]]
        }

    ch_donor_bams_united_sorted = ch_donor_bams_united
        .branch { meta_group, bams ->
            runnable: bams.size() > 1
            skip: true
                return [meta_group, bams[0]]
        }

    // Create process input channel for merge
    // channel: [ meta_merge, [bams, ...] ]
    ch_merge_tumor_inputs = WorkflowFastq2bam.restoreMeta(ch_tumor_bams_united_sorted.runnable, ch_inputs)
        .map { meta, bams ->
            def meta_merge = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
                sample_type: 'tumor'
            ]
            return [meta_merge, bams]
        }

    ch_merge_normal_inputs = WorkflowFastq2bam.restoreMeta(ch_normal_bams_united_sorted.runnable, ch_inputs)
        .map { meta, bams ->
            def meta_merge = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getNormalDnaSampleName(meta),
                sample_type: 'normal'
            ]
            return [meta_merge, bams]
        }

    ch_merge_donor_inputs = WorkflowFastq2bam.restoreMeta(ch_donor_bams_united_sorted.runnable, ch_inputs)
        .map { meta, bams ->
            def meta_merge = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getDonorDnaSampleName(meta),
                sample_type: 'donor'
            ]
            return [meta_merge, bams]
        }

    // Run merge processes
    SAMBAMBA_MERGE_TUMOR(ch_merge_tumor_inputs)
    SAMBAMBA_MERGE_NORMAL(ch_merge_normal_inputs)
    SAMBAMBA_MERGE_DONOR(ch_merge_donor_inputs)

    ch_versions = ch_versions.mix(
        SAMBAMBA_MERGE_TUMOR.out.versions.first(), 
        SAMBAMBA_MERGE_NORMAL.out.versions.first(),
        SAMBAMBA_MERGE_DONOR.out.versions.first()
    )

    //
    // MODULE: SAMtools sort
    //
    // Create process input channel for each sample type
    // channel: [ meta_sort, bam ]
    ch_sort_tumor_inputs = Channel.empty()
        .mix(
            WorkflowFastq2bam.restoreMeta(SAMBAMBA_MERGE_TUMOR.out.bam, ch_inputs),
            WorkflowFastq2bam.restoreMeta(ch_tumor_bams_united_sorted.skip, ch_inputs)
        )
        .map { meta, bam ->
            def meta_sort = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
                prefix: Utils.getTumorDnaSampleName(meta),
                sample_type: 'tumor'
            ]
            return [meta_sort, bam]
        }

    ch_sort_normal_inputs = Channel.empty()
        .mix(
            WorkflowFastq2bam.restoreMeta(SAMBAMBA_MERGE_NORMAL.out.bam, ch_inputs),
            WorkflowFastq2bam.restoreMeta(ch_normal_bams_united_sorted.skip, ch_inputs)
        )
        .map { meta, bam ->
            def meta_sort = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getNormalDnaSampleName(meta),
                prefix: Utils.getNormalDnaSampleName(meta),
                sample_type: 'normal'
            ]
            return [meta_sort, bam]
        }

    ch_sort_donor_inputs = Channel.empty()
        .mix(
            WorkflowFastq2bam.restoreMeta(SAMBAMBA_MERGE_DONOR.out.bam, ch_inputs),
            WorkflowFastq2bam.restoreMeta(ch_donor_bams_united_sorted.skip, ch_inputs)
        )
        .map { meta, bam ->
            def meta_sort = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getDonorDnaSampleName(meta),
                prefix: Utils.getDonorDnaSampleName(meta),
                sample_type: 'donor'
            ]
            return [meta_sort, bam]
        }
        

    // Run sort processes
    SAMTOOLS_SORT_TUMOR(ch_sort_tumor_inputs)
    SAMTOOLS_SORT_NORMAL(ch_sort_normal_inputs)
    SAMTOOLS_SORT_DONOR(ch_sort_donor_inputs)

    ch_versions = ch_versions.mix(
        SAMTOOLS_SORT_TUMOR.out.versions.first(),
        SAMTOOLS_SORT_NORMAL.out.versions.first(),
        SAMTOOLS_SORT_DONOR.out.versions.first()
    )

    //
    // MODULE: GATK4 markduplicates
    //
    // Create process input channel for tumor
    // channel: [ meta_markdups, bam ]
    ch_markdups_tumor_inputs = WorkflowFastq2bam.restoreMeta(SAMTOOLS_SORT_TUMOR.out.bam, ch_inputs)
        .map { meta, bam ->
            def meta_markdups = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
                sample_type: 'tumor'
            ]
            return [meta_markdups, bam]
        }

    ch_markdups_normal_inputs = WorkflowFastq2bam.restoreMeta(SAMTOOLS_SORT_NORMAL.out.bam, ch_inputs)
        .map { meta, bam ->
            def meta_markdups = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getNormalDnaSampleName(meta),
                sample_type: 'normal'
            ]
            return [meta_markdups, bam]
        }

    ch_markdups_donor_inputs = WorkflowFastq2bam.restoreMeta(SAMTOOLS_SORT_DONOR.out.bam, ch_inputs)
        .map { meta, bam ->
            def meta_markdups = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getDonorDnaSampleName(meta),
                sample_type: 'donor'
            ]
            return [meta_markdups, bam]
        }

    // Run markduplicates processes
    GATK4_MARKDUPLICATES_TUMOR(
        ch_markdups_tumor_inputs,
        [],
        [],
    )

    GATK4_MARKDUPLICATES_NORMAL(
        ch_markdups_normal_inputs,
        [],
        [],
    )

    GATK4_MARKDUPLICATES_DONOR(
        ch_markdups_donor_inputs,
        [],
        [],
    )

    ch_versions = ch_versions.mix(
        GATK4_MARKDUPLICATES_TUMOR.out.versions.first(),
        GATK4_MARKDUPLICATES_NORMAL.out.versions.first(),
        GATK4_MARKDUPLICATES_DONOR.out.versions.first()
    )

    // Combine BAMs and BAIs for tumor
    // channel: [ meta, bam, bai ]
    ch_tumor_bams_ready = WorkflowFastq2bam.groupByMeta(
        WorkflowFastq2bam.restoreMeta(GATK4_MARKDUPLICATES_TUMOR.out.bam, ch_inputs),
        WorkflowFastq2bam.restoreMeta(GATK4_MARKDUPLICATES_TUMOR.out.bai, ch_inputs),
    )

    // Combine BAMs and BAIs for normal
    ch_normal_bams_ready = WorkflowFastq2bam.groupByMeta(
        WorkflowFastq2bam.restoreMeta(GATK4_MARKDUPLICATES_NORMAL.out.bam, ch_inputs),
        WorkflowFastq2bam.restoreMeta(GATK4_MARKDUPLICATES_NORMAL.out.bai, ch_inputs),
    )

    // Combine BAMs and BAIs for donor
    ch_donor_bams_ready = WorkflowFastq2bam.groupByMeta(
        WorkflowFastq2bam.restoreMeta(GATK4_MARKDUPLICATES_DONOR.out.bam, ch_inputs),
        WorkflowFastq2bam.restoreMeta(GATK4_MARKDUPLICATES_DONOR.out.bai, ch_inputs),
    )

    // Set outputs
    // channel: [ meta, bam, bai ]
    ch_bam_tumor_out = Channel.empty()
        .mix(
            ch_tumor_bams_ready,
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_bam_normal_out = Channel.empty()
        .mix(
            ch_normal_bams_ready,
            ch_inputs_normal_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_bam_donor_out = Channel.empty()
        .mix(
            ch_donor_bams_ready,
            ch_inputs_donor_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    dna_tumor  = ch_bam_tumor_out  // channel: [ meta, bam, bai ]
    dna_normal = ch_bam_normal_out // channel: [ meta, bam, bai ]
    dna_donor  = ch_bam_donor_out  // channel: [ meta, bam, bai ]

    versions   = ch_versions       // channel: [ versions.yml ]
}