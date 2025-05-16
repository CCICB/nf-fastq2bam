import Constants
import Processes
import Utils


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse input samplesheet
// NOTE(SW): this is done early and outside of gpars so that we can access synchronously and prior to pipeline execution
inputs = Utils.parseInput(params.input, workflow.stubRun, log)

// Get run config
run_config = WorkflowMain.getRunConfig(params, inputs, log)

// Validate inputs
Utils.validateInput(inputs, run_config, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.isofox_counts,
    params.isofox_gc_ratios,
    params.isofox_gene_ids,
    params.isofox_tpm_norm,
]

if (run_config.stages.lilac) {
    if (params.genome_version.toString() == '38' && params.genome_type == 'alt' && params.containsKey('ref_data_hla_slice_bed')) {
        checkPathParamList.add(params.ref_data_hla_slice_bed)
    }
}

// TODO(SW): consider whether we should check for null entries here for errors to be more informative
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Used in Isofox subworkflow only
isofox_read_length = params.isofox_read_length !== null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH_TARGETED

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA    } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA    } from '../subworkflows/local/read_alignment_rna'
include { REDUX_PROCESSING      } from '../subworkflows/local/redux_processing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet = Utils.getFileObject(params.input)

workflow TARGETED {
    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = Channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data
    panel_data = PREPARE_REFERENCE.out.panel_data

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_align_dna_tumor_out = Channel.empty()
    ch_align_dna_normal_out = Channel.empty()
    ch_align_dna_donor_out = Channel.empty()
    ch_align_rna_tumor_out = Channel.empty()
    if (run_config.stages.alignment) {

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            params.max_fastq_records,
            params.fastp_umi,
            params.fastp_umi_location,
            params.fastp_umi_length,
            params.fastp_umi_skip,
        )

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ref_data.genome_star_index,
        )

        ch_versions = ch_versions.mix(
            READ_ALIGNMENT_DNA.out.versions,
            READ_ALIGNMENT_RNA.out.versions,
        )

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.dna_tumor)
        ch_align_dna_normal_out = ch_align_dna_normal_out.mix(READ_ALIGNMENT_DNA.out.dna_normal)
        ch_align_dna_donor_out = ch_align_dna_donor_out.mix(READ_ALIGNMENT_DNA.out.dna_donor)
        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)

    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_rna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }

    }

    //
    // SUBWORKFLOW: Run REDUX for DNA BAMs
    //
    // channel: [ meta, bam, bai ]
    ch_redux_dna_tumor_out = Channel.empty()
    ch_redux_dna_normal_out = Channel.empty()
    ch_redux_dna_donor_out = Channel.empty()

    // channel: [ meta, dup_freq_tsv, jitter_tsv, ms_tsv ]
    ch_redux_dna_tumor_tsv_out = Channel.empty()
    ch_redux_dna_normal_tsv_out = Channel.empty()
    ch_redux_dna_donor_tsv_out = Channel.empty()

    if (run_config.stages.redux) {

        REDUX_PROCESSING(
            ch_inputs,
            ch_align_dna_tumor_out,
            ch_align_dna_normal_out,
            ch_align_dna_donor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.unmap_regions,
            hmf_data.msi_jitter_sites,
            params.redux_umi,
            params.redux_umi_duplex_delim,
        )

        ch_versions = ch_versions.mix(REDUX_PROCESSING.out.versions)

        ch_redux_dna_tumor_out = ch_redux_dna_tumor_out.mix(REDUX_PROCESSING.out.dna_tumor)
        ch_redux_dna_normal_out = ch_redux_dna_normal_out.mix(REDUX_PROCESSING.out.dna_normal)
        ch_redux_dna_donor_out = ch_redux_dna_donor_out.mix(REDUX_PROCESSING.out.dna_donor)

        ch_redux_dna_tumor_tsv_out = ch_redux_dna_tumor_tsv_out.mix(REDUX_PROCESSING.out.dna_tumor_tsv)
        ch_redux_dna_normal_tsv_out = ch_redux_dna_normal_tsv_out.mix(REDUX_PROCESSING.out.dna_normal_tsv)
        ch_redux_dna_donor_tsv_out = ch_redux_dna_donor_tsv_out.mix(REDUX_PROCESSING.out.dna_donor_tsv)

    } else {

        ch_redux_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_redux_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_redux_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }

        ch_redux_dna_tumor_tsv_out = ch_inputs.map { meta -> [meta, [], [], [], []] }
        ch_redux_dna_normal_tsv_out = ch_inputs.map { meta -> [meta, [], [], [], []] }
        ch_redux_dna_donor_tsv_out = ch_inputs.map { meta -> [meta, [], [], [], []] }

    }

    //
    // TASK: Aggregate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true,
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
