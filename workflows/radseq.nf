/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRadseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config , params.genome, params.popmap]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                            } from '../subworkflows/local/input_check'
include { PROCESS_RAD                            } from '../subworkflows/local/fastq_fastp'
include { CDHIT_RAINBOW as DENOVO                } from '../subworkflows/local/fastq_cdhit_rainbow'
include { FASTQ_INDEX_ALIGN_BWA_MINIMAP as ALIGN } from '../subworkflows/local/fastq_index_align_bwa_minimap'
include { BAM_INTERVALS_BEDTOOLS                 } from '../subworkflows/local/bam_intervals_bedtools'
include { BAM_VARIANT_CALLING_FREEBAYES          } from '../subworkflows/local/bam_variant_calling_freebayes'
include { FASTQ_UNZIP                            } from '../subworkflows/local/fastq_unzip.nf'
include { VCF_BCFTOOLS_RADSEQ_FILTERS            } from '../subworkflows/local/vcf_bcftools_radseq_filters.nf'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RADSEQ {

    ch_versions = Channel.empty()
    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: remove/trim low quality reads, trim umi's
    //
    PROCESS_RAD (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(PROCESS_RAD.out.versions)

    // assign fasta channel based on method in config file
    switch ( params.method ) {
        // assign ch_reference (input for aligning subworkflow) to the reference in the params
        case 'reference':
            //
            // SUBWORKFLOW: unzips fasta if ends with .gz
            //
            ch_reference = FASTQ_UNZIP (params.genome).fasta
                .combine(PROCESS_RAD.out.trimmed_reads.map{it[0]}.first())
                .map{ WorkflowRadseq.addRefIdToChannels(params, it) }
            break
        case 'denovo':
            /* SUBWORKFLOW: Cluster READS after applying unique read thresholds within and among samples.
            *   option to provide a list of minimum depth thresholds. See nextflow.config for more details*/
            ch_reference = DENOVO (
                INPUT_CHECK.out.reads, 
                params.sequence_type // sequence type exe.: 'SE', 'PE', ''
            ).fasta
            ch_versions = ch_versions.mix(DENOVO.out.versions)
            break
        default:
            exit 1, "unknown method: ${method} \n supported options:" + params.method_options
    }

    // nf-core module index reference for bedtools + freebayes
    if (params.faidx) {
        ch_faidx = Channel.fromPath(params.faidx)
            .combine(PROCESS_RAD.out.trimmed_reads.map{it[0]}.first())
            .map{ WorkflowRadseq.addRefIdToChannels(params, it) }
    } else {
        ch_faidx = SAMTOOLS_FAIDX (
        ch_reference
        ).fai
            .combine(PROCESS_RAD.out.trimmed_reads.map{it[0]}.first())
            .map{ WorkflowRadseq.addRefIdToChannels(params, it) }
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    // embed function to have reference match bed and mbam_bai me

    //
    // SUBWORKFLOW: generate fasta indexes, align input files, dedup reads, index bam, calculate statistics
    //      if denovo and paired then pass length_stats to bwa mem
    ch_bam_bai = ALIGN (
        PROCESS_RAD.out.trimmed_reads, 
        ch_reference, 
        ch_faidx,
        params.sequence_type, 
        PROCESS_RAD.out.read_lengths
        ).bam_bai
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    //
    // SUBWORKFLOW: freebayes multithreading based on read coverage or genome samtools faidx index
    //
    ch_intervals = BAM_INTERVALS_BEDTOOLS (
        ALIGN.out.bam,
        ch_faidx,
        PROCESS_RAD.out.read_lengths,
        params.max_read_coverage_to_split
        ).intervals
    ch_versions = ch_versions.mix(BAM_INTERVALS_BEDTOOLS.out.versions)

    // expand channel across bed regions for variant calling multi-threading
    ch_bam_bai_bed_fasta_fai = ALIGN.out.mbam_bai
        .combine(ch_intervals, by: 0)
        .combine(ch_reference, by: [0])
        .combine(ch_faidx, by: [0])
        .map { meta, bam, bai, bed, fasta, fai -> 
            [[
                id:           meta.id,
                single_end:   meta.single_end,
                interval:     bed.getName().tokenize( '.' )[0],
                ref_id:       meta.ref_id
            ],
            bam, bai, file(bed), fasta, fai]
        }

    //
    // SUBWORKFLOW: freebayes parallel variant calling
    //
    ch_vcf = BAM_VARIANT_CALLING_FREEBAYES (
        ch_bam_bai_bed_fasta_fai,
        true
    ).vcf
    ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)

    ch_vcf_tbi_fasta = ch_vcf
        .join(BAM_VARIANT_CALLING_FREEBAYES.out.tbi,by:0)
        .join(ch_reference, by:0)

    //
    // SUBWORKFLOW: Apply RADseq Specific Filters (depth, missingness, allele counts)
    //
    VCF_BCFTOOLS_RADSEQ_FILTERS ( 
        ch_vcf_tbi_fasta
        )
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRadseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_RAD.out.fastp_json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.idxstats.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(BAM_VARIANT_CALLING_FREEBAYES.out.stats.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        [],
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
