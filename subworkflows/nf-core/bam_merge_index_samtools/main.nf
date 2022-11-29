//
// MERGE INDEX BAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_MERGE as MERGE_BAM       } from '../../../modules/nf-core/samtools/merge/main.nf'

workflow BAM_MERGE_INDEX_SAMTOOLS {
    take:
        bam // channel: [mandatory] meta, bam
        fasta
        fai

    main:
    ch_versions = Channel.empty()

    ch_bam_to_merge = bam.map {
        meta, bed -> 
            [['id':meta.id.split(/\d+/)[0]], bed ]
        }
        .groupTuple()

    MERGE_BAM(ch_bam_to_merge, fasta, fai)
    
    INDEX_MERGE_BAM(MERGE_BAM.out.bam)

    bam_bai = MERGE_BAM.out.bam
        .join(INDEX_MERGE_BAM.out.bai)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(INDEX_MERGE_BAM.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BAM.out.versions.first())

    emit:
        bam_bai = bam_bai
        versions = ch_versions
}