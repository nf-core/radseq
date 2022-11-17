include { SAMTOOLS_MERGE } from '../../modules/samtools/merge/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/samtools/index/main.nf'



workflow MERGE_INTERVALS_SAMTOOLS {

    take:
    bam_bai
    fasta // TODO: add [meta, fasta, fai]

    main:
    
    mbam = SAMTOOLS_MERGE (bam_bai).bam

    mbam_bai = SAMTOOLS_INDEX (mbam).bam_bai

    fai = SAMTOOLS_FAIDX (fasta).fai

    mbam_bai_bed = BAMTOBED (bam_bai)

    cov = COVERAGEBED ().cov

    COMBINE_FILTER_BED ()

    emit:
    mbam_bai

}