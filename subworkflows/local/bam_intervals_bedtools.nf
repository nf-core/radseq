include { BEDTOOLS_BAMTOBED } from '../../modules/nf-core/bedtools/bamtobed/main.nf'
include { BEDTOOLS_MERGE    } from '../../modules/nf-core/bedtools/merge/main.nf'

workflow BAM_INTERVALS_BEDTOOLS {

    take:
    bam
    fai // TODO: add [meta, fasta, fai]

    main:
    ch_versions = Channel.empty()

    bed = BEDTOOLS_BAMTOBED (bam).bed
    ch_versions = ch_versions.mix (BEDTOOLS_BAMTOBED.out.versions)
    
    bed.collect{meta, bed -> 
        def metaf = [:] 
        metaf.id = (meta.id =~ ).findAll().first()
        [metaf, bed]}.view()
    
    /*mbed = BEDTOOLS_MERGE (bed.collect{meta, bed -> 
        def metaf = [:] 
        metaf.id = (meta.id =~ /(?i)\b([a-z])[a-z]*\1\b/).findAll()*.first()) 
        [, bed]})
    ch_versions = ch_versions.mix (BEDTOOLS_MERGE.out.versions)*/

    //cov = COVERAGEBED ().cov

    //COMBINE_FILTER_BED ()

    //emit:
    //mbam_bai

}