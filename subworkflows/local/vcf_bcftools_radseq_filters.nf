include { BCFTOOLS_NORM            } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX     } from '../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_FILTERS_1 } from '../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_FITLERS_2 } from '../../modules/nf-core/tabix/tabix/main.nf'
include { BCFTOOLS_VIEW as FILTER1 } from '../../modules/nf-core/bcftools/view/main'
include { RADSEQ_FILTERS as RADSEQ_FILTERS_1 } from '../../modules/local/bcftools/main'
include { RADSEQ_FILTERS as RADSEQ_FILTERS_2 } from '../../modules/local/bcftools/main'

workflow VCF_BCFTOOLS_RADSEQ_FILTERS {

    take:
    vcf_tbi_fasta // [[id:],vcf,tbi,fasta]

    main:
    ch_versions = Channel.empty()
    
    ch_vcf = BCFTOOLS_NORM ( vcf_tbi_fasta ).vcf
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    TABIX ( ch_vcf )

    ch_vcf_tbi = ch_vcf.join(TABIX.out.tbi)

    ch_miss = Channel.fromList(params.fraction_missingness_list)
    ch_mac = Channel.fromList(params.minor_allele_count_list)

    FILTER1 ( ch_vcf_tbi, [], [], [], ch_miss, ch_mac)

    RADSEQ_FILTERS_1 ( FILTER1.out.vcf , [], 'first')
    
    RADSEQ_FILTERS_2 ( FILTER1.out.vcf, RADSEQ_FILTERS_1.out.txt, 'second')

}