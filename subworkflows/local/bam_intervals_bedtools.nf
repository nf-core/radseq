include { BEDTOOLS_BAMTOBED  } from '../../modules/nf-core/bedtools/bamtobed/main.nf'
include { BEDTOOLS_MERGE_BED } from '../../modules/nf-core/bedtools/merge/main.nf'
include { BEDTOOLS_COVERAGE  } from '../../modules/nf-core/bedtools/coverage/main.nf'
include { BEDTOOLS_MERGE_COV } from '../../modules/nf-core/bedtools/merge/main2.nf'
include { CREATE_INTERVALS     } from '../../modules/local/create_intervals.nf'
include { BEDTOOLS_MAKEWINDOWS } from '../../modules/nf-core/bedtools/makewindows/main.nf'
include { BEDTOOLS_INTERSECT   } from '../../modules/nf-core/bedtools/intersect/main.nf'

workflow BAM_INTERVALS_BEDTOOLS {

    take:
    bam
    faidx // TODO: add [meta, fasta, fai]
    read_lengths
    coverage_threshold

    main:
    ch_versions = Channel.empty()

    bed = BEDTOOLS_BAMTOBED (bam, faidx.first()).bed
    ch_versions = ch_versions.mix (BEDTOOLS_BAMTOBED.out.versions)

    // split based on number and return the first element and group bed files based on shared id string
    ch_bed_to_merge = bed.map {
        meta, bed -> 
            [['id':meta.id.split(/\d+/)[0]], bed ]
        }
        .groupTuple()

    mbed = BEDTOOLS_MERGE_BED (ch_bed_to_merge, faidx.first()).bed
    ch_versions = ch_versions.mix (BEDTOOLS_MERGE_BED.out.versions)

    // switch the order 
    cov = BEDTOOLS_COVERAGE (bed.combine(mbed.map{it[1]}).map{meta,bed,mbed -> [meta,mbed,bed]}, faidx.map{it[1]}.first()).cov
    ch_versions = ch_versions.mix (BEDTOOLS_COVERAGE.out.versions)

    ch_cov_to_merge = cov.map {
        meta, bed -> 
            [['id':meta.id.split(/\d+/)[0]], bed ]
        }
        .groupTuple()

    mcov = BEDTOOLS_MERGE_COV (ch_cov_to_merge).cov
    ch_versions = ch_versions.mix (BEDTOOLS_MERGE_COV.out.versions)

    ch_tab = BEDTOOLS_MAKEWINDOWS (mcov, true, read_lengths, params.splitByReadCoverage).tab
    ch_versions = ch_versions.mix (BEDTOOLS_MAKEWINDOWS.out.versions)

    ch_bedtointersect = mcov.join(ch_tab)

    ch_intersect = BEDTOOLS_INTERSECT (ch_bedtointersect, 'bed').intersect
    ch_versions = ch_versions.mix (BEDTOOLS_INTERSECT.out.versions)

    ch_intervals = CREATE_INTERVALS (ch_intersect, BEDTOOLS_MAKEWINDOWS.out.low_cov, read_lengths).intervals.transpose()
    //TODO: add channle versions

    emit:
    intervals = ch_intervals


    versions = ch_versions

}