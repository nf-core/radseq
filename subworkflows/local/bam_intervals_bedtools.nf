include { BEDTOOLS_BAMTOBED             } from '../../modules/nf-core/bedtools/bamtobed/main.nf'
include { BEDOPS_MERGE_BED              } from '../../modules/local/bedops/merge/main.nf'
include { BEDTOOLS_SORT                 } from '../../modules/nf-core/bedtools/sort/main.nf'
include { BEDTOOLS_COVERAGE             } from '../../modules/nf-core/bedtools/coverage/main.nf'
include { BEDTOOLS_MERGE_COV            } from '../../modules/local/bedtools/merge/main.nf'
include { CREATE_INTERVALS              } from '../../modules/local/create_intervals.nf'
include { BEDTOOLS_MAKEWINDOWS          } from '../../modules/local/bedtools/makewindows/main.nf'
include { BEDTOOLS_INTERSECT            } from '../../modules/nf-core/bedtools/intersect/main.nf'

workflow BAM_INTERVALS_BEDTOOLS {

    take:
    bam
    faidx
    read_lengths
    coverage_threshold

    main:
    ch_versions = Channel.empty()

    // Reduce the number of bam files to be passed in the subworkflow via parameters
        // Purpose: Memory constraints for large sample sizes at BEDTOOLS_MERGE_COV
        // .randomSample(# to subset to, random seed) random seed is critical to -resume functionality
    ch_bam = params.subset_intervals_channel ? bam.randomSample(params.subset_intervals_channel, 234) : bam

    ch_bed = BEDTOOLS_BAMTOBED (ch_bam).bed
    ch_versions = ch_versions.mix (BEDTOOLS_BAMTOBED.out.versions)

    ch_bed_to_merge = ch_bed.map {
        meta, bed -> 
            [['id':meta.id.split(/\d+/)[0]], bed ] // split based on number and return the first element and group bed files based on shared id string
        }
        .groupTuple()

    ch_mbed = BEDOPS_MERGE_BED (ch_bed_to_merge).bed
    ch_versions = ch_versions.mix(BEDOPS_MERGE_BED.out.versions)

    ch_sorted_mbed = BEDTOOLS_SORT (ch_mbed, faidx.first()).sorted
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    // Calculate read coverage across indv. samples 
    cov = BEDTOOLS_COVERAGE (ch_bed.combine(ch_sorted_mbed.map{it[1]}).map{meta,bed,mbed -> [meta,mbed,bed]}, faidx.first()).bed
    ch_versions = ch_versions.mix (BEDTOOLS_COVERAGE.out.versions)

    ch_cov_to_merge = cov.map {
        meta, bed -> 
            [['id':meta.id.split(/\d+/)[0]], bed ]
        }
        .groupTuple()

    // combines overlapping features into a single report
    ch_mcov = BEDTOOLS_MERGE_COV (ch_cov_to_merge, faidx.first()).cov
    ch_versions = ch_versions.mix (BEDTOOLS_MERGE_COV.out.versions)

    // split into 2 files: high and low then make intervals across the genome based
    ch_split_high_coverage = BEDTOOLS_MAKEWINDOWS (ch_mcov, true, read_lengths, coverage_threshold).tab
    ch_versions = ch_versions.mix (BEDTOOLS_MAKEWINDOWS.out.versions)

    // Writes overlapping regions into new bed file
    ch_intersect = BEDTOOLS_INTERSECT (ch_mcov.join(ch_split_high_coverage), 'bed').intersect
    ch_versions = ch_versions.mix (BEDTOOLS_INTERSECT.out.versions)

    ch_createintervals = ch_mcov.join(ch_intersect).join(BEDTOOLS_MAKEWINDOWS.out.low_cov)

    //TODO #2: Convert into Groovy function in nf-core-radseq/lib/WorkflowRadseq.groovy
    ch_intervals = CREATE_INTERVALS (ch_createintervals, read_lengths).intervals.transpose()
    ch_versions = ch_versions.mix (CREATE_INTERVALS.out.versions)
    
    emit:
    intervals = ch_intervals

    versions = ch_versions

}