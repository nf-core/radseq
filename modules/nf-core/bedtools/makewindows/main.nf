process BEDTOOLS_MAKEWINDOWS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1' }"

    input:
    tuple val(meta), path(regions)
    val(lengths)
    val(coverage_threshold)
    val(winsize)

    output:
    tuple val(meta), path("*.tab")           , optional: true, emit: tab
    tuple val(meta), path("*_cov.low.stats") , optional: true, emit: low_cov
    tuple val(meta), path("*_cov.high.stats"), optional: true, emit: high_cov
    tuple val(meta), path ("*.bed")          , optional: true, emit: bed
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_input = regions.extension in ["bed", "tab"] ? "-b ${regions}" : "-g ${regions}"
    if ("${regions}" == "${prefix}.bed") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if (regions.extension in ["bed", "tab", "cov"]) {
        """
        echo "${lengths.join("\n")}" > "${prefix}_lengths.txt"
        MaxLen=\$(awk '{ print length() | "sort -rn" }' "${prefix}_lengths.txt" | head -1)
        #split cov.stats file into high and low coverage intervals
        awk '\$4 > ${coverage_threshold}' ${regions} > ${prefix}_cov.high.stats
        awk '\$4 <= ${coverage_threshold}' ${regions} > ${prefix}_cov.low.stats
        MaxLen2=\$(("\$MaxLen" / 2))
        ML1=\$(("\$MaxLen2" + 1))
        bedtools \\
            makewindows \\
            -b ${prefix}_cov.high.stats \\
            -w \$MaxLen2 -s \$ML1 \\
            $args \\
            > ${prefix}.tab
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
        """
    } else {
        """
        slide=\$(expr ${winsize} - 100)
        bedtools \\
            makewindows \\
            ${arg_input} \\
            ${args} \\
            -w ${winsize} -s \$slide \\
            > ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
        """

    }
}