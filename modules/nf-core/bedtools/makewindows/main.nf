process BEDTOOLS_MAKEWINDOWS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1' }"

    input:
    tuple val(meta), path(regions)
    val(use_bed)
    val(lengths)
    val(coverage_threshold)

    output:
    tuple val(meta), path("*.tab")           , emit: tab
    tuple val(meta), path("*_cov.low.stats") , emit: low_cov
    tuple val(meta), path("*_cov.high.stats"), emit: high_cov
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_input = use_bed ? "-b $regions" : "-g $regions"
    """
    echo "${lengths.join("\n")}" > ${prefix}_lengths.txt
    MaxLen=\$(awk '{ print length() | "sort -rn" }' ${prefix}_lengths.txt| head -1)

    #split cov.stats file into high and low coverage intervals
	awk '\$4 > ${coverage_threshold}' $regions > ${prefix}_cov.high.stats
	awk '\$4 <= ${coverage_threshold}' $regions > ${prefix}_cov.low.stats
        
    MaxLen2=\$(( \$MaxLen / 2 ))
    ML1=\$(( \$MaxLen2 + 1 ))
    
    bedtools \\
        makewindows \\
        ${arg_input} \\
        -w \$MaxLen2 -s \$ML1 \\
        $args \\
        > ${prefix}.tab
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}