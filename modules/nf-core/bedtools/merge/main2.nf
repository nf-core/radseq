process BEDTOOLS_MERGE_COV {
    tag "$meta.id"
    label 'process_high_memory'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(cov)
    path(faidx)

    output:
    tuple val(meta), path('*.cov'), emit: cov
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$cov" == "${prefix}.cov") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    cat ${cov} > ${prefix}_tmp.cov
    bedtools sort -i ${prefix}_tmp.cov -faidx ${faidx} | \\
    bedtools \\
        merge \\
        -i - \\
        $args \\
        > ${prefix}.cov
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}