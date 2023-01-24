process RAINBOW_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::rainbow=2.0.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rainbow:2.0.4--hec16e2b_7' :
        'quay.io/biocontainers/rainbow:2.0.4--hec16e2b_7' }"

    input:
    tuple val (meta), path (rbdiv)
    val (type)
    val (save_assembly)

    output:
    tuple val (meta), path ("*_rbmerge.out"), emit: rbmerge
    tuple val (meta), path ('*_rbmerge.log'), emit: log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_assembly = save_assembly ? '-a' : ''
    def args = task.ext.args ?: ''
    """
    rainbow merge -i ${rbdiv} -o ${prefix}_rbmerge.out \\
    ${args} \\
    ${output_assembly} \\
    2> ${prefix}_rbmerge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rainbow: \$(rainbow | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}