process RAINBOW_DIV {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::rainbow=2.0.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rainbow:2.0.4--hec16e2b_7' :
        'quay.io/biocontainers/rainbow:2.0.4--hec16e2b_7' }"

    input:
    tuple val (meta), path (cluster)

    output:
    tuple val (meta), path ("*_rbdiv.out") , emit: rbdiv
    tuple val (meta), path ('*.log')               , emit: log
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-f 0.5 -K 10'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rainbow div -i ${cluster} -o ${prefix}_rbdiv.out \
    ${args} \
    2> rbdiv.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rainbow: \$(rainbow | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """


}