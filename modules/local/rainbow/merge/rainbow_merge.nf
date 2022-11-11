process RAINBOW_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::rainbow=2.0.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rainbow:2.0.4--hec16e2b_7' :
        'quay.io/biocontainers/rainbow:2.0.4--hec16e2b_7' }"

    input:
    tuple val (meta), path (rbdiv)
    val type
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
    def args = task.ext.args ?: '-r 2 -N10000 -R10000 -l 20 -f 0.75'
    // change default arguments based on data
    if (type == 'HYB') {
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
    } else {
        """
        #CLUSTER=(tail -1 ${rbdiv} | cut -f5)
        #CLUSTER1="\$((\$CLUSTER / 100 + 1))"
        #CLUSTER2="\$((\$CLUSTER1 + 100))"

        rainbow merge -i ${rbdiv} -o ${prefix}_rbmerge.out \
        ${args} \
        ${output_assembly} \
        2> ${prefix}_rbmerge.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rainbow: \$(rainbow | head -n 1 | cut -d ' ' -f 2)
        END_VERSIONS
        """

    }
}