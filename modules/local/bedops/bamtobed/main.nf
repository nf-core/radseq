process BEDOPS_BAMTOBED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bedops=2.4.41" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*_sed.bed'), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convert2bed --input=bam < ${bam} \\
    ${args} \\
    > ${prefix}.bed 
    sed -i 's,\\t*\$,,' ${prefix}.bed > ${prefix}_sed.bed 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedops: \$(bedops --version | sed -n '/version:/p' | cut -d' ' -f 5)
    END_VERSIONS
    """
}