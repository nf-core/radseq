process FREEBAYES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py38ha193a2f_3' :
        'quay.io/biocontainers/freebayes:1.3.5--py38ha193a2f_3' }"

    input:
    tuple val(meta), path(bam), path(bai), path(target_bed)
    val targets_file
    path fasta
    path fasta_fai
    path samples
    path populations
    path cnv

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def targets_file     = target_bed     ? "--target ${target_bed}"       : ""
    def samples_file     = samples        ? "--samples ${samples}"         : ""
    def populations_file = populations    ? "--populations ${populations}" : ""
    def cnv_file         = cnv            ? "--cnv-map ${cnv}"             : ""
    def interval         = meta.interval  ? '_' + meta.interval            : ''
    """
    freebayes \\
        -f $fasta \\
        $targets_file \\
        $samples_file \\
        $populations_file \\
        $cnv_file \\
        $args \\
        $bam | \\
    bgzip -c > ${prefix}${interval}.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}