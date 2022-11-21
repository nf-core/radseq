process UNZIP_FASTQ_LENGTHS {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    errorStrategy {task.exitStatus == 141 ? 'ignore' : 'terminate'}

    input:
    tuple val (meta), path (reads)

    output:
    path ('*_lengths.txt'), emit: lengths
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forward_read = reads[0]
    if ( meta.single_end ) {
        """
        gunzip -c ${forward_read} | head -2 | tail -1 > ${prefix}_lengths.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        END_VERSIONS
        """
    } else {
        def reverse_read = reads[1]
        """
        gunzip -c ${reverse_read} | head -2 | tail -1 > ${prefix}_lengths.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        END_VERSIONS
        """
    }
}