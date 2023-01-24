process CDHIT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::cd-hit=4.8.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--hdbcaa40_0' :
    'quay.io/biocontainers/cd-hit:4.8.1--hdbcaa40_0' }"

    input:
    tuple val (meta), path (fasta) // [[:], forward reads]
    tuple val (meta2), path (totaluniqseq)
    val type

    output:
    tuple val (meta), path ("*.clstr")       , emit: cdhit_cluster
    tuple val (meta), path ("uniq.F.fasta")  , emit: forward_uniq
    path 'versions.yml'                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    def args = task.ext.args ?: ''
    if (type == 'PE') {
        """
        sed -e 's/NNNNNNNNNN/	/g' ${fasta} | cut -f1 > uniq.F.fasta
                
        cd-hit-est -i uniq.F.fasta \\
        -o ${prefix} \\
        -M ${task.memory.mega} \\
        -T ${task.cpus} \\
        ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cdhit: \$(cd-hit -h | head -n 1 | sed 's/^.*====== CD-HIT version //;s/ (built on .*) ======//' )
        END_VERSIONS
        """
    } else { 
        // with random end assembly versions totaluniqseq is based only on the forward reads
        // uniq reduces down the dataset
        // Gets rid of any reduncy
        // cluster ids might match more than 1 reads
        // Random Shearing
        """
        sed -e 's/NNNNNNNNNN/   /g' ${totaluniqseq} | cut -f1 | sort | \\
        uniq | \\
        awk '{c= c + 1; print ">dDocent_Contig_" c "\\n" \$1}' > uniq.F.fasta
				
        cd-hit-est -i uniq.F.fasta \\
        -o ${prefix} \\
        -M ${task.memory.mega} \\
        -T ${task.cpus} \\
        ${args} 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cdhit: \$(cd-hit -h | head -n 1 | sed 's/^.*====== CD-HIT version //;s/ (built on .*) ======//' )
        END_VERSIONS
        """
    }
}