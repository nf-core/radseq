process CDHIT {
    tag "meta.id"
    label 'process_high'

    //publishDir '', contentType: 'test/html'

    conda (params.enable_conda ? 'bioconda::cd-hit=4.8.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--hdbcaa40_0' :
    'quay.io/biocontainers/cd-hit:4.8.1--hdbcaa40_0' }"


    input:
    tuple val (meta), path (fasta) // [[:], forward reads]
    tuple val (meta2), path (totaluniqseq)
    val type

    output:
    tuple val (meta), path ("*.clstr") , emit: cdhit_cluster
    tuple val (meta), path ('*_cdhit.log')   ,   emit: log
    //path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    def args = task.ext.args ?: '-M 0 -g 1 -d 100'
    // if parameter is provided, return the argument only if it's product minus 0.1 is above 0.8. Replaces dDocent python command w/ groovy
    def SEQUENCE_SIMILARITY = task.ext.simC ? println(Math.max(task.ext.simC.toFloat - 0.1, 0.8)) : '0.8'
    if (type == 'PE') {
        """
        sed -e 's/NNNNNNNNNN/	/g' ${fasta} | cut -f1 > uniq.F.fasta
                
        cd-hit-est -i uniq.F.fasta \\
        -o ${prefix} -c ${SEQUENCE_SIMILARITY} -T ${task.cpus} \\
        ${args} \\
        &> ${prefix}_cdhit.log
        """
    } else { 
        // with random end assembly versions totaluniqseq is based only on the forward reads
        // uniq reduces down the dataset
        // Gets rid of any reduncy
        // cluster ids might match more than 1 reads
        // Random Shearing 
        // Purple sea urchin is random sheared
        """
        sed -e 's/NNNNNNNNNN/	/g' ${totaluniqseq} | cut -f1 | sort -S 2G | \\
        uniq | \\
        awk '{c= c + 1; print ">dDocent_Contig_" c "\\n" \$1}' > ${prefix}.uniq.F.fasta
				
        cd-hit-est -i ${prefix}.uniq.F.fasta \\
        -o ${prefix} -c ${SEQUENCE_SIMILARITY} -T ${task.cpus} \\
        ${args} \\
        &> ${prefix}_cdhit.log
        """
    }
}