process CDHIT {
    tag "meta.id"
    label 'process_high'

    //publishDir '', contentType: 'test/html'

    conda (params.enable_conda ? 'bioconda::cd-hit=4.8.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5b5514e_7' :
    'https://quay.io/repository/biocontainers/cd-hit:4.8.1--h5b5514e_7' }"


    input:
    tuple val (meta), path (fasta) // [[:], forward reads]
    tuple val (meta), path (totaluniqseq)
    val type

    output:
    tuple val (meta), path ("*.clstr") , emit: cdhit_cluster
    tuple val (meta), path ("*.rclstr"), emit: rb_cluster
    tuple val (meta), path ('*.log')   ,   emit: log
    tuple val (meta), path ('*.contig.cluster.totaluniqseq'), emit: clstr_totaluniqseq 
    //path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    def args = task.ext.args ?: '-M 0 -g 1 -d 100'
    // if parameter is provided, return the argument only if it's product minus 0.1 is above 0.8. Replaces dDocent python command w/ groovy
    def SEQUENCE_SIMILARITY = task.ext.simC ? println(Math.max(task.ext.simC.toFloat - 0.1, 0.8)) : '0.8'
    if (type == 'PE'){
        """
        sed -e 's/NNNNNNNNNN/	/g' ${fasta} | cut -f1 > uniq.F.fasta
        
        # Sequence similarityffast
        
        cd-hit-est -i uniq.F.fasta \\
        -o ${prefix} -c ${SEQUENCE_SIMILARITY} -T ${task.cpus} \\
        ${args} \\
        &> ${prefix}.log
        
        awk '{if (\$1 ~ /Cl/) clus = clus + 1; else  print \$3 "\\t" clus}' ${prefix}.clstr | \\
        sed -e 's/[>dDocent_Contig_,...]//g' | \\
        sort -g -k1 -S 2G > ${prefix}.sort.contig.cluster.ids
        
        paste ${prefix}.sort.contig.cluster.ids ${totaluniqseq} > ${prefix}.contig.cluster.totaluniqseq

        # cd-hit TO rainbow cluster format
        sort -k2,2 -g ${prefix}.contig.cluster.totaluniqseq -S 2G | \\
        sed -e 's/NNNNNNNNNN/	/g' > ${prefix}.rclstr
        
        
        """
    } else { 
        // with random end assembly versions totaluniqseq is based only on the forward reads
        // uniq reduces down the dataset
        // Gets rid of any reduncy
        // cluster ids might match more than 1 reads
        // Random Shearing 
        // Purple sea urchin is random sheared
        """
        sed -e 's/NNNNNNNNNN/	/g' totaluniqseq | cut -f1 | sort -S 2G | \\
        uniq | \\
        awk '{c= c + 1; print ">dDocent_Contig_" c "\\n" \$1}' > ${prefix}.uniq.F.fasta
				
        cd-hit-est -i ${prefix}.uniq.F.fasta \\
        -o ${prefix} -c ${SEQUENCE_SIMILARITY} -T ${task.cpus} \\
        ${args} \\
        &> ${prefix}.log
		
        awk '{if (\$1 ~ /Cl/) clus = clus + 1; else  print \$3 "\\t" clus}' ${prefix}.clstr | \\
        sed -e 's/[>dDocent_Contig_,...]//g' | \\
        sort -g -k1 -S 2G > sort.contig.cluster.ids
		
        paste sort.contig.cluster.ids <(awk '!/>/' uniq.F.fasta) > ${prefix}.contig.cluster.Funiq
		
        sed -e 's/NNNNNNNNNN/	/g' totaluniqseq | \\
        sort -k1 -S 2G | \\
        awk '{print \$0 "\\t" NR}'  > ${prefix}.totaluniqseq.CN
		
        join -t \$'\\t' -1 3 -2 1 ${prefix}.contig.cluster.Funiq ${prefix}.totaluniqseq.CN -o 2.3,1.2,2.1,2.2 > ${prefix}.contig.cluster.totaluniqseq
        
        # cd-hit TO rainbow cluster format
        sort -k2,2 -g ${prefix}.contig.cluster.totaluniqseq -S 2G | \\
        sed -e 's/NNNNNNNNNN/	/g' > ${prefix}.rcluster
        """
    }
}