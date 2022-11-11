process CDHIT_TO_RBCLSTR {
    tag "meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::coreutils=8.25" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:8.25--0' :
        'quay.io/biocontainers/coreutils:8.25--0' }"
    
    input:
    tuple val(meta) , path(clstr)
    tuple val(meta2), path(totaluniqseq)
    val (type)

    output:
    tuple val(meta), path('*.rclstr'), emit: rbcluster
    tuple val (meta), path ('*.contig.cluster.totaluniqseq'), emit: clstr_totaluniqseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (type == 'PE'){
        """
        awk '{if (\$1 ~ /Cl/) clus = clus + 1; else  print \$3 "\\t" clus}' ${clstr} | \\
        sed -e 's/[>dDocent_Contig_,...]//g' | \\
        sort -g -k1 -S 2G > ${prefix}.sort.contig.cluster.ids
            
        paste ${prefix}.sort.contig.cluster.ids ${totaluniqseq} > ${prefix}.contig.cluster.totaluniqseq

        # cd-hit TO rainbow cluster format
        sort -k2,2 -g ${prefix}.contig.cluster.totaluniqseq -S 2G | \\
        sed -e 's/NNNNNNNNNN/	/g' > ${prefix}.rclstr    
        """
    } else {
        """
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