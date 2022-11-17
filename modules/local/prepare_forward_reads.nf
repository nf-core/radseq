/*
    # Take DNA characters out of fastq
    # Simplify data 
    # PE
    # RPE - use special unique function
    # SE
    # OL  
    # Using really long read technology
        # Short reads 
        # Long Reads
        # MiSEQ
    # HYB - Can leave out

    # Total up the data
*/

process PREPARE_FORWARD_READS {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::perl-ipc-sharelite=0.17' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b6fc09bed47d0dc4d8384ce9e04af5806f2cc91b:305092c6f8420acd17377d2cc8b96e1c3ccb7d26-0' :
        'quay.io/biocontainers/mulled-v2-b6fc09bed47d0dc4d8384ce9e04af5806f2cc91b:305092c6f8420acd17377d2cc8b96e1c3ccb7d26-0' }"

    when:
    task.ext.when == null || task.ext.when

    input:
    tuple val (meta), path (reads)
    val (type) // sequencing technology used. Changes how unique sequences are identified

    output:
    tuple val (meta), path ('*.uniq.seqs'), emit: indv_uniq_seqs


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (type == 'PE'){
        def forward_reads = "${reads[0]}"
        def reverse_reads = "${reads[1]}"
        """
        gunzip -c ${forward_reads} | \\
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | \\
        awk '!/>/' > ${prefix}.forward
		
        gunzip -c ${reverse_reads} | \\
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | \\
        awk '!/>/' > ${prefix}.reverse
        
        paste -d '-' ${prefix}.forward ${prefix}.reverse | \\
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | \\
        sed -e 's/-/NNNNNNNNNN/' | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\t\$k\n";}' > ${prefix}.uniq.seqs
        """
    } else if (type == 'RPE'){
        def forward_reads = "${reads[0]}"
        def reverse_reads = "${reads[1]}"
        """
        gunzip -c ${forward_reads} | \\
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | \\
        awk '!/>/' > ${prefix}.forward
		
        gunzip -c ${reverse_reads} | \\
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | \\
        awk '!/>/' > ${prefix}.reverse
        
        paste ${prefix}.forward ${prefix}.reverse | sort -k1 -S 200M > ${prefix}.fr
		cut -f1 ${prefix}.fr | uniq -c > ${prefix}.f.uniq && cut -f2 ${prefix}.fr > ${prefix}.r
		awk '{for(i=0;i<\$1;i++)print}' ${prefix}.f.uniq > ${prefix}.f.uniq.e
		
        paste -d '-' ${prefix}.f.uniq.e ${prefix}.r | \\
        awk '!/NNN/'| \\
        sed -e 's/-/NNNNNNNNNN/' | \\
        sed -e 's/^[ \\t]*//' | \\
        sed -e 's/\\s/\\t/g' > ${prefix}.uniq.seqs
		
        rm ${prefix}.f.uniq.e ${prefix}.f.uniq ${prefix}.r ${prefix}.fr
        """
    } else if (type == 'SE'){
        def forward_reads = "${reads[0]}"
        """
        gunzip -c ${forward_reads} | \\
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | \\
        awk '!/>/' | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' > ${prefix}.uniq.seqs
        """
    } else if (type == 'OL')  {
        """
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' ${assembled_fastq} | \\
        awk '!/>/' | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' > ${prefix}.uniq.seqs
        """
    } else {
        error "invalid sequence type specified or is not supported: ${type}"
    }
}