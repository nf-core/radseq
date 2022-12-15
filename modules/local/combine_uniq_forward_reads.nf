process COMBINE_UNIQUE_READS {
    tag "${meta.id}"
    label 'process_medium'

    // get a can't find conda dir. ? Check you have anaconda3 installed
    conda (params.enable_conda ? 'seqtk bioconda::fastp=0.23.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b6fc09bed47d0dc4d8384ce9e04af5806f2cc91b:305092c6f8420acd17377d2cc8b96e1c3ccb7d26-0' :
        'quay.io/biocontainers/mulled-v2-b6fc09bed47d0dc4d8384ce9e04af5806f2cc91b:305092c6f8420acd17377d2cc8b96e1c3ccb7d26-0' }"

    when:
    task.ext.when == null || task.ext.when

    input:
    tuple val (meta), path (reads) // loading all individual uniq sequence per collected
    val (type) // sequencing technology used. Changes how unique sequences are identified
    each withinIndv_MinDepth // within_individual
    each acrossIndv_MinDepth // number of unique individuals w/ reads

    output:
    tuple val (meta), path ('*_uniq.full.fasta'), emit: uniq_reads
    //tuple val (meta), path ('totaluniqseq'), emit: totaluniqseq

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (type == 'RPE' || type == 'ROL') {
        """
        awk -v x="${withinIndv_MinDepth}" '(\$1 >= x)' *.uniq.seqs | \\
        cut -f2 | \\
        sed -e 's/NNNNNNNNNN/-/' >  total.uniqs
        
        cut -f 1 -d "-" total.uniqs > total.u.F
        cut -f 2 -d "-" total.uniqs > total.u.R
        
        paste total.u.F total.u.R | \\
        sort -k1 -S 2G > total.fr

        awk -v x=${withinIndv_MinDepth} '\$1 >= x' *.uniq.seqs | \\
        cut -f2 | \\
        sed -e 's/NNNNNNNNNN/	/g' | \\
        cut -f1 | \\
        uniq | \\
        sort -S 2G | \\
        uniq -c > total.f.uniq

        join -1 2 -2 1 -o 1.1,1.2,2.2 total.f.uniq total.fr | \\
        awk '{print \$1 "\t" \$2 "NNNNNNNNNN" \$3}' | \\
        awk -v x=${acrossIndv_MinDepth} '\$1 >= x' > uniq.k.${withinIndv_MinDepth}.c.${acrossIndv_MinDepth}.seqs
        
        sort -k1 -r -n -S 2G uniq.k.${withinIndv_MinDepth}.c.${acrossIndv_MinDepth}.seqs | \\
        cut -f2 > totaluniqseq
        awk '{c= c + 1; print ">dDocent_Contig_" c "\\n" \$1}' totaluniqseq > ${prefix}_uniq.full.fasta
        """
    } else {
        """
        awk -v x=${withinIndv_MinDepth} '(\$1 >= x)' *.uniq.seqs | \\
        cut -f2 | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' | \\
        awk -v x=${acrossIndv_MinDepth} '(\$1 >= x)' > uniq.k.${withinIndv_MinDepth}.c.${acrossIndv_MinDepth}.seqs
        
        # order the sequences for reproducibility 
        sort -k1 -r -n -S 2G uniq.k.${withinIndv_MinDepth}.c.${acrossIndv_MinDepth}.seqs | \\
        cut -f2 > totaluniqseq 
        
        awk '{c= c + 1; print ">dDocent_Contig_" c "\\n" \$1}' totaluniqseq > ${prefix}_uniq.full.fasta
        """
    }
}