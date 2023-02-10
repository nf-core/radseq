process COMBINE_UNIQUE_READS {
    tag "${meta.id}"
    label 'process_medium'

    // get a can't find conda dir. ? Check you have anaconda3 installed
    conda (params.enable_conda ? 'bioconda::perl-sys-info-driver-linux=0.7905' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-sys-info-driver-linux:0.7905--pl5321hdfd78af_1' :
        'quay.io/upennlibraries/perl_apache' }"

    input:
    tuple val (meta), path (reads) // loading all individual uniq sequence per collected
    val (type) // sequencing technology used. Changes how unique sequences are identified
    each withinIndv_MinDepth // within_individual
    each acrossIndv_MinDepth // number of unique individuals w/ reads

    output:
    tuple val (meta), path ('*_uniq.full.fasta'), emit: uniq_reads
    tuple val (meta), path ('totaluniqseq')     , emit: totaluniqseq
    path 'versions.yml'                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
        END_VERSIONS
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
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
            perl: \$(perl --version | sed -n -E '/^This is/ s/.*\\(v([[:digit:].]+)\\).*/\\1/p')
        END_VERSIONS
        """
    }
}