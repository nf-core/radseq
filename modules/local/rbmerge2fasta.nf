process RBMERGE2FASTA {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

    when:
    task.ext.when == null || task.ext.when

    input:
    tuple val (meta), path (rbdiv)
    tuple val (meta), path (rbmerge)

    output:
    tuple val (meta), path ('*_rainbow.fasta'), emit: fasta
    path "versions.yml"                       , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	LENGTH1=\$(cut -f3 ${rbdiv} | awk '(NR==1||length<shortest){shortest=length} END {print shortest}')
	LENGTH=\$((\$LENGTH1 * 11 / 10))

	cat ${rbmerge} <(echo "E") | \\
    sed -e 's/[0-9]*:[0-9]*://g' | \\
    awk -v mlen=\$LENGTH '{
        if (NR == 1) e=\$2;
        else if (\$1 ~/E/ && lenp > len1) {
            c=c+1; print ">dDocent_A_Contig_" e "\\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=\$2;fclus=0;len1=0;freqp=0;lenf=0
        }
        else if (\$1 ~/E/ && lenp <= len1) {
            c=c+1; print ">dDocent_Contig_" e "\\n" seq1; seq1=0; seq2=0;lenp=0;e=\$2;fclus=0;len1=0;freqp=0;lenf=0
        }
        else if (\$1 ~/C/) clus=\$2;
        else if (\$1 ~/L/) len=\$2;
        else if (\$1 ~/S/) seq=\$2;
        else if (\$1 ~/N/) freq=\$2;
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 !~/1/ && len > lenf) {
            seq1 = seq; fclus=clus;lenf=len
        }
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 ~/1/ && \$0 ~/^R 0/ && len <= mlen) {
            seq1 = seq; fclus=clus;lenf=len
        }
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 ~/1/ && \$0 ~!/^R 0/ && len > mlen) {
            seq1 = seq; fclus=clus; len1=len
        }
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 ~/1/ && \$0 ~!/^R 0/ && len <= mlen) {
            seq1 = seq; fclus=clus; lenf=len
        }
        else if (\$1 ~/R/ && \$0 ~!/0/ && freq > freqp && len >= lenp || \$1 ~/R/ && \$0 ~!/0/ && freq == freqp && len > lenp) {
            seq2 = seq; lenp = len; freqp=freq
        }
    }' > ${prefix}_rainbow.fasta
    
	seqtk seq -r ${prefix}_rainbow.fasta > ${prefix}_rainbow.RC.fasta
	mv ${prefix}_rainbow.RC.fasta ${prefix}_rainbow.fasta

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
    seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
END_VERSIONS
    """
}