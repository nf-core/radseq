process READ_LENGTH_STATS {
    //tag
    //label

    conda (params.enable_conda ? "anaconda::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    val (lengths)

    output:
    //env(MLEN)
    //tuple env(MLEN2), env(ML1), emit: split_highcov_intervals
    //tuple env(INSERT), env(SD), env(INSERTH), env(INSERTL), emit: bwa_mem_denovo_param
    stdout emit: std
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """    
    echo "${lengths.join("\n")}" > lengths.txt
    
    # Interval Creation
    MLEN=\$(awk '{ print length() | "sort -rn" }' lengths.txt | head -1)
    MLEN2=\$(( \$MLEN / 2 ))
    ML1=\$(( \$MLEN2 + 1 ))

    # BWA MEM INSERT SIZE
    INSERT=\$(( \$MLEN * 2 ))
    INSERTH=\$(( \$INSERT + 100 ))
    INSERTL=\$(( \$INSERT - 100 ))
    SD=\$(( \$INSERT / 5 ))

    echo \$MLEN \$MLEN2 \$ML1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}