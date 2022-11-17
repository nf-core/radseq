process READ_LENGTH {
    tag
    label

    conda
    container

    input:
    tuple val(meta), path(lengths)

    output:
    path ('lengths.txt'), emit: read_length

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    """
}