process READ_LENGTH {
    //tag
    //label

    conda
    container

    input:
    path (reads)

    output:
    env(MLEN)
    tuple env(INSERT), env(INSERTH), env(INSERTL), env(SD)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: ${meta.id}
    def args = task.ext.args ?: ''
    if (meta.single_end && params.method == 'denovo') {
        """
        gunzip -c  | head -2 | tail -1 >> lengths.txt
        
        
        cat *_lengths.txt > lengths.txt

        MLEN=\$(awk '{ print length() | "sort -rn" }' lengths.txt | head -1)
        MLEN=\$(($MLEN / 2))
        
        INSERT=$(($MLEN * 2 ))
        INSERTH=$(($INSERT + 100 ))
        INSERTL=$(($INSERT - 100 ))
        SD=$(($INSERT / 5))
        """
    } 
}