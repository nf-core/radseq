process READ_LENGTH {

    //conda
    //container

    input:
    path (reads)

    output:
    env(MLEN), emit: mlen
    tuple env(INSERT), env(INSERTH), env(INSERTL), env(SD), emit: bwa_mem_params

    when:
    task.ext.when == null || task.ext.when

    script:
    //def forward_reads = reads.filter(x -> x.contains('.1.' || x.contains('.F.')))
    //def reverse_reads = meta.single_end ? reads.filter(x -> x.contains('.2.' || x.contains('.R.'))) : ''
    if (meta.single_end && params.method == 'denovo') {
        """
        gunzip -c  | head -2 | tail -1 >> lengths.txt
        
        
        #cat *_lengths.txt > lengths.txt

        MLEN=\$(awk '{ print length() | "sort -rn" }' lengths.txt | head -1)
        MLEN2=\$((\$MLEN / 2))
        MLEN1=\$((\$MLEN2 + 1))

        INSERT=$((\$MLEN * 2 ))
        INSERTH=$((\$INSERT + 100 ))
        INSERTL=$((\$INSERT - 100 ))
        SD=$((\$INSERT / 5))


        """
    } else {
        """

        """
    }
}