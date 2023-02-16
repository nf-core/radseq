process CREATE_INTERVALS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::perl-bioperl=1.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(cov), path(intersect), path(low_cov)
    val (lengths)

    output:
    tuple val(meta), path('mapped.*.bed'), emit: intervals
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    if (params.method == 'denovo') {
        """
        cat ${intersect} ${low_cov} > ${prefix}_cov.split.stats
        echo "${lengths.join("\n")}" > ${prefix}_lengths.txt
        MaxLen=\$(awk '{ print length() | "sort -rn" }' ${prefix}_lengths.txt| head -1)
        MaxLen2=\$(( \$MaxLen / 2 ))
        
        TT=\$(( \$MaxLen2 * 1000000 ))
        DP=\$(awk '{print \$4}' ${cov} | sort -rn | perl -e '\$d=.001;@l=<>;print \$l[int(\$d*@l)]')
        CC=\$( awk -v x=\$DP '\$4 < x' ${cov} | awk '{len=\$3-\$2;lc=len*\$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'${task.cpus}'"}')
        
        awk -v x=\$DP '\$4 < x' ${prefix}_cov.split.stats | sort -k1,1 -k2,2 | awk -v cutoff=\$CC -v tt=\$TT 'BEGIN{i=1}
            {len=\$3-\$2;lc=len*\$4;cov = cov + lc
            if (NR == 1 && lc > tt) {x="mapped."i".bed";print \$1"\\t"\$2"\\t"\$3 > x; i=i+1; e=1}
            else if ( cov < cutoff && lc < tt) {x="mapped."i".bed";print \$1"\\t"\$2"\\t"\$3 > x; e=0}
            else if (lc > tt && e > 0 ) {x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=0;i=i+1; e=1}
            else if (lc > tt && e < 1 ) {i=i+1; x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=0;i=i+1;e=1}
            else if (cov > cutoff && lc < tt ) {i=i+1; x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=lc;e=0}
            }'
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
            perl: \$(perl --version | sed -n -E '/^This is/ s/.*\\(v([[:digit:].]+)\\).*/\\1/p')
        END_VERSIONS
        """
    } else {
        """
        cat ${intersect} ${low_cov} > ${prefix}_cov.split.stats
        echo "${lengths.join("\n")}" > ${prefix}_lengths.txt
        MaxLen=\$(awk '{ print length() | "sort -rn" }' ${prefix}_lengths.txt| head -1)
        MaxLen2=\$(( \$MaxLen / 2 ))
        
        TT=\$(( \$MaxLen2 * 1000000 ))
        DP=\$(awk '{print \$4}' ${cov} | sort -rn | perl -e '\$d=.00005;@l=<>;print \$l[int(\$d*@l)]')
        CC=\$( awk -v x=\$DP '\$4 < x' ${cov} | awk '{len=\$3-\$2;lc=len*\$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'${task.cpus}'"}')
        
        awk -v x=\$DP '\$4 < x' ${prefix}_cov.split.stats | sort -k1,1 -k2,2 | awk -v cutoff=\$CC -v tt=\$TT 'BEGIN{i=1}
            {len=\$3-\$2;lc=len*\$4;cov = cov + lc
            if (NR == 1 && lc > tt) {x="mapped."i".bed";print \$1"\\t"\$2"\\t"\$3 > x; i=i+1; e=1}
            else if ( cov < cutoff && lc < tt) {x="mapped."i".bed";print \$1"\\t"\$2"\\t"\$3 > x; e=0}
            else if (lc > tt && e > 0 ) {x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=0;i=i+1; e=1}
            else if (lc > tt && e < 1 ) {i=i+1; x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=0;i=i+1;e=1}
            else if (cov > cutoff && lc < tt ) {i=i+1; x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=lc;e=0}
            }'
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
            perl: \$(perl --version | sed -n -E '/^This is/ s/.*\\(v([[:digit:].]+)\\).*/\\1/p')
        END_VERSIONS   
        """
    }
}