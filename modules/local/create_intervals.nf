process CREATE_INTERVALS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::perl=5.26.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(cov)
    tuple val(meta2), path(low_cov)
    val (lengths)

    output:
    tuple val(meta), path('mapped.*.bed'), emit: intervals

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    if (params.method == 'denovo') {
        """
        cat ${cov} ${low_cov} > ${prefix}_cov.split.stats
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
            BusyBox: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
            perl: \$()
        END_VERSIONS
        """
    } else {
        """
        cat ${cov} ${low_cov} > ${prefix}_cov.split.stats
        echo "${lengths.join("\n")}" > ${prefix}_lengths.txt
        MaxLen=\$(awk '{ print length() | "sort -rn" }' ${prefix}_lengths.txt| head -1)
        MaxLen2=\$(( \$MaxLen / 2 ))
        
        TT=\$(( \$MaxLen2 * 1000000 ))
        DP=\$(awk '{print \$4}' cov.stats | sort -rn | perl -e '\$d=.00005;@l=<>;print \$l[int(\$d*@l)]')
        CC=\$( awk -v x=\$DP '\$4 < x' cov.stats | awk '{len=\$3-\$2;lc=len*\$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'${task.cpus}'"}')
        
        awk -v x=\$DP '\$4 < x' ${prefix}_cov.split.stats | sort -k1,1 -k2,2 | awk -v cutoff=\$CC -v tt=\$TT 'BEGIN{i=1}
            {len=\$3-\$2;lc=len*\$4;cov = cov + lc
            if (NR == 1 && lc > tt) {x="mapped."i".bed";print \$1"\\t"\$2"\\t"\$3 > x; i=i+1; e=1}
            else if ( cov < cutoff && lc < tt) {x="mapped."i".bed";print \$1"\\t"\$2"\\t"\$3 > x; e=0}
            else if (lc > tt && e > 0 ) {x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=0;i=i+1; e=1}
            else if (lc > tt && e < 1 ) {i=i+1; x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=0;i=i+1;e=1}
            else if (cov > cutoff && lc < tt ) {i=i+1; x="mapped."i".bed"; print \$1"\\t"\$2"\\t"\$3 > x; cov=lc;e=0}
            }'
            
        """
    }
}