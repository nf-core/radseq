process RADSEQ_FILTERS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf)
    path(samples2remove)
    val round

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}") , emit: vcf
    path("*_RADseqPopFilters_indvFAIL.txt")            , emit: txt, optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ab       = task.ext.args_ab ? "bcftools filter ${task.ext.args_ab} | \\": '\\'
    def mq       = task.ext.args_mq ? "bcftools filter ${task.ext.args_mq} | \\": '\\'
    def pr       = task.ext.args_pr ? "bcftools filter ${task.ext.args_pr} | \\": '\\'
    def st       = task.ext.args_st ? "bcftools filter ${task.ext.args_st} | \\": '\\'
    def fmissing = task.ext.args_fmiss ?: ''
    def mac      = task.ext.args_mac ?: ''
    def mindp    = task.ext.args_mindp ?: ''
    def maxdp    = task.ext.args_maxdp ?: "\$MAX_DP"
    def args_indvfmiss = task.ext.args_indvfmiss ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def extension     = "vcf.gz"
    
    if (round == 'first') {
        """
        MAX_DP=\$(bcftools query -f '%DP\\n' ${vcf} | awk '{ sum += \$1; n++ } END { if (n > 0) print sum / n; }')
        cat $vcf | \
        $ab
        $mq
        $pr
        $st
        bcftools filter -i"F_MISSING<${fmissing} && MAC>${mac} && INFO/DP<${maxdp} && FORMAT/DP>${mindp}" -Oz -o ${prefix}_RADseqPopFilters.${extension}

        paste <(bcftools query -f '[%SAMPLE\\t]\\n' ${prefix}_RADseqPopFilters.${extension} | head -1 | tr '\\t' '\\n') \
        <(bcftools query -f '[%GT\\t]\\n' ${prefix}_RADseqPopFilters.${extension} | \
        awk -v OFS="\\t" '{for (i=1;i<=NF;i++) if (\$i == ".") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | \
        sort -k1,1n | \
        cut -f 2 ) | \
        awk -v prop="${args_indvfmiss}" '{ if (\$2 > prop) print \$1 }' > ${prefix}_RADseqPopFilters_indvFAIL.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        MAX_DP=\$(bcftools query -f '%DP\\n' ${vcf} | awk '{ sum += \$1; n++ } END { if (n > 0) print sum / n; }')        
        cat $vcf |\\
         $ab
         $mq
         $pr
         $st
         bcftools view -S ^${samples2remove} -i"F_MISSING<${fmissing} && MAC>${mac} && INFO/DP<${maxdp} && FORMAT/DP>${mindp}" -Oz -o ${prefix}_rmvindv_RADseqPopFilters.${extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """        
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}