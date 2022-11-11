/*
10 N's is a spacer and can be converted 

Combine uniq reads into two files
if [[ "$ATYPE" == "RPE" || "$ATYPE" == "ROL" ]]; then
  Could be optimized!!!! Index the list/array and join based on the index
  
  parallel --no-notice -j $NUMProc awk -v x=$withinIndv_MinDepth \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | sed -e 's/NNNNNNNNNN/-/' >  total.uniqs
  cut -f 1 -d "-" total.uniqs > total.u.F
  cut -f 2 -d "-" total.uniqs > total.u.R
  paste total.u.F total.u.R | $sort -k1 --parallel=$NUMProc -S 2G > total.fr
 
  parallel --no-notice --env special_uniq special_uniq $withinIndv_MinDepth {} ::: *.uniq.seqs  | $sort --parallel=$NUMProc -S 2G | uniq -c > total.f.uniq
  join -1 2 -2 1 -o 1.1,1.2,2.2 total.f.uniq total.fr | awk '{print $1 "\t" $2 "NNNNNNNNNN" $3}' | awk -v x=$acrossIndv_MinDepth '$1 >= x' > uniq.k.$withinIndv_MinDepth.c.$acrossIndv_MinDepth.seqs
  rm total.uniqs total.u.* total.fr total.f.uniq* 
  
else
    parallel --no-notice awk -v x=$withinIndv_MinDepth \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | awk -v x=$acrossIndv_MinDepth '$1 >= x' > uniq.k.$withinIndv_MinDepth.c.$acrossIndv_MinDepth.seqs
fi
*/
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
    tuple val (meta), path (reads)
    val (type) // sequencing technology used. Changes how unique sequences are identifiedF
    each withinIndv_MinDepth // within_individual
    each acrossIndv_MinDepth // number of unique individuals w/ reads

    output:
    tuple val (meta), path ('uniq.full.fasta'), emit: uniq_reads
    //tuple val (meta), path ('totaluniqseq'), emit: totaluniqseq

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    //def withinIndv_MinDepth = task.ext.withinIndv_MinDepth ?: ''
    /*if (type == 'RPE' || type == 'ROL') {
        """
        awk -v x="${withinIndv_MinDepth}" '(\$1 >= x)' *.uniq.seqs | cut -f2 | sed -e 's/NNNNNNNNNN/-/' >  total.uniqs
        cut -f 1 -d "-" total.uniqs > total.u.F
        cut -f 2 -d "-" total.uniqs > total.u.R
        paste total.u.F total.u.R | sort -k1 --parallel=${task.cpus} -S 2G > total.fr
 
        awk -v x=$1 '$1 >= x' $2  | cut -f2 | sed -e 's/NNNNNNNNNN/	/g' | cut -f1 | uniq ${withinIndv_MinDepth} *.uniq.seqs  | sort --parallel=${task.cpus} -S 2G | uniq -c > total.f.uniq
        join -1 2 -2 1 -o 1.1,1.2,2.2 total.f.uniq total.fr | awk '{print $1 "\t" $2 "NNNNNNNNNN" $3}' | awk -v x=${acrossIndv_MinDepth} '$1 >= x' > uniq.k.${withinIndv_MinDepth}.c.${acrossIndv_MinDepth}.seqs
        # Clean up the mess
        rm total.uniqs total.u.* total.fr total.f.uniq* 

        # for reproducibility 
        # read clustering is sensitive to the order of reads inputed, therefore sort reads to ensure reproducibility down the line 
        sort -k1 -r -n --parallel=${task.cpus} -S 2G uniq.k.${withinIndv_MinDepth}.c.${acrossIndv_MinDepth}.seqs | \\
        cut -f2 > ${prefix}_totaluniqseq
        
        awk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' ${prefix}_totaluniqseq > uniq.full.fasta
        
        LENGTH=$(awk '!/>/' uniq.full.fasta  | \\
        awk '(NR==1||length<shortest){shortest=length} END {print shortest}')
        
        LENGTH=$(($LENGTH * 3 / 4))
        
        # use to seqtk seq to Create DUMMY Quality scores for fastp
        seqtk seq -F I uniq.full.fasta > uniq.fq

        MaxLen=$(awk '!/>/' uniq.full.fasta  | awk '(NR==1||length<shortest){shortest=length} END {print shortest}')
        
        # Remove anything more than 25% adapter sequence
        fastp -i uniq.fq -o uniq.fq1 --thread ${task.cpus} -Q -l $MaxLen & fastp.log 
        
        # Fastq back to Fasta
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.fq1 | \\
        paste - - | \\
        sort -k1,1 -V | \\
        tr "\t" "\n" > ${prefix}_uniq.fasta
        awk '!/>/' ${prefix}_uniq.fasta > ${prefix}_totaluniqseq
        rm uniq.fq*
        """
    } else {*/
        """
        awk -v x=${withinIndv_MinDepth} '(\$1 >= x)' *.uniq.seqs | \\
        cut -f2 | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' | \\
        awk -v x=${acrossIndv_MinDepth} '(\$1 >= x)' > uniq.k."${withinIndv_MinDepth}".c."${acrossIndv_MinDepth}".seqs
        
        # for reproducibility 
        sort -k1 -r -n -S 2G uniq.k.${withinIndv_MinDepth}.c.${acrossIndv_MinDepth}.seqs | \\
        cut -f2 > totaluniqseq 
        
        awk '{c= c + 1; print ">dDocent_Contig_" c "\\n" \$1}' totaluniqseq > uniq.full.fasta
        """
    //}
} 