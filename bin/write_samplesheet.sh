#!/usr/bin/bash

# 
# Utility script for generating input csv file based on process_radtags output
#   NOTE. Make sure to remove *.rem.* files
#

echo "sample,fastq_1,fastq_2,umi_barcodes" > input.csv
paste -d',' <(for i in /mnt/d/nextflow_testing/radseq/simulated/Lpolyphemus/ddRAD_sbfi_mluci_chr26/*.1.fq.gz; do basename $i | cut -f1 -d'.' -; done)\
 <(ls /mnt/d/nextflow_testing/radseq/simulated/Lpolyphemus/ddRAD_sbfi_mluci_chr26/*.1.fq.gz) <(ls /mnt/d/nextflow_testing/radseq/simulated/Lpolyphemus/ddRAD_sbfi_mluci_chr26/*.2.fq.gz)\
 <(for i in /mnt/d/nextflow_testing/radseq/simulated/Lpolyphemus/ddRAD_sbfi_mluci_chr26/*.1.fq.gz; do if [[ "$i" =~ "Golden7".* ]]; then echo 'true'; else echo 'false'; fi; done)\
 >> input.csv
