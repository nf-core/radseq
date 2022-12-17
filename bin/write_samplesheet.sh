#!/usr/bin/bash


echo "sample,fastq_1,fastq_2,umi_barcodes" > input.csv
paste -d',' <(for i in $(pwd)/data/demult/*.1.fq.gz; do basename $i | cut -f1 -d'.' -; done)\
 <(ls $(pwd)/data/demult/*.1.fq.gz) <(ls $(pwd)/data/demult/*.2.fq.gz)\
 <(for i in $(pwd)/data/demult/*.1.fq.gz; do if [[ "$i" =~ "Golden7".* ]]; then echo 'true'; else echo 'false'; fi; done)\
 >> input.csv