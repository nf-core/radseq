# nf-core/radseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
    - [FastP](#fastp) - trim low-quality reads, umi-barcodes, adapters
- [Denovo Reference Construction](#denovo-reference-construction)
    - [Prepare forward reads](#prepare-forward-reads) - combine forward and reverse sequences separated by 'NNNNNNNNNN'
    - [Combine uniqe reads](#combine-uniq-reads) - retain reads present n number of times between and across individuals
    - [Seqtk](#seqtk-seq) - write dummy fasta file
    - [Denovo FastP](#denovo-fastp) - trim adapters
    - [CD-HIT-EST](#cdhit-est) - cluster similar sequences
    - [CD-HIT-EST to Rainbow div](#cdhit-to-rbdiv) - convert cdhit file output into rainbow div input file format
    - [Rainbow div](#rainbow-div) - distiguish sequence errors from heterozygote or variants between repetitive sequences
    - [Rainbow merge](#rainbow-merge) - merge potential heterozygous clusters
    - [Write Fasta](#write_fasta) - convert rainbow merge file output into fasta format
- [Alignment](#alignment)
    - [SAMtools](#samtools) - Sort, index and obtain alignment statistics
    - [BWA](#bwa) - short read aligner
    - [BWAMEM2](#bwa-mem2) - a faster short read aligner 
    - [UMI-tools dedup](#umi-tools-dedup) - UMI-based deduplication
- [Freebayes Intervals](#freebayes-intervals)
    - [BEDtools bamtobed](#bedtools-bamtobed): converts bam file into bed datastructure
    - [BEDOPS merge](#bedops-merge): merge indvidual bed files into a single bed file
    - [BEDtools sort](#bedtools-sort): sorts bed files
    - [BEDtools coverage](#bedtools-coverage): counts read depth
    - [BEDtools merge](#bedtools-merge): merges indv bed files and takes sum the read depths
    - [BEDtools makewindows](#bedtools-makewindows): split regions with coverage above `--max_read_coverage_to_split` to half read length
    - [BEDtools intersect](#bedtools-intersect): removes any overlapping regions between split reads and all merged reads
    - [Create intervals](#create-intervals): write regions for input to `freebayes`
- [Variant Calling](#variant-calling)
    - [FreeBayes](#freebayes) - a bayesian genotyper tool
- [Quality Control and Preprocessing](#qc-and-reporting)
    - [FastQC](#fastqc) - Raw read QC
    - [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
    - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

# Directory Structure

The default directory structure is as follows

```
{outdir}
├── denovo
│   ├── alignments
│   │   ├── samtools_index
│   │   ├── samtools_merge
│   │   ├── samtools_stats
│   │   └── umitools_dedup
│   │       └── stats
│   ├── reference
│   │   ├── cdhit
│   │   ├── cdhit_to_rbdiv
│   │   ├── rainbow_div
│   │   ├── rainbow_merge
│   │   └── write_fasta
│   └── variant_calling
├── fastp
├── fastqc
├── multiqc
│   └── multiqc_data
├── pipeline_info
└── reference
    ├── alignments
    │   ├── samtools_index
    │   ├── samtools_merge
    │   ├── samtools_stats
    │   └── umitools_dedup
    │       └── stats
    └── variant_calling
```

# Pre-Processing

Radseq pre-processes reads prior to the alignment step.

### FastP

[FastP](https://github.com/OpenGene/fastp) is a tool designed to be an all-in-one preprocessor for FastQ files. You can enable the saving of trimmed fq files in output directory through `--save_trimmed=true`.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/fastp/`
    * `*.fq.gz`: trimmed fq files

</details>

# Denovo Reference Construction

Radseq supports the construction of psuedoreference using a conglomerate of open source tools. By default the resulting fasta file is only outputed to enable the printing of intermediate files into `{outdir}/denovo/reference/` use `--denovo_intermediate_files=true`. 

### Prepare forward reads

Prior to clustering forward and reverse reads are joined into one sequence and seperated with a `NNNNNNNNNN`. Reads are reduced based on number of presences within an individual and among all individuals and combined into one file. 

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/unique_sequences`
    * `*.uniq.seqs`: all unique sequences in an individual 
    * `*_uniq.full.fasta`: all sequences

</details>

### Seqtk seq

[Seqtk](https://github.com/lh3/seqtk) is a tool for processing FASTA or FASTQ file formats.
<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/seqtk`
    * `*.seqtk-seq`: dummy fasta file

</details>

### Denovo FastP

Denovo Fastp uses [FastP](https://github.com/OpenGene/fastp) as a tool to trim any unwanted sequences prior to clustering like adapter content.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/fastp`
    * `*.uniq.fasta`: fasta format
    * `*.totaluniqseq`: all unique sequences remaining after data cutoffs and adapter trimming

</details>

### Cdhit est

[CD-HIT](https://sites.google.com/view/cd-hit) is used for clustering and comparing nucleotide sequences. 

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/cdhit`
    * `*_cdhit.logs`: log output from cdhit-est
    * `*.clstr`: clstr output used to convert into rainbow div format

</details>

### CD-HIT to Rainbow div

This module converts `CD-HIT` output into input for `Rainbow div`.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/cdhit_to_rbdiv`
    * `*.sort.contig.cluster.ids`: intermediate file used for conversion of `cd-hit` to `Rainbow` input and facilitates reproducibility
    * `*.contig.cluster.totaluniqseq`: used during assembly
    * `*.rcluster`: input for `Rainbow`

</details>

### Rainbow div

[Rainbow div](https://github.com/ChongLab/rainbow) is a tool used to divide heterozygote calls into putative haplotypes based on minimum thresholds passed as arguments. 

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/rainbow_div`
    * `*_rbdiv.out`: rainbow div output file
    * `*_rbdiv.log`: log file

</details>

### Rainbow merge

[Rainbow merge](https://github.com/ChongLab/rainbow) is a tool used to merge reads from `Rainbow div` and assemble into contigs based on minimum and maximum thresholds that are passed as arguments. 

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/rainbow_merge`
    * `*_rainbow.fasta`: final fasta file used in subsequent processes

</details>

### Write fasta

This module converts `Rainbow merge` into fasta format. These files are outputed by default to: 

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/write_fasta`
    * `*_rainbow.fasta`: denovo fasta file containing contig sequences

</details>

# Alignment

## Indices
enable the saving of reference indices with `--save_reference_indices true` generate from `samtools` and `bwa` for variant calling and short-read alignment respectiviely. 

### samtools faidx

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/reference/{samtools}/index/`
    * `*.fai`: samtools fai index

</details>

### bwa index

**Not Working**

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/reference/{aligner}/index/`
    * `*`: rainbow merge output file

</details>

## Aligners

To enable the output of bam files use `--save_bam_files=true`. Output of bam files from different alignment methods follow the same output structure. 

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/reference/{aligner}/bam/`
    * `*.bam`: sorted bam files

</details>

### BWA

[BWA](https://github.com/lh3/bwa) is a tool for mapping sequencies with low divergence against a reference genome. Aligned reads are then potentially filtered and coordinate-sorted using [samtools](https://github.com/samtools/samtools).

### BWA-mem2

[BWA-mem2](https://github.com/bwa-mem2/bwa-mem2)  is a tool next version of `bwa-mem` for mapping sequencies with low divergence against a reference genome with increased processing speed (~1.3-3.1x). Aligned reads are then potentially filtered and coordinate-sorted using [samtools](https://github.com/samtools/samtools).

### UMI-tools dedup

[UMI-tools dedup](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) is a tool for deduplicating reads and ensuring only a single representative read is retained in the bam file. Files are outputed using `--save_bam_files=true` as a result of many individual bed files being passed through.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/umitools_dedup`
    * `*.bam`: bam file
    * `*.tsv`: output statistics file

</details>

# FreeBayes Interval Construction

radseq supports multithreading of FreeBayes through the `bam_intervals_bedtools.nf` subworkflow that uses a collection of [BEDtools](https://bedtools.readthedocs.io/en/latest/index.html) software and [BEDOPS merge](https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#merge-m-merge).

To enable the saving of files generated for interval construction in the output-folder set `--save_bed_intervals=true`.

For large datasets, it may be useful to randomly subsample the number of individuals going into `bam_intervals_bedtools.nf` subworkflow using `--subset_intervals_channel=<integer>`. Particularly at the `BEDtools merge` module the memory footprint can be large. 

### BEDtools bamtobed

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedtools_bamtobed`
    * `*.bed`: individual bed file

</details>

### BEDOPS merge

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedops_merge`
    * `*.bed`: merged single bed file

</details>

### BEDtools sort

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedtools_sort`
    * `*.bed`: single sorted bed file

</details>

### BEDtools coverage

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedtools_coverage`
    * `*.cov`: Individual bed files with fourth column describing number of reads in a particular region

</details>

### BEDtools merge cov

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedtools_merge`
    * `*.cov`: single bed file with fourth column describing number of reads in a particular region

</details>

### BEDtools makewindows

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedtools_makewindows`
    * `*_cov.low.stats`: regions less than `max_read_coverage_to_split`
    * `*_cov.high.stats`: regions equal to or greater than `--max_read_coverage_to_split`
    * `*.tab`: regions from `*_cov.high.stats` split to half their read length

</details>

### BEDtools intersect

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedtools_intersect`
    * `*.bed`: Final single bed file to be subsequently split

</details>

### Create intervals

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/bedtools_intersect`
    * `mapped.*.bed`: Split bed files passed into FreeBayes for multithreading

</details>


# Variant Calling

### FreeBayes

[FreeBayes](https://github.com/freebayes/freebayes) is a bayesian genetic variant detector capable of genotyping SNPs, indels, MNPs, and complex events smaller than the length of short-read sequencing alignment. 

By default radseq will output a single vcf joined from independent runs. To enable the outputting of VCF files based on regions passed to FreeBayes use `--save_freebayes_intervals=true`.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/variant_calling/`
    * `*.vcf.gz`: Final VCF

* `{outdir}/{method}/variant_calling/intervals`
    * `*_*.vcf.gz`: unsorted VCF intervals
    * `*_sort_*.vcf.gz`: sorted VCF intervals
    * `*_sort_*.vcf.gz.tbi`: sorted VCF intervals tabix index

</details>

# Quality Control and Visualization

### bcftools stats

[bcftools stats](https://samtools.github.io/bcftools/bcftools.html) is a tool for collecting statistics from a VCF or BCF file that can be interpretted by MultiQC.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/variant_calling/`
    * `*_stats.txt`: 

</details>

### FastP

[FastP](https://github.com/OpenGene/fastp) is a tool designed to be an all-in-one preprocessor for FastQ files. Statistics are passed to MultiQC.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/fastp/`
    * `*_fastp.html`: Fastp report containing quality metrics
    * `*_fastp.log`: Log output containing statistics

</details>

### samtools

radseq supports generation of statistics from alignment files with [samtools stats](http://www.htslib.org/doc/samtools-stats.html), [samtools flagstat](http://www.htslib.org/doc/samtools-flagstat.html), [samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html).

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/{method}/alignments/{aligner}/samtools_stats`
    * `*.flagstat`: 
    * `*.stats`:
    * `*.idxstats`:

</details>

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
