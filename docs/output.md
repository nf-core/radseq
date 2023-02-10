# nf-core/radseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
    - [FastP](#fastp) - trim low-quality reads, umi-barcodes, adapters
- [Denovo Reference Construction](#denovo-reference-construction)
    - [Prepare Forward Reads](#prepare-forward-reads) - combine forward and reverse sequences separated by 'NNNNNNNNNN'
    - [Combine Uniqe Reads](#combine-uniq-reads) - retain reads present n number of times between and across individuals
    - [Seqtk](#seqtk-seq) - write dummy fasta file
    - [Denovo Fastp](#denovo-fastp) - trim adapters
    - [CDHIT-est](#cdhit-est) - cluster similar sequences
    - [cdhit_to_rbdiv](#cdhit-to-rbdiv) - convert cdhit file output into rainbow div input file format
    - [Rainbow div](#rainbow-div) - distiguish sequence errors from heterozygote or variants between repetitive sequences
    - [Rainbow merge](#rainbow-merge) - merge potential heterozygous clusters
    - [write_fasta](#write_fasta) - convert rainbow merge file output into fasta format
- [Alignment](#alignment)
    - [SAMtools](#samtools) - Sort, index and obtain alignment statistics
    - [BWA](#bwa) - short read aligner
    - [BWA-mem2](#bwa-mem2) - a faster short read aligner 
    - [UMI-tools dedup](#umi-tools-dedup) - UMI-based deduplication
- [Freebayes Intervals](#freebayes-intervals)
    - [bedtools bamtobed](#bedtools-bamtobed): converts bam file into bed datastructure
    - [bedops merge](#bedops-merge): merge indvidual bed files into a single bed file
    - [bedtools sort](#bedtools-sort): sorts bed files
    - [bedtools coverage](#bedtools-coverage): counts read depth
    - [bedtools merge](#bedtools-merge): merges indv bed files and takes sum the read depths
    - [bedtools makewindows](#bedtools-makewindows): split regions with coverage above `--max_read_coverage_to_split` to half read length
    - [bedtools intersect](#bedtools-intersect): removes any overlapping regions between new bed file from `bedtools makewindows` and from `bedtools merge`
    - [create intervals](#create-intervals): write regions for input to `freebayes`
- [Variant Calling](#variant-calling)
    - [Freebayes](#freebayes) - a bayesian genotyper tool
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

# Preprocessing

Radseq pre-processes reads prior to the alignment step.

### FastP

[FastP](https://github.com/OpenGene/fastp) is a tool designed to be an all-in-one preprocessor for FastQ files. You can enable the saving of trimmed fq files in output directory through `--save_trimmed=true`.

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/fastp/`
    * `*_fastp.html`: Fastp report containing quality metrics.
    * `*_fastp.log`: Log output containing statistics
    * `*.fq.gz`: trimmed fq files

</details>

# Denovo Reference Construction

Radseq supports the construction of psuedoreference using a conglomerate of open source tools. By default the resulting fasta file is only output to enable the output of intermediate files into `{outdir}/denovo/reference/` use `--denovo_intermediate_files=true`. 

### Prepare forward reads

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference`
    * `*.uniq.seqs`: all unique sequences

</details>

### Seqtk seq

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/seqtk`
    * `*.seqtk-seq`: dummy fasta file

</details>

### Denovo FastP

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/fastp`
    * `*.uniq.fasta`: fasta format
    * `*.totaluniqseq`: all unique sequences remaining after data cutoffs and adapter trimming

</details>

### Cdhit est

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/cdhit`
    * `*_cdhit.logs`: log output from cdhit-est
    * `*.clstr`: clstr output used to convert into rainbow div format

</details>

### Cdhit to rbdiv

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/cdhit_to_rbdiv`
    * `*.sort.contig.cluster.ids`: file used for conversion of `cd-hit` to `Rainbow` input
    * `*.contig.cluster.totaluniqseq`: used during assembly
    * `*.rcluster`: input for `Rainbow`

</details>

### Rainbow div

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/rainbow_div`
    * `*_rbdiv.out`: rainbow div output file
    * `*.log`: log file

</details>


### Rainbow merge

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/rainbow_div`
    * `*_rainbow.fasta`: final fasta file used in subsequent processes

</details>

### Write fasta

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/rainbow_div`
    * `*_rbmerge.out`: rainbow merge output file
    * `*_rbmerge.log`: log file


</details>

# Freebayes Intervals

# Alignment

## Indices
enable the saving of reference indices with `--save_reference_indices true` generate from `samtools` and `bwa` for variant calling and short-read alignment respectiviely. 

### samtools faidx

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/index`
    * `*.fai`: samtools fai index

</details>

### bwa index

<details markdown="1">
<summary>Output files</summary>

* `{outdir}/denovo/reference/index`
    * `*_rbmerge.out`: rainbow merge output file

</details>

# Variant Calling

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
