# LRPKD

## Overview

LRPKD is an analysis pipeline for haplotype-based variant detection with three independent modules including read alignment with Single Nucleotide Variation (SNV) and Structural Variation (SV) calling, downstream analysis, and the read evidence validation from linked-read sequencing data. 

LRPKD can identify germline variants with high accuracy and confidence in complex loci such as PKD1 by effectively eliminating the false-positives and capturing the true-positives. With phasing the variants by long-range phase blocks that most other methods lack, it provides the haplotype information across genes such as TSC2 and PKD1 associated with the TSC2/PKD1 contiguous gene syndrome.

## Modules

* Read alignment & SNV, SV calling: processing the 10x Genomics Chromium sequencing data by aligning reads, de-duplication, applying filters and using the Chromium molecular barcodes to call and phase SNPs, indels, and structural variants in 10x Genomics Sequencing (Long Ranger, Zheng et al., 2016).

* Downstream analysis: false-positive call investigation using built-in filters, VEP annotation (McLaren et al., 2016), formatting (vcf2maf), read count (bam-readcount), haplotype collection, and VAF calculation, and then returns a summary. 

* Variant validation: checking the mapping quality, read level mismatch rate (Rausch et al., 2018), and barcode-phaseblock matching of the supporting reads.

LRPKD generates a variant calling report, validation report with figures and summary statistics. Github: https://github.com/ding-lab/LRPKD

## Installation

### Prerequisition
* Long Ranger and reference: please follow the instruction at https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation
* git
* conda

### Installation steps
* Download the LRPKD source.
```sh
  git clone https://github.com/ding-lab/LRPKD.git
```
* Download and install Anaconda®.
If conda is not installed, download the installer for Anaconda® for python 2.7.x here (https://www.anaconda.com/download/) or download and install via command line as follows (please copy link for the .sh file appropriate for your operational system in the webpage provided):
```sh
  wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh
```

Install anaconda:
```sh
  bash Anaconda2-5.2.0-Linux-x86_64.sh
```

* Prepare for LRPKD
Create a virtual environment:
```sh
  conda create --name LRPKD bam-readcount ensembl-vep blast alfred samtools
```

## Run
* Users need to get into the folder with the LRPKD source and define config.sh

* bash LRPKD_main.sh &

 









