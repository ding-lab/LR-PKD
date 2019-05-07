# LRPKD

## Overview

LRPKD is an analysis pipeline for haplotype-based variant detection with three independent modules including read alignment with Single Nucleotide Variation (SNV) and Structural Variation (SV) calling, downstream analysis, and the read evidence validation from linked-read sequencing data. LRPKD can identify germline variants with high accuracy and confidence in complex loci such as PKD1 by effectively eliminating the false-positives and capturing the true-positives. With phasing the variants by long-range phase blocks that most other methods lack, it provides the haplotype information across genes such as TSC2 and PKD1 associated with the TSC2/PKD1 contiguous gene syndrome.

## Modules

* 1.Read alignment & SNV, SV calling: processing the 10x Genomics Chromium sequencing data by aligning reads, de-duplication, applying filters and using the Chromium molecular barcodes to call and phase SNPs, indels, and structural variants in 10x Genomics Sequencing (Long Ranger, Zheng et al., 2016).

* 2.Downstream analysis: false-positive call investigation using built-in filters, VEP annotation (McLaren et al., 2016), formatting (vcf2maf), read count (bam-readcount), haplotype collection, and VAF calculation, and then returns a summary. 

* 3.Variant validation: checking the mapping quality, read level mismatch rate (Rausch et al., 2018), and barcode-phaseblock matching of the supporting reads.

LRPKD generates a variant calling report, validation report with figures and summary statistics. Github: https://github.com/ding-lab/LRPKD

## Installation

### Prerequisition
* Long Ranger and reference: please follow the instruction at https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation
* git
* conda

### Installation steps
* (1) &nbsp; &nbsp; &nbsp; &nbsp; Download the LRPKD source.
```sh
  git clone https://github.com/ding-lab/LRPKD.git
```
* 



samtools: conda install samtools
## bam-readcount: conda install bam-readcount
## VEP: conda install ensembl-vep or follow the instruction at https://github.com/Ensembl/ensembl-vep
## BLAST: conda install blast
## Alfred: conda install -c bioconda alfred









