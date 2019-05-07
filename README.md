# LRPKD

## Overview

LRPKD is an analysis pipeline for haplotype-based variant detection with three independent modules including read alignment with Single Nucleotide Variation (SNV) and Structural Variation (SV) calling, downstream analysis, and the read evidence validation from linked-read sequencing data. LRPKD can identify germline variants with high accuracy and confidence in complex loci such as PKD1 by effectively eliminating the false-positives and capturing the true-positives. With phasing the variants by long-range phase blocks that most other methods lack, it provides the haplotype information across genes such as TSC2 and PKD1 associated with theTSC2/PKD1 contiguous gene syndrome.

## Modules

* Read alignment and SNV, SV calling: processing the 10x Genomics Chromium sequencing data by aligning reads, de-duplication, applying filters and using the Chromium molecular barcodes to call and phase SNPs, indels, and structural variants in 10x Genomics Sequencing (Long Ranger, Zheng et al., 2016).

* Downstream analysis: false-positive call investigation using built-in filters, VEP annotation (McLaren et al., 2016), formatting (vcf2maf), read count (bam-readcount), haplotype collection, and VAF calculation, and then returns a summary. 

* Variant validation: checking the mapping quality, read level mismatch rate (Rausch et al., 2018), and barcode-phaseblock matching of the supporting reads.

LRPKD generates a validation report with figures and summary statistics. Github: https://github.com/ding-lab/LRPKD









