# LR-PKD

## Overview

LR-PKD is an analysis pipeline for haplotype-based variant detection with three independent modules including read alignment with Single Nucleotide Variation (SNV) and Structural Variation (SV) calling, downstream analysis, and the read evidence validation from linked-read sequencing data. 

LR-PKD can identify germline variants with high accuracy and confidence in complex loci such as PKD1 by effectively eliminating the false-positives and capturing the true-positives. With phasing the variants by long-range phase blocks that most other methods lack, it provides the haplotype information across genes such as TSC2 and PKD1 associated with the TSC2/PKD1 contiguous gene syndrome.

## Modules

* Read alignment & SNV, SV calling: processing the 10x Genomics Chromium sequencing data by aligning reads, de-duplication, applying filters and using the Chromium molecular barcodes to call and phase SNPs, indels, and structural variants in 10x Genomics Sequencing (Long Ranger, Zheng et al., 2016).

* Downstream analysis: false-positive call investigation using built-in filters, VEP annotation (McLaren et al., 2016), formatting (vcf2maf), read count (bam-readcount), haplotype collection, and VAF calculation, and then returns a summary. 

* Variant validation: checking the mapping quality, read level mismatch rate (Rausch et al., 2018), and barcode-phaseblock matching of the supporting reads.

LR-PKD generates a variant calling report, validation report with figures and summary statistics. Github: https://github.com/ding-lab/LR-PKD

## Installation

### Prerequisition
* Long Ranger and reference: please follow the instruction at https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation
* git
* conda

### Installation steps
* Download the LR-PKD source.
```sh
  git clone https://github.com/ding-lab/LR-PKD.git
```
* Download Anaconda®.

If conda is not installed, download the installer for Anaconda® for python 2.7.x here (https://www.anaconda.com/download/) or download and install via command line as follows (please copy link for the .sh file appropriate for your operational system in the webpage provided):
```sh
  wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh
```

* Install anaconda.
```sh
  bash Anaconda2-5.2.0-Linux-x86_64.sh
```

* Create a virtual environment for LR-PKD.
```sh
  conda create --name LR-PKD bam-readcount ensembl-vep blast alfred samtools
```

* Intall VEP data libraries.

The ensembl-vep installs only the variant effect predictor (VEP) library code. To install data libraries, use the 'vep_install' command installed along with it. For example, to install the VEP library for human GRCh38 to a directory.
```sh
  vep_install -a cf -s homo_sapiens -y GRCh38 -c /output/path/to/GRCh38/vep --CONVERT
```
/output/path/to/GRCh38/vep should be identical to vep_cache_path in the config.sh.

## Run
* Users need to get into the folder with the LR-PKD source to run the pipeline.
```sh
  # To activate this environment, use:
  conda activate LR-PKD
  # To deactivate an active environment, use:
  conda deactivate
```

* Customize the arguments in config.sh as the inputs. If needed, users can update their list of genes with interest in Gene_interest.txt.

* Run the pipeline.
```sh
  bash LR-PKD_main.sh &
```
 
## Contact
Yize Li yize.li@wustl.edu









