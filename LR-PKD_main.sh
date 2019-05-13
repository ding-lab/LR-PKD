#!/bin/bash

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# LR-PKD
# Author: Yize Li yize.li@wustl.edu
# Updated 05/12/2019

## Overview
### LR-PKD is an analysis pipeline for haplotype-based variant detection with three independent modules including read alignment with Single Nucleotide Variation (SNV) and Structural Variation (SV) calling, downstream analysis, and the read evidence validation from linked-read sequencing data.
### LR-PKD can identify germline variants with high accuracy and confidence in complex loci such as PKD1 by effectively eliminating the false-positives and capturing the true-positives. With phasing the variants by long-range phase blocks that most other methods lack, it provides the haplotype information across genes such as TSC2 and PKD1 associated with the TSC2/PKD1 contiguous gene syndrome.
### LR-PKD generates a variant calling report, validation report with figures and summary statistics. Github: https://github.com/ding-lab/LR-PKD

## Modules
### 1.Read alignment & SNV, SV calling: processing the 10x Genomics Chromium sequencing data by aligning reads, de-duplication, applying filters and using the Chromium molecular barcodes to call and phase SNPs, indels, and structural variants in 10x Genomics Sequencing (Long Ranger, Zheng et al., 2016).
### 2.Downstream analysis: false-positive call investigation using built-in filters, VEP annotation (McLaren et al., 2016), formatting (vcf2maf), read count (bam-readcount), haplotype collection, and VAF calculation, and then returns a summary.
### 3.Variant validation: checking the mapping quality, read level mismatch rate (Rausch et al., 2018), and barcode-phaseblock matching of the supporting reads.

## Installation
### Prerequisition
#### Long Ranger and reference: please follow the instruction at https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation
#### git
#### conda

### Installation steps
#### 1.Download the LR-PKD source: git clone https://github.com/ding-lab/LR-PKD.git
#### 2.Download Anaconda®: wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh
#### 3.Install anaconda: bash Anaconda2-5.2.0-Linux-x86_64.sh
#### 4.Prepare for LR-PKD Create a virtual environment: conda create --name LR-PKD bam-readcount ensembl-vep blast alfred samtools

## Run
### Users need to get into the folder with the LR-PKD source to run the pipeline. Customize the arguments in config.sh as the inputs. If needed, users can update their list of genes with interest in Gene_interest.txt.
### Run the pipeline: bash LR-PKD_main.sh &
########################################################################################################################################################################################################################################################################################################################################################################################################################################################

# Pipeline
source activate LR-PKD # Switch to virtual environment
source config.sh # Argumments set by the users in config.sh

# Module 1: Long Ranger is called to process the 10x Genomics Chromium se-quencing data to align reads, de-duplication, filtering and using the Chromium molecular barcodes to call and phase SNPs, indels, and structural variants in 10X genome sequencing (longranger wgs) and 10X exome sequencing (longranger target)Long Ranger is called to process the 10x Genomics Chromium se-quencing data to align reads, de-duplication, filtering and using the Chromium molecular barcodes to call and phase SNPs, indels, and structural variants in 10X genome sequencing (longranger wgs) and 10X exome sequencing (longranger target).
## Input
### If WES, BED file associated with the pulldown used for this Chromium library and a BED file indicating poorly performing targets or problematic genomic regions need to be provided (can be downloaded at https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/targeted).
#### The target BED file supplied to the pipeline via the --targets is used in computing metrics such as on-target coverage. The SV and CNV calling algorithms in Long Ranger also use the target BED file to define regions of interest. For Agilent SureSelect Human All Exon V6, we strongly recommend using the latest BED file released by Agilent. The BED file is available as SureSelect Human All Exon V6 r2 from Agilent.
#### Long Ranger includes a CNV caller that detects exon-scale deletions in targeted regions. It is important to supply the --cnvfilter argument a BED file which masks problematic regions where baits perform poorly, to prevent false-positive calls. We have created a cnvfilter file tailored to the SureSelect Human All Exon V6 r2 BED file, which is available for download: Agilent Exome V6 r2 CNV Filter BED.

## Calling Long Ranger
if [ $input_datatype = "WGS" ]; then
	echo "Th WGS sample of ${input_sample_name} was provided to LR-PKD."
	echo "The variant calling results will be saved in the .${input_sample_name}/outs under the current working path ${working_path}."
	for i in ${input_sample_name}
	do
	nohup longranger wgs --id=${i} --reference=${input_ref} --fastqs=${input_fastq} --vcmode=freebayes --somatic --localcores=60
	done
elif [ $input_datatype = "WES" ]; then
	echo "Th WES sample of ${input_sample_name} was provided to LR-PKD."
	echo "The variant calling results will be saved in the .${input_sample_name}/outs under the current working path ${working_path}."
	for i in ${input_sample_name}
	do
	nohup longranger targeted --id=${i} --reference=${input_ref} --fastqs=${input_fastq} --targets=${input_target_bed} --cnvfilter=${input_CNV_bed} --vcmode=freebayes --somatic --localcores=30
	done
else
	echo "The data type should be either wgs or wes. Wrong data type was provided by the user. Please double check."
fi

# Module 2: Based on the variant calling results, LR-PKD provides the downstream analysis including VEP annotation,  driver gene extraction, read count and haplotype information collection, VAF calculation and return a summary table with all the information.
## VEP annotation
### Path of the Long Ranger phased VCF in the current working folder
path_variant=${working_path}/${input_sample_name}/outs

### Extract the information from the Long Ranger phased VCF
zcat ${path_variant}/phased_variants.vcf.gz > ${working_path}/input.vcf # The VCF file should contain the column with the name of input_sample_name

### VEP annotation argument
reffasta="${input_ref}/fasta/genome.fa" # The reference genome should be identical to the one downloaded for running Long Ranger from 10x
in="${working_path}/input.vep.vcf"
out="${working_path}/input.vep.maf" # The output is a VEP annotated VCF file.
vep_path=$(which vep)
vep_only_path=${vep_path%/vep}

### Docker image
 #docker run -t -i --rm -u $(id -u):$(id -g)       \
 docker run --rm -u $(id -u):$(id -g)       \
            -v ${vep_cache_path}:/home/vep/.vep    \
            -v $PWD:/preprocess_root                \
            ensemblorg/ensembl-vep:release_91.3      \
                vep --cache --offline                 \
                    --assembly ${input_build}                  \
                    --cache_version 91                  \
                    --fork 4                             \
                    --no_escape --hgvs --vcf              \
                    -i /preprocess_root/input.vcf          \
                    -o /preprocess_root/input.vep.vcf

## Formatting: VCF to MAF; VEP annotation
mkdir temp
mkdir fasta
cp ${input_ref}/fasta/* fasta

### Change the format using vcf2maf.pl which was adjusted from the version at https://github.com/mskcc/vcf2maf
if [ ! -f ${working_path}/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz ]; then
	echo "ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz not found, skip ExAC annotation."
	perl ${working_path}/vcf2maf.v1.6.16.no_vep_no_ref.pl --input-vcf ${in} --output-maf ${out} --tumor-id PKD_K --normal-id PKD_N --ref-fasta ${reffasta} --filter-vcf 0 --ncbi-build GRCh38 #--retain-info GT,DP,RO,QR,AO,QA,GL,BX,PS
else
	perl ${working_path}/vcf2maf.v1.6.16.no_vep_no_ref.pl --input-vcf ${in} --output-maf ${out} --tumor-id PKD_K --normal-id PKD_N --ref-fasta ${reffasta} --filter-vcf ${working_path}/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --ncbi-build GRCh38 #--retain-info GT,DP,RO,QR,AO,QA,GL,BX,PS
fi

sed '/^#/ d' ${out} > input_d1.vep.maf # Remove the first row of output

#### Extract the phasing genotype, barcode and phase block information to return file input.vep.GTPSBX.maf
python GT_PS_BX_extraction.py input.vep.vcf # This requires pysam and pandas
cut -f2- input.vep.GTPSBX.temp.maf > input.vep.GTPSBX.maf

# Extract the results of certain genes with interest (User can provide a list of genes with interest; if not provided, all genes in the results will be kept)
if [ ! -f ${input_gene_interest} ]; then
	echo "Gene list with interest is not provided by the user. All the genes in the results will be kept."
	cut -f 1 input.vep.GTPSBX.maf | uniq > Gene_all.txt
	perl Gene_noextraction.pl # This will return two files Genes_interest_all.results and Genes_interest_coding.results.
else
	perl Gene_extraction.pl # This will return two files Genes_interest_all.results and Genes_interest_coding.results.
fi

# Replace empty space with NA in .results file
awk -F"\t" -v OFS="\t" '{
	for (i=1;i<=NF;i++) {
		if ($i == "") $i="NA"}
	print $0}' Genes_interest_all.results > Genes_interest_all.results.new && mv Genes_interest_all.results.new Genes_interest_all.results

awk -F"\t" -v OFS="\t" '{
	for (i=1;i<=NF;i++) {
		if ($i == "") $i="NA"}
	print $0}' Genes_interest_coding.results > Genes_interest_coding.results.new && mv Genes_interest_coding.results.new Genes_interest_coding.results

# Reformat
perl -e 'while(<>){chomp; @l = split(/\t/,); print "$l[2]\t$l[3]\t$l[3]\t\.\t".join("\_\_",@l)."\t\+\n";}' ${working_path}/Genes_interest_all.results > ${working_path}/Genes_interest_all.bed
perl -e 'while(<>){chomp; @l = split(/\t/,); print "$l[2]\t$l[3]\t$l[3]\t\.\t".join("\_\_",@l)."\t\+\n";}' ${working_path}/Genes_interest_coding.results > ${working_path}/Genes_interest_coding.bed

# Get the bam-readcount
bc_path=$(which bam-readcount)
${bc_path} -l ${working_path}/Genes_interest_all.bed -f ${reffasta} ${path_variant}/phased_possorted_bam.bam  > Genes_interest_all.readcount
${bc_path} -l ${working_path}/Genes_interest_coding.bed -f ${reffasta} ${path_variant}/phased_possorted_bam.bam  > Genes_interest_coding.readcount

# Combine all information into a final report and add the built-in LR-PKD filter
perl final_report.pl Genes_interest_all
perl final_report.pl Genes_interest_coding

# Add column Chr:Pos
awk -v OFS="\t" 'NR==1{$23="Chr:Pos"}NR>1{$23=$3":"$4"-"$5}1' Genes_interest_all.output.tsv > Genes_interest_all.final.tsv
awk -v OFS="\t" 'NR==1{$23="Chr:Pos"}NR>1{$23=$3":"$4"-"$5}1' Genes_interest_coding.output.tsv > Genes_interest_coding.final.tsv

# Module 3: LR-PKD validates the PKD1 read evidence by checking the mapping quality, read level mismatch rate, cross-checking barcodes and the corresponding phase block IDs and aligning reads supporting the muta-tions to the reference genome to make sure the reads truly align to PKD1 over its six pseudogenes.

## Read collection, mapping quality, read level mismatch rate etc (Alfred Alignment Quality Control)
### Agilent BED file (hg38) in PKD1 has been downloaded to LR-PKD folder. If needed, users can customize the file or download the hg19 version with naming it as PKD1.bed.
samtools_path=$(which samtools)
alfred_path=$(which alfred)
${samtools_path} view -b ${path_variant}/phased_possorted_bam.bam -L PKD1.bed -o PKD1_10x${input_datatype}.bam
${samtools_path} index PKD1_10x$input_datatype.bam
${alfred_path} qc -r ${reffasta} -o PKD1_10x${input_datatype}.tsv.gz PKD1_10x${input_datatype}.bam
git clone https://github.com/tobiasrausch/alfred.git
Rscript alfred/scripts/stats.R PKD1_10x${input_datatype}.tsv.gz

## Read alignment (BLAST)
makeblastdb -in ${reffasta} -dbtype nucl -parse_seqids # Building a BLAST database with local sequences

cat Genes_interest_coding.final.tsv | awk '{ print $23}' | grep chr > variant_list # Getting the chr:start-end
idx=0
while read chrpos; do
  idx=$((idx+1))
  samtools view -b PKD1_10x${input_datatype}.bam ${chrpos} > temp.bam
  samtools bam2fq temp.bam | seqtk seq -A > temp.fa
  blastn -db ${reffasta} -query temp.fa -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > blast_${idx}.txt
  rm temp.bam
  rm temp.fa
done <variant_list
mv variant_list temp

## MAPQ (ambiguity checking)
while read CHR START END; do
  samtools view -b PKD1_10x${input_datatype}.bam ${CHR}:${START}-${END} -o MAPQ_${START}.bam
  samtools view MAPQ_${START}.bam | cut -f5 > MAPQ_${START}.out
  rm MAPQ_${START}.bam
done <PKD1.bed

cat MAPQ_*.out > MAPQ_PKD1_exon.txt
rm MAPQ_*.out

## MAPQ density plot
Rscript MAPQ_stat_histogram.R ${working_path}/MAPQ_PKD1_exon.txt > ${working_path}/MAPQ_PKD1_exon.txt.stat.txt

## Barcode and phase-block matching is addressed in the final output in module 2 (Genes_interest_coding.final.tsv and Genes_interest_all.final.tsv).

# Module 4: Clean-up
## Create the output folders
mkdir outputs
mkdir outputs/1.Variant_calling
mkdir outputs/2.Downstream_analysis
mkdir outputs/3_1.Validation_summary_stat
mkdir outputs/3_2.Validation_BLAST
mkdir outputs/3_3.Validation_MAPQ

## Move outputs
mv ${working_path}/input.vcf ${working_path}/outputs/1.Variant_calling
mv Genes_interest_*.final.tsv ${working_path}/outputs/2.Downstream_analysis
mv PKD1_*.tsv.gz.* ${working_path}/outputs/3_1.Validation_summary_stat
mv blast_* ${working_path}/outputs/3_2.Validation_BLAST
mv MAPQ_PKD1_exon.txt* ${working_path}/outputs/3_3.Validation_MAPQ

## Move all intermediate files to temp folder.
mv ${working_path}/Genes_interest_* ${working_path}/temp
mv ${working_path}/input.vep.* ${working_path}/temp
mv ${working_path}/input_d1.vep.maf ${working_path}/temp
mv ${working_path}/PKD1_* ${working_path}/temp
mv ${working_path}/fasta ${working_path}/temp
mv ${working_path}/alfred ${working_path}/temp

## Users can decide whether to delete the temp folder manually to save space.
### rm -r ${working_path}/temp
