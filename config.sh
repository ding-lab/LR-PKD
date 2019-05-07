
  input_sample_name=" " # Sample name
  input_ref=" " # Path saving the reference (e.g. path/name, no need to include the *.fa file name) 
  input_fastq=" " # Path saving the FASTQ files (e.g. path/name, no need to include the *.fastq.gz file name)
  input_datatype=" " # WES or WGS
  input_build="GRCh38" # Reference built
  input_target_bed="/diskmnt/Projects/Users/yize.li/Reference/S31285117_hs_hg38/S31285117_Regions.bed"
  input_CNV_bed="/diskmnt/Projects/Users/yize.li/PKD/Reference/agilent_v6r2_cnvfilter.bed"
  working_path=$(pwd -P) # Current working path
  input_gene_interest=${working_path}/Gene_interest.txt
  vep_cache_path="/diskmnt/Datasets/VEP"
