## Arguments
input_sample_name=" " # Sample name
input_ref=" " # Path saving the reference (e.g. path/name, no need to include the *.fa file name) 
input_fastq=" " # Path saving the FASTQ files (e.g. path/name, no need to include the *.fastq.gz file name)
input_datatype=" " # WES or WGS
input_build="GRCh38" # Reference built GRCh38 or GRCh37
input_target_bed=" " # Target BED files. Users can choose to use their own BED or use the one provided in LR-PKD source
input_CNV_bed=" " # CNV BED files. Users can choose to use their own BED or use the one provided in LR-PKD source
working_path=$(pwd -P) # Current working path
input_gene_interest=${working_path}/Gene_interest.txt # Users can modified the list of genes with interest
vep_cache_path=" " # The path of VEP cache 
