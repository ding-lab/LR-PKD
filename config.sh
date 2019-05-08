## Arguments
input_sample_name=" " # Sample name
input_ref=" " # Path saving the reference (Users should follow the reference folder structure used in Long Ranger. If /reference/fasta/ref.fa, set /reference here) 
input_fastq=" " # Path saving the FASTQ files (Users should follow the data folder structure used in Long Ranger. If /sample/1.fastq.gz, set /sample here)
input_datatype=" " # WES or WGS
input_build="GRCh38" # Reference built GRCh38 or GRCh37
input_target_bed=" " # Target BED files. Users can choose to use their own BED or use the one provided in LR-PKD source
input_CNV_bed=" " # CNV BED files. Users can choose to use their own BED or use the one provided in LR-PKD source
working_path=$(pwd -P) # Current working path
input_gene_interest=${working_path}/Gene_interest.txt # Users can modified the list of genes with interest
vep_cache_path=" " # The path of VEP cache

### In order to directly call the tools, make sure the paths are saved in ./~bashrc. For instance,
### export PATH="/path/longranger-2.2.2/longranger-cs/2.2.2/bin:$PATH"
### export PATH="/path/ensembl-vep/vep:$PATH" 
