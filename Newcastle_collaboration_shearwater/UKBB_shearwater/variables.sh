#!/bin/bash

##### before running the scripts:
# 1. make directory called "resources" inside working_dir
# 2. list names of desired cram files in resources/all_queried_crams.txt
# 3. place sitesList file in resources/
#####

working_dir="/lustre/scratch114/projects/ukbb_team152/ar31/"

# script dir on NFS
script_dir="/nfs/users/nfs_a/ar31/shearwater/shearwater_UKBB/"

# output dir for bams and counts
output_dir="${working_dir}bams/"

# directory for -o and -e
clusterOutput="${working_dir}clusterOutputFiles_shearwater/"

# paths to software - in team152 anaconda environment
samtools_path="/software/hgi/installs/anaconda3/envs/team152/bin/samtools" # version 1.9
bam_readcount_path="/software/hgi/installs/anaconda3/envs/team152/bin/bam-readcount" # version 0.8

# paths to cram files
cram_fileList="${working_dir}resources/UKBB_cram_list_shearwater.txt"

#names of failed crams for re-running
#cram_fileList="${working_dir}resources/failed_cram_files.txt"

# Reference genome in  FASTA format required by samtools view for bam to cram conversion
referenceFile="/lustre/scratch115/resources/ref/Homo_sapiens/GRCh38_15_no_chr/Homo_sapiens.GRCh38_15.no_chr.fa"

# site list for converting bam to cram - 10bp window around the CHIP site
convert_siteList="${working_dir}resources/CHIP_sites_UKBB_for_conversion.txt"

# site list for bam-readcount and summarising results - exact location in hg38 (no chr notation for UKBB)
siteList="${working_dir}resources/unaffected_CHIP_sites_UKBB_and_CHIP_no_chr.txt"
