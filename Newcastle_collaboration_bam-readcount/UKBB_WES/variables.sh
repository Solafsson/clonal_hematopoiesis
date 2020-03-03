#!/bin/bash

##### before running the scripts:
# 1. make directory called "resources" inside working_dir
# 2. list names of desired cram files in resources/all_queried_crams.txt
# 3. place sitesList file in resources/
#####

working_dir="/lustre/scratch114/projects/ukbb_team152/ar31/"

# script dir on NFS
script_dir="/nfs/users/nfs_a/ar31/UKBB_WES_scripts/"

# output dir for bams and counts
output_dir="${working_dir}results/"

# directory for -o and -e
clusterOutput="${working_dir}clusterOutputFiles/"

# paths to software - in team152 anaconda environment
samtools_path="/software/hgi/installs/anaconda3/envs/team152/bin/samtools" # version 1.9
bam_readcount_path="/software/hgi/installs/anaconda3/envs/team152/bin/bam-readcount" # version 0.8

# paths to cram files
cram_fileList="${working_dir}resources/all_queried_crams_test.txt"

#names of failed crams for re-running
#cram_fileList="${working_dir}resources/failed_cram_files.txt"

# Reference genome in  FASTA format required by samtools view for bam to cram conversion
referenceFile="/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/UKBB_control_WES/resources/Homo_sapiens_assembly38_chr.fasta"

# variables for bam-readcount and summarising results
siteList="${working_dir}resources/unaffected_CHIP_sites_UKBB_and_CHIP_chr.txt"
#siteList="/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/Interval_WES/resources/SNV_site_list.txt"
