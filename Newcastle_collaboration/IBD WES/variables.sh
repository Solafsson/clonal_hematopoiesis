#!/bin/bash

##### before running the scripts:
# 1. make directory called "resources" inside working_dir
# 2. list names of desired cram files in resources/cram_files.txt
# 3. place sitesList and phaseKey files in resources/
#####

working_dir="/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/IBD_WES/"

# script dir on NFS
script_dir="/nfs/users/nfs_a/ar31/IBD_WES_scripts/"

# output dir for bams and counts
output_dir="${working_dir}results/"

# directory for -o and -e
clusterOutput="${working_dir}clusterOutputFiles/"

# paths to software - in team152 anaconda environment
samtools_path="/software/hgi/installs/anaconda3/envs/team152/bin/samtools" # version 1.9
bam_readcount_path="/software/hgi/installs/anaconda3/envs/team152/bin/bam-readcount" # version 0.8

# location of cram files
genomes_dir="/lustre/scratch118/humgen/hgi/projects/ibdx10/crams_20200107/"

# names of crams from genomes_dir
cram_fileList="${working_dir}resources/cram_files.txt"

#cram_fileList="${working_dir}resources/failed_cram_files.txt"

# Reference genome in  FASTA format required by samtools view for bam to cram conversion
referenceFile="/lustre/scratch119/humgen/teams/anderson/users/so11/SVs_in_ibd/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta"

# variables for bam-readcount and summarising results
siteList="$working_dir/resources/SNV_site_list.txt"
#siteList_sum="$working_dir/resources/SNV_site_list_for_summarising_results.txt"

