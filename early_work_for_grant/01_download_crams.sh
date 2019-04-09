#!/bin/bash

## Usage: This script takes in a 'batch-file' and downloads all the data from irods
##        You don't need to submit this script to the farm. Run it locally.
## Input:
## Output: The cram files for the batch have been downloaded and indexed. Data is stored in output_dir


download_fileList=$1

working_dir=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/
script_dir=/nfs/users/nfs_s/so11/phd/clonal_hematopoiesis/
samtools_path=/software/hgi/pkglocal/samtools-1.3.1/bin/samtools
## Reference genome (hg38 for ibd data)
referenceFile=/lustre/scratch119/humgen/teams/anderson/users/so11/SVs_in_ibd/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta


output_dir=${working_dir}genomes/
mkdir -p ${output_dir}

mkdir -p ${working_dir}clusterOutputFiles/
mkdir -p ${working_dir}clusterOutputFiles/download_and_index/
mkdir -p ${working_dir}clusterOutputFiles/cram2bam/

clusterOutput=${working_dir}clusterOutputFiles/


nrGenomes=$( wc -l ${download_fileList} | awk '{print $1}' )

## Only allows 20 jobs running at a time.
# I got spanked by John Constable for overloading irods.
bsub -J "download_and_index[1-${nrGenomes}]%15" -M100 -R'span[hosts=1] select[mem>100] rusage[mem=100]' \
-e ${clusterOutput}download_and_index/download_and_index_errors.%J.%I  -o ${clusterOutput}download_and_index/download_and_index_output.%J.%I  \
bash ${script_dir}01_download_crams_wrapper.sh ${script_dir}01_download_crams_worker.sh ${output_dir} \
${download_fileList} ${samtools_path} ${referenceFile} ${script_dir} ${clusterOutput}cram2bam/


exit $?
