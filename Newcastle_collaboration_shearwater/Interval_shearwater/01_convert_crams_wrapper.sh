#!/bin/bash

script=$1
output_dir=$2
cram_fileList=$3
samtools_path=$4
referenceFile=$5
genomes_dir=$6
script_dir=$7

cramFile=$( awk -v INDEX="${LSB_JOBINDEX}" 'NR==INDEX {print}' < "${cram_fileList}" )

${script} ${output_dir} ${cramFile} ${samtools_path} ${referenceFile} ${genomes_dir} ${script_dir} 

exit $?
