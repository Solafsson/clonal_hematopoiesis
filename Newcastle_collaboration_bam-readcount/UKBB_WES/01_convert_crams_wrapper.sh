#!/bin/bash

script=$1
output_dir=$2
cram_fileList=$3
samtools_path=$4
referenceFile=$5
script_dir=$6

cramFile=$( awk -v INDEX="${LSB_JOBINDEX}" 'NR==INDEX {print}' < "${cram_fileList}" )

${script} ${output_dir} ${cramFile} ${samtools_path} ${referenceFile} ${script_dir} 

exit $?
