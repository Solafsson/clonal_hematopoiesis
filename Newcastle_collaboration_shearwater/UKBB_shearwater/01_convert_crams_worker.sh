#!/bin/bash

output_dir=$1
cramFile=$2
samtools_path=$3
referenceFile=$4
script_dir=$5
convert_siteList=$6

# removing .cram extension for bam file name
sampleName=$( echo ${cramFile} | awk 'BEGIN {FS="/"} {print substr($NF, 1, length($NF)-5)}' )

echo "Converting cram to bam"

## converting cram to bam
# options:
# -b Output in the BAM format
# -M removes duplicate sequences
# -T FASTA format reference FILE
# -L list of sites

${samtools_path} view -b -M -T ${referenceFile} ${cramFile} -L ${convert_siteList} > ${output_dir}${sampleName}.bam

#indexing BAM
${samtools_path} index ${output_dir}${sampleName}.bam

exit $?
