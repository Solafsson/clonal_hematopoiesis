#!/bin/bash

output_dir=$1
cramFile=$2
samtools_path=$3
referenceFile=$4
genomes_dir=$5
script_dir=$6

# removing .cram extension for bam file name
sampleName=$( echo ${cramFile} | awk 'BEGIN {FS="/"} {print substr($NF, 1, length($NF)-5)}' )

echo "Converting cram to bam"

## converting cram to bam
# options:
# -@ Number of BAM compression threads to use in addition to main thread
# -b Output in the BAM format
# -T A FASTA format reference FILE
${samtools_path} view -b -T ${referenceFile} ${genomes_dir}${cramFile} > ${output_dir}${sampleName}.bam

#indexing BAM
${samtools_path} index ${output_dir}${sampleName}.bam

## calling count_reads script here so that runs in the same job as the cram2bam conversion
echo "Counting the reads"
echo "${script_dir}02_count_reads.sh ${referenceFile} ${output_dir}${sampleName}.bam ${sampleName}"

${script_dir}02_count_reads.sh ${referenceFile} ${output_dir}${sampleName}.bam ${sampleName}

exit $?
