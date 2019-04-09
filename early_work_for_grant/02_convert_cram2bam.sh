#!/bin/bash

genomes_dir=$1
cramFile=$2
sammtools_path=$3
referenceFile=$4
script_dir=$5

sampleName=$( echo ${cramFile} | awk 'BEGIN {FS="/"} {print substr($NF, 1, length($NF)-5)}' )

echo "Converting cram to bam"
/software/hgi/pkglocal/samtools-1.3.1/bin/samtools view -@ 8 -b -T ${referenceFile} ${cramFile} > ${genomes_dir}${sampleName}.bam
/software/hgi/pkglocal/samtools-1.3.1/bin/samtools index ${genomes_dir}${sampleName}.bam

## Remove the cram file to save space.
echo "Removing the cramFile: ${cramFile}"
rm -f ${cramFile}


## I'm going to call the count_reads script here so that runs in the same job as the cram2bam conversion
echo "Counting the reads"
echo "${script_dir}03_count_reads.sh ${referenceFile} ${genomes_dir}${sampleName}.bam ${sampleName}"
${script_dir}03_count_reads.sh ${referenceFile} ${genomes_dir}${sampleName}.bam ${sampleName}


exit $?