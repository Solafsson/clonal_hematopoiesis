#!/bin/bash

## https://github.com/genome/bam-readcount

referenceFile=$1
bamFile=$2
sampleName=$3

siteList=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/resources/snv_siteList_hg38.txt
output_dir=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/results/


/software/hgi/pkglocal/bam-readcount-0.8.0/bin/bam-readcount -f ${referenceFile} -l ${siteList} ${bamFile} | \
awk -v SAM="$sampleName" '$1~/chr.*/ {print SAM, $0}' > ${output_dir}${sampleName}_snv_results.txt

## rm bam file
rm ${bamFile}
rm ${bamFile}.bai

exit $?