#!/bin/bash

source variables.sh

referenceFile=$1
bamFile=$2
sampleName=$3

## bam-redcount options:
# -b - minimum base quality at a position to use the read for counting
# -q - minimum mapping quality of reads used for counting
# -w - maximum number of warnings of each type to emit (set to avoid getting repeated "Couldn't find single-end mapping quality" WARNINGS)

${bam_readcount_path} -b 30 -q 30 -w 0 -f ${referenceFile} -l ${siteList} ${bamFile} | \
awk -v SAM="$sampleName" '{print SAM, $0}' > ${output_dir}${sampleName}_snv_results.txt


## rm bam file
#rm ${bamFile}
#rm ${bamFile}.bai

exit $?
