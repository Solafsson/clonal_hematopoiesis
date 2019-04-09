#!/bin/bash

script=$1
bam_dir=$2
jobFile=$3
samtools_path=$4
reference_genome=$5
script_dir=$6
clusterOutput=$7



${script} ${bam_dir} ${jobFile} ${samtools_path} ${script_dir} ${clusterOutput} ${reference_genome}  ${LSB_JOBINDEX}

exit $?

script_dir=$5
clusterOutput=$6
referenceFile=$7