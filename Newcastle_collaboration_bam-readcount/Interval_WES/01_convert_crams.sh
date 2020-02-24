#!/bin/bash

# run script locally: ./01_convert_crams.sh

source variables.sh

mkdir -p ${output_dir}
mkdir -p ${clusterOutput}

nrGenomes=$(wc -l ${cram_fileList} | awk '{print $1}' )

# runnng as singe-threaded as it's more core effiecient
bsub -J "cram2bam[1-${nrGenomes}]%1000" -M600 -R'span[hosts=1] select[mem>600] rusage[mem=600]' \
-e ${clusterOutput}cram2bam_errors.%J.%I  -o ${clusterOutput}cram2bam_output.%J.%I  \
bash ${script_dir}01_convert_crams_wrapper.sh ${script_dir}01_convert_crams_worker.sh ${output_dir} \
${cram_fileList} ${samtools_path} ${referenceFile} ${genomes_dir} ${script_dir}

exit $?
