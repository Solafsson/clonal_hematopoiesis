#!/bin/bash
source variables.sh

cat ${output_dir}*_snv_results.txt | awk -F":" '{print $1, $15, $28, $41, $54, $67}' | awk -F"=" '{print $1, $2}' > ${output_dir}all_results_combined.txt

rm ${output_dir}*.bam
rm ${output_dir}*.bam.bai

exit $?
