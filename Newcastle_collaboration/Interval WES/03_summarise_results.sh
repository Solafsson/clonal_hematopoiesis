#!/bin/bash
source variables.sh

#the operations below are done one at a time because piping did not work - was creating a 100+G file
#possibly sample names (supplier id) are to blame

cat ${output_dir}*_bqsr_snv_results.txt > "${output_dir}all_results_combined_1.txt"
awk -F":" '{print $1, $15, $28, $41, $54, $67}'  "${output_dir}all_results_combined_1.txt" > "${output_dir}all_results_combined_2.txt"
awk -F"=" '{print $1, $2}'  "${output_dir}all_results_combined_2.txt" > "${output_dir}all_results_combined.txt"

rm "${output_dir}all_results_combined_1.txt" "${output_dir}all_results_combined_2.txt"
exit $?
