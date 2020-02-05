#!/bin/bash
source variables.sh

cat ${output_dir}*txt | awk -F":" '{print $1, $15, $28, $41, $54, $67}' | \
awk -F"=" '{print $1, $2}' > "${output_dir}/all_results_combined.txt"

exit $?
