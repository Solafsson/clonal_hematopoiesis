#!/bin/bash


results_dir=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/results/
sitesList=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/resources/snv_sites_Gene_Change.txt
phaseKey=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/resources/EGAN_phase_key.txt

while read chr bp bp2 hg37 annot; do
    cat ${results_dir}*txt | awk -v CHR="$chr" -v BP="$bp" '$2==CHR && $3==BP {print}' | \
    awk 'BEGIN {OFS=":"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' | \
    awk -v ANN="$annot" 'BEGIN {FS=":"} {print $1, $2"_"$3, $4, ANN, $5, $7, $21, $35, $49, $63, $77}' > ${results_dir}${chr}_${bp}_summary

    join -1 1 -2 1 ${phaseKey} <(sort -k1,1 ${results_dir}${chr}_${bp}_summary) | awk 'BEGIN {print "Sample", "Phase", "Pos", "REF","Annotation", "Depth", "=", "A", "C", "G", "T", "N"} {print}' > ${results_dir}${chr}_${bp}_summary_PhaseIncl
done < ${sitesList}


cat ${results_dir}*summary_PhaseIncl | awk 'NR==1 {print} $1!~/Sample/ {print}' > ~/phd/clonal_hematopoiesis/results/all_results_combined.txt



exit $?