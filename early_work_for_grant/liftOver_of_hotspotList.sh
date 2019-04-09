#!/bin/bash

## This script will not run multiple times. It's just a documentation of how the liftOver was performed.

## I got sent a list of myeloid hotspots from Peter Campbell which I assume are in build hg37 (like everything in CASM).
## I want to convert this into hg38 since all our data is in hg38.

liftOver_exec=/software/pathogen/external/apps/usr/bin/liftOver
chain_file=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/resources/hg19ToHg38.over.chain
hotSpot_list=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/resources/myeloid_hotspots_all_fromPeter_hg37.txt

output_dir=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/resources/

## Bed format assumes 0-based indexing but Peter has used 1-based indexing.
## For now I want only snvs, thus $2==$3
## Some position appear multiple times with different alternative allele. Just look up once
## $4 is the reference allele. Write that out for sanity checking.
awk 'NR>1 && $2==$3 {print "chr"$1, $2-1, $3, $4}' < ${hotSpot_list} | sort -u > ${output_dir}snv_hotspot_pos_hg37.bed

${liftOver_exec} ${output_dir}snv_hotspot_pos_hg37.bed ${chain_file} ${output_dir}snv_hotspot_pos_lifted_hg38.bed ${output_dir}snv_hotspot_pos_unlifted_hg38.bed

awk '{print $1, $2+1, $3}' < ${output_dir}snv_hotspot_pos_lifted_hg38.bed > ${output_dir}snv_siteList_hg38.txt

exit $?



