#!/bin/bash

# the script looks up the 395 CHIP sites from literature and the surronding 10 bp (location +/- 5bp) in gVCF files
# the input bed file contains 3 columns: chr, start_pos, end_pos - end pos must be diffent than start
# two files are produced containing:
#1. all identified sites that have counts for both a ref and an alt allele;
#2. all identified sites with a homozygous ref genotype

while IFS= read -r gvcf; do
/software/hgi/installs/anaconda3/envs/team152/bin/bcftools query -f'%CHROM\t%POS\t%REF\t%ALT[%GT\t%AD\t%DP]\n' /lustre/scratch118/humgen/hgi/projects/interval_wes/crams/$gvcf \
-R /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater_Interval/resources/CHIP_loci_shearwater.bed \
| grep "," | uniq \
| awk -v gVCF="$gvcf" '{print $0, "\t", gVCF}'\
>> /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/gVCF_pileup_Interval/results/mut_in_CHIP_sites.txt
done < /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/gVCF_pileup_Interval/resources/gVCF_files.txt

grep '0/0\|0|0' /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/gVCF_pileup_Interval/results/mut_in_CHIP_sites.txt \
> /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/gVCF_pileup_Interval/results/somatic_mut_in_CHIP_sites.txt
