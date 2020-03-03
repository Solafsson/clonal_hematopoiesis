#!/bin/bash

cat /lustre/scratch119/humgen/teams/anderson/users/so11/SVs_in_ibd/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta | sed 's/>chr/>/g' > /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/UKBB_control_WES/resources/Homo_sapiens_assembly38_chr.fasta
