#!/bin/bash

# script for checking overlap between the CHIP-associated sites and incorrectly mapped UKBB data
# http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1911

working_dir="/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/UKBB_control_WES/resources/"
bedtools_dir="/software/hgi/installs/anaconda3/envs/team152/bin/bedtools"

# -v flag in bedtools -  only reports those entries in A that have no overlap in B
${bedtools_dir} intersect -v -a "${working_dir}SNV_site_list.bed"  -b "${working_dir}xgen_plus_spikein_b38_alt_affected.bed" |uniq > "${working_dir}unaffected_CHIP_sites_UKBB_and_CHIP.bed"

exit $?

