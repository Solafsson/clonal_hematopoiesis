#!/bin/bash


sampleList=$1

genomes_dir=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/genomes/
script_dir=/nfs/users/nfs_s/so11/phd/clonal_hematopoiesis/
samtools_path=/software/hgi/pkglocal/samtools-1.3.1/bin/samtools
## Reference genome (hg38 for ibd data)
referenceFile=/lustre/scratch119/humgen/teams/anderson/users/so11/SVs_in_ibd/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta
clusterOutput=/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/clusterOutputFiles/cram2bam/

while read sampleName; do
    ${script_dir}03_count_reads.sh ${referenceFile} ${genomes_dir}${sampleName}.bam ${sampleName}
done < ${sampleList}