library(tidyverse)
setwd("/Users/ar31/Documents/early_work_for_grant_anas_copy/ana")

all_res_IBD <- read.csv("all_res_with_age_IBD.csv", row.names = 1)
sample_age = distinct(all_res_IBD[,c(1,11)])

#slect the 12 samples below 18
young_as_normal = sample_age[sample_age$age_sampling <18,]
young_as_normal$cram = paste0(young_as_normal$sample,".cram")
young_as_normal$bam = paste0(young_as_normal$sample,".bam")
#select 50 random samples avove 60
random_old_as_case = sample_age[sample_age$age_sampling > 60,] %>% {.[sample(1:nrow(.),50,replace=F),]}
random_old_as_case$cram = paste0(random_old_as_case$sample,".cram")
random_old_as_case$bam = paste0(random_old_as_case$sample,".bam")

write(as.character(young_as_normal$cram), "./IDB_WES_shearwater/cram_list_normals.txt")
write(as.character(random_old_as_case$cram), "./IDB_WES_shearwater/cram_list_cases.txt")

write(as.character(young_as_normal$bam), "./IDB_WES_shearwater/bam_list_normals.txt")
write(as.character(random_old_as_case$bam), "./IDB_WES_shearwater/bam_list_cases.txt")
