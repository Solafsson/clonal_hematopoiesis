library(magrittr)
library(stringr)

# read phenotype data linked to anderson ID, remove withdrawn samples
withdrawn_samples = read.csv("/lustre/scratch114/projects/ukbb_team152/45669/withdrawn_samples_20200224/w45669_20200204.csv", header = F)
phenotype_anderson = read.table("/lustre/scratch114/projects/ukbb_team152/ar31/resources/ukbb_phenotype_20191030.tsv", sep = "\t", header =T) %>% .[!.$f.eid %in% withdrawn_samples$V1,]

#read file linking anderson and hurles IDs, include hurles ID in phenotype data
map_hurles = read.table("/lustre/scratch114/projects/ukbb_team152/ar31/resources/linked_fam_files.stripped.txt", header = T)
phenotype_merged = merge(phenotype_anderson, map_hurles, by.x = "f.eid", by.y = "anderson")

#read map linking UKBB ID and path to gVCF files, get first 7 didgits of file name, equivalent to hurles ID
#then merge file path and hurles id with phenotype
map_ukbb = read.table("/lustre/scratch114/projects/ukbb_team152/ar31/resources/ukbb_exomes_gvcf_map.txt")
colnames(map_ukbb) = c("UKBB_id", "gvcf_path")
map_ukbb$hurles_ID = str_sub(as.character(map_ukbb$gvcf_path) ,-26,-20)

#create map that includes all ID and file paths to WES data
map_merged = merge(map_ukbb, phenotype_merged , by.x = "hurles_ID", by.y = "hurles")

#replace gvcf path in map with path to crams
map_merged$cram_path = gsub("gVCF", "cram", map_merged$gvcf_path) %>% gsub(".g.vcf.gz" , ".cram", .)
map_merged$sample = str_sub(as.character(map_merged$gvcf_path) ,-26,-10)

#select samples for bam-readcount analysis: old as cases, younger as controls - around 4k each
summary(map_merged$age)
old_UKBB = map_merged[!is.na(map_merged$CTR),] %>% .[.$CTR == 1,]   %>% .[.$age > 65 ,] 
young_UKBB_as_control = map_merged[!is.na(map_merged$CTR),] %>% .[.$CTR == 1,]   %>% .[.$age < 45,] 

old_UKBB$group = "Oldest"
young_UKBB_as_control$group = "Youngest"

HC_phenotypes = map_merged[!is.na(map_merged$CTR),] %>% .[.$CTR == 1,] 

#write files for bam-readcount
write.table(map_merged,"/lustre/scratch114/projects/ukbb_team152/ar31/resources/complete_map_UKBB_WES.txt", row.names = F, quote = F)
write.table(old_UKBB$cram_path,"/lustre/scratch114/projects/ukbb_team152/ar31/resources/old_UKBB_cram_list.txt", row.names = F, col.names = F, quote = F)
write.table(young_UKBB_as_control$cram_path,"/lustre/scratch114/projects/ukbb_team152/ar31/resources/young_UKBB_cram_list.txt", row.names = F, col.names = F, quote = F)
write.table(HC_phenotypes[,c("sample", "age")], "/lustre/scratch114/projects/ukbb_team152/ar31/resources/HC_phenotypes.txt")

#write file with 10bp window around CHIP sites - used for cram2bam conversion
setwd("/lustre/scratch114/projects/ukbb_team152/ar31/resources")
CHIP_sites_UKBB = read.table("unaffected_CHIP_sites_UKBB_and_CHIP_no_chr.txt")
CHIP_sites_UKBB$start = CHIP_sites_UKBB$V2 - 5
CHIP_sites_UKBB$end = CHIP_sites_UKBB$V2 + 5
write.table(CHIP_sites_UKBB[,c(1,4,5)], "CHIP_sites_UKBB_for_conversion.txt", row.names = F, col.names = F, quote = F)

#write files for shearwater
#select samples for shearwater analysis: old as cases, younger as controls - 10k old, 500 young
summary(map_merged$age)
old_UKBB_shearwater = map_merged[!is.na(map_merged$CTR),] %>% .[.$CTR == 1,]   %>% .[.$age > 62 ,] 
young_UKBB_as_control_shearwater = map_merged[!is.na(map_merged$CTR),] %>% .[.$CTR == 1,]   %>% .[.$age < 41,] 

write.table(rbind(old_UKBB_shearwater,young_UKBB_as_control_shearwater)$cram_path,"/lustre/scratch114/projects/ukbb_team152/ar31/resources/UKBB_cram_list_shearwater.txt", row.names = F, col.names = F, quote = F)

#HC and case bam list for shearwater -first remove bams that failed due to lack of indexing

old_UKBB_shearwater_rm = old_UKBB_shearwater[!old_UKBB_shearwater$cram_path %in% c(rbind(old_UKBB_shearwater,young_UKBB_as_control_shearwater)$cram_path[1326], rbind(old_UKBB_shearwater,young_UKBB_as_control_shearwater)$cram_path[4815]),]
old_UKBB_shearwater_rm$bam_path = paste0("/lustre/scratch114/projects/ukbb_team152/ar31/bams/",old_UKBB_shearwater_rm$sample,".bam")
  
young_UKBB_as_control_shearwater_rm = young_UKBB_as_control_shearwater[!young_UKBB_as_control_shearwater$cram_path %in% c(rbind(young_UKBB_as_control_shearwater,young_UKBB_as_control_shearwater)$cram_path[1326], rbind(young_UKBB_as_control_shearwater,young_UKBB_as_control_shearwater)$cram_path[4815]),]
young_UKBB_as_control_shearwater_rm$bam_path = paste0("/lustre/scratch114/projects/ukbb_team152/ar31/bams/",young_UKBB_as_control_shearwater_rm$sample,".bam")

write(as.character(old_UKBB_shearwater_rm$bam_path), "/lustre/scratch114/projects/ukbb_team152/ar31/resources/bam_list_cases.txt")
write(as.character(young_UKBB_as_control_shearwater_rm$bam_path), "/lustre/scratch114/projects/ukbb_team152/ar31/resources/bam_list_normals.txt")
