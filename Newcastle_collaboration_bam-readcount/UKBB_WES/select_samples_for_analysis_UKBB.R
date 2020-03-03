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

#select sample for analysis: old as cases, younger as controls
summary(map_merged$age)
old_UKBB = map_merged[!is.na(map_merged$CTR),] %>% .[.$CTR == 1,]   %>% .[.$age > 55 ,]
young_UKBB_as_control = map_merged[!is.na(map_merged$CTR),] %>% .[.$CTR == 1,]   %>% .[.$age < 43,]

length(unique(old_UKBB$hurles_ID))
length(unique(young_UKBB_as_control$hurles_ID))

#write to results to file
write.table(map_merged,"/lustre/scratch114/projects/ukbb_team152/ar31/resources/complete_map_UKBB_WES.txt", row.names = F, quote = F)
write.table(old_UKBB$cram_path,"/lustre/scratch114/projects/ukbb_team152/ar31/resources/old_UKBB_cram_list.txt", row.names = F, col.names = F, quote = F)
write.table(young_UKBB_as_control$cram_path,"/lustre/scratch114/projects/ukbb_team152/ar31/resources/young_UKBB_cram_list.txt", row.names = F, col.names = F, quote = F)

