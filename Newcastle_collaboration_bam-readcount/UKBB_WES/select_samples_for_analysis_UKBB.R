library(dplyr)
phenotype = read.table("/lustre/scratch114/projects/ukbb_team152/ar31/ukbb_phenotype_20191030.tsv", sep = "\t", header =T)
withdrawn_samples = read.csv("/lustre/scratch114/projects/ukbb_team152/45669/withdrawn_samples_20200224/w45669_20200204.csv", header = F)
phenotype_old_ctrl = phenotype[phenotype$CTR == 1,] %>% .[.$age > 50 ,] %>% .[!.$f.eid %in% withdrawn_samples$V1,]
system( )


write.table(phenotype_old_ctrl$f.eid, "/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/UKBB_control_WES/resources/old_ctrl_UKBB.txt")