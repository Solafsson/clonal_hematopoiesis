library(reshape2)
library(ggplot2)
#library(lubridate)
library(stringr)
library(dplyr)

setwd("/Users/ar31/Documents/early_work_for_grant_anas_copy/ana")
all_res <- read.table("./Interval_WES/all_results_combined.txt", header = F)
colnames(all_res) = c("sample", "chr", "pos", "ref_allele", "depth", "A", "C", "G", "T", "N")
#remove rows with zero coverage
all_res = all_res[all_res$depth>0,]

phenotype = read.table("./Interval_WES/INTERVAL_WES_phenotypes_basic.txt", header = T)

map = read.table("./Interval_WES/interval_wes_IRODS_EGA.txt", header = F)
map$V6 = gsub(".cram","_bqsr", map$V6)

siteList <- read.csv("combined_list_of_SNVs.csv")

#remove introns form siteList
siteList = siteList[!grepl("c.e", siteList$optional_id),]

#merge phenotype and map info
merged_phenotype = merge(as.data.frame(map[,c(1,6)]), phenotype, by.x = c("V1"), by.y =c("ID"))
colnames(merged_phenotype)[1] = c("accession_number")
colnames(merged_phenotype)[2] = c("supplier_id")


#keep only samples present in all_res
merged_phenotype = merged_phenotype[merged_phenotype$supplier_id %in% all_res$sample,]
#remove counts for samples with unknown phenotype
all_res = all_res[all_res$sample %in% merged_phenotype$supplier_id,]


#merge all_res with phenotype, change sample names of all_res from supplier_id to accession_number
all_res = distinct(merge(merged_phenotype[,c("accession_number","supplier_id")], all_res, by.x = "supplier_id", by.y = "sample"))
all_res = all_res[,2:ncol(all_res)]

colnames(all_res)[1] = c("sample")

#round age 
merged_phenotype$age_sampling = round(as.numeric(as.vector(merged_phenotype$age)), digits = 0)

####get top 1000 youngest and top 1000 oldest
merged_ord_y = merged_phenotype[order(merged_phenotype$age_sampling),][1:1000,]
all_res_y = all_res[all_res$sample %in% merged_ord_y$accession_number,]
all_res_y$group = "Youngest"
all_res_y$Pos = paste0(all_res_y$chr,"_",all_res_y$pos)

merged_ord_o = merged_phenotype[order(-merged_phenotype$age_sampling),][1:1000,]
all_res_o = all_res[all_res$sample %in% merged_ord_o$accession_number,]
all_res_o$group = "Oldest"
all_res_o$Pos = paste0(all_res_o$chr,"_",all_res_o$pos)

all_res_comb = rbind(all_res_y, all_res_o)

###get format for pileup
#merge all_res with siteList to get type of mutation
all_res_comb = distinct( merge(all_res_comb, siteList[,c("chr","pos_in_hg38","CHIP_driver_or_control")], by.x = c("chr","pos"), by.y = c("chr","pos_in_hg38")))
all_res_comb$Pos = paste0(all_res_comb$chr,"_",all_res_comb$pos)

#rearrange columns 
all_res_comb = all_res_comb[, c("sample","group","Pos","CHIP_driver_or_control","ref_allele","depth","A","C","G","T","N")]

#get only sites that were picked up in both gropus 
siteList$Pos = paste0(siteList$chr,"_",siteList$pos)
siteList_use = distinct(siteList[(siteList$Pos %in% intersect(all_res_o$Pos, all_res_y$Pos)),c("chr", "ref_allele","var_allele","pos_in_hg38","Pos")])

### get pileup
CC_pileup <- function(res, site="chr1_114713908") {
  tmp <- subset(res, Pos==site)
  
  if(tmp$ref_allele[1]=="A") {
    ref_col_nr <- 7
    alt_col_nrs <- c(8:11)
  } else if(tmp$ref_allele[1]=="C") {
    ref_col_nr <- 8
    alt_col_nrs <- c(7,9:11)
  } else if(tmp$ref_allele[1]=="G") {
    ref_col_nr <- 9
    alt_col_nrs <- c(7,8,10,11)
  } else if(tmp$ref_allele[1]=="T") {
    ref_col_nr <- 10
    alt_col_nrs <- c(7:9,11)
  }
  
  alt_reads_total <- tapply(tmp[, alt_col_nrs[1]], tmp$group, sum) + tapply(tmp[, alt_col_nrs[2]], tmp$group, sum) + tapply(tmp[, alt_col_nrs[3]], tmp$group, sum) + tapply(tmp[, alt_col_nrs[4]], tmp$group, sum) 
  total_reads <- tapply(tmp$depth, tmp$group, sum)
  pileup <- t(data.frame(c(as.numeric(alt_reads_total/total_reads), site)))
  return(pileup)
}

pileups <- data.frame()
for(j in 1:nrow(siteList_use[,c("chr","pos_in_hg38")])) {
  pileups <- rbind(pileups, CC_pileup(res=all_res_comb, paste(siteList_use$chr[j], siteList_use$pos_in_hg38[j], sep="_")))
}

colnames(pileups) <- c(levels(as.factor(all_res_comb$group)), "Site")
pileups = pileups[order(pileups$Oldest),]
m <- melt(pileups, id="Site")
m$value <- as.numeric(m$value)
m2 <- distinct(merge(m, unique(all_res_comb[, c("Pos", "CHIP_driver_or_control")]), by.x="Site", by.y="Pos"))

  
#### plot bar charts for all genes
ggplot(m2[m2$CHIP_driver_or_control == "CHIP_driver",], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site",y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


ggplot(m2[m2$CHIP_driver_or_control == "negative_control",], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site",y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


#### bar charts for top20 CHIP sites accoridng to Watson et al doi: http://dx.doi.org/10.1101/569566. 
top_chip = read.csv("Top20_CHIP_sites_Watson_et_al.csv", header = F)
top_chip_pos = merge(top_chip, siteList, by.x ="V1", by.y = "optional_id")

ggplot(m2[m2$Site %in% top_chip_pos$Pos,], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site",y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


## bar charts for top 25 CHIP sites from Peter Cambell
PC_sitelist = read.table("snv_siteList_hg38.txt")
PC_chip_pos = merge(PC_sitelist, siteList, by.x =c("V1","V2"), by.y = c("chr", "pos_in_hg38"))

ggplot(m2[m2$Site %in% PC_chip_pos$Pos,], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site",y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


##### write all_res to file along with age
all_res_with_age = merge(all_res, merged_phenotype[,c("accession_number", "age_sampling")], by.x = "sample", by.y = "accession_number")
write.csv(all_res_with_age, "all_res_with_age_interval.csv")
