#young IBD vs young Interval

library(reshape2)
library(ggplot2)
#library(lubridate)
library(stringr)
library(dplyr)

setwd("/Users/ar31/Documents/early_work_for_grant_anas_copy/ana")
all_res_IBD <- read.csv("all_res_with_age_IBD.csv", row.names = 1)
all_res_Interval <- read.csv("all_res_with_age_Interval.csv", row.names = 1)

#inspect age range in the cohorts
summary(all_res_IBD$age_sampling)
summary(all_res_Interval$age_sampling)

#remove donors from the IBD cohort that can't be age-matched
all_res_IBD  = all_res_IBD[(all_res_IBD$age_sampling >=18) & (all_res_IBD$age_sampling <= 76),]
summary(all_res_IBD$age_sampling)

#load sitelist and remove introns
siteList <- read.csv("combined_list_of_SNVs.csv")
siteList = siteList[!grepl("c.e", siteList$optional_id),]

###
#age <35 on one hand and >55 on the other


####get top 1000 youngest in both datasets
IBD_y = all_res_IBD[all_res_IBD$age_sampling <=35,]
IBD_y$group = "Young_IBD"
IBD_y$Pos = paste0(IBD_y$chr,"_",IBD_y$pos)
summary(IBD_y$age_sampling)
length(unique(IBD_y$sample))

Interval_y = all_res_Interval[all_res_Interval$age_sampling <=35,]
Interval_y$group = "Young_Interval"
Interval_y$Pos = paste0(Interval_y$chr,"_",Interval_y$pos)
summary(Interval_y$age_sampling)
length(unique(Interval_y$sample))

all_res_comb = rbind(IBD_y, Interval_y)



###get format for pileup
#merge all_res with siteList to get type of mutation
all_res_comb = distinct( merge(all_res_comb, siteList[,c("chr","pos_in_hg38","CHIP_driver_or_control")], by.x = c("chr","pos"), by.y = c("chr","pos_in_hg38")))

#rearrange columns 
all_res_comb = all_res_comb[, c("sample","group","Pos","CHIP_driver_or_control","ref_allele","depth","A","C","G","T","N")]

#get only sites that were picked up in both gropus 
siteList$Pos = paste0(siteList$chr,"_",siteList$pos)
siteList_use = distinct(siteList[(siteList$Pos %in% intersect(IBD_y$Pos, Interval_y$Pos)),c("chr", "ref_allele","var_allele","pos_in_hg38","Pos")])


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

m <- melt(pileups, id="Site")
m$value <- as.numeric(m$value)
m2 <- distinct(merge(m, unique(all_res_comb[, c("Pos", "CHIP_driver_or_control")]), by.x="Site", by.y="Pos"))


#plot bar charts with all sites
ggplot(m2[m2$CHIP_driver_or_control == "CHIP_driver",], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site",y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


ggplot(m2[m2$CHIP_driver_or_control == "negative_control",], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site", y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


#### bar charts for top20 CHIP sites accoridng to Watson et al doi: http://dx.doi.org/10.1101/569566. 
top_chip = read.csv("Top20_CHIP_sites_Watson_et_al.csv", header = F)
top_chip_pos = merge(top_chip, siteList, by.x ="V1", by.y = "optional_id")

ggplot(m2[m2$Site %in% top_chip_pos$Pos,], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site",y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


### bar charts for top 25 CHIP sites from Peter Cambell
PC_sitelist = read.table("snv_siteList_hg38.txt")
PC_chip_pos = merge(PC_sitelist, siteList, by.x =c("V1","V2"), by.y = c("chr", "pos_in_hg38"))

ggplot(m2[m2$Site %in% PC_chip_pos$Pos,], aes(x=reorder(Site, value), y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) + labs( x= "Site",y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())

