
library(reshape2)
library(ggplot2)

siteList <- read.table("/Users/so11/phd/so11_nfs/clonal_hematopoiesis/results/snv_siteList_hg38.txt")

all_res <- read.table("/Users/so11/phd/so11_nfs/clonal_hematopoiesis/results/all_results_combined.txt", h=T)
all_res$Case_Ctrl <- "CTRL"
all_res$Case_Ctrl[grep("IBD", all_res$Phase)] <- "Case"


CC_pileup <- function(res, site="chr1_114713908") {
  tmp <- subset(res, Pos==site)
  
  if(tmp$REF[1]=="A") {
    ref_col_nr <- 8
    alt_col_nrs <- c(7,12,9,10,11)
  } else if(tmp$REF[1]=="C") {
    ref_col_nr <- 9
    alt_col_nrs <- c(8,7,12,10,11)
  } else if(tmp$REF[1]=="G") {
    ref_col_nr <- 10
    alt_col_nrs <- c(7,9,8,12,11)
  } else if(tmp$REF[1]=="T") {
    ref_col_nr <- 11
    alt_col_nrs <- c(10,7,8,9,12)
  }
  
  alt_reads_total <- tapply(tmp[, alt_col_nrs[1]], tmp$Case_Ctrl, sum) + tapply(tmp[, alt_col_nrs[2]], tmp$Case_Ctrl, sum) + tapply(tmp[, alt_col_nrs[3]], tmp$Case_Ctrl, sum) + tapply(tmp[, alt_col_nrs[4]], tmp$Case_Ctrl, sum) + tapply(tmp[, alt_col_nrs[5]], tmp$Case_Ctrl, sum)
  total_reads <- tapply(tmp$Depth, tmp$Case_Ctrl, sum)
  pileup <- t(data.frame(c(as.numeric(alt_reads_total/total_reads), site)))
  return(pileup)
}

pileups <- data.frame()
for(j in 1:nrow(siteList)) {
  pileups <- rbind(pileups, CC_pileup(res=all_res, paste(siteList$V1[j], siteList$V2[j], sep="_")))
}
colnames(pileups) <- c("Case", "CTRL", "Site")
m <- melt(pileups, id="Site")
m$value <- as.numeric(m$value)
m2 <- merge(m, unique(all_res[, c("Pos", "Annotation")]), by.x="Site", by.y="Pos")
m2$Site <- paste(m2$Site, m2$Annotation, sep="_")

ggplot(m2, aes(x=Site, y=value)) + geom_bar(aes(fill=variable), position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle=90)) + labs(y="Fraction of reads supporting \nany alternative allele") +
  theme(legend.title = element_blank())


## Do a sign-test to test for differences:
## See http://rcompanion.org/handbook/F_07.html

library(BSDA)
cases <- m2$value[m2$variable=="Case"]
ctrls <- m2$value[m2$variable=="CTRL"]

SIGN.test(x=cases, y=ctrls)







