# RUN like this:
# cd /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater/scripts
# nohup /software/R-3.6.1/bin/Rscript annotate_mutations.R >& annotate_mutations.log


# Inigo Martincorena - Feb 2017
# Annotating mutations from the output of esopipe_perpatient.R
#
# fa8: I've done several modifications:
#       1. After FDR correction I do not remove all non-significant p-values, I set q<0.1
#          This is done to allow for the posterior rescue of variants: 
#             when a variant is reliable in a sample we can relax the threshold in another sample
#          In this version, variants that have q>=0.01 and q <0.1 are not removed but labelled as no-fdr
#          This rescuing worked very well for some samples but not for others.
#          All the variants "rescued" are labelled as OK-rescued.
#       2. I have included a piece of code to estimate contamination based on alternative homozygotes (AA)
#          Any read matching the reference is considered a contaminant (not necessarily true but works)
#          The code also prints the list of AA sites and read counts, allowing to identify sites that are
#          not good contamination estimators (i.e. sites that show a higher "contamination" proportion for
#          all the samples)
#       3. As with FDR, variants filtered based on indels, strandness, etc are not filtered now, they
#          are just labelled accordingly. In the end, we will only be interested in OK variants or OK + OK-rescued
#          (depending on whether we want to include the rescued variants or not)
#       4. There is a initial step, before FDR and any filtering in which consecutive indels (- or INS, 
#          but not mixed) are merged together if they have consistent VAFs (consistent meaning not considered
#          different in a Fisher's exact test). Later, if one of the indels is considered reliable (it passed
#          all the filters), all the associated indels will be included in the output.
#          ** Possible improvement: instead of comparing the VAF of indels, we could check that the indels
#             happen in the same reads
#       5. There are some variables that have to be hard coded at present (not provided as arguments), like the 
#          sample_table and the path to the bam files.
        
# Input arguments:
# 1. Dataset name (e.g. PD prefix for the patient)
#
# Example:
# nohup Rscript annotate_mutations.R dataset_name [bait_regions.bed] [bam_path] &



####################################################################################################
## 1. Environment
####args = commandArgs(TRUE)
####dataset_name = args[1]; #"18003";
####dataset_name = gsub(" ","",dataset_name);
####if (length(args)>1) { 
####	baits_bed = args[2];
####} else { 
####	baits_bed = "/lustre/scratch116/casm/cgp/users/fa8/Katerina/v2_panel_new.bed";
####}
####outdir = "Mutation_calls"; system(sprintf("mkdir %s",outdir));
#####bams_path = "/lustre/scratch116/casm/cgp/users/fa8/BLADDER_SHEARWATER/BAMS_REALIGNED/"; #change accordingly
#####bams_path = "/lustre/scratch112/sanger/casm/cgp/im3/Projects/Oesophagus/20170814/BAM_files/"
####bams_path = "/lustre/scratch116/casm/cgp/users/fa8/Katerina/BAM_files/"; #link files in 1 directory 
#####sample_table = read.table(sprintf("samples_%s.txt",dataset_name), header=1, sep="\t", stringsAsFactors=F) # changed by fede to:
####sample_table = read.table("/lustre/scratch116/casm/cgp/users/fa8/Katerina/1484_samples.txt", header=1, sep="\t", stringsAsFactors=F) #list of sample ids
#####sample_table = read.table("/lustre/scratch116/casm/cgp/users/fa8/BLADDER_SHEARWATER/1385.list", header=1, sep="\t", stringsAsFactors=F)
#####The include_REF_REF variable can be modified below to change the way contamination is estimated
#####sample_table = read.table("/lustre/scratch112/sanger/casm/cgp/im3/Projects/Oesophagus/20170814/all_samples.txt", header=0, sep="\t", stringsAsFactors=F)
#####colnames(sample_table) = c("pid","sampleID")
#####sample_table = sample_table[substr(sample_table$sampleID,1,7)==dataset_name,]


# Libraries
library("GenomicRanges")
library("rtracklayer")
library("deepSNV")
library("data.table")
library("ggplot2")
library("stringr")
#install.packages("dplyr", repos = "http://cran.r-project.org")
library("dplyr") #added by ar31
setwd("/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater/results")

#MODIFY THE FOLLOWING DEPENDING ON THE EXPERIMENT!
dataset_name <- "CHIP_050220"
baits_bed = "/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater/resources/GENCODE_v33_genes_100L.bed"; #same bed used to run shearwater (used to calculate number of sites for  multiple hypothesis correction)
outdir = "Mutation_calls"; system(sprintf("mkdir %s",outdir));
bams_path = "/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater/results/bams"; #need linked folder containing all bam files
sample_table = read.table("/lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater/resources/bam_list_cases.txt", header=F, sep="\t", stringsAsFactors=F)
if(colnames(sample_table)[1] != "sampleID") {
	colnames(sample_table) <- c("sampleID")
}

####################################################################################################
## 2. Calling mutations from the Shearwater output files

baits = read.table(baits_bed, header=0, sep="\t", stringsAsFactors=F)
baits = GRanges(baits[,1], IRanges(baits[,2],baits[,3]))
numsegments_per_job = 20
entry_start = seq(from=1, to=length(baits), by=numsegments_per_job)
entry_end = pmin(entry_start+numsegments_per_job-1, length(baits))

####################################################################################################
# a. Loading the table of putative mutations from each patient
mutations = NULL
for (h in 1:length(entry_start)) {
	if(file.exists(sprintf("CHIP_050220/shearwater_temp_%s/mismatches_%s_%s.txt", dataset_name, entry_start[h], entry_end[h]))) {
	    m = read.table(file=sprintf("CHIP_050220/shearwater_temp_%s/mismatches_%s_%s.txt", dataset_name, entry_start[h], entry_end[h]), header=1, sep="\t", stringsAsFactors=F)
	    mutations = rbind(mutations,m)
	}
}

mutations$ref[which(mutations$ref == "TRUE")] <- "T"
mutations$mut[which(mutations$mut == "TRUE")] <- "T"

indels_f <- length(grep("[-I]",mutations[,"mut"]));
subs_f   <- nrow(mutations)-indels_f;
cat("#INITIAL_NUMBER_OF_MUTATIONS\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
mutations$mut_site <- paste(mutations$chr,mutations$pos,mutations$ref,mutations$mut);

L = sum(end(reduce(baits))-start(reduce(baits))+1) # Bait footprint
s = unique(sample_table$sampleID) # List of samples

mutations$vaf <- (mutations$xfw + mutations$xbw) / (mutations$xfw + mutations$xbw + mutations$nfw + mutations$nbw)

#added by ar31 to remove duplicate rows in "mutations" data frame
mutations = distinct(mutations)

####################################################################################################
# fa8: Let's merge neighbor indels and check they have consistent VAFs
rownames(mutations) <- paste(mutations$sampleID,mutations$chr,mutations$pos,mutations$mut,sep="|")
indels         <- mutations[which(mutations$mut=="-"),];
indels         <- indels[order(indels$sampleID, indels$chr, indels$pos),];
i              <- 1;
group_counter  <- 1;
while(i <= nrow(indels)) {
	indel_from  <- indels[i,"pos"];
	indel_to    <- indel_from;
	from_index  <- i;
	to_index    <- i;
	deleted_seq <- indels[i,"ref"];
	sample_f    <- indels[i,"sampleID"];
	if(i>nrow(indels)) {
		break;
	}	
	i <- i+1;
	for(j in c(i:nrow(indels))) {
		if(j>nrow(indels)) {
			i <- j 
			break; #finish
		}
		else if(indels[j,"pos"]>(indels[(j-1),"pos"]+1)) {
			i <- j; #start a new indel
			break;
		} else {
			if(indels[j,"chr"] != indels[(j-1),"chr"]) {
				i <- j; #start a new indel
				break;
			} else {
				#This is a candidate, but check their VAFs are compatible with a Fishers exact test:
				mat <- matrix(nrow=2,ncol=2,0);
				#mat[1,] <- c(indels[max(i,j-1),  "xfw"]+indels[max(i,j-1),  "xbw"], indels[max(i,j-1),  "nfw"]+indels[max(i,j-1),  "nbw"])
				mat[1,] <- c(indels[j-1,  "xfw"]+indels[j-1,  "xbw"], indels[j-1,  "nfw"]+indels[j-1,  "nbw"])
				mat[2,] <- c(indels[j,    "xfw"]+indels[j,    "xbw"], indels[j,    "nfw"]+indels[j,    "nbw"])
				mat[1,2] <- mat[1,2]-mat[1,1];
				mat[2,2] <- mat[2,2]-mat[2,1];
				pvalue <- fisher.test(mat)$p.value;
				if(pvalue < 0.01) {
					cat(" Breaking up indel because VAFs do not match\n");
					cat("             ",j-1, " vs ", j, ": pval=", pvalue, " [",indels[j-1,"pos"],"-",indels[j,"pos"],"]",sep="");
					cat("    (mat=", mat[1,1],",",mat[1,2],",",mat[2,1],",",mat[2,2],")\n",sep="");
					i <- j;
					break;
				}		
				indel_to <- indels[j,"pos"];
				to_index <- j;
				i <- j;
				deleted_seq <- paste(deleted_seq,indels[j,"ref"],sep="");
			}
		}
	}
	#cat("   [Sample=",sample_f,"] Indel goes from=",indel_from,", to=", indel_to," [",deleted_seq,">-]\n",sep="");
	indels[c(from_index:to_index),"groupID"    ] <- group_counter;
	indels[c(from_index:to_index),"deleted_seq"] <- deleted_seq;
	group_counter <- group_counter + 1;
}
mutations$indel_group <- NA;
mutations$deleted_seq <- NA;

mutations[rownames(indels),"indel_group"] <- indels$groupID
mutations[rownames(indels),"deleted_seq"] <- indels$deleted_seq


####################################################################################################
# b. Identifying putative germline or somatic indels to flag mutations near them
putative_indelsites = mutations[mutations$mut %in% c("-","INS"),]
s = unique(sample_table$sampleID) # List of samples 

putative_indelsites$qval = p.adjust(putative_indelsites$pval, method="BH", n=L*length(s)*2)
putative_indelsites = unique(putative_indelsites[putative_indelsites$qval<0.20, c("sampleID","chr","pos")])
indel_flank = 10
putative_indelsites_gr = GRanges(putative_indelsites$chr, IRanges(putative_indelsites$pos-indel_flank, putative_indelsites$pos+indel_flank))



####################################################################################################
# d. FDR calculation: significant mutations (after removing SNPs, to avoid inflating the FDR adjustment) 
mutations$label <- "";

L = sum(end(reduce(baits))-start(reduce(baits))+1) # Bait footprint
s = unique(sample_table$sampleID) # List of samples 

mutations$qval = p.adjust(mutations$pval, method="BH", n=L*length(s)*5)
mutations[which(mutations$qval>=0.01),"label"] = "no-fdr;";
prefdr.mutations <- mutations;                        # This will be the matrix used for the rescuing
mutations <- mutations[which(mutations$qval < 0.1),]; # To make the matrix smaller 
mutations = mutations[order(mutations$chr,mutations$pos),]
mutations$vaf = (mutations$xfw+mutations$xbw)/(mutations$nfw+mutations$nbw)
indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
cat("#AFTER_FDR\t",nrow(mutations[which(mutations$label == ""),]),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");


####################################################################################################
# f. Requesting at least 1 supporting read from both strands and annotating substitutions near indels (somatic or germline)

rmpos = (mutations$xfw < 2 | mutations$xbw < 2) # Asking for 2 supporting read in both strands for all mutations
filt2 = mutations[rmpos,]; 
if(nrow(filt2)>0) {
	filt2$filter = "Strandness"
}
mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"strandness;",sep="");
#####mutations = mutations[!rmpos,]
samples = unique(mutations$sampleID)
rmpos = NULL
for (h in 1:length(samples)) {
    m = mutations[mutations$sampleID==samples[h] & !(mutations$mut %in% c("-","INS")),]
    m_gr = GRanges(m$chr, IRanges(m$pos,m$pos))
    i_gr = putative_indelsites_gr[putative_indelsites$sampleID==samples[h]]
    ol = as.matrix(suppressWarnings(findOverlaps(m_gr, i_gr, type="any", select="all")))
    rmpos = c(rmpos, rownames(m)[unique(ol[,1])])
}
filt3 = mutations[rmpos,]; 
if(nrow(filt3) > 0) {
	filt3$filter = "Near_indel"
}
mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"near_indel;",sep="");
###########mutations = mutations[!(rownames(mutations) %in% rmpos),]
write.table(mutations, file=sprintf("%s/mutations_q01_%s.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
#write.table(rbind(filt1,filt2,filt3), file=sprintf("%s/filteredout_q01_%s.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
cat("#AFTER_BOTH_STRANDS_FILT\t",nrow(mutations[which(mutations$label == ""),]),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");


####################################################################################################
# fa8: Write table before merging consecutive subs/indels 
#      Interesting to save this because after the merging all error labels will be lost (only OK and OK-rescued 
#      will pass)
mutations[which(mutations$label == ""),"label"] <- "OK;";
write.table(mutations, file=sprintf("%s/mutations_including_failed_ones_%s.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)


####################################################################################################
# g. Annotating possible dinucleotides or runs of changes [WARNING: This annotation assumes that consecutive changes in a sample belong to a complex event. This may not be the case in 100% of the cases]
# fa8: We could compare VAFs and make sure they are consistent (Fisher's exact test) [not done yet, only for deletions at the beginning]
#    : Now I process subs, del, ins separately (they were being mixed sometimes
#    : For indels I use the group information defined at the beginning of the script

#mutations[which(mutations$label == ""),"label"] <- "OK";

ok_muts <- mutations[grep("OK",mutations$label),         ]
sub     <- ok_muts[grep("[-I]",ok_muts[,"mut"],invert=T),]
ins     <- ok_muts[grep("I",   ok_muts[,"mut"]),         ]
del     <- ok_muts[grep("-",   ok_muts[,"mut"]),     ]
sub     <- sub[order(sub$sampleID, sub$chr, sub$pos),    ]
ins     <- ins[order(ins$sampleID, ins$chr, ins$pos),    ]
del     <- del[order(del$sampleID, del$chr, del$pos),    ]

#To store the new data
new_mutations <- mutations[0,]

# Deletions (defined in mutations$indel_group):
# For every "OK" deletion, get its del-groupID and find all the other deletions belonging
# to that group. Merge them and create a new entry in mutations: combine pvalues, vaf, etc
# For every "OK" deletion first check it hasn't been already merged
for(j in 1:nrow(del)) {
	indel_group       <- del[j,"indel_group"]
	if(nrow(new_mutations[which(new_mutations$indel_group==indel_group),]) > 0) {
		next; #we already have one from the group of indels
	}
	indels_from_group                                  <- mutations[which(mutations$indel_group==indel_group),]
	new_mutations                                      <- rbind(new_mutations,indels_from_group[1,])
	new_mutations[nrow(new_mutations),"pos"          ] <- min  (indels_from_group$pos              )
	new_mutations[nrow(new_mutations),"vaf"          ] <- mean (indels_from_group$vaf              )
	new_mutations[nrow(new_mutations),"tum_globalvaf"] <- mean (indels_from_group$tum_globalvaf    )
	new_mutations[nrow(new_mutations),"pval"         ] <- min  (indels_from_group$pval             )
	new_mutations[nrow(new_mutations),"qval"         ] <- min  (indels_from_group$qval             )
	new_mutations[nrow(new_mutations),"label"        ] <- paste(indels_from_group$label,collapse="")
	new_mutations[nrow(new_mutations),"ref"          ] <- indels_from_group[1,"deleted_seq"]
}

# Insertions. No need to look for consecutive INS. Just add them to new_mutations
new_mutations <- rbind(new_mutations,ins);


# Substitutions: merge consecutive... [using Iñigo's code]
d = sub$pos-(1:nrow(sub))
runs = rle(d)
rmpos = rep(0,nrow(sub))
runstarts = cumsum(runs$length)-runs$length+1
for (h in 1:length(runs$length)) {
    if (runs$length[h]>1) { # Adjacent mutations
        mutcluster                         = runstarts[h]:(runstarts[h]+runs$lengths[h]-1)
        rmpos[mutcluster[-1]             ] = 1 # Removing all the affected rows except the first one (which we will edit to capture the complex event)
        sub[mutcluster[1],"ref"          ] = paste(sub[mutcluster,"ref"          ],collapse="")
        sub[mutcluster[1],"mut"          ] = paste(sub[mutcluster,"mut"          ],collapse="")
        sub[mutcluster[1],"mu"           ] = mean (sub[mutcluster,"mu"           ]            )
        sub[mutcluster[1],"tum_globalvaf"] = mean (sub[mutcluster,"tum_globalvaf"]            )
        sub[mutcluster[1],"vaf"          ] = mean (sub[mutcluster,"vaf"          ]            )
        sub[mutcluster[1],"pval"         ] = min  (sub[mutcluster,"pval"         ]            )
        sub[mutcluster[1],"qval"         ] = min  (sub[mutcluster,"qval"         ]            )
    }
}
sub = sub[!rmpos,]        
new_mutations <- rbind(new_mutations,sub);

mutations.old <- mutations
mutations     <- new_mutations
mutations[which(mutations$label == ""),"label"] <- "OK;";
write.table(mutations, file=sprintf("%s/mutations_q01_%s.annotruns.1.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
indels_f <- length(grep("[-I]",mutations[,"mut"]));
subs_f   <- nrow(mutations)-indels_f;
cat("#AFTER_MERGING_RUNS\t",nrow(mutations[which(mutations$label == ""),]),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");

####################################################################################################
# h. Filtering out mutations if reads carry indels within 10 bp of the mutation (allowing for mapQ<30)
#    Using mapQ>=30 in the calling of indels tends to remove reads with genuine indels and
#    so misses potential false positives associated with the indel

flank = 10 # bp around the mutation where we will look for indels
rmpos = NULL
for (h in 1:length(samples)) {
	cat("sample=",h," out of=", length(samples),"\n")
    m = mutations[mutations$sampleID==samples[h] & !(mutations$mut %in% c("-","INS")),]
    s <- samples[h]
    if (nrow(m)>0) { 
        for (r in 1:nrow(m)) {
        	f <- NULL;
        	#if(length(grep("out",m$sampleID[r]))>=1) {
        	#	f = sprintf("/lustre/scratch116/casm/cgp/users/fa8/BLADDER_SHEARWATER/BAMS_REALIGNED/%s.bam",m$sampleID[r])
				f = sprintf("%s/%s.bam",bams_path,s)
				#cat("Going to read: ", f, "\n");
        	#} else {
	        #    f = sprintf("/lustre/scratch116/casm/cgp/users/fa8/BLADDER_SHEARWATER/BAMS/%s.bam",m$sampleID[r])
	        #}
			f <- gsub(" ","",f);
			#cat("     ", m$chr[r], m$pos[r]-flank, m$pos[r]+flank," f=",f,"\n",sep=" : ")
            n_indels = sum(bam2R(f, m$chr[r], m$pos[r]-flank, m$pos[r]+flank, q=30, mask=3844, mq=10)[,c("-","INS","_","ins")]) # Number of reads with an indel around the mutation
            cat("n_indels=",n_indels," for mutation=",r, " of sample=",s,", xfw+xbw=",m$xfw[r]+m$xbw[r],"\n")
            #if (n_indels > 5*(m$xfw[r]+m$xbw[r])) { 
            if (n_indels > 3*(m$xfw[r]+m$xbw[r])) {  # Modified by Fede
                rmpos = c(rmpos, rownames(m)[r])
            }
        }
    }
}

if (length(rmpos)>0) {
    filt4 = mutations[rmpos,]; filt4$filter = "Near_hidden_indel"
	mutations[which(mutations$label == "OK;"),"label"] <- "";
	mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"near_hidden_indel;",sep="");
	mutations[which(mutations$label == ""),"label"] <- "OK;";
    write.table(mutations, file=sprintf("%s/mutations_q01_%s.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
    #write.table(rbind(filt1,filt2,filt3,filt4), file=sprintf("%s/filteredout_q01_%s.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
}
write.table(mutations, file=sprintf("%s/mutations_q01_%s.annotruns.2.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
indels_f <- length(grep("[-I]",mutations[grep("OK",mutations$label),"mut"]));
subs_f   <- nrow(mutations[grep("OK",mutations$label),])-indels_f;
cat("#AFTER_NEARBY_INDELS\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="")

############################################################################################################
#Create vcf file and run vagrent
# 
# mutations <- read.table(paste0("Mutation_calls/mutations_q01_", dataset_name,".annotruns.2.txt"), stringsAsFactors = F, header = T, sep = "\t")
# 
# vcf_file = mutations[,c(2:5)]
# names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
# 
# #Add the required null columns
# vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
# 
# #Rearrange into correct order
# vcf_file = vcf_file[,c(1,2,5,3,4,6,7,8)]
# 
# #Write files
# write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0(dataset_name,"_mutations_all.vcf"), row.names = FALSE)
# 
# 
# vcf_header_path = "/lustre/scratch119/casm/team154pc/em16/CHIP/shearwater/VCF_header_for_VaGrent.txt"
# vcf_path = paste0("/lustre/scratch119/casm/team154pc/em16/CHIP/shearwater/",dataset_name,"_mutations_all.vcf")
# vagrent_input_path = paste0("/lustre/scratch119/casm/team154pc/em16/CHIP/shearwater/",dataset_name, "_mutations_all_header.vcf")
# vagrent_output_path = paste0(vagrent_input_path,".annot")
#  
# #1. paste vcf file to a dummy header file
# system(paste0("cat ",vcf_header_path," ",vcf_path," > ", vagrent_input_path))
# #2. commands to run vagrent
# system(paste0("perl-5.16.3 -I /software/CGP/canpipe/live/lib/perl5 /software/CGP/canpipe/live/bin/AnnotateVcf.pl -i ",vagrent_input_path," -o ",vagrent_output_path," -sp Human -as NCBI37 -c /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75/vagrent.cache.gz"))
# #3. import vagrent output
# vagrent_output = fread(vagrent_output_path,skip = "#CHROM")
# annot_info = as.data.frame(str_split(vagrent_output$INFO, pattern = ";",simplify = TRUE), stringsAsFactors = FALSE)
# colnames(annot_info) <- c("VT","VD","VC","VW")
#  
# annot_info$VC <- gsub(x=annot_info$VC, pattern = "VC=", replacement = "")
# annot_info$VT <- gsub(x=annot_info$VT, pattern = "VT=", replacement = "")
# annot_info$VW <- gsub(x=annot_info$VW, pattern = "VW=", replacement = "")
# annot_info$VD <- gsub(x=annot_info$VD, pattern = "VD=", replacement = "")
#  
#  
# #Attempt to functionalize neatly
# split_vagrent_output = function(df,split_col,col_IDs = c("Gene","Transcript","RNA","CDS","Protein","Type","SO_codes")) {
#   col = df[[split_col]]
#   output = matrix(nrow = nrow(df), ncol = length(col_IDs))
#   for(i in 1:length(col_IDs)) {
#     output[,i] = str_split(col, pattern = "\\|", simplify = TRUE)[,i]
#   }
#   colnames(output) = col_IDs
#   return(as.data.frame(output))
# }
# 
# split_output  <- split_vagrent_output(df = annot_info,split_col = "VD")
#  
#  write.table(split_output, file = paste0(dataset_name,"vagrent_output.txt"), sep = "\t", quote = FALSE)
# 
# ############################################################################################################
# #Merge shearwater and vagrent outputs
# 
# 
# 
# ############################################################################################################
# #Remove germline based on VAF
# #>0.44
# 
# 
# ############################################################################################################
# 
# #location of perl lib and vagrent installation - ar31
# #system(paste0("perl5.26.2 -I /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater/anaconda_varcall_env/lib/perl5 /lustre/scratch119/humgen/projects/cnv_15x/clonal_hematopoiesis_ibd/newcastle_collaboration/shearwater/anaconda_varcall_env/bin/AnnotateVcf.pl -i ",vagrent_input_path," -o ",vagrent_output_path," -sp Human -as NCBI37 -c /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75/vagrent.cache.gz"))
