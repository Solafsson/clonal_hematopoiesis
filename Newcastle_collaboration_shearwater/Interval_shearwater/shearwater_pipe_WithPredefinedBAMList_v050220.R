# Inigo Martincorena - Feb 2017
# Modified by Federico Abascal - April 2017, to:
#   -make the preparation of bam files faster
#   -combine the preparation of bams and the running of shearwater into one single script
# Post-processing and filtering of the potential mutations is to be done by a different
#   script


# Input arguments:
# 1. Analysis name
# 2. List of case bams
# 3. List of normal bams
# 4. Bait capture bed file
# 5. Output directory

###########################################################################################
# Modify in non-Sanger environments:
R_program = "/software/R-3.6.1/bin/Rscript"
#
# You will need to modify the line containing the call to the queueing system:
#   sprintf("bsub -J %s -q normal -R 'select[me...
#
###########################################################################################


## 1. Environment

args = commandArgs(TRUE)
analysis_name = args[1];                # e.g. "PD12345_oeso_kk"
samples1_file = args[2];                # "case" list of samples
samples2_file = args[3];                # "normal" list of samples
baits_bed     = args[4];				# Bedfile with the list of regions captured
outdir        = args[5];
numsegments_per_job = 100
shearwater_fun = "/nfs/users/nfs_a/ar31/shearwater/shearwater_Interval/wrapper_shearwaterML_multiplebams_v050220.R"; # Add path to this script if necessary

cat("Analysis_name = ", analysis_name,"\n")
cat("samples1_file = ", samples1_file,"\n")
cat("samples2_file = ", samples2_file,"\n")


# Libraries
library("GenomicRanges")
library("rtracklayer")


## 2. Loading and preparing input files

# Loading input files
cases = read.table(samples1_file, header=0, sep="\t", stringsAsFactors=F);
normals = read.table(samples2_file, header=0, sep="\t", stringsAsFactors=F);
b = read.table(baits_bed, header=0, sep="\t", stringsAsFactors=F)
baits = GRanges(b[,1], IRanges(b[,2],b[,3]))

# Prepare lists of bam files to run shearwater
system(paste("mkdir -p",outdir))
setwd(outdir);
write.table(cases,"case.list.tsv", col.names=F,row.names=F,sep="\t",quote=F)
write.table(normals,"normal.list.tsv",col.names=F,row.names=F,sep="\t",quote=F)


## 3. Running Shearwater

entry_start = seq(from=1, to=length(baits), by=numsegments_per_job)
entry_end = pmin(entry_start+numsegments_per_job-1, length(baits))
cmds = NULL
job_prefix = paste(sample(letters,5,TRUE),collapse="")
job_count  = 1
job_vec    = vector(length=length(entry_start))
cat("Preparing to run ", length(job_vec), " shearwater jobs...\n")
for (h in 1:length(entry_start)) {
	job_id   <- paste(job_prefix,job_count,sep="")
	# Change the following line if necessary:
    cmds[length(cmds)+1] = sprintf("bsub -J %s -q long -n 10 -R 'select[mem>=3000] rusage[mem=3000]' -M3000 -e bsub.err -o bsub.out %s %s %s %0.0f %0.0f %s %s %s/shearwater_temp_%s", job_id, R_program, shearwater_fun, baits_bed, entry_start[h], entry_end[h], "case.list.tsv", "normal.list.tsv", outdir, analysis_name)
	cat("   Submitting job: ", cmds[length(cmds)], "\n");
	system(cmds[length(cmds)])
	job_vec[job_count] <- job_id
	job_count <- job_count + 1
}
