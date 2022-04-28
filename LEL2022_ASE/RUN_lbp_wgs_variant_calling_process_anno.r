#!/usr/bin/Rscript
 
##
## Libraries ---------------------------
##

suppressMessages(library(data.table))
 
##
## Arguments ---------------------------
##

args <- commandArgs(trailingOnly=TRUE)
lof <- args[[1]]
hlp <- args[[2]]
out <- args[[3]]
###lof <- "/sc/arion/projects/psychgen/lbp/scratch/wgs_annotation/chr22/TMP_sites_only.noinfo.loftee.vcf.gz.rinput"
###hlp <- "/sc/arion/projects/psychgen/lbp/scratch/wgs_annotation/chr22/TMP_sites_only.noinfo.loftee.vcf.gz"
###out <- "/sc/arion/projects/psychgen/lbp/scratch/wgs_annotation/chr22/TMP_loftee_clean.tsv"
hp1 <- system(paste("zgrep -m1 CSQ", hlp), intern=TRUE)
hp2 <- system(paste("zgrep -m1 CHROM", hlp), intern=TRUE)
hp1 <- strsplit(hp1, split=" ")[[1]][7]
hp1 <- strsplit(hp1, split="\"")[[1]][1]
hp2 <- strsplit(hp2, split="#")[[1]][2]
hp2 <- gsub("\t", "|", hp2)

##
## Data ---------------------------
##

f1 <- fread(lof, na=c("NA", "."), sep='\t')
ar1 <- paste0("LOFTEE.", unlist(strsplit(hp1, split="|", fixed=T)))
ar2 <- unlist(strsplit(hp2, split="|", fixed=T))

##
## Format loftee ---------------------------
##

colnames(f1) <- ar2
f1 <- f1[, list(CSQ = unlist(strsplit(INFO, ",", fixed=T))), by=list(CHROM, POS, REF, ALT)]
f1[,CSQ:=gsub("CSQ=", "", CSQ)]
f1[,c(ar1):=tstrsplit(CSQ,split="|",fixed=T)]
for(i in ar1) f1[get(i)=="", try(i):=NA]
f1[,CSQ:=NULL]

##
## Write --------------------------------
##

fwrite(f1, row=F, quo=F, sep='\t', file=out)

