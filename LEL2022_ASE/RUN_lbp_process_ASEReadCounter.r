#!/usr/bin/Rscript

## FILTERS APPLIED HERE:
## 1)  WGS GT of 0/1 1/0 0|1 1|0
## 2)  Depth >= 10 in WGS and RNAseq

##
## Libraries ---------------------------
##

suppressMessages(library(data.table))
 
##
## Arguments ---------------------------
##

args <- commandArgs(trailingOnly=TRUE)
sid <- args[[1]] ###sid <- "LBPSEMA4BRAIN695"

##
## VARIANTS ---------------------------
##

vcf <- fread(paste0("/sc/arion/projects/psychgen/lbp/data/dna/wgs_ASEReadCounterInput/", sid, ".vcf.gz.hetsites"))
colnames(vcf) <- c("CHROM", "POS", "REF", "ALT", "SID", "GT", "AD")
vcf <- vcf[GT %in% c("0/1", "1/0", "0|1", "1|0")]
vcf[,c("REF_AD", "ALT_AD"):=tstrsplit(AD, split=",")]
vcf[,REFCOUNT:=as.integer(REF_AD)]
vcf[,ALTCOUNT:=as.integer(ALT_AD)]
vcf[,TOTCOUNT:=ALTCOUNT+REFCOUNT]
vcf[,AAF:=ALTCOUNT/TOTCOUNT]
vcf[,AD:=NULL]
vcf[,GT:=NULL]

##
## ASEReadCounter ---------------------------
##

ase <- fread(paste0("/sc/arion/projects/psychgen/lbp/results/ASEReadCounter/", sid, ".output.table"))
ase[,variantID:=NULL]
colnames(ase) <- c("CHROM", "POS", "REF", "ALT",
                  "ASE_REFCOUNT", "ASE_ALTCOUNT", "ASE_TOTCOUNT", "ASE_lowMAPQDepth", "ASE_lowBaseQDepth",
                  "ASE_rawDepth", "ASE_otherBases", "ASE_improperPairs")
add <- ase[,.(CHROM, POS, REF, ALT, ASE_REFCOUNT, ASE_ALTCOUNT, ASE_TOTCOUNT)]
add[,ASE_AAF:=ASE_ALTCOUNT/ASE_TOTCOUNT]
clean <- merge(vcf, add)[TOTCOUNT>=10 & ASE_TOTCOUNT>=10]

##
## LOFTEE ---------------------------
##

lof <- fread("/sc/arion/projects/psychgen/lbp/data/dna/wgs_GenotypeRefinement/strictAnnoFiltPassLoftee.tsv", na=c(".", "", "NA"), sep='\t')
lof[,LOFTEE.inLoftee:=1]
clean <- merge(clean, lof, all.x=T, by=c("CHROM", "POS", "REF", "ALT"))
clean[is.na(LOFTEE.inLoftee), LOFTEE.inLoftee:=0]
rm(lof)

##
## GNOMAD --------------------------------
##

gno <- fread("/sc/arion/projects/psychgen/lbp/data/dna/wgs_GenotypeRefinement/strictAnnoFiltPassGnomad.tsv", na=c(".", "", "NA"), sep='\t')
colnames(gno)[5:ncol(gno)] <- paste0("GNOMAD.", colnames(gno)[5:ncol(gno)])
clean <- merge(clean, gno, all.x=T, by=c("CHROM", "POS", "REF", "ALT"))
clean[is.na(GNOMAD.inGnomad), GNOMAD.inGnomad:=0]
rm(gno)
 
##
## SAVE  --------------------------------
##

fout1 <- paste0("/sc/arion/projects/psychgen/lbp/results/ASEReadCounter/", sid, ".output.table.gt10x.anno")
fout2 <- paste0("/sc/arion/projects/psychgen/lbp/results/ASEReadCounter/", sid, ".output.table.gt10x.anno.success")
fwrite(clean, row=F, quo=F, sep='\t', file=fout1)
system(paste0("touch ", fout2))

