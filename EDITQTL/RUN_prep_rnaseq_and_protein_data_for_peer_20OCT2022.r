#!/usr/bin/Rscript

##SETUP | load libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
suppressMessages(library(data.table))
suppressMessages(library(peer))
suppressMessages(library(limma))

##SETUP | parse command line argumentss ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 1: reading arguments\n")
args <- commandArgs(trailingOnly=TRUE)
nsv <- as.integer(args[[1]]) ## nsv <- 2
ome <- args[[2]] ## ome <- "rna"

##LOAD DATAs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 2: loading data\n")
data <- readRDS("/sc/arion/projects/psychgen/lbp/data/PEER/LIVPM_20OCT2022_RNA_PRO_RNAEDITING.RDS")[[ome]]
qcMtx <- data$peerMatrices$QC2CALLSET
qsMtx <- data$peerMatrices$QC2STRICTCALLSET
qcDes <- data$peerDesigns$QC2CALLSET
qsDes <- data$peerDesigns$QC2STRICTCALLSET
qcCov <- data$peerSampleMetadata$QC2CALLSET
qsCov <- data$peerSampleMetadata$QC2STRICTCALLSET
qcVob <- data$peerVobjects$QC2CALLSET
qsVob <- data$peerVobjects$QC2STRICTCALLSET
myFrm <- data$peermyFrm
myCov <- unlist(strsplit(as.character(myFrm)[2], split=" "))
myCov <- c(myCov[myCov != "+"], paste0("PEER",1:nsv))
myCv2 <- myCov[myCov!="mymet_postmortem"]
myFrm <- as.formula(paste("~", paste(myCov,  collapse=" + " )))
myFr2 <- as.formula(paste("~", paste(myCv2,  collapse=" + " )))

##PEER  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 3: running peer\n")
qcMod <- PEER()
qsMod <- PEER()
PEER_setPhenoMean(qcMod,qcMtx)
PEER_setPhenoMean(qsMod,qsMtx)
if (ome == "rna") {
    PEER_setPhenoVar(qcMod,t(1/qcVob$weights))
    PEER_setPhenoVar(qsMod,t(1/qsVob$weights))
}
PEER_setCovariates(qcMod, qcDes)
PEER_setCovariates(qsMod, qsDes)
PEER_setNk(qcMod, nsv)
PEER_setNk(qsMod, nsv)
PEER_update(qcMod)
PEER_update(qsMod)
qcK <- ncol(qcDes)+1
qsK <- ncol(qsDes)+1
qcSV <- PEER_getX(qcMod)[,qcK:(qcK+nsv-1),drop=F]
qsSV <- PEER_getX(qsMod)[,qsK:(qsK+nsv-1),drop=F]
colnames(qcSV) <- paste0("PEER",1:nsv)
colnames(qsSV) <- paste0("PEER",1:nsv)
qcSV <- data.table( IID = rownames(qcMtx), qcSV )
qsSV <- data.table( IID = rownames(qsMtx), qsSV )

##RESIDUALS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 4: calculating residuals\n")
qcCov <- as.data.table(qcCov, keep.rownames="IID_ISMMS")
qsCov <- as.data.table(qsCov, keep.rownames="IID_ISMMS")
qcCov <- merge(qcCov, qcSV, by.x="IID_ISMMS", by.y="IID")
qsCov <- merge(qsCov, qsSV, by.x="IID_ISMMS", by.y="IID")
qcMet <- qcCov[, c("IID_ISMMS",myCov), with=F]
qsMet <- qsCov[, c("IID_ISMMS",myCov), with=F]
qcMet <- as.data.frame(qcMet)
qsMet <- as.data.frame(qsMet)
rownames(qcMet) <- qcMet$IID_ISMMS
rownames(qsMet) <- qsMet$IID_ISMMS
qcMet$IID_ISMMS <- NULL
qsMet$IID_ISMMS <- NULL
qcMet$mymet_postmortem <- NULL #so we dont regress out the effect of LIVPM status
qsMet$mymet_postmortem <- NULL
qcDes <- model.matrix(myFr2, qcMet)
qsDes <- model.matrix(myFr2, qsMet)
qcDes <- qcDes[colnames(qcVob),]
qsDes <- qsDes[colnames(qsVob),]
qcFit  <- lmFit(qcVob, qcDes)
qsFit  <- lmFit(qsVob, qsDes)
qcRes <- residuals(qcFit, qcVob,)
qsRes <- residuals(qsFit, qsVob,)

## BED FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 4: formatting data\n")
if (ome=="rna" | ome=="protein"){
    setwd("/sc/arion/projects/psychgen2/lbp/data/RAW/rna/bulk/fromSema4/Merged_Batches/201112_A00734_0092_AHLJH3DSXY.201112_A00734_0093_BHLGNFDSXY.Merged/LBPSEMA4BLOOD018/RAPiD")
    bed <- fread("featureCounts/LBPSEMA4BLOOD018.primary.txt")[,.(Geneid, Chr, Strand, Length, Start, End)]
    Starts <- sapply(bed$Start,function(x){unlist(strsplit(x,";",fixed=T))})
    Ends <- sapply(bed$End,function(x){unlist(strsplit(x,";",fixed=T))})
    Strands <- sapply(bed$Strand,function(x){unique(unlist(strsplit(x,";",fixed=T)))})
    Chrs <- sapply(bed$Chr,function(x){unique(unlist(strsplit(x,";",fixed=T)))})
    TSS <- matrix(data=NA,nrow=length(Starts),ncol=1)
    Gids <- bed$Geneid
    for(i in 1:nrow(TSS)){
        if(length(Strands[i])>1){
            cat("error on index",i,"\n")
        }
        if(length(Chrs[i])>1){
            cat("error on index",i,"\n")
        }
        if(Strands[i]=="+"){
            TSS[i]=min(as.numeric(Starts[[i]]))
        }else if(Strands[i]=="-"){
            TSS[i]=max(as.numeric(Ends[[i]]))
        }
    }
    bed <- data.table(Chr=Chrs, Start=TSS[,1], End=TSS[,1]+1, Geneid=Gids, Strand=Strands)
    bed[,Chr:=tstrsplit(Chr, split=";", fixed=T, keep=1L)]
    bed[,Chr:=gsub("chr", "", Chr)]
    bed[,Strand:=tstrsplit(Strand, split=";", fixed=T, keep=1L)]
    bed <- bed[,.(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand)]
    bed[,pid:=tstrsplit(pid, split=".", fixed=T, keep=1L)]
    bed[,gid:=tstrsplit(gid, split=".", fixed=T, keep=1L)]
    bed <- bed[gid %in% intersect(bed$gid, rownames(qcRes))]
    bed <- bed[order(as.integer(`#Chr`), start)]
} else if (ome=="rnaediting") {
    bed <- data$peerRedBed
    colnames(bed)[ncol(bed)] <- "strand"
}
qcRes <- qcRes[bed$gid,]
qsRes <- qsRes[bed$gid,]
qcRes <- cbind(bed, qcRes)
qsRes <- cbind(bed, qsRes)

## SAVE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 5: saving\n")
odr <- paste0("/sc/arion/projects/psychgen/lbp/data/runEQTL/20OCT2022/", toupper(ome), "/")
qcOut <- paste0(odr, "QC2CALLSET_nsv", nsv, ".bed")
qsOut <- paste0(odr, "QC2STRICTCALLSET_nsv", nsv, ".bed")
fwrite(qcRes, sep='\t', row=F, quo=F, file=qcOut)
fwrite(qsRes, sep='\t', row=F, quo=F, file=qsOut)
#system(paste0("ml tabix;bgzip ", lvOut, " && tabix -p bed ", lvOut, ".gz"))
#system(paste0("ml tabix;bgzip ", pmOut, " && tabix -p bed ", pmOut, ".gz"))
cat("DONE!")
