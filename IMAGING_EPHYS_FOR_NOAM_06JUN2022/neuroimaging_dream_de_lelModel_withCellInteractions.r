   
## Load Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
options(stringsAsFactors=F)
library(edgeR)
library(batchtools)
library(devtools)
library(withr)
library(limma)
library(variancePartition)
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(BiocParallel)
library(sp)
library(gsubfn)
library(Matrix)
library(Biobase)
library(qvalue)
Sys.setenv(OMP_NUM_THREADS = 6)
setwd("/sc/arion/projects/psychgen/lbp/data/neuroimaging/LEO/GMWM19JUL2021/dataForDream19JUL2021")

## Read Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
args <- commandArgs(trailingOnly=TRUE)
myFeature <- args[[1]]
myOutDir <- args[[2]]
####for testing:
##myFeature <- "WMF3734" ## has lots of degs
##myFeature <- "GMF001" ##has 0 degs
##myOutDir <- "/sc/arion/projects/psychgen/lbp/results/neuroimaging/awcDreamTest/features/CellInteractions/WMF3734"
##myOutDir <- "/sc/arion/projects/psychgen/lbp/results/neuroimaging/awcDreamTest/features/CellInteractions/GMF001"
 
## Read Imaging Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gmdata <- readRDS("dataForDream19JUL2021_imagingData_gm_data.RDS")
wmdata <- readRDS("dataForDream19JUL2021_imagingData_wm_data.RDS")
cvdat2 <- readRDS("dataForDream19JUL2021_covDay.RDS")
cvdat3 <- readRDS("dataForDream19JUL2021_pfcPC.RDS")
gFeatures <- rownames(gmdata)
wFeatures <- rownames(wmdata)
isGray <- myFeature %in% gFeatures
if (isGray){
    imagingData <- as.data.frame(t(gmdata[myFeature,]))
} else {
    imagingData <- as.data.frame(t(wmdata[myFeature,]))
}

## Read Expression Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

exdata <- readRDS("dataForDream19JUL2021_expressionData.RDS")
exdata$covariates <- merge(exdata$covariates, cvdat2[,.(IID_ISMMS=iid, lbpday)], by="IID_ISMMS", all.x=T)
exdata$covariates <- merge(exdata$covariates, cvdat3[,.(IID_ISMMS=iid, PFCPC)], by="IID_ISMMS", all.x=T)
vobjDream <- exdata$vobjDream


## Make Formulas  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

formulas <- list()
cells <- c("AST", "ODC", "MG", "GABA", "GLU")
catVar <- paste(paste0("(1|", c("IID_ISMMS", "mymet_sex", "mymet_depletionbatch", "mymet_phe"), ")"), collapse=" + ")
numVar <- paste(c("mymet_age", "lbpday", "mymet_rin", "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
              "RNASeqMetrics_PCT_MRNA_BASES", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
              "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR"), collapse=" + ")
####formulas$full <- as.formula(paste0("~ imagingData + ", catVar, " + ", numVar))
for (i in cells){
    not_i <- cells[cells != i]
    not_i <- not_i[2:length(not_i)] ##to deal with fact that when all 5 cell types are in model it throws error (because they add to 1?)
    celVar <- paste("imagingData", i, sep="*")
    frcVar <- paste(not_i, collapse=" + ")
    formulas[[i]] <- as.formula(paste0("~", celVar, "+", catVar, "+", numVar, "+", frcVar))
}
    
## Make Covariates Tables  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

allVar <- c("AST", "ODC", "MG", "GABA", "GLU", "IID_ISMMS", "mymet_sex", "mymet_depletionbatch", "mymet_phe",
           "mymet_age", "lbpday", "mymet_rin", "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
           "RNASeqMetrics_PCT_MRNA_BASES", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
           "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR")
numVar <- c("AST", "ODC", "MG", "GABA", "GLU", "mymet_age", "lbpday", "mymet_rin", "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
           "RNASeqMetrics_PCT_MRNA_BASES", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
           "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR")
met <- as.data.frame(exdata$covariates[,c(allVar),with=F])
for (i in numVar) {met[,i] <- scale(met[,i])}
rownames(met) <- exdata$covariates$SAMPLE_ISMMS
   
## Harmonize Samples  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

iSam <- rownames(imagingData)
eSam <- rownames(met)
cSam <- colnames(vobjDream)
kSam <- iSam[iSam %in% eSam & iSam %in% cSam]
imagingData <- imagingData[kSam,,drop=F]
met <- met[kSam,]
vobjDream <- vobjDream[,kSam]
colnames(imagingData) <- "imagingData"
met <- cbind(met, imagingData)

## Run Dream   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

piVals <- c()
for (i in names(formulas)){
    myOutFil <- paste0(myOutDir, "/", i, ".tsv")
    dfit <- dream( vobjDream, formulas[[i]], met, BPPARAM = MulticoreParam(5))
    mycoef <- grep("imagingData:", colnames(dfit$coef), value=T) 
    de <- as.data.table(topTable(dfit, coef=mycoef, number=nrow(vobjDream)), keep.rownames="gene")
    fwrite(de, sep='\t', quo=F, row=F, file=myOutFil)
    mypi1 <- 1 - qvalue(de$P.Value)$pi0   
    add <- data.table(cell=i, pi1=mypi1)
    piVals <- rbind(piVals, add)
}
myOutFil <- paste0(myOutDir, "/PI1.tsv")
fwrite(piVals, sep='\t', quo=F, row=F, file=myOutFil)


