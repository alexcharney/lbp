## setup
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
Sys.setenv(OMP_NUM_THREADS = 20)
library(sp)
library(gsubfn)
library(Matrix)
setwd("/sc/arion/projects/psychgen/lbp/data/neuroimaging/GMWM19JUL2021")
args <- commandArgs(trailingOnly=TRUE)
myFeature <- args[[1]]
myOutput <- args[[2]]
##myFeature <- "GMF199"

## read in data
data <- readRDS("/sc/arion/projects/psychgen/lbp/data/neuroimaging/GMWM19JUL2021/dataForDream19JUL2021.RDS")

## expression data
vobjDream <- data$expressionData$vobjDream

## get imaging data
gFeatures <- rownames(data$gm$data)
wFeatures <- rownames(data$wm$data)
isGray <- myFeature %in% gFeatures
if (isGray){
    imagingData <- as.data.frame(t(data$gm$data[myFeature,]))
} else {
    imagingData <- as.data.frame(t(data$wm$data[myFeature,]))
}

## make formulas
lelVar <- data$covariates$lel
leoVar <- data$covariates$leo
lelCatVar <- c("IID_ISMMS", "mymet_sex", "mymet_depletionbatch")
lelNumVar <- c("mymet_rin", "neuronal", "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
               "RNASeqMetrics_PCT_MRNA_BASES", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
              "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR")
leoCatVar <- c("IID_ISMMS", "mymet_sex", "mymet_bank")
leoNumVar <- c("mymet_rin", "mymet_rna_conc_ngul", "RNASeqMetrics_PCT_INTRONIC_BASES", "RNASeqMetrics_PCT_CODING_BASES", "STAR_pct_of_reads_unmapped_other")
lelMod0 <- paste(paste0("(1|", lelCatVar, ")"), collapse=" + ")
leoMod0 <- paste(paste0("(1|", leoCatVar, ")"), collapse=" + ")
lelMod1 <- paste(lelNumVar, collapse=" + ")
leoMod1 <- paste(leoNumVar, collapse=" + ")
lelForm <- as.formula(paste0("~ imagingData + ",lelMod0, " + ", lelMod1))
leoForm <- as.formula(paste0("~ imagingData + ",leoMod0, " + ", leoMod1))

## covariates
lelMet <- as.data.frame(data$expressionData$Covariates[,lelVar,with=F])
leoMet <- as.data.frame(data$expressionData$Covariates[,leoVar,with=F])
rownames(lelMet) <- data$expressionData$Covariates$SAMPLE_ISMMS
rownames(leoMet) <- data$expressionData$Covariates$SAMPLE_ISMMS
lelMet <- lelMet[colnames(vobjDream$E),]
leoMet <- leoMet[colnames(vobjDream$E),]

## subset for samples with imaging
iSam <- rownames(imagingData)
eSam <- rownames(lelMet)
cSam <- colnames(vobjDream)
kSam <- iSam[iSam %in% eSam & iSam %in% cSam]
imagingData <- imagingData[kSam,,drop=F]
lelMet <- lelMet[kSam,]
leoMet <- leoMet[kSam,]
vobjDream <- vobjDream[,kSam]
colnames(imagingData) <- "imagingData"
lelMet <- cbind(lelMet, imagingData)
leoMet <- cbind(leoMet, imagingData)
for (i in lelNumVar) {lelMet[,i] <- scale(lelMet[,i])}
for (i in leoNumVar) {leoMet[,i] <- scale(leoMet[,i])}

## dream
coefcol <- "imagingData"
lelFit <- dream( vobjDream, lelForm, lelMet, BPPARAM = MulticoreParam(5))
lelGroupDE <- as.data.table(topTable(lelFit, coef=coefcol, number=nrow(vobjDream)), keep.rownames=T)
colnames(lelGroupDE)[1] <- "gene"
leoFit <- dream( vobjDream, leoForm, leoMet, BPPARAM = MulticoreParam(5))
leoGroupDE <- as.data.table(topTable(leoFit, coef=coefcol, number=nrow(vobjDream)), keep.rownames=T)
colnames(leoGroupDE)[1] <- "gene"
x <- merge(leoGroupDE, lelGroupDE, by="gene", suffixes=c(".leo", ".lel"))
fwrite(x, sep='\t', quo=F, row=F, file=myOutput)


