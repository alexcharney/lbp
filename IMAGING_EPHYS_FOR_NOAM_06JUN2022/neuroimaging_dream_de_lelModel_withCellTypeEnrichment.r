
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
Sys.setenv(OMP_NUM_THREADS = 6)
setwd("/sc/arion/projects/psychgen/lbp/data/neuroimaging/LEO/GMWM19JUL2021/dataForDream19JUL2021")

## Read Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
args <- commandArgs(trailingOnly=TRUE)
myFeature <- args[[1]]
runDe <- args[[2]]
myOutpt1 <- args[[3]]
myOutpt2 <- gsub(".tsv", ".noNeu.tsv", myOutpt1)
myOutpt3 <- gsub(".tsv", ".noNeu.pfcCorr.tsv", myOutpt1)
myOutpt4 <- gsub(".tsv", ".CellTypeEnrichment.tsv", myOutpt1)
####for testing:
##myFeature <- "GMF000"
##runDe <- "TRUE"
##myOutpt1 <- "/sc/arion/projects/psychgen/lbp/results/neuroimaging/awcDreamTest/features/GMF000.tsv"
##myOutpt2 <- gsub(".tsv", ".noNeu.tsv", myOutpt1)
##myOutpt3 <- gsub(".tsv", ".noNeu.pfcCorr.tsv", myOutpt1)

## Run Dream ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (runDe == "TRUE"){
    
    ## read imaging data
    exdata <- readRDS("dataForDream19JUL2021_expressionData.RDS")
    gmdata <- readRDS("dataForDream19JUL2021_imagingData_gm_data.RDS")
    wmdata <- readRDS("dataForDream19JUL2021_imagingData_wm_data.RDS")
    cvdat2 <- readRDS("dataForDream19JUL2021_covDay.RDS")
    cvdat3 <- readRDS("dataForDream19JUL2021_pfcPC.RDS")
    exdata$covariates <- merge(exdata$covariates, cvdat2[,.(IID_ISMMS=iid, lbpday)], by="IID_ISMMS", all.x=T)
    exdata$covariates <- merge(exdata$covariates, cvdat3[,.(IID_ISMMS=iid, PFCPC)], by="IID_ISMMS", all.x=T)
    
    ## expression data
    vobjDream <- exdata$vobjDream
    
    ## get imaging data
    gFeatures <- rownames(gmdata)
    wFeatures <- rownames(wmdata)
    isGray <- myFeature %in% gFeatures
    if (isGray){
        imagingData <- as.data.frame(t(gmdata[myFeature,]))
    } else {
        imagingData <- as.data.frame(t(wmdata[myFeature,]))
    }
    
    ## make formulas
    lelCatVar <- c("IID_ISMMS", "mymet_sex", "mymet_depletionbatch", "mymet_phe")
    lelNumVr1 <- c("mymet_age", "lbpday", "mymet_rin", "neuronal", "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
                  "RNASeqMetrics_PCT_MRNA_BASES", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
                  "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR")
    lelNumVr2 <- c("mymet_age", "lbpday", "mymet_rin", "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
                  "RNASeqMetrics_PCT_MRNA_BASES", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
                  "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR")
    lelNumVr3 <- c("mymet_age", "lbpday", "mymet_rin", "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
                  "RNASeqMetrics_PCT_MRNA_BASES", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
                  "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR", "PFCPC")
    lelMod0 <- paste(paste0("(1|", lelCatVar, ")"), collapse=" + ")
    lelMod1 <- paste(lelNumVr1, collapse=" + ")
    lelMod2 <- paste(lelNumVr2, collapse=" + ")
    lelMod3 <- paste(lelNumVr3, collapse=" + ")
    lelFrm1 <- as.formula(paste0("~ imagingData + ",lelMod0, " + ", lelMod1))
    lelFrm2 <- as.formula(paste0("~ imagingData + ",lelMod0, " + ", lelMod2))
    lelFrm3 <- as.formula(paste0("~ imagingData + ",lelMod0, " + ", lelMod3))
    
    ## covariates
    lelMt1 <- as.data.frame(exdata$covariates[,c(lelCatVar,lelNumVr1),with=F])
    lelMt2 <- as.data.frame(exdata$covariates[,c(lelCatVar,lelNumVr2),with=F])
    lelMt3 <- as.data.frame(exdata$covariates[,c(lelCatVar,lelNumVr3),with=F])
    rownames(lelMt1) <- exdata$covariates$SAMPLE_ISMMS
    rownames(lelMt2) <- exdata$covariates$SAMPLE_ISMMS
    rownames(lelMt3) <- exdata$covariates$SAMPLE_ISMMS
    lelMt1 <- lelMt1[colnames(vobjDream$E),]
    lelMt2 <- lelMt2[colnames(vobjDream$E),]
    lelMt3 <- lelMt3[colnames(vobjDream$E),]
    
    ## subset for samples with imaging
    iSam <- rownames(imagingData)
    eSam <- rownames(lelMt1)
    cSam <- colnames(vobjDream)
    kSam <- iSam[iSam %in% eSam & iSam %in% cSam]
    imagingData <- imagingData[kSam,,drop=F]
    lelMt1 <- lelMt1[kSam,]
    lelMt2 <- lelMt2[kSam,]
    lelMt3 <- lelMt3[kSam,]
    vobjDream <- vobjDream[,kSam]
    colnames(imagingData) <- "imagingData"
    lelMt1 <- cbind(lelMt1, imagingData)
    lelMt2 <- cbind(lelMt2, imagingData)
    lelMt3 <- cbind(lelMt3, imagingData)
    for (i in lelNumVr1) {lelMt1[,i] <- scale(lelMt1[,i])}
    for (i in lelNumVr2) {lelMt2[,i] <- scale(lelMt2[,i])}
    for (i in lelNumVr3) {lelMt3[,i] <- scale(lelMt3[,i])}
    
    ## dream
    lelFit <- dream( vobjDream, lelFrm1, lelMt1, BPPARAM = MulticoreParam(5))
    de1 <- as.data.table(topTable(lelFit, coef="imagingData", number=nrow(vobjDream)), keep.rownames="gene")
    fwrite(de1, sep='\t', quo=F, row=F, file=myOutpt1)
    
    ## dream - noNeu
    lelFt2 <- dream( vobjDream, lelFrm2, lelMt2, BPPARAM = MulticoreParam(5))
    de2 <- as.data.table(topTable(lelFt2, coef="imagingData", number=nrow(vobjDream)), keep.rownames="gene")
    fwrite(de2, sep='\t', quo=F, row=F, file=myOutpt2)
    
    ## dream - no neu, with pfccorr 
    lelFt3 <- dream( vobjDream, lelFrm3, lelMt3, BPPARAM = MulticoreParam(5))
    de3 <- as.data.table(topTable(lelFt3, coef="imagingData", number=nrow(vobjDream)), keep.rownames="gene")
    fwrite(de3, sep='\t', quo=F, row=F, file=myOutpt3)

} else {
    de1 <- fread(myOutpt1)
    de2 <- fread(myOutpt2)
    de3 <- fread(myOutpt3)
}

## Run CT Enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##format de results
mycen1 <- de1[,.(gene, p=P.Value, padj=adj.P.Val, logFC)]
mycen2 <- de2[,.(gene, p=P.Value, padj=adj.P.Val, logFC)]
mycen3 <- de3[,.(gene, p=P.Value, padj=adj.P.Val, logFC)]
mycen1[,gene:=tstrsplit(gene, split=".", fixed=T, keep=1L)]
mycen2[,gene:=tstrsplit(gene, split=".", fixed=T, keep=1L)]
mycen3[,gene:=tstrsplit(gene, split=".", fixed=T, keep=1L)]
mycen1[, nomUpDEG:=0]
mycen1[, nomDwDEG:=0]
mycen1[, sigUpDEG:=0]
mycen1[, sigDwDEG:=0]
mycen2[, nomUpDEG:=0]
mycen2[, nomDwDEG:=0]
mycen2[, sigUpDEG:=0]
mycen2[, sigDwDEG:=0]
mycen3[, nomUpDEG:=0]
mycen3[, nomDwDEG:=0]
mycen3[, sigUpDEG:=0]
mycen3[, sigDwDEG:=0]
mycen1[p<0.05 & logFC>0, nomUpDEG:=1]
mycen1[p<0.05 & logFC<0, nomDwDEG:=1]
mycen1[padj<0.05 & logFC>0, sigUpDEG:=1]
mycen1[padj<0.05 & logFC<0, sigDwDEG:=1]
mycen2[p<0.05 & logFC>0, nomUpDEG:=1]
mycen2[p<0.05 & logFC<0, nomDwDEG:=1]
mycen2[padj<0.05 & logFC>0, sigUpDEG:=1]
mycen2[padj<0.05 & logFC<0, sigDwDEG:=1]
mycen3[p<0.05 & logFC>0, nomUpDEG:=1]
mycen3[p<0.05 & logFC<0, nomDwDEG:=1]
mycen3[padj<0.05 & logFC>0, sigUpDEG:=1]
mycen3[padj<0.05 & logFC<0, sigDwDEG:=1]

##cell type markers
DECONVOLUTION_REFERENCE <- "/sc/arion/projects/psychgen/lbp/files/lake_for_cibersort_3.Rdata"
reference <- new.env()
load(DECONVOLUTION_REFERENCE, env=reference)
sce2 <- exprs(reference$eset)
all_cell_type <- reference$eset@phenoData@data$cluster_name2
dacList <- reference$dacListAll1
gluGenes <- dacList$GLU
gabGenes <- dacList$GABA
astGenes <- dacList$AST
odcGenes <- dacList$ODC
micGenes <- dacList$MG
neuGenes <- unique(c(gluGenes, gabGenes))
gliGenes <- unique(c(astGenes, odcGenes, micGenes))
rm(reference)
rm(dacList)
mycen1[, GLU:=0]
mycen1[, GABA:=0]
mycen1[, AST:=0]
mycen1[, ODC:=0]
mycen1[, MG:=0]
mycen1[, NEU:=0]
mycen1[, GLI:=0]
mycen2[, GLU:=0]
mycen2[, GABA:=0]
mycen2[, AST:=0]
mycen2[, ODC:=0]
mycen2[, MG:=0]
mycen2[, NEU:=0]
mycen2[, GLI:=0]
mycen3[, GLU:=0]
mycen3[, GABA:=0]
mycen3[, AST:=0]
mycen3[, ODC:=0]
mycen3[, MG:=0]
mycen3[, NEU:=0]
mycen3[, GLI:=0]
mycen1[gene %in% gluGenes, GLU:=1]
mycen1[gene %in% gabGenes, GABA:=1]
mycen1[gene %in% astGenes, AST:=1]
mycen1[gene %in% odcGenes, ODC:=1]
mycen1[gene %in% micGenes, MG:=1]
mycen1[gene %in% neuGenes, NEU:=1]
mycen1[gene %in% gliGenes, GLI:=1]
mycen2[gene %in% gluGenes, GLU:=1]
mycen2[gene %in% gabGenes, GABA:=1]
mycen2[gene %in% astGenes, AST:=1]
mycen2[gene %in% odcGenes, ODC:=1]
mycen2[gene %in% micGenes, MG:=1]
mycen2[gene %in% neuGenes, NEU:=1]
mycen2[gene %in% gliGenes, GLI:=1]
mycen3[gene %in% gluGenes, GLU:=1]
mycen3[gene %in% gabGenes, GABA:=1]
mycen3[gene %in% astGenes, AST:=1]
mycen3[gene %in% odcGenes, ODC:=1]
mycen3[gene %in% micGenes, MG:=1]
mycen3[gene %in% neuGenes, NEU:=1]
mycen3[gene %in% gliGenes, GLI:=1]
modresCT <- c()
for (i in c("nomUpDEG", "nomDwDEG", "sigUpDEG", "sigDwDEG")){
    for (j in c("GLU", "GABA", "AST", "ODC", "MG", "NEU", "GLI")){
        check1 <- uniqueN(mycen1[[i]]) == 2 & uniqueN(mycen1[[j]]) == 2
        check2 <- uniqueN(mycen2[[i]]) == 2 & uniqueN(mycen2[[j]]) == 2
        check3 <- uniqueN(mycen3[[i]]) == 2 & uniqueN(mycen3[[j]]) == 2        
        if (check1 == TRUE){
            curRes <- fisher.test(table(mycen1[[i]], mycen1[[j]]))
            nDeg <- sum(mycen1[[i]]==1)
            nCel <- sum(mycen1[[j]]==1)
            nBot <- sum(mycen1[[i]]==1 & mycen1[[j]]==1)
            or <- curRes$estimate
            pv <- curRes$p.value
            add <- data.table(deRun="full", degStatus=i, cell=j, fisherOR=or, fisherP=pv, ndeg=nDeg, nmarkers=nCel, nshared=nBot)
            modresCT <- rbind(modresCT, add)
        }
        if (check2 == TRUE){
            curRes <- fisher.test(table(mycen2[[i]], mycen2[[j]]))
            nDeg <- sum(mycen2[[i]]==1)
            nCel <- sum(mycen2[[j]]==1)
            nBot <- sum(mycen2[[i]]==1 & mycen2[[j]]==1)
            or <- curRes$estimate
            pv <- curRes$p.value
            add <- data.table(deRun="noNeu", degStatus=i, cell=j, fisherOR=or, fisherP=pv, ndeg=nDeg, nmarkers=nCel, nshared=nBot)
            modresCT <- rbind(modresCT, add)
        }
        if (check3 == TRUE){
            curRes <- fisher.test(table(mycen3[[i]], mycen3[[j]]))
            nDeg <- sum(mycen3[[i]]==1)
            nCel <- sum(mycen3[[j]]==1)
            nBot <- sum(mycen3[[i]]==1 & mycen3[[j]]==1)
            or <- curRes$estimate
            pv <- curRes$p.value
            add <- data.table(deRun="noNeuPfcCorr", degStatus=i, cell=j, fisherOR=or, fisherP=pv, ndeg=nDeg, nmarkers=nCel, nshared=nBot)
            modresCT <- rbind(modresCT, add)
        }
    }
}
modresCT <- data.table( "imagingFeatureIndex" = myFeature, modresCT )
fwrite(modresCT, sep='\t', quo=F, row=F, file=myOutpt4)
