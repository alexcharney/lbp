
## SETUP
rm(list=ls())
options(stringsAsFactors=F)
suppressMessages(library(edgeR))
suppressMessages(library(withr))
suppressMessages(library(limma))
suppressMessages(library(variancePartition))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(BiocParallel))
suppressMessages(library(qvalue))
Sys.setenv(OMP_NUM_THREADS = 6)
args <- commandArgs(trailingOnly=TRUE)
idx <- args[[1]] #idx <- 1
setwd("/sc/arion/projects/psychgen/lbp/data/proteomics")

## READ - rnaseq data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rna <- readRDS("/sc/arion/work/charna02/symlinks/lbp/liharska2021/final.everything.RDS")
vob <- rna$vobjDream
met <- rna$covariates

## READ - protein data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pro <- readRDS("./LBP_Proteomics_fromWeipingMa_toLbpTeam_01APR2022/Data/20220214_MSSM_Living_Brain_Tissue_Global_34MPs_Combined_awcParsed.RDS")$imputed
pnm <- rownames(pro)[idx]
pro <- data.table( SAMPLE_ISMMS=colnames(pro), PROTEIN=pro[pnm,] )
pro <- pro[grep("_1", SAMPLE_ISMMS)] #keep one replicate per sample
pro[,SAMPLE_ISMMS:=gsub("_1", "", SAMPLE_ISMMS)]

## SUBSET - samples in both rna and protein data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sam <- intersect(pro$SAMPLE_ISMMS, met$SAMPLE_ISMMS)
vob <- vob[,sam]
met <- met[SAMPLE_ISMMS %in% sam]
pro <- pro[SAMPLE_ISMMS %in% sam]

## MAKE - de variable ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
met <- merge(met, pro, by="SAMPLE_ISMMS")
met[mymet_postmortem==1, DE:="PM"]
met[mymet_postmortem==0, DE:="LIV"]

## FORMAT - metadata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
met <- as.data.frame(met)
rownames(met) <- met$SAMPLE_ISMMS
met <- met[colnames(vob),]

## MAKE - de formulas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
form1 <- ~ DE + DE:PROTEIN - PROTEIN + (1|mymet_phe) + (1|mymet_sex) + (1|mymet_bank) +
            mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS +
            RNASeqMetrics_PCT_MRNA_BASES + (1|IID_ISMMS) + (1|mymet_depletionbatch) +
            InsertSizeMetrics_MEDIAN_INSERT_SIZE + AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR
form2 <- ~0 + PROTEIN*mymet_postmortem + (1|mymet_phe) + (1|mymet_sex) + (1|mymet_bank) +
            mymet_rin + neuronal + RNASeqMetrics_MEDIAN_3PRIME_BIAS +            
            RNASeqMetrics_PCT_MRNA_BASES + (1|IID_ISMMS) + (1|mymet_depletionbatch) +
            InsertSizeMetrics_MEDIAN_INSERT_SIZE + AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR

## RUN - de of protein level ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dfit <- dream( vob, form1, met, BPPARAM = MulticoreParam(5) )
lvProDe <- topTable( dfit, coef="DELIV:PROTEIN", number=nrow(vob) )
pmProDe <- topTable( dfit, coef="DEPM:PROTEIN", number=nrow(vob) )
lvProDe <- data.table( gene = rownames(lvProDe), lvProDe )
pmProDe <- data.table( gene = rownames(pmProDe), pmProDe )
lvPi1 <- 1-pi0est(lvProDe$P.Value)$pi0
pmPi1 <- 1-pi0est(pmProDe$P.Value)$pi0
lpMer <- merge(lvProDe, pmProDe, by="gene")
lpCor <- cor(lpMer$logFC.x, lpMer$logFC.y, method="spearman")
lvDeg <- nrow(lvProDe[])

## RUN - de of protein level interaction with livpm status ~~~~~~~~~~~~~~~~~~~~~
dfit <- dream( vob, form2, met, BPPARAM = MulticoreParam(5) )
inProDe <- topTable( dfit, coef="PROTEIN:mymet_postmortem", number=nrow(vob) )
inAgeDe <- data.table( gene = rownames(intAgeDe), intAgeDe )
inPi1 <- 1-pi0est(inProDe$P.Value)$pi0

## SAVE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pn2 <- gsub(".", "_", pnm, fixed=T)
odr <- "/sc/arion/projects/psychgen/lbp/results/rna_protein_de/"
out <- paste0(odr, pn2, ".RDS")
saveRDS(list("LIV_DE"=lvProDe,
             "PM_DE"=pmProDe,
             "INTERACTION_DE"=intProDe,
             "LIV_PI1"=lvPi1,
             "PM_PI1"=pmPi1,
             "INTERACTION_PI1"=inPi1,
             "LIVPM_RHO"=lpCor), file=out)
