### downsampling LIV to one samnple per person and downsampling PM to be same sample size as LIV

## SETUP
cat("MESSAGE 1 | Setting up environment\n")
rm(list=ls())
options(stringsAsFactors=F)
suppressMessages(library(edgeR))
suppressMessages(library(withr))
suppressMessages(library(limma))
suppressMessages(library(variancePartition))
suppressMessages(library(data.table))
suppressMessages(library(qvalue))
args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[[1]]) #idx <- 1
setwd("/sc/arion/projects/psychgen/lbp/data/proteomics")
source("/sc/arion/work/charna02/scripts/lbp/pi0estFunctionUnbounded.r")

## READ - rnaseq data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("MESSAGE 2 | reading in rnaseq data\n")
rna <- readRDS("/sc/arion/work/charna02/symlinks/lbp/liharska2021/final.everything.RDS")
met <- rna$covariates
vob <- rna$vobjDream ##verified elsewhere using vobject from full model vs variation on model doesnt matter
met[,mymet_bank:=as.character(mymet_bank)]

## READ - protein data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("MESSAGE 3 | reading in protein data\n")
prt <- readRDS("./LBP_Proteomics_fromWeipingMa_toLbpTeam_01APR2022/Data/20220214_MSSM_Living_Brain_Tissue_Global_34MPs_Combined_awcParsed_update26JUN2022.RDS")
pro <- prt$imputed_26JUN2022
pcv <- prt$covariates
pnm <- rownames(pro)[idx]
pro <- data.table( SAMPLE_ISMMS=colnames(pro), PROTEIN=pro[pnm,] )
pro <- pro[grep("_1", SAMPLE_ISMMS)] #keep one replicate per sample
pro[,SAMPLE_ISMMS:=gsub("_1", "", SAMPLE_ISMMS)]
pn2 <- gsub(".", "_", pnm, fixed=T)
odr <- "/sc/arion/projects/psychgen/lbp/scratch/TMP2/"
out <- paste0(odr, pn2, ".RDS")

## DOWNSAMPLE - pre-determined set of LIV and PM samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("MESSAGE 4 | downsampling to pre-determined set of LIV and PM samples\n")
sam <- readRDS("/sc/arion/projects/psychgen/lbp/scratch/TMP_n318_for_rnaprotein_associativity_troubleshooting.RDS")$SAMPLE_ISMMS
met <- met[SAMPLE_ISMMS %in% sam]
vob <- vob[,sam]
pro <- pro[SAMPLE_ISMMS %in% sam]

## MAKE - livpm variable ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("MESSAGE 5 | making livpm variable\n")
met <- merge(met, pro, by="SAMPLE_ISMMS")
met[mymet_postmortem==1, DE:="PM"]
met[mymet_postmortem==0, DE:="LIV"]

## FORMAT - metadata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("MESSAGE 6 | formatting metadata\n")
met <- as.data.frame(met)
rownames(met) <- met$SAMPLE_ISMMS
met <- met[colnames(vob),]

## MAKE - de formulas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("MESSAGE 7 | making de formulas\n")
allCovariates <- c("mymet_phe", "mymet_sex", "mymet_bank","mymet_rin","neuronal",
                  "RNASeqMetrics_MEDIAN_3PRIME_BIAS", "RNASeqMetrics_PCT_MRNA_BASES",
                  "mymet_depletionbatch", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
                  "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR")
lelCovariates <- c("mymet_sex", "mymet_rin", "neuronal",
                  "RNASeqMetrics_MEDIAN_3PRIME_BIAS", "RNASeqMetrics_PCT_MRNA_BASES",
                  "mymet_depletionbatch", "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
                  "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR")
myCovariates <- list("none" = "", "all" = allCovariates, "lel" = lelCovariates,
                  "LOI.mymet_phe" = "mymet_phe",
                  "LOI.mymet_sex" = "mymet_sex",
                  "LOI.mymet_bank" = "mymet_bank",
                  "LOI.mymet_rin" = "mymet_rin",
                  "LOI.neuronal" = "neuronal",
                  "LOI.RNASeqMetrics_MEDIAN_3PRIME_BIAS" = "RNASeqMetrics_MEDIAN_3PRIME_BIAS",
                  "LOI.RNASeqMetrics_PCT_MRNA_BASES" = "RNASeqMetrics_PCT_MRNA_BASES",
                  "LOI.mymet_depletionbatch" = "mymet_depletionbatch",
                  "LOI.InsertSizeMetrics_MEDIAN_INSERT_SIZE" = "InsertSizeMetrics_MEDIAN_INSERT_SIZE",
                  "LOI.AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR" = "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR",
                  "LOO.mymet_phe" = setdiff(allCovariates, "mymet_phe"),
                  "LOO.mymet_sex" = setdiff(allCovariates, "mymet_sex"),
                  "LOO.mymet_bank" = setdiff(allCovariates, "mymet_bank"),
                  "LOO.mymet_rin" = setdiff(allCovariates, "mymet_rin"),
                  "LOO.neuronal" = setdiff(allCovariates, "neuronal"),
                  "LOO.RNASeqMetrics_MEDIAN_3PRIME_BIAS" = setdiff(allCovariates, "RNASeqMetrics_MEDIAN_3PRIME_BIAS"),
                  "LOO.RNASeqMetrics_PCT_MRNA_BASES" = setdiff(allCovariates, "RNASeqMetrics_PCT_MRNA_BASES"),
                  "LOO.mymet_depletionbatch" = setdiff(allCovariates, "mymet_depletionbatch"), 
                  "LOO.InsertSizeMetrics_MEDIAN_INSERT_SIZE" = setdiff(allCovariates, "InsertSizeMetrics_MEDIAN_INSERT_SIZE"),
                  "LOO.AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR" = setdiff(allCovariates, "AlignmentSummaryMetrics_STRAND_BALANCE_FIRST_OF_PAIR"))

## RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 8 | running DE of protein levels\n")
myOut <- c()
for (i in names(myCovariates)){
    cat("             ", i, "\n")
    mycov.covlist = myCovariates[[i]]
    if (i=="none"){
        mycov.formula <- as.formula("~DE + DE:PROTEIN - PROTEIN")
        mycov.alllist <- c("DE", "PROTEIN")
    } else {        
        mycov.formula <- as.formula(paste("~DE + DE:PROTEIN - PROTEIN +",paste(mycov.covlist,collapse="+")))
        mycov.alllist <- c("DE", "PROTEIN", mycov.covlist)
    }
    mycov <- met[,mycov.alllist, drop=F]
    myexp <- vob[,rownames(mycov)]
    mydesign <- model.matrix(mycov.formula, mycov)
    colnames(mydesign) <- make.names(colnames(mydesign))
    mylmgroup <- lmFit(myexp, mydesign)
    mylmgroup_DE <- eBayes(mylmgroup)
    livDe <- topTable(mylmgroup_DE, coef="DELIV.PROTEIN", number=nrow(mylmgroup_DE))
    pmtDe <- topTable(mylmgroup_DE, coef="DEPM.PROTEIN", number=nrow(mylmgroup_DE))   
    livDe <- data.table( gene = rownames(livDe), livDe )
    pmtDe <- data.table( gene = rownames(pmtDe), pmtDe )
    lvPi1 <- 1-pi0est(livDe$P.Value)$pi0
    pmPi1 <- 1-pi0est(pmtDe$P.Value)$pi0
    lvPi1ub <- 1-pi0estUnbounded(livDe$P.Value)$pi0
    pmPi1ub <- 1-pi0estUnbounded(pmtDe$P.Value)$pi0
    lpMer <- merge(livDe, pmtDe, by="gene")
    lpCor <- cor(lpMer$logFC.x, lpMer$logFC.y, method="spearman")    
    add <- data.table("protein"=pnm, "model"=i,
                     "LIV_PI1"=lvPi1, "LIV_PI1UB"=lvPi1ub,
                     "PM_PI1"=pmPi1, "PM_PI1UB"=pmPi1ub, "LIVPM_RHO"=lpCor)
    myOut <- rbind(myOut, add)
}
    
## SAVE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("MESSAGE 11 | saving results to", out, "\n")
saveRDS(myOut, file=out)
cat("SUCCESS\n")
