## note: I sanity checked that the "_" in the eqtl sample IDs
##       here are only the ones plink added (ie, miami dont have "_")

## SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(variancePartition))
setwd("/sc/arion/projects/psychgen/lbp/results/runEQTL/20OCT2022")

## ARGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
args <- commandArgs(trailingOnly=TRUE)
gIx1 <- as.integer(args[[1]]) ##gIx1 <- 1
gIx2 <- as.integer(args[[2]]) ##gIx2 <- 10
fOut <- args[[3]]

## MAIN PROTEIN DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pro <- readRDS("../../LBP_LIVPM_PROTEIN_DE_INPUT_DATA_AND_RESULTS_USING_LELQC_19OCT2022.RDS")
met <- as.data.table(pro$meta)
vob <- pro$expression
met[mymet_postmortem==1, IID_ISMMS:=gsub("_", "-", IID_ISMMS)]
rownames(vob) <- unlist(tstrsplit(rownames(vob), split=".", fixed=T, keep=1L))
vob <- vob[pro$livpmDE[!is.na(gene)]$refseq,]
rownames(vob) <- pro$livpmDE[!is.na(gene)]$gene

## SUBSET PROTEIN DATA FOR SAMPLES ANALYZED IN EQTL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
erd <- readRDS("/sc/arion/projects/psychgen/lbp/data/PEER/LIVPM_20OCT2022_RNA_PRO_RNAEDITING.RDS")
emt <- erd$protein$peerSampleList$sid
kpMet <- met[SAMPLE_ISMMS %in% emt]
kpVob <- vob[,emt]

## PQTL RESULTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qtl <- readRDS("RESULTS_UPDATED24OCT2022.RDS")
eql <- qtl$fullResults$PROTEIN$nsv30$QC2CALLSET
eql <- eql[!is.na(IDTopVar)]
gLst <- eql$gene[gIx1:gIx2]

## PQTL GENOTYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gty <- fread("RESULTS_UPDATED24OCT2022_PROTEIN_QC2CALLSET_nsv30_IDTopVar.traw", na=c("NA", ".", "", "-9"))
gty[,`(C)M`:=NULL]
gty <- merge(eql[,.(gene, SNP=IDTopVar)], gty)

## DEFINE VARIABLES FOR VP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frm <- ~1 + EQTL + mymet_postmortem
kpCol <- c("EQTL", "mymet_postmortem")
nwCol <- c("gene", "eqtl_ve", "livpm_ve")

## RUN VP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output <- c()
for (myGen in gLst){

    print(myGen)

    ##genotypes for the eqtl for this gene
    lvSnp <- gty[gene==myGen]
    lvVar <- lvSnp$SNP
    lvAdd <- data.table( IID_ISMMS = colnames(lvSnp)[7:ncol(lvSnp)], EQTL=unlist(lvSnp[,7:ncol(lvSnp)]) )
    lvAdd[,IID_ISMMS:=tstrsplit(IID_ISMMS, split="_", fixed=T, keep=1L)]
    lvMet <- merge(kpMet, lvAdd, by="IID_ISMMS")[!is.na(EQTL)]
    lvMet <- as.data.frame(lvMet)
    rownames(lvMet) <- lvMet$SAMPLE_ISMMS
    check <- nrow(lvMet)

    if (check>0){

        ##subset vobject for individuals with nonmissing genotypes
        lvVob <- kpVob[myGen,rownames(lvMet),drop=F]

        ##run variancePartition
        vpLiv <- fitExtractVarPartModel( lvVob, frm, lvMet)
        vpLvd <- as.data.table(vpLiv[,kpCol], keep.rownames="gene")
        colnames(vpLvd) <- nwCol
        vpLvd$nLIV <- sum(lvMet$mymet_postmortem==0)
        vpLvd$nPM <- sum(lvMet$mymet_postmortem==1)

        ##format metadata on the eqtl for this gene from the eqtl analysis
        lvAdd <- lvSnp[,.(gene, eqtl_source="LIVPM", eqtl_snp=SNP, eqtl_counted=COUNTED, eqtl_notcounted=ALT)]
        lvAd2 <- eql[IDTopVar==lvVar,.(gene, eqtl_snp=IDTopVar, eqtl_chr=chr, eqtl_tss=tss, eqtl_strand=strand,
                                      eqtl_nvar=nvar, eqtl_distToTopVar=distToTopVar,
                                      eqtl_nominalPVal=nominalPVal, eqtl_regressionSlope=regressionSlope,
                                      eqtl_storey=st)]
        lvAdd <- merge(lvAdd, lvAd2, by=c("gene", "eqtl_snp"))

        ##add eqtl metadata to variancePartition
        vpLvd <- merge(vpLvd, lvAdd, by="gene")

        ##add indicator to show what samples are subject of the result
        vpLvd <- data.table(sampleset="LIVPM", vpLvd)

        ##add to output
        output <- rbind(output, vpLvd)
        
    }
}

## ADD DE SUMSTATS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add <- pro$livpmDE
add[,LIVPM:="NONE"]
add[DEG=="DEG" & LFC=="NEGLFC",LIVPM:="LIV"]
add[DEG=="DEG" & LFC=="POSLFC",LIVPM:="PM"]
add <- add[,.(gene, livpmDE_DEG=LIVPM, livpmDE_logFC=logFC, livpmDE_padj=adj.P.Val,
              livpmDE_AveExpr=AveExpr, livpmDE_computedAvgExp=computedAvgExp,
              livpmDE_computedAvgExpLiv=mean_liv, livpmDE_computedAvgExpPm=mean_pm)]
add <- add[!is.na(gene)]
output <- merge(output, add, by="gene")

## REORDER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
myOrd <- c("gene", "sampleset", "nLIV", "nPM",
           "eqtl_source", "eqtl_snp", "eqtl_counted", "eqtl_notcounted", "eqtl_chr",
           "eqtl_tss", "eqtl_strand", "eqtl_nvar", "eqtl_distToTopVar",
           "eqtl_nominalPVal", "eqtl_regressionSlope", "eqtl_storey",
           "livpmDE_DEG", "livpmDE_logFC", "livpmDE_padj", "livpmDE_AveExpr",
           "livpmDE_computedAvgExp", "livpmDE_computedAvgExpLiv", "livpmDE_computedAvgExpPm",
          "eqtl_ve", "livpm_ve")
output <- output[,c(myOrd),with=F]

## WRITE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(output, file=fOut)
