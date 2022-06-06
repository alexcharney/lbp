   
## Load Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
options(stringsAsFactors=F)
library(data.table)
library(IHW)

## Read Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
args <- commandArgs(trailingOnly=TRUE)
pvl <- readRDS(args[[1]])
pad <- readRDS(args[[2]])
out <- args[[3]]

## Function for IHW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ihwWrapper <- function(dreamPvalues, dreamPvaluesAdj) {
    cat(date(), " | Making vectors for IHW\n")
    pvec <- c(dreamPvalues)
    apvec <- c(dreamPvaluesAdj)
    nvec <- unlist(lapply(colnames(dreamPvalues),rep,nrow(dreamPvalues)))
    cat(date(), " | Running IHW\n")
    myihw <- ihw(pvalues=pvec, covariates=as.factor(nvec), alpha=0.05, covariate_type="nominal")
    cat(date(), " | Formatting output of IHW\n")
    xx <- data.table(name=nvec, dreamP=pvec, dreamPADJ=apvec)
    yy <- as.data.table(myihw@df)
    colnames(yy)[1] <- "dreamP"
    yy[,dreamPADJ:=xx$dreamPADJ]
    yy[,covariate:=NULL]
    cat("   number of dream P.Value below 0.05 = ",  nrow(yy[dreamP<0.05]), "\n")
    cat("   number of dream adj.P.Val below 0.05 = ",  nrow(yy[dreamPADJ<0.05]), "\n")
    cat("   number of IHW adj_pvalue below 0.05 = ",  nrow(yy[adj_pvalue<0.05]), "\n")
    cat("   number of IHW adj_pvalue and dream adj.P.Val below 0.05 = ",  nrow(yy[adj_pvalue<0.05 & dreamPADJ<0.05]), "\n")
    cat(date(), " | Running sanity check for reformatting output of IHW into matrix\n")
    sanityMtx <- matrix(yy$dreamP,nrow=nrow(dreamPvalues)) #sanity check
    colnames(sanityMtx) <- colnames(dreamPvalues)
    rownames(sanityMtx) <- rownames(dreamPvalues)
    if (identical(sanityMtx, dreamPvalues)){
        cat("   SUCCESS\n")
    } else {
        stop("   ERROR\n")
    }
    cat(date(), " | Reformatting output of IHW into matrix\n")
    ihwMtx <- matrix(yy$adj_pvalue,nrow=nrow(pvl))
    colnames(ihwMtx) <- colnames(pvl)
    rownames(ihwMtx) <- rownames(pvl)
    ihwMtx
}

## Run IHW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ihwOut <- ihwWrapper(pvl, pad)
saveRDS(ihwOut, file=out)
 
