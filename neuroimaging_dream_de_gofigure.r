## Performs GO annotation of the genes of interest
## inputs are mel, a 2 column data.frame with number of rows corresponding
## to your gene sample size, where the first column is gene name and the second
## column is a binary vector of 0s and 1s, with 0s for genes in background
## (not significantly differentially expressed genes for example) and 1s
## for genes of interest (significantly differentially expressed genes for
## example); sample_name is a string containing the name of the file that
## will be saved with the output of the analysis; and date, a string of the
## date of the experiment that will also be used in the name of the output
## file. The function has no output, output is directly saved as a table.

## setup
options(stringsAsFactors = FALSE)
library(ggplot2)
library(data.table)
library("R.matlab")
library(goseq)
library(topGO)
library(org.Hs.eg.db)
library(Rgraphviz)
setwd("/sc/arion/projects/psychgen/lbp/data/neuroimaging/GMWM19JUL2021")
args <- commandArgs(trailingOnly=TRUE)
myFeature <- args[[1]]
myInput <- args[[2]]
myOutput <- args[[3]]
##myFeature <- "GMF199"
##myInput <- "/sc/arion/projects/psychgen2/lbp/results/neuroimaging/awcDreamTest/GMF199.tsv"
##myOutput <- "/sc/arion/projects/psychgen2/lbp/results/neuroimaging/awcDreamTest/GMF199.GO.tsv"

## read in data
dt <- fread(myInput)
dt[,gene:=tstrsplit(gene, split=".", fixed=T, keep=1L)]
    
## make "mel" objects
netlelUpSig <- dt[,.(gene,deg=0)]
netlelDwSig <- dt[,.(gene,deg=0)]
netlelUpNom <- dt[,.(gene,deg=0)]
netlelDwNom <- dt[,.(gene,deg=0)]
netleoUpSig <- dt[,.(gene,deg=0)]
netleoDwSig <- dt[,.(gene,deg=0)]
netleoUpNom <- dt[,.(gene,deg=0)]
netleoDwNom <- dt[,.(gene,deg=0)]
lelUpSigG <- dt[logFC.lel>0 & adj.P.Val.lel<0.05]$gene
lelDwSigG <- dt[logFC.lel<0 & adj.P.Val.lel<0.05]$gene
leoUpSigG <- dt[logFC.leo>0 & adj.P.Val.leo<0.05]$gene
leoDwSigG <- dt[logFC.leo<0 & adj.P.Val.leo<0.05]$gene
lelUpNomG <- dt[logFC.lel>0 & P.Value.lel<0.05]$gene
lelDwNomG <- dt[logFC.lel<0 & P.Value.lel<0.05]$gene
leoUpNomG <- dt[logFC.leo>0 & P.Value.leo<0.05]$gene
leoDwNomG <- dt[logFC.leo<0 & P.Value.leo<0.05]$gene
netlelUpSig[gene %in% lelUpSigG,deg:=1]
netlelDwSig[gene %in% lelDwSigG,deg:=1]
netlelUpNom[gene %in% lelUpNomG,deg:=1]
netlelDwNom[gene %in% lelDwNomG,deg:=1]
netleoUpSig[gene %in% leoUpSigG,deg:=1]
netleoDwSig[gene %in% leoDwSigG,deg:=1]
netleoUpNom[gene %in% leoUpNomG,deg:=1]
netleoDwNom[gene %in% leoDwNomG,deg:=1]
netlelUpSig <- as.data.frame(netlelUpSig)
netlelDwSig <- as.data.frame(netlelDwSig)
netlelUpNom <- as.data.frame(netlelUpNom)
netlelDwNom <- as.data.frame(netlelDwNom)
netleoUpSig <- as.data.frame(netleoUpSig)
netleoDwSig <- as.data.frame(netleoDwSig)
netleoUpNom <- as.data.frame(netleoUpNom)
netleoDwNom <- as.data.frame(netleoDwNom)
myList <- list("lel" = list("sig"=list("up"=netlelUpSig,"down"=netlelDwSig),
                      "nom"=list("up"=netlelUpNom,"down"=netlelDwNom)),
              "leo" = list("sig"=list("up"=netleoUpSig,"down"=netleoDwSig),
                           "nom"=list("up"=netleoUpNom,"down"=netleoDwNom)))

## gene map
gene.map  <- getgo(netlelUpSig[,1],'hg38','ensGene')
gene.map <- gene.map[!is.na(names(gene.map))] # removes genes in our data with no GO terms
number_of_GO <- length(unique(unlist(gene.map)))
counter <- 0
counter <- counter + 1
flush.console()

## run go enrichment
myResults <- c()
for (i in c("lel", "leo")){
    for (j in c("sig", "nom")){
        for (k in c("up", "down")){
            cur <- myList[[i]][[j]][[k]]
            a <- cur$deg
            if (length(table(a))>1){
                names(a) <- cur$gene
                bp <- new("topGOdata", description = "GO_BP", ontology = c("BP"),
                         allGenes = as.factor(a), geneSel = names(a[a==1]),
                         nodeSize = 1, annot = annFUN.gene2GO, gene2GO = gene.map)
                mf <- new("topGOdata", description = "GO_MF", ontology = c("MF"),
                         allGenes = as.factor(a), geneSel = names(a[a==1]),
                         nodeSize = 1, annot = annFUN.gene2GO, gene2GO = gene.map)
                cc <- new("topGOdata", description = "GO_CC", ontology = c("CC"),
                         allGenes = as.factor(a), geneSel = names(a[a==1]),
                         nodeSize = 1, annot = annFUN.gene2GO, gene2GO = gene.map)
                bp.res <- runTest(bp, algorithm = "classic", statistic = "fisher")
                mf.res <- runTest(mf, algorithm = "classic", statistic = "fisher")
                cc.res <- runTest(cc, algorithm = "classic", statistic = "fisher")
                bp.fin <- data.table(imagingFeature=myFeature, ontology="goBP", recipe=i, degDef=j, degDir=k, GenTable(bp, classic=bp.res, topNodes=length(bp.res@score)))
                mf.fin <- data.table(imagingFeature=myFeature, ontology="goMF", recipe=i, degDef=j, degDir=k, GenTable(mf, classic=mf.res, topNodes=length(mf.res@score)))          
                cc.fin <- data.table(imagingFeature=myFeature, ontology="goCC", recipe=i, degDef=j, degDir=k, GenTable(cc, classic=cc.res, topNodes=length(cc.res@score)))
                add <- rbind(bp.fin, mf.fin, cc.fin)
                add[,fold_enrichment:=Significant/Expected]
                add[,BH:=p.adjust(classic,method="BH", n=number_of_GO)]
                myResults <- rbind(myResults, add)
            }
        }
    }
}
fwrite(myResults, sep='\t', quo=F, row=F, file=myOutput)
