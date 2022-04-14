## R/4.0.3
##
## README: performs DE using dream on pseudobulk data from scrnaseq experiment

## setup
options(stringsAsFactors=F)
suppressMessages(library(Seurat))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(data.table))
suppressMessages(library(edgeR))
suppressMessages(library(variancePartition))
suppressMessages(library(qvalue))
Sys.setenv(OMP_NUM_THREADS = 6)
args <- commandArgs(trailingOnly=TRUE)
myDat <- args[[1]]
myHlp <- args[[2]]
myIdx <- as.integer(as.character(args[[3]]))
myOut <- args[[4]]
print(myDat)
print(myHlp)
print(myIdx)
print(myOut)
##myDat <- "/sc/arion/projects/psychgen/lbp/data/rna/brbl_scrnaseq_olink_pseudobulk_de_data_CellType.RDS"
##myHlp <- "/sc/arion/projects/psychgen/lbp/data/rna/brbl_scrnaseq_olink_pseudobulk_de_data_CellType_HELPER.tsv"
##myIdx <- "1"
##myOut <- "/sc/arion/projects/psychgen/lbp/results/brbl_scrnaseq_olink_pseudobulk_de/tracker_index_1"

## read in helper
myCel <- fread(myHlp, sep="\t", header=F)[myIdx]$V1
myPro <- fread(myHlp, sep="\t", header=F)[myIdx]$V2

## read in data
myDat <- readRDS(myDat)
myExp <- myDat$data
myOlk <- myDat$olink

## formula
form <- ~ protein + (1|iid)

## process data
y <- DGEList(assay(myExp[[myCel]]), samples=colData(myExp[[myCel]]), remove.zeros=TRUE)
y <- calcNormFactors(y, method="TMM")
y$samples$iid <- tstrsplit(rownames(y$samples), split="|", fixed=T, keep=1L)[[1]]
colnames(y) <- gsub("|brain", "", fixed=TRUE, colnames(y))
colnames(y) <- gsub("|blood", "", fixed=TRUE, colnames(y))
add <- myOlk[colnames(y),myPro,drop=F]
y <- y[,rownames(add)]
y$samples$protein <- add[,1]
keep <- filterByExpr(y, min.count=10)
vobjDream <- voomWithDreamWeights(y[keep,], form, y$samples, BPPARAM = MulticoreParam(5) )

## run DREAM
dreamResults <- dream(vobjDream, form, y$samples, BPPARAM = MulticoreParam(5))

## format de results
lmgroup_DE <- dreamResults
coefcol <- "protein"
Group_DE_tab <- topTable(lmgroup_DE, coef=coefcol, number=nrow(vobjDream))
de <- data.table( gene = rownames(Group_DE_tab), Group_DE_tab)
de <- de[order(logFC)]
de[adj.P.Val<0.05, DEG:="DEG"]
de[adj.P.Val>0.05, DEG:="NOTDEG"]
de[P.Value<0.05, NOMDEG:="NOMDEG"]
de[P.Value>0.05, NOMDEG:="NOTNOMDEG"]
de[logFC<0, LFC:="NEGLFC"]
de[logFC>0, LFC:="POSLFC"]

## calculate pi1 statistic
pi1 <- 1 - qvalue(de$P.Value)$pi0

## write
fout1 <- paste0(myOut, "_de.tsv") 
fout2 <- paste0(myOut, "_pi1.tsv")
fwrite(de, sep='\t', row=F, quo=F, na="NA", file=fout1)
fwrite(data.table(cell=myCel, protein=myPro, value=pi1), sep='\t', row=F, quo=F, col=F, na="NA", file=fout2)
