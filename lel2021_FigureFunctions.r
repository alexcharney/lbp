library(data.table)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(scales)
library(rlang)
library(assertthat)
library(dplyr)
library(magrittr)
library(stringr)
library(grid)
library(gridExtra)
library(egg)
library(cowplot)
library(patchwork)
library(seriation)

## fishers function for module-module correlations between networks (not used here, but its what lora used to make cen_subsets_overlap.RDS)
fishModSim <- function(cen.list, cen1, cen2) {
    out <- c()
    c1 <- cen.list[[cen1]]
    c2 <- cen.list[[cen2]]
    c1[,gene:=tstrsplit(gene, split=".", fixed=T, keep=1L)]
    c2[,gene:=tstrsplit(gene, split=".", fixed=T, keep=1L)]
    c3 <- merge(c1, c2, by="gene", suffixes=c(".cen1", ".cen2"))
    mod1 <- unique(c1$module)
    mod2 <- unique(c2$module)
    mod.com <- expand.grid( mod1, mod2 )
    for (i in 1:nrow(mod.com)){
        cur <- c3[,.(gene, CEN1=FALSE, CEN2=FALSE)]
        m1 <- as.character(mod.com[i,1])
        m2 <- as.character(mod.com[i,2])
        g1 <- c1[module==m1]$gene
        g2 <- c2[module==m2]$gene
        cur[gene %in% g1, CEN1:=TRUE]
        cur[gene %in% g2, CEN2:=TRUE]
        fish <- fisher.test(table(cur$CEN1, cur$CEN2), alternative="greater")
        res <- data.table("cen1"=cen1, "cen2"=cen2, "mod1"=m1, "mod2"=m2, "mod1.size"=length(g1), "mod2.size"=length(g2),
                         "overlap"=length(intersect(g1, g2)), "fisher.estimate"=fish$estimate, "fisher.p"=fish$p.value)
        out <- rbind(out, res)
    }
    out
}

## function to find consensus modules from a list of networks
findConsensusModules <- function(stats, refnet, nets, cenlist, N=ceiling(length(nets)/2)){
    out <- c()
    refmod <- unique(stats[indexCEN==refnet]$indexMod)
    othnets <- nets[nets != refnet]
    for (mod in refmod){
        cur <- stats[indexCEN==refnet & indexMod==mod & compareCEN %in% othnets & fisher.p.adjust<0.05]
        if (nrow(cur) == length(othnets)){
            refvec <- c(mod, cur$compareMod)
            names(refvec) <- c(refnet, cur$compareCEN)
            refvec <- refvec[order(names(refvec))]
            check <- 1
            for (othmod in 1:nrow(cur)){
                oNet <- cur[othmod]$compareCEN
                oMod <- cur[othmod]$compareMod
                othnets2 <- nets[nets != oNet]
                oCur <- stats[indexCEN==oNet & indexMod==oMod & compareCEN %in% othnets2 & fisher.p.adjust<0.05]
                othvec <- c(oMod, oCur$compareMod)
                names(othvec) <- c(oNet, oCur$compareCEN)
                othvec <- othvec[order(names(othvec))]
                if (identical(refvec, othvec)) check <- check+1
            }
            if (check==length(nets)) {
                add <- as.data.table(t(as.data.frame(refvec)))
                out <- rbind(out, add)
            }
        }
    }
    out <- out[rowSums(out == "turquoise" | out == "blue") %in% c(0,ncol(out))]
    out2 <- copy(out)
    for (i in 1:ncol(out)){
        set <- colnames(out)[i]
        for(j in 1:nrow(out)){
            mod <- out[[set]][j]
            siz <- unique(stats[indexCEN==set & indexMod==mod]$indexModSize)
            new <- paste(mod,siz,sep="|")
            out2[[set]][j] <- new
        }
    }
    list("meta"=out, "sizes"=out2)
}

## function to find consensus modules in Group A network(s) that overlap with grey modules in Group B network(s)
findConsensusModulesVsGrey <- function(stats, colorCEN, greyCEN, cenlist){
    statsNoGrey <- stats[indexMod!="grey" & compareMod!="grey"]
    if (length(colorCEN) > 1 & length(greyCEN)==1){
        colorMods <- findConsensusModules(stats = statsNoGrey, cenlist=cenlist, nets = colorCEN, refnet=colorCEN[1])[[1]]
        aMod <- tstrsplit(colorMods[[colorCEN[1]]], split="|", fixed=T, keep=1L)[[1]]
        aCen <- colorCEN[1]
        gMod <- stats[indexCEN==aCen & indexMod %in% aMod & compareCEN==greyCEN & compareMod=="grey"]$indexMod
        ret <- colorMods[ get(aCen) %in% gMod ]
    }
    if (length(colorCEN) == 1 & length(greyCEN) > 1){
        ret <- stats[indexCEN==colorCEN & indexMod!="grey" & compareCEN %in% greyCEN & compareMod=="grey"]
        ret <- ret[,.N,list(indexMod)][N==length(greyCEN)]$indexMod
    }
    ret
}

## function to find intersection of genes in conserved modules
findConsensusModuleGenes <- function(conmod, cenlist){
    glist <- c()
    nnet <- ncol(conmod)
    out <- c()
    for (i in 1:nrow(conmod)){
        nametracker <- c()
        for (j in 1:ncol(conmod)){
            net <- colnames(conmod)[j]
            mod <- conmod[[net]][i]
            name <- paste(net,mod,sep=".")
            if (j==1) nametracker <- name
            if (j!=1) nametracker <- paste(nametracker, name, sep="|")
            add <- cenlist[[net]][module==mod,.(gene,network=net)]
            glist <- rbind(glist,add)
        }
        out <- rbind(out, glist[,.N,gene][N==nnet][,.(gene,module=paste0("module",i),moduleOrigins=nametracker)])
    }
    out
}

## function for kegg
runKegg <- function(dt, dename, klist){
    out <- c()
    count <- 1
    for (i in names(klist)){
        pct <- round(100*(count/length(klist)), 0)
        if (count %% 100 == 0 ) cat(dename, pct,"%\n")
        count <- count+1
        ##cat(i,'\n')
        dt[,kegg:=0]
        dt[gene %in% klist[[i]],kegg:=1]
        n1 <- nrow(dt[deg==1])
        n2 <- nrow(dt[kegg==1])
        if (n2>0){
            n3 <- nrow(dt[deg==1 & kegg==1])
            ft <- fisher.test(table(dt$deg, dt$kegg), alternative="greater")
            or <- ft$estimate
            pv <- ft$p.value
            add <- data.table(awcid=i, ndeg=n1, nkegg=n2, nintersect=n3, or=or, pval=pv)
            out <- rbind(out, add)
        }
    }
    out
}

## ma plot function
do_maplot <- function(tt, fdr_threshold = 0.05, num_top_genes = 40, label_col = "symbol") {
    assert_that(is_string(label_col))
    tt <- tt %>% select(all_of(label_col), gene, logFC, z.std, AveExpr, adj.P.Val, de.status)
    names(tt)[[1]] <- "Gene_Label"
    y_limit <- max(abs(tt$logFC)) * 1.0 #highest logFC by abs value
    top_genes <- rbind( head(tt[!is.na(Gene_Label)][logFC<0][order(adj.P.Val)], num_top_genes),
                       head(tt[!is.na(Gene_Label)][logFC>0][order(adj.P.Val)], num_top_genes))  %>%
        {
            bind_rows(head(., num_top_genes/2),
                      tail(., num_top_genes/2))
        } %>%
        distinct %>%
        mutate(
            nudge_y = sign(logFC) * median(abs(logFC)) * 0.2,
            nudge_x = AveExpr %>% subtract(mean(.)) %>% multiply_by(0.2)
        )
    top_up_genes <- top_genes %>% filter(logFC > 0)
    top_down_genes <- top_genes %>% filter(logFC < 0)
    all_genes <- tt %>%
        select(-Gene_Label) %>%
        arrange(AveExpr, abs(logFC))
    ggplot(all_genes) +
        aes(x = AveExpr, y = logFC) +
        geom_point(pch=1, size=3, alpha=0.4, aes(col=de.status)) +
        scale_color_manual(values=c("#d80f8c", "#2cace2", "#999999")) +
    geom_density2d(show.legend = FALSE, alpha = 0.7) + #the "density map" aka the squiggles
        geom_smooth(show.legend = FALSE) +
        scale_size(
            name = "FDR (-log10)",
            limits = c(0, 10), range = c(0, 0.5)) +
        coord_cartesian(ylim = c(-1,1) * y_limit) +
        theme_bw() +
        theme(legend.position = "none")
}
