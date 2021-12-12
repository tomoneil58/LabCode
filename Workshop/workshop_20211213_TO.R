
library(Seurat)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)

download.file("https://github.com/tomoneil58/LabCode/raw/main/Workshop/listForHarmanBulk.rds", "ListForHarmanBulk.rds")
harmanList <- readRDS("ListForHarmanBulk.rds")

cpm <- harmanList[[4]]
lcpm <- harmanList[[5]]
meta <- harmanList[[6]]

cpmBox <- function(gene, heading = "Bulk RNA Expression (cpm) - N=4") {
  temp <- meta
  temp$gene <- cpm[gene,]
  ggplot(temp, aes(x=group, y=gene, colour=group)) + ggtitle(heading)+geom_boxplot()+geom_point()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()
}

lcpmBox <- function(gene, heading = "Bulk RNA Expression (lcpm) - N=4") {
  temp <- meta
  temp$gene <- lcpm[gene,]
  ggplot(temp, aes(x=group, y=gene, colour=group)) + ggtitle(heading)+geom_boxplot()+geom_point()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()
}

topDiff <- function(ident1, ident2, numgenes = 25, heading = "Average Gene Expressions", xthresh=50000, ythresh=50000, lowxthresh = 1, lowythresh=1) {
  ref <- grep(ident1, meta$group)
  group1 <- cpm[,ref]
  group1m <- meta[ref,]
  xlab <- ident1
  
  ref <- grep(ident2, meta$group)
  group2 <- cpm[,ref]
  group2m <- meta[ref,]
  ylab <- ident2
  
  means <- data.frame(Genes = rownames(cpm), ident.1 = rowMeans(group1), ident.2 = rowMeans(group2))
  
  #get Log2 and Log2FC
  for (i in 1:nrow(means)) {
    
    a <- log2(means[i,"ident.1"])
    b <- log2(means[i, "ident.2"])
    means$ident.1.log2[i] <-a 
    means$ident.2.log2[i] <-b
    means$log2fc[i] <- b-a
  }
  
  #filter shit genes - Inf and -Inf from 0s
  for (i in 1:nrow(means)) {
    if (means$log2fc[i] == 'Inf' | means$log2fc[i] == '-Inf') {
      means$log2fc[i] <- NA
    } 
  }
  means <- na.omit(means)
  
  #filter threshold
  means <- means[means$ident.1<xthresh,]
  means <- means[means$ident.2<ythresh,]
  
  ## using top by log2fc
  tempmeans <- means
  tempmeans <- tempmeans[tempmeans$ident.1>lowxthresh & tempmeans$ident.2>lowythresh,]
  
  topgroup1 <- tempmeans %>%  top_n(n = numgenes, wt = log2fc)
  topgroup1 <- topgroup1[order(- topgroup1$log2fc),]
  
  topgroup2 <- tempmeans %>%  top_n(n = -numgenes, wt = log2fc)
  topgroup2 <- topgroup2[order(topgroup2$log2fc),]
  p2<-ggplot(means, aes(x=ident.1, y=ident.2)) + geom_point() + ggtitle(heading)+xlab(xlab) + ylab(ylab)+geom_abline(slope=1)+theme_classic()
  LabelPoints(plot = p2, points = c(topgroup1$Genes, topgroup2$Genes), repel = TRUE)
  
}




