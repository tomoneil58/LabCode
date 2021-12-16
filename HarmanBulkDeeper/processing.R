
setwd("~/Desktop/DC bulk_new/")

meta <- read.csv('meta.csv')
data <- read.csv('data.csv')
genes <- read.csv('genes.csv')

counts <- data[,-c(1:4)]
rownames(counts) <- make.names(data$Gene.Symbol, unique=T)

ENTREZID <- mapIds(Homo.sapiens, keys = as.character(rownames(counts)), column = c("ENTREZID"), 
                    keytype = "SYMBOL", multiVals = "first")
genes <- data.frame(ENTREZID, SYMBOL = rownames(counts))

chooseGenes <- which(!duplicated(genes$SYMBOL) & !duplicated(genes$ENTREZID))

genes <- genes[chooseGenes, ]
counts <- counts[chooseGenes, ]

rownames(counts) <- genes$SYMBOL

group <- meta$group %>% factor()
individual <- meta$ident%>% factor()

samples <- data.frame(group, individual)

dgeDCs <- DGEList(counts = counts, samples = samples, genes = genes)
dgeDCs <- calcNormFactors(dgeDCs, method = "TMM")
dgeDCs$samples$norm.factors

cpm <- cpm(dgeDCs)
lcpm <- cpm(dgeDCs, log=T)


bulk1 <- meta
gene = 'FOSB'
bulk1$gene <- lcpm[gene,]
ggplot(bulk1, aes(x=group, y=gene, colour=tissue)) + geom_boxplot() +ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()

description <- c("This is the most recent bulkRNA data given to me by Andrew. First this description, then the raw counts, dge, cpm, lcpm and meta")
list <- list(description,counts, dgeDCs,cpm, lcpm, meta)
saveRDS(list, "ListofObject_BulkNewDCData.rds")


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
  LabelPoints(plot = p2, points = c(topgroup1$Genes, topgroup2$Genes), repel = TRUE, xnudge=0,ynudge=0)
  
}

meta <- list[[6]]
counts <- list[[2]]

ref <- grep("cDC2", meta$group)
counts <- counts[,ref]
meta <- meta[ref,]
group <- meta$group %>% factor()

dgeDCs <- DGEList(counts = counts, samples = meta, genes = genes)
dgeDCs <- calcNormFactors(dgeDCs, method = "TMM")
dgeDCs$samples$norm.factors
cpm <- cpm(dgeDCs)
lcpm <- cpm(dgeDCs, log = TRUE)


# Perform PCA
pca <- prcomp(t(lcpm), scale = TRUE)
# Plot results
autoplot(pca, data = meta, colour = "group")+theme_classic()

# Another method for dimension reduction is called MDS.
plotMDS(lcpm, labels = group)
plotMDS(lcpm, labels = individual)

# I also like this cool interactive version of a MDS plot which can be generated
# using Glimma
glMDSPlot(lcpm, labels = paste(group, individual, sep = "_"), groups = BLCL$samples[, 
                                                                                    c("group", "individual")], launch = TRUE)

design <- model.matrix(~ group)
design

v <- voom(dgeDCs, design, plot = TRUE)

fit <- lmFit(v, design)
efit <- eBayes(fit)

### EXPLORE ####

dt_fdr <- decideTests(efit, adjust.method = "fdr", p.value = 0.05)
summary(dt_fdr)

# The coefficient we are interested in is 'groupLCL'.  Does the number of DE
# genes change if we instead use the more stringent bonferroni correction?
dt_bonf <- decideTests(efit, adjust.method = "bonferroni", p.value = 0.05)
summary(dt_bonf)

# histogram of p value distribution for groupLCL colnames(efit)
hist(efit$p.value[, "groupCDC1"])


# volcano plot of results for groupLCL
volcanoplot(efit, coef = "groupcDC2LangP", highlight = 100, names = efit$genes$SYMBOL)


# genes for LangP v Langn

ref <- grep("cDC2", meta$group)
cpm <- cpm[,ref]
meta <- meta[ref,]

pdf("cDC2Genes.pdf")
cpmBox("CD207")
cpmBox("LAMP3")
cpmBox("BIRC3")
cpmBox("IRF4")
cpmBox("POLR2A")
cpmBox("STAT5A")
cpmBox("TNFRSF9")
cpmBox("ITGA1")
cpmBox("IL1B")
cpmBox("IL10")
cpmBox("IER3")
cpmBox("ISG20")
cpmBox("CCL3")
cpmBox("CXCL2")
cpmBox("C3AR1")
dev.off()


