BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
BiocManager::install("edgeR")
library(edgeR)
BiocManager::install("Homo.sapiens")
library(Homo.sapiens)

dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()                                                         # Tests that this is correct
PrimDir <- getwd()

data <- read.delim(list.files()[3])
meta <- read.csv(list.files()[6])
description <-  c("Data was normalised using TMM and cpm calculated. Bulk data function uses bulk data in list.")

data <- list(description, cpm #calculated from below
             , meta)

bulk1 <- meta
gene = 'CD40LG'
bulk1$gene <- cpm[gene,]
ggplot(bulk1, aes(x=resident, y=gene, colour=tissue)) + geom_boxplot() +ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()+facet_wrap(~tcell)

description <-  c("Data was normalised using TMM and cpm calculated. Bulk data function uses bulk data in list.")

# storing data
data <- list(description, cpm, meta)

ENTREZID2 <- mapIds(Homo.sapiens, keys = as.character(rownames(data)), column = c("ENTREZID"), 
                   keytype = "SYMBOL", multiVals = "first")
genes2 <- data.frame(ENTREZID2, SYMBOL = rownames(data))

chooseGenes2 <- which(!duplicated(genes2$SYMBOL) & !duplicated(genes2$ENTREZID))

genes2 <- genes2[chooseGenes2, ]
Data2 <- data[chooseGenes2, ]

counts2 <- Data2
rownames(counts2) <- genes2$SYMBOL

group <- meta$group[1:30] %>% factor()
individual <- meta$donor [1:30]%>% factor()

samples2 <- data.frame(group, individual)

dgekumar <- DGEList(counts = counts2, samples = samples2, genes = genes2)
dgekumar <- calcNormFactors(dgekumar, method = "TMM")
dgekumar$samples$norm.factors
keep <- filterByExpr(dgekumar)
dgekumar <- dgekumar[keep, ]

cpm <- cpm(dgekumar) # use this for function above
# now calculate the log of these CPM values
lcpm <- cpm(dgekumar, log = TRUE)

design <- model.matrix(~ group)
design

v <- voom(dgekumar, design, plot = TRUE)

fit <- lmFit(v, design)
efit <- eBayes(fit)

dt_fdr <- decideTests(efit, adjust.method = "fdr", p.value = 0.05)
summary(dt_fdr)
dt_bonf <- decideTests(efit, adjust.method = "bonferroni", p.value = 0.05)
summary(dt_bonf)

hist(efit$p.value[, "groupCD69pos"])
volcanoplot(efit, coef = "groupCD69pos", highlight = 100, names = efit$genes$SYMBOL)

## Log2FC ####

### GROUP 1 ####
#tissue 
## tissue and blood
ref <- grep("Spleen", meta$tissue)
ref <- c(ref, grep("Lung", meta$tissue))
ref <- c(ref, grep("Blood", meta$tissue))
## tissue only
ref <- grep("Spleen", meta$tissue)
ref <- c(ref, grep("Lung", meta$tissue))
## lung
ref <- grep("Lung", meta$tissue)

## spleen
ref <- grep("Spleen", meta$tissue)

## blood
ref <- grep("Blood", meta$tissue)

#temp 
temp <- cpm[,ref]
tempmeta <- meta[ref,]

#resident
## CD69-
ref <- grep("CD69neg", tempmeta$resident)

## CD69+
ref <- grep("CD69pos", tempmeta$resident)

## both
ref <- grep("CD69neg", tempmeta$resident)
ref <- c(ref, grep("CD69pos", tempmeta$resident))

#temp
temp <- temp[,ref]
tempmeta <- tempmeta[ref,]

#Tcell subset
## CD4
ref <- grep("CD4", tempmeta$tcell)

## CD8
ref <- grep("CD8", tempmeta$tcell)

## Both 
ref <- grep("CD4", tempmeta$tcell)
ref <- c(ref, grep("CD8", tempmeta$tcell))

temp <- temp[,ref]
tempmeta <- tempmeta[ref,]
#
group1 <- temp
group1m <- tempmeta

### GROUP 2 ####
#tissue 
## tissue and blood
ref <- grep("Spleen", meta$tissue)
ref <- c(ref, grep("Lung", meta$tissue))
ref <- c(ref, grep("Blood", meta$tissue))
## tissue only
ref <- grep("Spleen", meta$tissue)
ref <- c(ref, grep("Lung", meta$tissue))
## lung
ref <- grep("Lung", meta$tissue)

## spleen
ref <- grep("Spleen", meta$tissue)

## blood
ref <- grep("Blood", meta$tissue)

#temp 
temp <- cpm[,ref]
tempmeta <- meta[ref,]

#resident
## CD69-
ref <- grep("CD69neg", tempmeta$resident)

## CD69+
ref <- grep("CD69pos", tempmeta$resident)

## both
ref <- grep("CD69neg", tempmeta$resident)
ref <- c(ref, grep("CD69pos", tempmeta$resident))

#temp
temp <- temp[,ref]
tempmeta <- tempmeta[ref,]

#Tcell subset
## CD4
ref <- grep("CD4", tempmeta$tcell)

## CD8
ref <- grep("CD8", tempmeta$tcell)

## Both 
ref <- grep("CD4", tempmeta$tcell)
ref <- c(ref, grep("CD8", tempmeta$tcell))

temp <- temp[,ref]
tempmeta <- tempmeta[ref,]
#
group2 <- temp
group2m <- tempmeta

#means
means <- data.frame(Genes = rownames(cpm), ident.1 = rowMeans(group1), ident.2 = rowMeans(group2))
for (i in 1:nrow(means)) {
  
  a <- log2(means[i,"ident.1"])
  b <- log2(means[i, "ident.2"])
  means$ident.1.log2[i] <-a 
  means$ident.2.log2[i] <-b
  means$log2fc[i] <- b-a
}


for (i in 1:nrow(means)) {
  if (means$log2fc[i] == 'Inf' | means$log2fc[i] == '-Inf') {
    means$log2fc[i] <- NA
  } 
}
means <- na.omit(means)
#filter out ribosomal genes
means<-means[-grep("RPL", rownames(means)),]
means<-means[-grep("RPS", rownames(means)),]
means[grep("^RPS", rownames(means)),]



## using top by log2fc
num = 20
threshx = 10000
threshy = 10000

tempmeans <- means
tempmeans <- tempmeans[tempmeans$ident.1>threshx & tempmeans$ident.2>threshy,]
topgroup1 <- tempmeans %>%  top_n(n = num, wt = log2fc)
topgroup1 <- topgroup1[order(- topgroup1$log2fc),]

topgroup2 <- tempmeans %>%  top_n(n = -num, wt = log2fc)
topgroup2 <- topgroup2[order(topgroup2$log2fc),]

xlab = 'CD4+ CD69-'
ylab = 'CD4+CD69+'

p2<-ggplot(means, aes(x=ident.1, y=ident.2)) + geom_point() + ggtitle("Log2FC")+xlab(xlab) + ylab(ylab)+geom_abline(slope=1)+theme_classic()


LabelPoints(plot = p2, points = c(topgroup1$Genes, topgroup2$Genes), repel = TRUE)

## using top by log2fc
num = 20
threshx = 20000
threshy = 20000

tempmeans <- means
tempmeans <- tempmeans[tempmeans$ident.1>threshx & tempmeans$ident.2>threshy,]


xlab = 'CD4+ CD69-'
ylab = 'CD4+CD69+'

p2<-ggplot(means, aes(x=ident.1, y=ident.2)) + geom_point() + ggtitle("Log2FC")+xlab(xlab) + ylab(ylab)+geom_abline(slope=1)+theme_classic()

LabelPoints(plot = p2, points = c(tempmeans$Genes), repel = TRUE)


names <-c("CD4+CD69- Spleen", #1
          "CD4+CD69+ Spleen", 
          "CD8+CD69- Spleen", 
          "CD8+CD69+ Spleen", 
          "CD4+CD69- Lung", #5
          "CD4+CD69+ Lung", 
          "CD8+CD69- Lung", 
          "CD8+CD69+ Lung", 
          "CD4+CD69- Spleen2", 
          "CD4+CD69+ Spleen2", #10
          "CD8+CD69- Spleen2", 
          "CD8+CD69+ Spleen2", 
          "CD4+CD69- Lung2", 
          "CD4+CD69+ Lung2", 
          "CD8+CD69- Lung2", #15
          "CD8+CD69+ Lung2",  
          "CD4+CD69- Spleen3", #1
          "CD4+CD69+ Spleen3", 
          "CD8+CD69- Spleen3", 
          "CD8+CD69+ Spleen3", 
          "CD4+CD69- Lung3", #5
          "CD4+CD69+ Lung3", 
          "CD8+CD69- Lung3", 
          "CD8+CD69+ Lung3", 
          "CD4+CD69- Blood1", #25
          "CD8+CD69- Blood1", 
          "CD4+CD69- Blood2", 
          "CD8+CD69- Blood2", 
          "CD4+CD69- Blood3", 
          "CD8+CD69- Blood3" #30
)


