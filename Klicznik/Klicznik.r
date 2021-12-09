
#load normalised counts
klic <- as.data.frame(read.delim(list.files()[2]))

rownames(klic) <- make.names( klic[,1], unique=T)
klic <- klic[,-1]
klic <- klic[-1,]

for (i in 1:ncol(klic)){
klic[,i] <- as.numeric(unlist(klic[,i]))
}

meta <- read.csv("meta.csv")
meta <- meta[,-1]



ENTREZID <- mapIds(Homo.sapiens, keys = as.character(rownames(klic)), column = c("ENTREZID"), 
                   keytype = "SYMBOL", multiVals = "first")
genes <- data.frame(ENTREZID, SYMBOL =rownames(klic))
genes <- na.omit(genes)

chooseGenes <- which(!duplicated(genes$SYMBOL) & !duplicated(genes$ENTREZID))

genes <- genes[chooseGenes, ]
klic <- klic[chooseGenes, ]

samples <- data.frame(samples = meta$sample)
samples$group <- meta$group
samples$individual <- meta$donor

klicdge <- DGEList(counts = klic, samples = samples, genes = genes)
klicdge <- calcNormFactors(klicdge, method = "TMM")
klicdge$samples$norm.factors

cpm <- cpm(klicdge)
lcpm <- cpm(klicdge, log=T)



bulk1 <- meta
#lungCD69+vsCD69- CXCR4+TXNIP+TSC22D3+DUSP1+KLF6+FOS+  TPT1-AHNAK-EEF1A1-EEF2- 
gene = c('GATA3') #LGALS3
bulk1$gene <- lcpm[gene,]
ggplot(bulk1, aes(x=group, y=gene, colour=group)) + geom_boxplot() +ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()

skinref <- grep("skin", meta$tissue)

cpmskin <- cpm[,skinref]
lcpmskin <- lcpm[,skinref]
skinmeta <- meta[skinref,]

bulkskin <- skinmeta
#lungCD69+vsCD69- CXCR4+TXNIP+TSC22D3+DUSP1+KLF6+FOS+  TPT1-AHNAK-EEF1A1-EEF2- 
gene = c('TBXAS1') #LGALS3
bulkskin$gene <- cpmskin[gene,]
ggplot(bulkskin, aes(x=group, y=gene, colour=group)) + geom_boxplot() +ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()
description <- c("This list has this description, raw counts, cpm, lcpm, dge, and meta objects.")
klicznik <- c(description, klic, cpm, lcpm, klicdge, meta)
saveRDS(klicznik, "20211210_klicznikDataList.rds")

### LOG2FC ####

### GROUP 1 ####
#tissue 
## skin and blood
ref <- grep("blood", meta$tissue)
ref <- c(ref, grep("skin", meta$tissue))
## skin only
ref <- grep("skin", meta$tissue)

## blood
ref <- grep("blood", meta$tissue)

#temp 
temp <- cpm[,ref]
tempmeta <- meta[ref,]

#resident
## CD103+
ref <- grep("cd103p", tempmeta$cd103)

## CD103-
ref <- grep("cd103n", tempmeta$cd103)

## both
ref <- grep("cd103p", tempmeta$cd103)
ref <- c(ref, grep("cd103n", tempmeta$cd103))

#temp
temp <- temp[,ref]
tempmeta <- tempmeta[ref,]

#ccr7 subset
## ccr7+
ref <- grep("ccr7p", tempmeta$ccr7)

## ccr7-
ref <- grep("ccr7n", tempmeta$ccr7)

## Both 
ref <- grep("ccr7p", tempmeta$ccr7)
ref <- c(ref, grep("ccr7n", tempmeta$ccr7))

temp <- temp[,ref]
tempmeta <- tempmeta[ref,]
#
group1 <- temp
group1m <- tempmeta
xlab = 'skin CD103+'

### GROUP 2 ####
#tissue 
## skin and blood
ref <- grep("blood", meta$tissue)
ref <- c(ref, grep("skin", meta$tissue))
## skin only
ref <- grep("skin", meta$tissue)

## blood
ref <- grep("blood", meta$tissue)

#temp 
temp <- cpm[,ref]
tempmeta <- meta[ref,]

#resident
## CD103+
ref <- grep("cd103p", tempmeta$cd103)

## CD103-
ref <- grep("cd103n", tempmeta$cd103)

## both
ref <- grep("cd103p", tempmeta$cd103)
ref <- c(ref, grep("cd103n", tempmeta$cd103))

#temp
temp <- temp[,ref]
tempmeta <- tempmeta[ref,]

#ccr7 subset
## ccr7+
ref <- grep("ccr7p", tempmeta$ccr7)

## ccr7-
ref <- grep("ccr7n", tempmeta$ccr7)

## Both 
ref <- grep("ccr7p", tempmeta$ccr7)
ref <- c(ref, grep("ccr7n", tempmeta$ccr7))

temp <- temp[,ref]
tempmeta <- tempmeta[ref,]
#
group2 <- temp
group2m <- tempmeta
ylab = 'skin CD103-'

#### means ####
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
means<-means[-grep("MT.CO", rownames(means)),]
means<-means[-grep("MT.ND", rownames(means)),]



## using top by log2fc
num = 20
threshx =10
threshy =10

tempmeans <- means
tempmeans <- tempmeans[tempmeans$ident.1>threshx,]
tempmeans <- tempmeans[tempmeans$ident.2>threshy,]
tempmeans$diff <- tempmeans$ident.1-tempmeans$ident.2

topgroup1 <- tempmeans %>%  top_n(n = num, wt = log2fc)
topgroup1 <- topgroup1[order(- topgroup1$log2fc),]

topgroup2 <- tempmeans %>%  top_n(n = -num, wt = log2fc)
topgroup2 <- topgroup2[order(topgroup2$log2fc),]


p2<-ggplot(means, aes(x=ident.1.log2, y=ident.2.log2)) + geom_point() + ggtitle("Log2FC")+xlab(xlab) + ylab(ylab)+geom_abline(slope=1)+theme_classic()
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




