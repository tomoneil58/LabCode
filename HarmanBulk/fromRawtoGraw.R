data <- read.csv("GSE130804_AHcounts.csv.gz")

rownames(data) <- make.names(data$X, unique=T)
data <- data[, -1]
data <-data[,order(colnames(data))]
colnames(data)
datatest <- data.frame(X = rownames(data))

for (i in 1:nrow(data)) {
  datatest$AH1[i] <- sum(data[i, 1:4]) 
}

for (i in 1:nrow(data)) {
  datatest$AH10[i] <- sum(data[i, 5:8]) 
}

for (i in 1:nrow(data)) {
  datatest$AH11[i] <- sum(data[i, 9:12]) 
}

for (i in 1:nrow(data)) {
  datatest$AH12[i] <- sum(data[i, 13:16]) 
}

for (i in 1:nrow(data)) {
  datatest$AH13[i] <- sum(data[i, 17:20]) 
}

for (i in 1:nrow(data)) {
  datatest$AH14[i] <- sum(data[i, 21:24]) 
}
for (i in 1:nrow(data)) {
  datatest$AH15[i] <- sum(data[i, 25:28]) 
}
for (i in 1:nrow(data)) {
  datatest$AH16[i] <- sum(data[i, 29:32]) 
}
for (i in 1:nrow(data)) {
  datatest$AH17[i] <- sum(data[i, 33:36]) 
}
for (i in 1:nrow(data)) {
  datatest$AH18[i] <- sum(data[i, 37:40]) 
}
for (i in 1:nrow(data)) {
  datatest$AH19[i] <- sum(data[i, 41:44]) 
}
for (i in 1:nrow(data)) {
  datatest$AH2[i] <- sum(data[i, 45:48]) 
}
for (i in 1:nrow(data)) {
  datatest$AH20[i] <- sum(data[i, 49:52]) 
}
for (i in 1:nrow(data)) {
  datatest$AH21[i] <- sum(data[i, 53:56]) 
}
for (i in 1:nrow(data)) {
  datatest$AH22[i] <- sum(data[i, 57:60]) 
}
for (i in 1:nrow(data)) {
  datatest$AH23[i] <- sum(data[i, 61:64]) 
}
for (i in 1:nrow(data)) {
  datatest$AH24[i] <- sum(data[i, 65:68]) 
}
for (i in 1:nrow(data)) {
  datatest$AH25[i] <- sum(data[i, 69:72]) 
}
for (i in 1:nrow(data)) {
  datatest$AH26[i] <- sum(data[i, 73:76]) 
}
for (i in 1:nrow(data)) {
  datatest$AH27[i] <- sum(data[i, 77:80]) 
}
for (i in 1:nrow(data)) {
  datatest$AH28[i] <- sum(data[i, 81:84]) 
}
for (i in 1:nrow(data)) {
  datatest$AH29[i] <- sum(data[i, 85:88]) 
}
for (i in 1:nrow(data)) {
  datatest$AH3[i] <- sum(data[i, 89:92]) 
}
for (i in 1:nrow(data)) {
  datatest$AH30[i] <- sum(data[i, 93:96]) 
}
for (i in 1:nrow(data)) {
  datatest$AH31[i] <- sum(data[i, 97:100]) 
}
for (i in 1:nrow(data)) {
  datatest$AH4[i] <- sum(data[i, 101:104]) 
}
for (i in 1:nrow(data)) {
  datatest$AH5[i] <- sum(data[i, 105:108]) 
}
for (i in 1:nrow(data)) {
  datatest$AH6[i] <- sum(data[i, 109:112]) 
}
for (i in 1:nrow(data)) {
  datatest$AH7[i] <- sum(data[i, 113:116]) 
}
for (i in 1:nrow(data)) {
  datatest$AH8[i] <- sum(data[i, 117:120]) 
}
for (i in 1:nrow(data)) {
  datatest$AH9[i] <- sum(data[i, 121:124]) 
}

#test1
num <- c()
for (i in 1:(ncol(data)/4)) {
  a = i*4
  b = a-3
  s = sum(data[1, b:a])
  num <- c(num, s)
}
num 
datatest[1,]

#test2
num <- c()
for (i in 1:(ncol(data)/4)) {
  a = i*4
  b = a-3
  num <- c(num, sum(data[345, b:a]))
}
num 
datatest[345,]

datasummed <-datatest
rownames(datasummed) <- datasummed$X
datasummed <- datasummed[,-1]

samples <- data.frame(label=colnames(datasummed))

samples$individual <- c("A", "B", "B", "B", "B", #AH1, 10, 11, 12, 13
                        "B", "B", "C", "C", "C", #AH14, 15, 16, 17, 18
                        "C", "A", "C", "C", "C", #AH19, 2, 20, 21, 22
                        "C", "D", "D", "D", "D", #AH23, 24, 25, 26, 27
                        "D", "D", "A", "D", "D", #AH28, 29, 3, 30, 31
                        "A", "A", "A", "A", "B", "B" #AH4,5, 6, 7, 8, 9
)
samples$group <- c("Mac", "CD33Low", "CDC2", "Mac", "MDM", 
                        "CDC1", "Lang", "LC", "CD11c", "CD33Low", 
                        "CDC2", "MDM", "Mac", "MDM", "Lang", 
                        "CDC1", "LC", "CD11c", "CD33Low", "CDC2", 
                        "Mac", "MDM", "Lang", "CDC1", "Lang", 
                        "CDC1", "CD33Low", "CD11c", "LC", "LC", "CD11c"
)

samples$label2 <- paste0(samples$group, "_", samples$individual)

colnames(datasummed) <- samples$label2  

saveRDS(datasummed, "counts_lanesSummed.rds" )
saveRDS(samples, "metadata_samples.rds" )


ENTREZID <- mapIds(Homo.sapiens, keys = as.character(rownames(datasummed)), column = c("ENTREZID"), 
                   keytype = "SYMBOL", multiVals = "first")
genes <- data.frame(ENTREZID, SYMBOL =rownames(datasummed))
genes <- na.omit(genes)

chooseGenes <- which(!duplicated(genes$SYMBOL) & !duplicated(genes$ENTREZID))

genes <- genes[chooseGenes, ]
datasummed <- datasummed[chooseGenes, ]

samplesx <- samples %>% arrange(group)


BLCL <- DGEList(counts = datasummed, samples = samples, genes = genes)
BLCL <- calcNormFactors(BLCL, method = "TMM")
BLCL$samples$norm.factors

cpm <- cpm(BLCL)

bulk1 <- samples
gene = 'ITGAX'
bulk1$gene <- cpm[gene,]
ggplot(bulk1, aes(x=group, y=gene, colour=group)) + geom_boxplot()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()

#group 1 
ref <- grep("CDC1", samples$group)

#temp

group1 <- cpm[,ref]
group1m <- samples[ref,]

#group 1 
ref <- grep("CDC2", samples$group)

#temp
group2 <- cpm[,ref]
group2m <- samples[ref,]

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
means[grep("^MT.ND", rownames(means)),]



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

xlab = 'CDC1'
ylab = 'CDC2'

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