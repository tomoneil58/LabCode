
data <- read.delim(list.files()[7])
meta <- read.csv(list.files()[12])

## VOLCANO PLOT ####
ENTREZID <- mapIds(Homo.sapiens, keys = as.character(rownames(data)), column = c("ENTREZID"), 
                   keytype = "SYMBOL", multiVals = "first")
genes <- data.frame(ENTREZID, SYMBOL = rownames(data))

chooseGenes <- which(!duplicated(genes$SYMBOL) & !duplicated(genes$ENTREZID))

genes <- genes[chooseGenes, ]
Data <- data[chooseGenes, ]

group <- meta$group %>% factor()
individual <- meta$donor%>% factor()

samples <- data.frame(group, individual)

dgekumar <- DGEList(counts = Data, samples = samples, genes = genes)
dgekumar <- calcNormFactors(dgekumar, method = "TMM")
dgekumar$samples$norm.factors

cpm <- cpm(dgekumar)
# now calculate the log of these CPM values
lcpm <- cpm(dgekumar, log = TRUE)

list <- list(description, dgekumar, Data, cpm, lcpm, meta)
saveRDS(list, "KumarListofObjects.rds")

#func
cpmBox <- function(gene, heading = "Bulk RNA Expression (cpm)") {
  temp <- meta
  temp$gene <- cpm[gene,]
  ggplot(temp, aes(x=group, y=gene, colour=tcell)) + ggtitle(heading)+geom_boxplot()+geom_point()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()+facet_wrap(~tissue)
}

lcpmBox <- function(gene, heading = "Bulk RNA Expression (lcpm) - N=4") {
  temp <- meta
  temp$gene <- lcpm[gene,]
  ggplot(temp, aes(x=group, y=gene, colour=tcell)) + ggtitle(heading)+geom_boxplot()+geom_point()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()+facet_wrap(~tissue)
}

getCSV <- function(data, geneList = c(), name="") {
  geneDF <- meta
  cpm <- data
  plots <- list()
  for (i in 1:length(geneList)) { 
    geneDF[[geneList[i]]] <- cpm[geneList[[i]],]
  }
  date <- date()
  mon <- substr(date, 5,7)
  day <- substr(date, 9,10)
  year <- substr(date, 21,24)
  write.csv(geneDF, paste0("selectedGenes_", name,"_", year,mon,day, ".csv"))
  pdf(paste0("Plots_", year,mon,day, ".pdf"))
  for (i in 1:length(geneList)) {
    plot <- cpmBox(geneList[i])
    print(plot)
  }
  dev.off()
}

### volcano ####
design <- model.matrix(~ group)
design

v <- voom(dgekumar, design, plot = TRUE)

fit <- lmFit(v, design)
efit <- eBayes(fit)

dt_fdr <- decideTests(efit, adjust.method = "fdr", p.value = 0.05)
summary(dt_fdr)
dt_bonf <- decideTests(efit, adjust.method = "bonferroni", p.value = 0.05)
summary(dt_bonf)

hist(efit$p.value[, "groupCD8CD69pos"])
volcanoplot(efit, coef = "groupCD8CD69pos", highlight = 200, names = efit$genes$SYMBOL)

test <- data.frame(Genes = efit$genes$SYMBOL)
test$pvalue <- efit$p.value
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



## using top by log2fc
num = 20
threshx =10
threshy =100

tempmeans <- means
tempmeans <- tempmeans[tempmeans$ident.1<threshx,]
tempmeans <- tempmeans[tempmeans$ident.2>threshy,]
tempmeans$diff <- tempmeans$ident.1-tempmeans$ident.2

topgroup1 <- tempmeans %>%  top_n(n = num, wt = log2fc)
topgroup1 <- topgroup1[order(- topgroup1$log2fc),]

topgroup2 <- tempmeans %>%  top_n(n = -num, wt = log2fc)
topgroup2 <- topgroup2[order(topgroup2$log2fc),]

xlab = 'Lung CD4+ CD69+'
ylab = 'Lung CD4+ CD69-'

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

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[meta$tcell]

gene = c("KLF2", "KLF6")
hlcpm <- lcpm[gene,]
meta$heatmap <- c("S4N", 
                  "S4P", 
                  "S8N", 
                  "S8P",
                  "L4N", 
                  "L4P",
                  "L8N", 
                  "L8P", 
                  "S4N", 
                  "S4P",
                  "S8N", 
                  "S8P",
                  "L4N", 
                  "L4P", 
                  "L8N", 
                  "L8P",
                  "S4N", 
                  "S4P",
                  "S8N", 
                  "S8P", 
                  "L4N", 
                  "L4P",
                  "L8N", 
                  "L8P",
                  "B4N", 
                  "B8P", 
                  "B4N", 
                  "B8P",
                  "B4N", 
                  "B8P")
colnames(hlcpm) <- meta$heatmap

# Plot the heatmap
heatmap.2(hlcpm,col=rev(morecols(50)),trace="none", main="Genes title",ColSideColors=col.cell,scale="row", labCol=group, annCol=meta[, 2:4])

