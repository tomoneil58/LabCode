
library(dplyr)
library(ggplot2)
library(Seurat)

proportions <- function(data, ident.1, ident.2, position) {
  x<- FetchData(data,c(ident.1,ident.2))
  colnames(x) <- c('ident.2', 'ident.1')
  x%>% group_by(ident.1) %>%
    mutate(prop=1/length(ident.2)) %>%
    ungroup() %>%
    group_by(ident.2,ident.1) %>%
    summarise(totprop=sum(prop)) %>%
    ggplot(aes(x=ident.2,fill=ident.1,y=totprop)) +  
    geom_bar(position=position, stat='identity') + theme(axis.text.x = 
                                                           element_text(angle = 45,hjust=1))+scale_y_continuous(name="Cluster
    Proportion")+ theme_classic()
}


setwd("C:/Users/thomas.oneil/Desktop/Luoma/raw")

#ignoring the csv TCR annotations

# took a bit of work to untar each folder - had to do most of it manually and havent included it here.
# first untarred the GSE144469_RAW, which made a huge list of tar files - had to then untar each of those
# the csv were just moved to another folder. 
# then each donor had a folder which had barcodes, features and the matrix -
# below is the data wrangling to read in each automatically, set rows as 
#  genes and columns as barcode (cell name) and load it into a list 


donor <- c()

setwd(prim)
donor <- c(donor, list.files()[1])
setwd(list.files()[1])
barcode <- read.table(list.files()[1])

features <- read.table(list.files()[2])
matrix <- as.data.frame(readMM(list.files()[3]))

rownames(matrix) <- make.names(features[,2], unique=T)
colnames(matrix) <- barcode[,1]

list <- list()
list[[1]] <- matrix
names(list)[1] <- donor[1]

#loop for the remaining folders! 
#do not repeat! It creates a 38+GB object ####
for(i in 2:38) {
  setwd(prim)
  donor <- c(donor, list.files()[i])
  setwd(list.files()[i])
  barcode <- read.table(list.files()[1])
  features <- read.table(list.files()[2])
  matrix <- as.data.frame(readMM(list.files()[3]))
  
  rownames(matrix) <- make.names(features[,2], unique=T)
  colnames(matrix) <- barcode[,1]
  list[[i]] <- matrix
  names(list)[i] <- donor[i]
  
}

#making the metadata for donor
BATCH <- c(rep(donor[1], dim(list[[1]])[2]),
           rep(donor[2], dim(list[[2]])[2]),
           rep(donor[3], dim(list[[3]])[2]),
           rep(donor[4], dim(list[[4]])[2]),
           rep(donor[5], dim(list[[5]])[2]),
           rep(donor[6], dim(list[[6]])[2]),
           rep(donor[7], dim(list[[7]])[2]),
           rep(donor[8], dim(list[[8]])[2]),
           rep(donor[9], dim(list[[9]])[2]),
           rep(donor[10], dim(list[[10]])[2]),
           rep(donor[11], dim(list[[11]])[2]),
           rep(donor[12], dim(list[[12]])[2]),
           rep(donor[13], dim(list[[13]])[2]),
           rep(donor[14], dim(list[[14]])[2]),
           rep(donor[15], dim(list[[15]])[2]),
           rep(donor[16], dim(list[[16]])[2]),
           rep(donor[17], dim(list[[17]])[2]),
           rep(donor[18], dim(list[[18]])[2]),
           rep(donor[19], dim(list[[19]])[2]),
           rep(donor[20], dim(list[[20]])[2]),
           rep(donor[21], dim(list[[21]])[2]),
           rep(donor[22], dim(list[[22]])[2]),
           rep(donor[23], dim(list[[23]])[2]),
           rep(donor[24], dim(list[[24]])[2]),
           rep(donor[25], dim(list[[25]])[2]),
           rep(donor[26], dim(list[[26]])[2]),
           rep(donor[27], dim(list[[27]])[2]),
           rep(donor[28], dim(list[[28]])[2]),
           rep(donor[29], dim(list[[29]])[2]),
           rep(donor[30], dim(list[[30]])[2]),
           rep(donor[31], dim(list[[31]])[2]),
           rep(donor[32], dim(list[[32]])[2]),
           rep(donor[33], dim(list[[33]])[2]),
           rep(donor[34], dim(list[[34]])[2]),
           rep(donor[35], dim(list[[35]])[2]),
           rep(donor[36], dim(list[[36]])[2]),
           rep(donor[37], dim(list[[37]])[2]),
           rep(donor[38], dim(list[[38]])[2])
)
names <- c(colnames(list[[1]]),
           colnames(list[[2]]),
           colnames(list[[3]]),
           colnames(list[[4]]),
           colnames(list[[5]]),
           colnames(list[[6]]),
           colnames(list[[7]]),
           colnames(list[[8]]),
           colnames(list[[9]]),
           colnames(list[[10]]),
           colnames(list[[11]]),
           colnames(list[[12]]),
           colnames(list[[13]]),
           colnames(list[[14]]),
           colnames(list[[15]]),
           colnames(list[[16]]),
           colnames(list[[17]]),
           colnames(list[[18]]),
           colnames(list[[19]]),
           colnames(list[[20]]),
           colnames(list[[21]]),
           colnames(list[[22]]),
           colnames(list[[23]]),
           colnames(list[[24]]),
           colnames(list[[25]]),
           colnames(list[[26]]),
           colnames(list[[27]]),
           colnames(list[[28]]),
           colnames(list[[29]]),
           colnames(list[[30]]),
           colnames(list[[31]]),
           colnames(list[[32]]),
           colnames(list[[33]]),
           colnames(list[[34]]),
           colnames(list[[35]]),
           colnames(list[[36]]),
           colnames(list[[37]]),
           colnames(list[[38]]))
meta <- data.frame(names, BATCH)

setwd("C:/Users/thomas.oneil/Desktop/Luoma/")
saveRDS(meta, '20210824_donormetadata.rds')
saveRDS(list, '20210824_listeddata.rds')

rm(matrix);rm(barcode);rm(features);rm(BATCH);rm(donor);rm(i); rm(names)

### Create seurat objects within the listed data #### 

for (i in 1:length(list)) {
  list[[i]] <- CreateSeuratObject(list[[i]])
}
saveRDS(list, '20210824_listasseurat.rds')

### QC: removing mt high ####

# JUST THE STANDARD REMOVING MT HIGH - hashed out is an example
# typically genes are labelled as MT-, but these are MT. - just clarifying, so we dont pick up unwanted genes here
genes <- rownames(list[[1]])[grep("^MT.ND", rownames(list[[1]]))]
genes <- c(genes, rownames(list[[1]])[grep("^MT.CO", rownames(list[[1]]))])
genes <- c(genes, rownames(list[[1]])[grep("^MT.CYB", rownames(list[[1]]))])
genes <- c(genes, rownames(list[[1]])[grep("^MT.ATP", rownames(list[[1]]))])

for (i in 1:length(list)) {
  #list[[i]][['percent.mt']] <- PercentageFeatureSet(list[[i]], features=genes)
  print(table(list[[i]]$percent.mt <15))
}
## Around no more than 400 cells per sample were >15% MT - will remove these

## Also checked the nFeature_RNA and want to remove those that are too low
for(i in 1:length(list)) {
  print(table(list[[i]]$nFeature_RNA <200))
}

## SUBSET THESE CELLS
for (i in 1:length(list)) {
  list[[i]]<- subset(list[[i]], subset = nFeature_RNA >200 & percent.mt<15)
  print(table(list[[i]]$percent.mt <15))
  print(table(list[[i]]$nFeature_RNA <200))
}

#healthy
ncx <- list[28:38]
#no-colitis
ctx <- list[15:27]
#colitis 
cx <- list[1:14]

list <- cx

### INTEGRATE DATA VIA rPCA #### 
#individual donors categorised above were combined because integrating all together was too difficult
# the aim is to cluster in smaller groups, pull out cells of interest and repeat this entirely
for (i in 1:length(list)){
  list[[i]]$donor <- names(list[i])
}

list <- lapply(X=list, FUN = function(x){
  x<-NormalizeData(x, verbose=F)
  x <- FindVariableFeatures(x, verbose=F)
})
features <- SelectIntegrationFeatures(object.list=list)
list <-lapply(X=list, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=F)
  x<-RunPCA(x, features=features, verbose=F)
})
list <- FindIntegrationAnchors(object.list=list, reduction='rpca', dims=1:30)
list <- IntegrateData(anchorset=list, dims=1:30)


tisCCA <- ScaleData(tisCCA, verbose=F)
tisCCA <- RunPCA(tisCCA, verbose=F)
data <- RunUMAP(data, dims=1:12)
data <- RunTSNE(data, dims=1:12)

data <- FindNeighbors(data, dims=1:12)
tisCCA <- FindClusters(tisCCA, resolution=0.5)

TSNEPlot(tisCCA, group.by='donor')

DefaultAssay(tisCCA) <- 'RNA'
tisCCA <- NormalizeData(tisCCA)
FeaturePlot(tisCCA, features=c("CD3E", "CD19", "FCER1A", "CD3D"))


#noncolitis 
saveRDS(tisCCA, "20210824_NonColDonors_full_scaled.rds")
#healthy
saveRDS(tisCCA, "20210824_HealthyDonors_full_scaled.rds")
#colitis 
saveRDS(tisCCA, "20210824_colitisDonors_full_scaled.rds")


#2022Oct
setwd("~/Desktop/Luoma/")
list.files()
data <- readRDS(list.files()[2])

DefaultAssay(data)<- 'integrated'

data <- RunUMAP(data, dims=1:12)
data <- RunTSNE(data, dims=1:12)

data <- FindNeighbors(data, dims=1:12)
data <- FindClusters(data, resolution=0.4)

DefaultAssay(data) <- 'RNA'

FeaturePlot(data, 'CD207')


data[['percent.mt']]<- PercentageFeatureSet(data, pattern='^MT.')
VlnPlot(data, 'percent.mt') + geom_hline(yintercept=15)

UMAPPlot(data,label=T) ->p1
FeaturePlot(data, label=T, "CD3E")

saveRDS(data, "20221030_Reclustered.rds")
pdf("selectionOfTcells.pdf")
FeaturePlot(data, label=T, "CD3E")
FeaturePlot(data, label=T, "percent.mt")
FeaturePlot(data, label=T, "MS4A1")
FeaturePlot(data, label=T, "JCHAIN")
FeaturePlot(data, label=T, "MKI67")
dev.off()

data <- subset(data, idents=c("2", "3", "6", "5", "4", "0", "9", "7", "14"))


# forget this - the data above is only of healthy donors. Instead, I'm going back and I'll run ALL donors from healthy > colitis. 
# start again at line 145
# nevermind R crashed

saveRDS(data, "20221030_Reclustered.rds")

#colitis data
data <- readRDS(list.files()[1])

DefaultAssay(data)<- 'integrated'

data <- RunUMAP(data, dims=1:12)
#data <- RunTSNE(data, dims=1:12)

data <- FindNeighbors(data, dims=1:12)
data <- FindClusters(data, resolution=0.3)

DefaultAssay(data) <- 'RNA'

FeaturePlot(data, 'CD207')

data <- subset(data, idents=c("4", "8", "6", "0", "2", "5", "1", "7"))
saveRDS(data, "20221030_ReclusteredColitis.rds")

# this wont work - instead, I'll subset the FOXP3 subsets + adjacent. 
#data <- merge(data, readRDS("20221030_Reclustered.rds"))

#read healthy
data1 <- subset(data, idents=c("5", "4")) #by FOXP3
data1$disease <- 'Healthy'
#read colitis
data2 <- subset(data, idents=c("5", "2")) #by FOXP3
data2$disease <- 'Colitis'

data <-merge(data1, data2); rm(data1); rm(data2)

data <- SplitObject(data, split.by='donor')

data <- lapply(X=data, FUN = function(x){
  DefaultAssay(x) <- 'RNA'
  x<-NormalizeData(x, verbose=F)
  x <- FindVariableFeatures(x, verbose=F)
})
features <- SelectIntegrationFeatures(object.list=data)
data <-lapply(X=data, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=F)
  x<-RunPCA(x, features=features, verbose=F)
})
data <- FindIntegrationAnchors(object.list=data, reduction='rpca', dims=1:30)
data <- IntegrateData(anchorset=data, dims=1:30)
saveRDS(data, 
        "colitishealthycombinedTreg.rds")

data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)
data <- RunUMAP(data, dims=1:15)
#data <- RunTSNE(data, dims=1:12)
saveRDS(data, 
        "colitishealthycombinedTreg.rds")
data <- FindNeighbors(data, dims=1:15)
data <- FindClusters(data, resolution=0.5)

DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
UMAPPlot(data, label=T)+UMAPPlot(data, label=T, group.by='disease')+UMAPPlot(data, label=T, group.by='donor')+FeaturePlot(data, "FOXP3")
UMAPPlot(data, label=T, pt.size=1.5)

DotPlot(data, features=c("KLRB1","CCR5", "IL10","CCR7","KLF6","LAG3","IL6R", "DUSP1", "CXCR6","CCR4","CCR6","RORC",
                         "IL7R", "TIGIT", "ENTPD1", "PDCD1", "ITGAE","TNFRSF18", "TNFRSF9", "CCR8", 
                         "CD27",  "AHR","RORA","CXCR4", "CXCR3", "KLF2", "HLA.DRB1", "HLA.DQB1", "TGFB1" ), idents=c("3", "4", "8"))+RotatedAxis()




markers <- FindAllMarkers(subset(data, idents=c("3", "4", "8")), only.pos=T)
markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top20
saveRDS(data, "Tregs.rds")

UMAPPlot(data, label=T, pt.size=1.5, split.by='disease')


data <- subset(data, idents=c("3", "4", "8"))

DefaultAssay(data) <- 'integrated'

data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)
data <- RunUMAP(data, dims=1:8)
data <- RunTSNE(data, dims=1:8)

UMAPPlot(data, label=T, pt.size=1.5)

data <- FindNeighbors(data, dims=1:8)
data <- FindClusters(data, resolution=0.15)

DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)

markers <- FindAllMarkers(data, only.pos=T)
markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top20

data <- RenameIdents(data, 
                     "0" = "0",
                     "1" = "1", 
                     "2" = "0", 
                     "3" = "2"
                     )

DotPlot(data[,data$disease == 'Healthy'], features=c("KLRB1","CCR5", "IL10","CCR7","KLF6","LAG3","IL6R", "DUSP1", "CXCR6","CCR4","CCR6","RORC",
                         "IL7R", "TIGIT", "ENTPD1", "PDCD1", "ITGAE","TNFRSF18", "TNFRSF9", "CCR8", 
                         "CD27",  "AHR","RORA","CXCR4", "CXCR3", "KLF2", "HLA.DRB1", "HLA.DQB1", "TGFB1" ))+RotatedAxis()

DotPlot(data, features=c(grep("^ITGB", rownames(data), value=T)))+RotatedAxis()

# add noncolitis
data <- readRDS("colitishealthycombinedTreg.rds")
x <- readRDS(list.files()[3]) #all full scaled
DefaultAssay(x) <- 'RNA'

FeaturePlot(x, "CD3E", label=T)
FeaturePlot(x, "CD8A", label=T)
FeaturePlot(x, "MKI67", label=T)

x <- subset(x, idents=c("9", "19", "12", "13", "3", "24", "3"))
x$disease <- 'Non-Colitis'
data <- merge(x,data)


data <- SplitObject(data, split.by='donor')

data <- lapply(X=data, FUN = function(x){
  DefaultAssay(x) <- 'RNA'
  x<-NormalizeData(x, verbose=F)
  x <- FindVariableFeatures(x, verbose=F)
})
features <- SelectIntegrationFeatures(object.list=data)
data <-lapply(X=data, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=F)
  x<-RunPCA(x, features=features, verbose=F)
})
data <- FindIntegrationAnchors(object.list=data, reduction='rpca', dims=1:30)
data <- IntegrateData(anchorset=data, dims=1:30, k.weight=90)
saveRDS(data, 
        "colitisNChealthycombinedTreg.rds")

data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)
data <- RunUMAP(data, dims=1:10)
UMAPPlot(data, label=T, group.by='donor')+NoLegend()+UMAPPlot(data, label=T, group.by='disease')

DefaultAssay(data) <- 'integrated'
data <- FindNeighbors(data, dims=1:10)
data <- FindClusters(data, resolution=0.3)
DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
FeaturePlot(data, "FOXP3", label=T)

genes <- rownames(data)[grep("^MT.ND", rownames(data))]
genes <- c(genes, rownames(data)[grep("^MT.CO", rownames(data))])
genes <- c(genes, rownames(data)[grep("^MT.CYB", rownames(data))])
genes <- c(genes, rownames(data)[grep("^MT.ATP", rownames(data))])

data[['percent.mt']] <- PercentageFeatureSet(data, features=genes)

saveRDS(data, 
        "colitisNChealthycombinedTreg.rds")

VlnPlot(data, "FOXP3")

data <-subset(data, idents=c("2", "3", "10"))

DefaultAssay(data) <- 'integrated'
data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)
data <- RunUMAP(data, dims=1:10)
data <- RunTSNE(data, dims=1:10)

UMAPPlot(data, label=T, group.by='donor')+NoLegend()+UMAPPlot(data, label=T, group.by='disease')
TSNEPlot(data, label=T, group.by='donor')+NoLegend()+TSNEPlot(data, label=T, group.by='disease')

DefaultAssay(data) <- 'integrated'
data <- FindNeighbors(data, dims=1:10)
data <- FindClusters(data, resolution=0.15)

UMAPPlot(data, label=T, pt.size=1.5)

DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
FeaturePlot(data, "FOXP3", label=T)
VlnPlot(data, "FOXP3")

saveRDS(data, "NovTregColNCHealthy.rds") ######***#####


datax <- subset(data, idents=c("0", "1", "2","4"))
### new from here - 
DefaultAssay(datax) <- 'integrated'
datax <- ScaleData(datax, verbose=F)
datax <- RunPCA(datax, verbose=F)
ElbowPlot(datax)
datax <- RunUMAP(datax, dims=1:8)
datax <- RunTSNE(datax, dims=1:8)
UMAPPlot(datax, label=T)
datax <- FindNeighbors(datax, dims=1:8)
datax <- FindClusters(datax, resolution=0.2)
UMAPPlot(datax, label=T, pt.size=1.5)

DefaultAssay(datax) <- 'RNA'
datax <- NormalizeData(datax)
FeaturePlot(datax, "FOXP3", label=T)
VlnPlot(datax, "FOXP3")

datax <- subset(datax, idents=c("0","1", '2', "3"))


proportions(datax, "disease", 'ident', "fill")
data@active.ident <- factor(data@active.ident, 
                            levels=c("3", "1", "0", "2"))
DotPlot(datax[,datax$disease=='Healthy'], features=c("KLRB1","IKZF3", "DUSP6", "CCR5", "IL10","CCR7","KLF6","LAG3","IL6R", "DUSP1", "CXCR6","CCR4","CCR6","RORC",
                                                    "IKZF2", "IL7R", "TIGIT", "ENTPD1", "PDCD1", "ITGAE","TNFRSF18", "TNFRSF9", "CCR8", 
                                                     "CD27",  "AHR","RORA","CXCR4", "CXCR3", "KLF2", "HLA.DRB1", "HLA.DQB1", "TGFB1" ))+RotatedAxis()
UMAPPlot(data, label=T, pt.size=1.5, split.by='disease' )

markers <- FindAllMarkers(datax, only.pos=T)
markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top20
datax <- RenameIdents(datax, 
                      "0" = "0", "1"= "1", "2"="1", "4" ="2")
cluster.averages <- AverageExpression(data, return.seurat = TRUE)

DoHeatmap(cluster.averages, features=c(
  #0
  "LAG3", "IL10", "IKZF3", "KLRB1", "CD2", "CD6", "IL6R", "TNFRSF1B", "DUSP6", "CCL20", "RGS1", "LGALS3", "DUSP16","AQP3",
  "CCR4", "CXCR6",
  #1
  "IL32", "HLA.DRB1", "HLA.DQB1", "LAYN", "CCR8", "IKZF2", "TNFRSF18", "CARD16", "HLA.DQA1", "FAM129A","TNFRSF9", "TNFRSF4", "STAT4", "DUSP4",
  "LGALS1", "IL7R", "ENTPD1", "CD27", "TIGIT", "LAIR2", "CASP1", "RORA", "CCR6", "GATA3", "CXCR3", "STAT1", 
  #2
  "LEF1", "KLF2", "TCF7", "S1PR1", "SELL", "CCR7", "CXCR4", 'KLF3', "ITGA1", "ID2", "EGR3", "IRF4", "IL21", "IFNG", 'NME1', "TCF7", "NR4A1"
  #4
  
), draw.lines = F, size=3)+ scale_fill_gradientn(colors = c("blue", "white", "red"))

DotPlot(datax[,datax$disease=='Healthy'], features=c("LAG3", "IL10", "IKZF3", "KLRB1", "CD2", "CD6", "IL6R", "TNFRSF1B", "DUSP6", "CCL20", "RGS1", "LGALS3", "DUSP16","AQP3",
                            "CCR4", "CXCR6", "PDCD1", "TIGIT", "IL7R", "TGFB1", "IFNG", "TBX21", "IL21")) + RotatedAxis()
                            #1
                            #"IL32", "HLA.DRB1", "HLA.DQB1", "LAYN", "CCR8", "IKZF2", "TNFRSF18", "CARD16", "HLA.DQA1", "FAM129A","TNFRSF9", "TNFRSF4", "STAT4", "DUSP4",
                            #"LGALS1", "IL7R", "ENTPD1", "CD27", "TIGIT", "LAIR2", "CASP1", "RORA", "CCR6", "GATA3", "CXCR3", "STAT1", 
                            #2
                            #"LEF1", "KLF2", "TCF7", "S1PR1", "SELL", "CXCR4", 'KLF3', "ITGA1", "ID2", "EGR3", "IRF4", "IL21", "IFNG", 'NME1', "TCF7", "NR4A1"
#))
data$fraction <- paste0(data@active.ident, "_", data$disease)
data$FinalIdents <- Idents(data)
Idents(data) <- data$fraction
data@active.ident <- factor(data@active.ident, 
                            levels=c("3_Healthy", "3_Non-Colitis","3_Colitis", 
                                     "1_Healthy", "1_Non-Colitis","1_Colitis", 
                                     "0_Healthy", "0_Non-Colitis","0_Colitis", 
                                     "2_Healthy", "2_Non-Colitis","2_Colitis"))

pdf("dotplot1.pdf", width=11, height=3)
DotPlot(data[,data$FinalIdents =='3'], features=c("TGFB1", "IFNG", "TBX21", "IL21", "IL17F", "CD69","TNFRSF9", "TNFRSF18","PDCD1", 
                                                                            "LAG3", "IL10", "KLRB1", "CCR6","CCR7", "CXCR6","TNFRSF1B", "IL6R", "CCR4", 
                                                                            "TIGIT", "IL7R", "ENTPD1", "CTLA4", "CD27","CARD16", "HLA.DRB1", "HLA.DQB1", 
                                                                            "SELL", "S1PR1", "KLF2", "KLF3", "ITGA4", "ITGB1"), cols=c("blue", "red"), scale.max=90) +RotatedAxis()+xlab("")+ylab("")+NoLegend()
dev.off()
pdf("dotplot2.pdf", width=11, height=3)

DotPlot(data[,data$FinalIdents =='1'], features=c("TGFB1", "IFNG", "TBX21", "IL21", "IL17F", "CD69","TNFRSF9", "TNFRSF18","PDCD1", 
                                                  "LAG3", "IL10", "KLRB1", "CCR6","CCR7", "CXCR6","TNFRSF1B", "IL6R", "CCR4", 
                                                  "TIGIT", "IL7R", "ENTPD1", "CTLA4", "CD27","CARD16", "HLA.DRB1", "HLA.DQB1", 
                                                  "SELL", "S1PR1", "KLF2", "KLF3", "ITGA4", "ITGB1"), cols=c("blue", "red"), scale.max=90) +RotatedAxis()+xlab("")+ylab("")+NoLegend()
dev.off()
pdf("dotplot3.pdf", width=11, height=3)


DotPlot(data[,data$FinalIdents =='0'], features=c("TGFB1", "IFNG", "TBX21", "IL21", "IL17F", "CD69","TNFRSF9", "TNFRSF18","PDCD1", 
                                                  "LAG3", "IL10", "KLRB1", "CCR6","CCR7", "CXCR6","TNFRSF1B", "IL6R", "CCR4", 
                                                  "TIGIT", "IL7R", "ENTPD1", "CTLA4", "CD27","CARD16", "HLA.DRB1", "HLA.DQB1", 
                                                  "SELL", "S1PR1", "KLF2", "KLF3", "ITGA4", "ITGB1"), cols=c("blue", "red"), scale.max=90) +RotatedAxis()+xlab("")+ylab("")+NoLegend()
dev.off()
pdf("dotplot4.pdf", width=11, height=3)

DotPlot(data[,data$FinalIdents =='2'], features=c("TGFB1", "IFNG", "TBX21", "IL21", "IL17F", "CD69","TNFRSF9", "TNFRSF18","PDCD1", 
                                                  "LAG3", "IL10", "KLRB1", "CCR6","CCR7", "CXCR6","TNFRSF1B", "IL6R", "CCR4", 
                                                  "TIGIT", "IL7R", "ENTPD1", "CTLA4", "CD27","CARD16", "HLA.DRB1", "HLA.DQB1", 
                                                  "SELL", "S1PR1", "KLF2", "KLF3", "ITGA4", "ITGB1"), cols=c("blue", "red"), scale.max=90) +RotatedAxis()+xlab("")+ylab("")+NoLegend()
dev.off()
pdf("alldotplot.pdf", width=11, height=6)
DotPlot(data, features=c("TGFB1", "IFNG", "TBX21", "IL21", "IL17F", "CD69","TNFRSF9", "TNFRSF18","PDCD1", 
                                                  "LAG3", "IL10", "KLRB1", "CCR6","CCR7", "CXCR6","TNFRSF1B", "IL6R", "CCR4", 
                                                  "TIGIT", "IL7R", "ENTPD1", "CTLA4", "CD27","CARD16", "HLA.DRB1", "HLA.DQB1", 
                                                  "SELL", "S1PR1", "KLF2", "KLF3", "ITGA4", "ITGB1"), cols=c("blue", "red"), scale.max=90) +RotatedAxis()+xlab("")+ylab("")
dev.off()
pdf("alldotplot2.pdf", width=11, height=6)
DotPlot(data, features=c("TGFB1", "IFNG", "TBX21", "IL21", "IL17F", "CD69","TNFRSF9", "TNFRSF18","PDCD1", 
                         "LAG3", "IL10", "KLRB1", "CCR6","CCR7", "CXCR6","TNFRSF1B", "IL6R", "CCR4", 
                         "TIGIT", "IL7R", "ENTPD1", "CTLA4", "CD27","CARD16", "HLA.DRB1", "HLA.DQB1", 
                         "SELL", "S1PR1", "KLF2", "KLF3", "ITGA4", "ITGB1"), cols=c("blue", "red"), scale.max=90) +RotatedAxis()+xlab("")+ylab("")+NoLegend()


dev.off()

dev.off()






### all T cells ######

healthy <- readRDS(list.files()[2])

FeaturePlot(healthy, "CD3E", label=T)+NoLegend()
#rough outline of Cd3E +
healthy <- subset(healthy, idents=c("10", "14", "3", "19", "1", "5", "22", "29", "8", "21", "7", "13", "12", "4","0", "9", "15", "25"))
saveRDS(healthy, "Nov04HealthyCD3.rds")

colitis <- readRDS(list.files()[1])
DefaultAssay(colitis) <- 'RNA'
FeaturePlot(colitis, "CD3E", label=T, order=T, pt.size=2)+NoLegend()
#rough outline of Cd3E +
colitis <- subset(colitis, idents=c("8", "32", "14", "10", "18", 
                                    "33", "23", "27", "1", "5", "13", "2", "3", "35", 
                                    "12", "16", "15", "6", "7", "30", "11", "0", "4", "21"))
saveRDS(colitis, "Nov04ColitisCD3.rds")

nc <- readRDS(list.files()[5])
DefaultAssay(nc) <- 'RNA'
FeaturePlot(nc, "MKI67", label=T, order=F, pt.size=2)+NoLegend()
#rough outline of Cd3E +
nc <- subset(nc, idents=c("21", "2", "1", "27", "5", "0", "4", "3", "24", "9", "19", "12", "13", "10", "16", "23"))
saveRDS(nc, "Nov04NonColitisCD3.rds")


#healthy November 7 2022
data$disease <- 'Healthy'
data <- SplitObject(data, split.by='donor')

data <- lapply(X=data, FUN = function(x){
  DefaultAssay(x) <- 'RNA'
  x<-NormalizeData(x, verbose=F)
  x <- FindVariableFeatures(x, verbose=F, nfeatures=3000)
})
features <- SelectIntegrationFeatures(object.list=data, nfeatures=3000)
data <-lapply(X=data, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})
gc();data <- FindIntegrationAnchors(object.list=data, reduction='rpca', dims=1:30);gc()
data <- IntegrateData(anchorset=data, dims=1:30, k.weight=100)


DefaultAssay(data) <- 'integrated'
data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)
data <- RunUMAP(data, dims=1:10)
#data <- RunTSNE(data, dims=1:10)
UMAPPlot(data, label=T)
data <- FindNeighbors(data, dims=1:10)
data <- FindClusters(data, resolution=0.3)
UMAPPlot(data, label=T, pt.size=1.5)
DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)
FeaturePlot(data, "TBX21", label=T)

FeaturePlot(data, c("CD4", "CD8A"), label=T)

DotPlot(data, idents=c("0", "3", "5", "6", "7", "8", "9", "10", "11", "12", "14"), 
        features=c("FOXP3", "IL17A", "CCR6", "TOX2", "CXCR5" ))
UMAPPlot(data, label=T, pt.size=1.5)

genes <- rownames(data)[grep("^MT.ND", rownames(data))]
genes <- c(genes, rownames(data)[grep("^MT.CO", rownames(data))])
genes <- c(genes, rownames(data)[grep("^MT.CYB", rownames(data))])
genes <- c(genes, rownames(data)[grep("^MT.ATP", rownames(data))])

data[['percent.mt']] <- PercentageFeatureSet(data, features=genes)

Idents(data) <- 'integrated_snn_res.0.3'
data <- RenameIdents(data, 
                     "0" = "Th1",
                     "1" = "CD8a",
                     "2" = "CD8b",
                     "3" = "Tfh",
                     "4" = "CD8c",
                     "5" = "Th17",
                     "6" = "MtHigh",
                     "7" = "Treg",
                     "8" = "Recirc1",
                     "9" = "Recirc2",
                     "10" = "MtHigh2/Bcell",
                     "11" = "NK Cells", #	NCR3
                     
                     "12" = "MtHigh3",
                     "13" = "CD8d",
                     "14" = "Contam. Monoyctes")

markers <- FindAllMarkers(subset(data, idents=c("0", "3", "5", "6", "7", "8", "9", "10", "11", "12", "14")), only.pos=T)
markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top20

data <- subset(data, idents=c("Th1", "CD8a", "CD8b", "Tfh", "CD8c", "Th17", "Treg", "Recirc1", "Recirc2", "NK Cells", "CD8d"))
##recluster
data <- SplitObject(data, split.by='donor')

data <- lapply(X=data, FUN = function(x){
  DefaultAssay(x) <- 'RNA'
  x<-NormalizeData(x, verbose=F)
  x <- FindVariableFeatures(x, verbose=F, nfeatures=3000)
})
features <- SelectIntegrationFeatures(object.list=data, nfeatures=3000)
data <-lapply(X=data, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})
gc();data <- FindIntegrationAnchors(object.list=data, reduction='rpca', dims=1:30);gc()
data <- IntegrateData(anchorset=data, dims=1:30, k.weight=100)


DefaultAssay(data) <- 'integrated'
data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)
data <- RunUMAP(data, dims=1:10)
#data <- RunTSNE(data, dims=1:10)
UMAPPlot(data, label=T, pt.size=2)
data <- FindNeighbors(data, dims=1:10)
data$lowresclust <- Idents(data)
data <- FindClusters(data, resolution=0.6)
UMAPPlot(data, label=T, pt.size=1.5)
DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)

markers <- FindAllMarkers(subset(data, idents=c("5", "6", "7", "12")), only.pos=T)
markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top20

data <- RenameIdents(data, 
                     "0" = "Th1",
                     "1" = "CD8a",
                     "2" = "CD8b",
                     "3" = "CD8c",
                     "4" = "Th17",
                     "5" = "Central Memory",
                     "6" = "Naive",
                     "7" = "Non-TRM",
                     "8" = "Treg1",
                     "9" = "Treg2",
                     "10" = "Tfh",
                     "11" = "NK cells",
                     "12" = "Cytotoxic",
                     "13" = "CD8d")
data$FinalIdents <- Idents(data)
saveRDS(data,"20221107_HealthyTcellsFullAnno.rds")



data <- CellCycleScoring(
       object = data,
       g2m.features = cc.genes$g2m.genes,
       s.features = cc.genes$s.genes
   )        
proportions(data, "ident", "Phase", "fill")





#Colitis November 7 2022
list.files()
data <- readRDS("Nov04ColitisCD3.rds")
data$disease <- 'Colitis'

genes <- rownames(data)[grep("^MT.ND", rownames(data))]
genes <- c(genes, rownames(data)[grep("^MT.CO", rownames(data))])
genes <- c(genes, rownames(data)[grep("^MT.CYB", rownames(data))])
genes <- c(genes, rownames(data)[grep("^MT.ATP", rownames(data))])
data[['percent.mt']] <- PercentageFeatureSet(data, features=genes)
data <- subset(data, subset = percent.mt<15)

Idents(data) <- 'donor'
data <- subset(data, idents=c("C1-CD3",
                              "C2-CD3",
                              "C3-CD3",
                              "C4-CD3",
                              "C5-CD3",
                              "C6-CD3",
                              "C7-CD3",
                              "C8-CD3"))


data <- SplitObject(data, split.by='donor')

data <- lapply(X=data, FUN = function(x){
  DefaultAssay(x) <- 'RNA'
  x<-NormalizeData(x, verbose=F)
  x <- FindVariableFeatures(x, verbose=F, nfeatures=3000)
})

features <- SelectIntegrationFeatures(object.list=data, nfeatures=3000)
data <-lapply(X=data, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})
gc();data <- FindIntegrationAnchors(object.list=data, reduction='rpca', dims=1:30);gc()
data <- IntegrateData(anchorset=data, dims=1:30, k.weight=100)

DefaultAssay(data) <- 'integrated'
data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)

data <- RunUMAP(data, dims=1:20)
data <- RunTSNE(data, dims=1:20)
TSNEPlot(data, label=T, pt.size=1.5)

DefaultAssay(data) <- 'integrated'
data <- FindNeighbors(data, dims=1:20)
data <- FindClusters(data, resolution=0.3)
UMAPPlot(data, label=T, pt.size=1.5)

DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)

FeaturePlot(data, c("CD4", "CD8A", "MKI67"), label=T)
DotPlot(data,  
        features=c("CD4", "CD8A", "MKI67", "TYMS" ))+
  UMAPPlot(data, label=T, pt.size=1.2)



#Idents(data) <- 'integrated_snn_res.0.X'###
data <- RenameIdents(data, 
                     "0"= "CD8a", 
                     "1" = "Th1",
                     "2"= "Th17", 
                     "3" = "Circulating",
                     "4"= "CD8b", 
                     "5" = "Treg",
                     "6"= "Tfh", 
                     "7" = "Cycling",
                     "8"= "Cycling", 
                     "9" = "NK cells", #NCR1 KLRA8
                     "10"= "CD8c", 
                     "11" = "CD8d",
                     "12"= "MitoEnriche", 
                     "13" = "ILC", #LST1+IL7R+
                     "14"= "CD8a"
                     )

markers <- FindAllMarkers(subset(data, idents=c("1","2", "3","4", "5", "12", "13", "14")), only.pos=T)
markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top20

saveRDS(data,"20221107_ColitisTcellsFullAnno.rds")


### combining Healthy+colitis #####
data <- readRDS(list.files()[6])
data2 <- readRDS(list.files()[5])
data2$disease <- 'Colitis'
Idents(data) <- 'donor'
data <- subset(data, idents=c("CT1-CD3","CT2-CD3","CT3-CD3","CT4-CD3",
                              "CT5-CD3","CT6-CD3","CT7-CD3","CT8-CD3"))
Idents(data) <- 'FinalIdents'
#just CD4
table(Idents(data))
data <- subset(data, idents=c("Th1", "Th17", "Central Memory", "Naive", "Non-TRM", "Treg1", "Treg2", "Tfh", "Cytotoxic" ))
table(Idents(data2))
data2 <- subset(data2, idents=c("Th1", "Th17", "Circulating", "Treg", "Tfh", "MitoEnriche"))
data <- merge(data, data2, add.cell.ids=c("Healthy", "Colitis")); rm(data2)

data <- SplitObject(data, split.by='donor')

data <- lapply(X=data, FUN = function(x){
  DefaultAssay(x) <- 'RNA'
  x<-NormalizeData(x, verbose=F)
  x <- FindVariableFeatures(x, verbose=F, nfeatures=2000)
})

features <- SelectIntegrationFeatures(object.list=data, nfeatures=2000)
data <-lapply(X=data, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})
gc();data <- FindIntegrationAnchors(object.list=data, reduction='rpca', dims=1:30);gc()
data <- IntegrateData(anchorset=data, dims=1:30, k.weight=100)


DefaultAssay(data) <- 'integrated'
data <- ScaleData(data, verbose=F)
data <- RunPCA(data, verbose=F)
ElbowPlot(data)

data <- RunUMAP(data, dims=1:15)
UMAPPlot(data, label=T, pt.size=1, split.by='disease')

data <- RunTSNE(data, dims=1:15)
TSNEPlot(data, label=T, pt.size=1.5, split.by='disease')

data$SepIdents <- Idents(data)

DefaultAssay(data) <- 'integrated'
data <- FindNeighbors(data, dims=1:15)
data <- FindClusters(data, resolution=0.4)
TSNEPlot(data, label=T, pt.size=1.5)+TSNEPlot(data, label=T, group.by='SepIdents')

DefaultAssay(data) <- 'RNA'
data <- NormalizeData(data)

DotPlot(data, features=c("FOXP3", "TOX2", "PDCD1", "TOX", "CXCR5", "BCL6", "TIGIT", "CTLA4"))
data <-RenameIdents(data, 
                    "0"= "Th1",
                    "1"= "Non-TRM",
                    "2"= "Th17",
                    "3"= "Central Memory",
                    "4"= "Treg1",
                    "5"= "Unhealthy Cells",
                    "6"= "Functional Th17",
                    "7"= "Treg2",
                    "8"= "CD8 T cells",
                    "9"= "Naive",
                    "10"= "TFH",
                    "11"= "Cytotoxic",
                    "12"= "NK Cells" )

UMAPPlot(data, pt.size=1.5, label=T, label.box=T)



data$clusters <- paste0(data$disease,"_", Idents(data))
data$FinalIdents <- Idents(data)

saveRDS(data, "20221109_CD4ColAndHealt.rds")





### TESTING ####
### Kegg #####
Idents(data) <- 'clusters'
m <- FindAllMarkers(data, only.pos=T, logfc.threshold = 0.15)
m <- m[m$p_val_adj<0.05,]
#install.packages('clusterProfiler')
#BiocManager::install("clusterProfiler")
library(clusterProfiler)

ENT<-bitr(unique(rownames(data@assays$RNA@data)), fromType="SYMBOL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")

m$ENTREZID<-ENT[match(m$gene,ENT$SYMBOL),"ENTREZID"]
keggEnrich<-compareCluster(ENTREZID~cluster,data=m, fun="enrichKEGG", universe=ENT$ENTREZID)
dotplot(keggEnrich)+RotatedAxis()

sa <- data.frame(keggEnrich)
sax <- data.frame(Cluster=sa$Cluster, 
                  Description = sa$Description, 
                  Ratio = sa$GeneRatio, 
                  pvalue= sa$p.adjust, 
                  ratio1 = substr(sa$GeneRatio, 1,unlist(gregexpr('/', sa$GeneRatio))-1),
                  ratio2 =substr(sa$GeneRatio, unlist(gregexpr('/', sa$GeneRatio))+1, nchar(sa$GeneRatio))
)
sax$actualratio <- as.numeric(sax$ratio1)/as.numeric(sax$ratio2)


keggBar <-function(grep) {
  ggplot(sax[grep(grep, sax$Description),], aes(Cluster, actualratio*100))+
    geom_col(aes(fill=pvalue))+
    scale_fill_gradient2(low = "red",mid = 'yellow', high='light yellow' ,midpoint = mean(sax$pvalue[grep(grep, sax$Description)]))+
    theme_classic()+coord_flip()+ylab("%")+xlab("")+facet_wrap(~Description)
};
keggBar("Inflam")
dotplot(keggEnrich, showCategory=1)+RotatedAxis()
unique(sax$Description)


data <- CellCycleScoring(
  object = data,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)        
proportions(data, "ident", "Phase", "fill")

DiffPlots <- function(data, ident, nlabel=50) {
  cells <- subset(data, idents = ident)
  Idents(cells) <- "disease"
  markstomap <- FindAllMarkers(cells, only.pos=T)
  markstomap %>%
    group_by(cluster) %>%
    top_n(n =nlabel, wt = avg_log2FC) -> top50
  cells <- as.data.frame(log1p(AverageExpression(cells, verbose = FALSE)$RNA))
  cells$gene <- rownames(cells)
  cells$color <- NA
   for(i in 1:nrow(top50)){
    for(a in 1:nrow(cells)){
      if(top50$gene[i] == cells$gene[a] & cells$Healthy[a]>cells$Colitis[a]) {
        cells$color[a] <- 'Healthy'
      } else if(top50$gene[i] == cells$gene[a] & cells$Healthy[a]<cells$Colitis[a]) {
        cells$color[a] <- 'Colitis'
      }
    }
   }
  
  for (i in 1:nrow(cells)){
    if(cells$gene[i] %in% top50$gene) {
    } else {
      cells$color[i] <- 'NA'
    }
  }
  p1 <- ggplot(cells, aes(Healthy, Colitis, color=color )) + geom_point(alpha=0.5) + ggtitle(ident)
  p1 <- LabelPoints(plot = p1, points = top50$gene, repel = TRUE)
  p1+theme_classic() + NoLegend()
}

DiffPlots(data, names(table(Idents(data))), nlabel=50)

### Kegg #####

m <- FindAllMarkers(data, only.pos=T, logfc.threshold = 0.15)
m <- m[m$p_val_adj<0.05,]
#install.packages('clusterProfiler')
#BiocManager::install("clusterProfiler")
#library(clusterProfiler)

ENT<-bitr(unique(rownames(data@assays$RNA@data)), fromType="SYMBOL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")

m$ENTREZID<-ENT[match(m$gene,ENT$SYMBOL),"ENTREZID"]
keggEnrich<-compareCluster(ENTREZID~cluster,data=m, fun="enrichKEGG", universe=ENT$ENTREZID)
dotplot(keggEnrich, )
sa <- data.frame(keggEnrich)
sax <- data.frame(Cluster=sa$Cluster, 
                  Description = sa$Description, 
                  Ratio = sa$GeneRatio, 
                  pvalue= sa$p.adjust, 
                  ratio1 = substr(sa$GeneRatio, 1,unlist(gregexpr('/', sa$GeneRatio))-1),
                  ratio2 =substr(sa$GeneRatio, unlist(gregexpr('/', sa$GeneRatio))+1, nchar(sa$GeneRatio))
)
sax$actualratio <- as.numeric(sax$ratio1)/as.numeric(sax$ratio2)


keggBar <-function(grep) {
  ggplot(sax[grep(grep, sax$Description),], aes(Cluster, actualratio*100))+
    geom_col(aes(fill=pvalue))+
    scale_fill_gradient2(low = "red",mid = 'yellow', high='light yellow' ,midpoint = mean(sax$pvalue[grep(grep, sax$Description)]))+
    theme_classic()+coord_flip()+ylab("%")+xlab("")+facet_wrap(~Description)
};keggBar("IL-17")
dotplot(keggEnrich, showCategory=1)+RotatedAxis()
unique(sax$Description)




