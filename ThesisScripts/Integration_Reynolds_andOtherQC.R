source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

library(Seurat)
library(ggplot2)
library(dplyr)
set.seed(34325)


# cell type labels
monocytes <- c("DC1", "DC2", "Inf_mono","Macro_1", "Macro_2","LC_1", "LC_2", "LC_3", "LC_4", "MigDC", "Mono", 'moDC_1', 'moDC_2', 'moDC_3')
tcells <- c("ILC1_3", "ILC1_NK","ILC2", "NK", "Tc", "Th", "Treg")
bcells <- c("Plasma")
nonimmune <-c("Differentiated_KC", "Differentiated_KC*", "F1", "F2", "F3", "LE1", "LE2", "Mast_cell", "Melanocyte", "Pericyte_1_non_inflamm", "Pericyte_2_inflamm", "Proliferating_KC", "Schwann1", "Schwann2", "Undifferentiated_KC*", "VE1", "VE2", "VE3")

### FIRST UMAP OF WHOLE TISSUE ####

tissue <- readRDS("20210903_inflamed_unprocessed.rds")

tissue[['percent.mt']] <- PercentageFeatureSet(tissue, pattern="^MT-")
Idents(tissue)<- 'Site'
# Get just lesion
tissue <- subset(tissue, idents = c('lesion')) 

Idents(tissue) <- 'final_clustering'
# Next line is to remove nan cells
tissue <- subset(tissue, idents = c(tcells, bcells, monocytes, nonimmune)) #remove nan

# Split data into donors (in a list)
tissue <- SplitObject(tissue, split.by='donor_id')

# run reciprocal PCA - we can discuss the different methods later! 
# for each donor in the list, normalize the data and identify the variable features
tissue <- lapply(X=tissue, FUN = function(x){
  x<-NormalizeData(x, verbose=F)
  x<-FindVariableFeatures(x, verbose=F)
})
features <- SelectIntegrationFeatures(object.list=tissue)

# Scale and run PCA on each donor using the features identified above
tissue <-lapply(X=tissue, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})

tissue <- FindIntegrationAnchors(object.list=tissue, reduction='rpca', dims=1:50)

tissue <- IntegrateData(anchorset=tissue, dims=1:50)

tissue <- ScaleData(tissue, verbose=F)
tissue <- RunPCA(tissue, verbose=F)
tissue <- RunUMAP(tissue, dims=1:50)
tissue <- RunTSNE(tissue, dims=1:50)
tissue <- FindNeighbors(tissue, dims=1:50)
tissue <- FindClusters()

pdf("Plot1.pdf", 
    width=, #enter a width and height that you want for the pdf OR ignore this and use Export
    height= ) 

UMAPPlot(tissue, label=T, group.by='final_clustering') # For reference 
UMAPPlot(tissue, label=F, group.by='final_clustering')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='final_clustering') + NoAxes() # can add your own axes in powerpoint
UMAPPlot(tissue, label=F, group.by='final_clustering') + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
TSNEPlot(tissue, label=T, group.by='final_clustering')
TSNEPlot(tissue, label=F, group.by='final_clustering')
TSNEPlot(tissue, label=F, group.by='final_clustering') + NoAxes()
TSNEPlot(tissue, label=F, group.by='final_clustering') + NoAxes()+NoLegend()

dev.off() #Closes off the pdf - ignore this if you dont use the pdf() function

### Plot 2: Clustering #### ONLY USING LESION (NO HEALTHY) BECAUSE R KEEPS CRASHING

Idents(tissue) <- 'final_clustering'
reference <- subset(tissue, idents=monocytes)

# for merged, we can restart the process with no integration methods

tissue <- reference
DefaultAssay(tissue) <- 'RNA'

tissue<- NormalizeData(tissue)
tissue <- FindVariableFeatures(tissue)
tissue <- ScaleData(tissue, features=VariableFeatures(tissue))
tissue <- RunPCA(tissue, dims=1:10)
tissue <- RunUMAP(tissue, dims=1:10)
tissue <- RunTSNE(tissue, dims=1:15)
UMAPPlot(tissue, label=F, pt.size=.8, group.by='donor_id') + NoAxes()+NoLegend()+ggtitle("")
UMAPPlot(tissue, label=T, label.box=T, pt.size=.8,group.by="final_clustering") + NoAxes()+NoLegend()+ggtitle("")

tissue <- FindNeighbors(tissue, dims=1:15)
tissue <- FindClusters(tissue, resolution = 0.3)
UMAPPlot(tissue, label=F) + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
proportions(tissue, 'final_clustering', 'donor_id', position='fill')

pdf("Plot2_Merged.pdf", 
    width=, #enter a width and height that you want for the pdf OR ignore this and use Export
    height= ) 
 
UMAPPlot(tissue, label=F, group.by='donor_id')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='donor_id') + NoAxes() # can add your own axes in powerpoint
UMAPPlot(tissue, label=F, group.by='donor_id') + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
TSNEPlot(tissue, label=F, group.by='donor_id')
TSNEPlot(tissue, label=F, group.by='donor_id') + NoAxes()
TSNEPlot(tissue, label=F, group.by='donor_id') + NoAxes()+NoLegend()

dev.off() #Closes off the pdf - ignore this if you dont use the pdf() function


## BEER ##
# Needs a bit more R work - need to get it down to raw matricies

list <- SplitObject(reference, split.by='donor_id')

P1 <- list[[5]]@assays$RNA@counts
P2 <- list[[6]]@assays$RNA@counts
P3 <- list[[7]]@assays$RNA@counts
E1 <- list[[4]]@assays$RNA@counts
E2 <- list[[1]]@assays$RNA@counts
E3 <- list[[2]]@assays$RNA@counts
E4 <- list[[3]]@assays$RNA@counts


BATCH=c(rep('P1',ncol(P1)),
        rep('P2',ncol(P2)),
        rep('P3',ncol(P3)),
        rep('E1',ncol(E1)),
        rep('E2',ncol(E2)),
        rep('E3',ncol(E3)),
        rep('E4',ncol(E4))
)

Ss =.simple_combine(P1, P2, FILL=TRUE)$combine; rm(P1); rm(P2)
Ss =.simple_combine(Ss, P3, FILL=TRUE)$combine; rm(P3);
Ss =.simple_combine(Ss, E1, FILL=TRUE)$combine; rm(E1);
Ss =.simple_combine(Ss, E2, FILL=TRUE)$combine; rm(E2);
Ss =.simple_combine(Ss, E3, FILL=TRUE)$combine; rm(E3);
DATA=.simple_combine(Ss, E4, FILL=TRUE)$combine; rm(E4)

rm(list); rm(Ss)

beer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE )

beerseurat <- beer$seurat
beerseurat <- RunUMAP(object = beerseurat, reduction='pca',dims = beer$select, check_duplicates=FALSE)
beerseurat <- RunTSNE(object = beerseurat, reduction='pca',dims = beer$select, check_duplicates=FALSE)
beerseurat <- AddMetaData(beerseurat, 
                          metadata = reference@meta.data) 
UMAPPlot(beerseurat, pt.size=.8, label=F, group.by='donor_id') + NoAxes()+NoLegend() +ggtitle("")#gets rid of legend - all these help create the exact 
UMAPPlot(beerseurat, pt.size=.8, label=T, label.box=T, group.by='final_clustering') + NoAxes()+NoLegend() +ggtitle("")#gets rid of legend - all these help create the exact 

beerseurat <- FindNeighbors(beerseurat, dims=1:15)
beerseurat <- FindClusters(beerseurat, resolution = 0.3)
UMAPPlot(beerseurat, label=F) + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
proportions(beerseurat, 'ident', 'donor_id', position='fill')


pdf("Plot2_BEER.pdf", 
    width=, #enter a width and height that you want for the pdf OR ignore this and use Export
    height= ) 

UMAPPlot(tissue, label=F, group.by='donor_id')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='donor_id') + NoAxes() # can add your own axes in powerpoint
UMAPPlot(tissue, label=F, group.by='donor_id') + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
TSNEPlot(tissue, label=F, group.by='donor_id')
TSNEPlot(tissue, label=F, group.by='donor_id') + NoAxes()
TSNEPlot(tissue, label=F, group.by='donor_id') + NoAxes()+NoLegend()

dev.off() #Closes off the pdf - ignore this if you dont use the pdf() function

rm(beerseurat);rm(DATA);rm(beer);rm(BATCH)

### rPCA ###

tissue <- SplitObject(reference, split.by='donor_id'); rm(reference)
tissue <- lapply(X=tissue, FUN = function(x){
  x<-NormalizeData(x, verbose=F)
  x<-FindVariableFeatures(x, verbose=F)
})
features <- SelectIntegrationFeatures(object.list=tissue)

tissue <-lapply(X=tissue, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})

tissue <- FindIntegrationAnchors(object.list=tissue, reduction='rpca', dims=1:30)
tissue <- IntegrateData(anchorset=tissue, dims=1:30)

tissue <- ScaleData(tissue, verbose=F)
tissue <- RunPCA(tissue, verbose=F)
tissue <- RunUMAP(tissue, dims=1:30)
tissue <- RunTSNE(tissue, dims=1:30)

UMAPPlot(tissue, pt.size=.8, label=F, group.by='donor_id') + NoAxes()+NoLegend() +ggtitle("")#gets rid of legend - all these help create the exact 
UMAPPlot(tissue, pt.size=.8, label=T, label.box=T, group.by='final_clustering') + NoAxes()+NoLegend() +ggtitle("")#gets rid of legend - all these help create the exact 

pdf("Plot2_rPCA.pdf", 
    width=, #enter a width and height that you want for the pdf OR ignore this and use Export
    height= ) 
UMAPPlot(tissue, label=F, group.by='donor_id')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='donor_id') + NoAxes() # can add your own axes in powerpoint
UMAPPlot(tissue, label=F, group.by='donor_id') + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
TSNEPlot(tissue, label=F, group.by='donor_id')
TSNEPlot(tissue, label=F, group.by='donor_id') + NoAxes()
TSNEPlot(tissue, label=F, group.by='donor_id') + NoAxes()+NoLegend()

dev.off() #Closes off the pdf - ignore this if you dont use the pdf() function


## CCA 

DefaultAssay(tissue) <- 'RNA'
tissue <- SplitObject(tissue, split.by='donor_id'); 
tissue <- lapply(X=tissue, FUN = function(x){
  x<-NormalizeData(x, verbose=F)
  x<-FindVariableFeatures(x, verbose=F)
})
features <- SelectIntegrationFeatures(object.list=tissue)

tissue <-lapply(X=tissue, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})

tissue <-  FindIntegrationAnchors(object.list = tissue, anchor.features = features)
tissue <- IntegrateData(anchorset=tissue, dims=1:30)

tissue <- ScaleData(tissue, verbose=F)
tissue <- RunPCA(tissue, verbose=F)
tissue <- RunUMAP(tissue, dims=1:30)
tissue <- RunTSNE(tissue, dims=1:30)

UMAPPlot(tissue, pt.size=.8, label=F, group.by='donor_id') + NoAxes()+NoLegend() +ggtitle("")#gets rid of legend - all these help create the exact 
UMAPPlot(tissue, pt.size=.8, label=T, label.box=T, group.by='final_clustering') + NoAxes()+NoLegend() +ggtitle("")#gets rid of legend - all these help create the exact 


### QC: removing mitochondrial ####
tissue[['percent.mt']] <- PercentageFeatureSet(tissue, pattern='^MT-')

tissue$mito <- ifelse(tissue$percent.mt <5, "<5", ">5")

pdf("Plot3_QC.pdf", 
    width=, #enter a width and height that you want for the pdf OR ignore this and use Export
    height= ) 
VlnPlot(tissue, 'percent.mt', group.by='donor_id')

UMAPPlot(tissue, label=T, group.by='final_clustering')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='mito')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='mito') + NoAxes() # can add your own axes in powerpoint
UMAPPlot(tissue, label=F, group.by='mito') + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
TSNEPlot(tissue, label=F, group.by='mito')
TSNEPlot(tissue, label=F, group.by='mito') + NoAxes()
TSNEPlot(tissue, label=F, group.by='mito') + NoAxes()+NoLegend()

UMAPPlot(tissue, label=F, group.by='final_clustering')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='final_clustering') + NoAxes() # can add your own axes in powerpoint
UMAPPlot(tissue, label=F, group.by='final_clustering') + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
TSNEPlot(tissue, label=F, group.by='final_clustering')
TSNEPlot(tissue, label=F, group.by='final_clustering') + NoAxes()
TSNEPlot(tissue, label=F, group.by='final_clustering') + NoAxes()+NoLegend()


dev.off() #Closes off the pdf - ignore this if you dont use the pdf() function

tis2 <- subset(tissue, subset = percent.mt <5)

tissue <- SplitObject(tis2, split.by='donor_id')
tissue <- lapply(X=tissue, FUN = function(x){
  x<-NormalizeData(x, verbose=F)
  x<-FindVariableFeatures(x, verbose=F)
})
features <- SelectIntegrationFeatures(object.list=tissue)

tissue <-lapply(X=tissue, FUN=function(x){
  x<-ScaleData(x, features=features, verbose=T)
  x<-RunPCA(x, features=features, verbose=F)
})

tissue <- FindIntegrationAnchors(object.list=tissue, reduction='rpca', dims=1:50)
tissue <- IntegrateData(anchorset=tissue, dims=1:50)

tissue <- ScaleData(tissue, verbose=F)
tissue <- RunPCA(tissue, verbose=F)
tissue <- RunUMAP(tissue, dims=1:50)
tissue <- RunTSNE(tissue, dims=1:50)

pdf("Plot3_QC2.pdf", 
    width=, #enter a width and height that you want for the pdf OR ignore this and use Export
    height= ) 
VlnPlot(tissue, 'percent.mt', group.by='donor_id')

UMAPPlot(tissue, label=T, group.by='final_clustering')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='final_clustering')# Better to draw your own shapes around the clusters 
UMAPPlot(tissue, label=F, group.by='final_clustering') + NoAxes() # can add your own axes in powerpoint
UMAPPlot(tissue, label=F, group.by='final_clustering') + NoAxes()+NoLegend() #gets rid of legend - all these help create the exact 
TSNEPlot(tissue, label=F, group.by='final_clustering')
TSNEPlot(tissue, label=F, group.by='final_clustering') + NoAxes()
TSNEPlot(tissue, label=F, group.by='final_clustering') + NoAxes()+NoLegend()

dev.off() #Closes off the pdf - ignore this if you dont use the pdf() function



#### functions #####


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








