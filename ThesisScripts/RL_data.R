#Applying to the skin data

#healthy only

skin #<- #tcells
skindc #<- dcs

skin$idents.forRL <- Idents(skin)
skindc$idents.forRL <- Idents(skindc)

Idents(skindc) <- 'Status'
skin$celltype <- "Tcells"
skindc$celltype <- 'DCs'
skindc <- subset(skindc, idents=c("Healthy"))

Idents(skin) <- 'idents.forRL'
Idents(skindc) <- 'idents.forRL'

skins <- merge(skin, skindc)
skins <- SplitObject(skins, split.by='donor_id')

skins <- lapply(X = skins, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})

features <- SelectIntegrationFeatures(object.list = skins, nfeatures = 4000)

skins <- lapply(X = skins, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

skins <- FindIntegrationAnchors(object.list = skins, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
skins <- IntegrateData(anchorset = skins)


skins <- ScaleData(skins)
skins <- RunPCA(skins)
ElbowPlot(skins)
skins <- RunUMAP(skins, dims=1:30)

UMAPPlot(skins, group.by='donor_id')

Idents(skins) <- 'idents.forRL'
skins <- RenameIdents(skins,
                      "cDC2/MDDC1" = 'cDC2/MDDC', 
                      "cDC2/MDDC2" = 'cDC2/MDDC',
                      "cDC2/MDDC3" = 'cDC2/MDDC',
                      "cDC2/MDDC4" = 'cDC2/MDDC',
                      "Cont.T cell" = 'Junk',
                      "Inf.cDC2-1" = 'Inf.cDC2',
                      "Inf.cDC2-2" = 'Inf.cDC2',
                      "Inf.cDC2-3" = 'Inf.cDC2',
                      "Inf.cDC2-4" = 'Inf.cDC2',
                      "MDM1" = 'MDM',
                      "MDM2" = 'MDM',
                      "Mig.cDC2/MDDC1" = "Mig.cDC2/MDDC",
                      "Mig.cDC2/MDDC2" = "Mig.cDC2/MDDC",
                      "Mig.cDC2/MDDC3" = "Mig.cDC2/MDDC",
                      "Mig.cDC2/MDDC4" = "Mig.cDC2/MDDC"
                      
)
skins$summarisedDCs <- Idents(skins)
skins <- RenameIdents(skins, 
                      "1" = "T cells", 
                      "2" = "T cells", 
                      "3" = "T cells", 
                      "4" = "T cells", 
                      "5" = "T cells", 
                      "6" = "T cells", 
                      "7" = "T cells", 
                      "8" = "T cells", 
                      "Treg1" = "Treg", 
                      "Treg2" = "Treg"
                      
)
skins$summarisedDCsandT <- Idents(skins)
DefaultAssay(skins) <- 'RNA'
skins <- NormalizeData(skins)
saveRDS(skins, "Skin/DCandT.rds")
skins <- subset(skins, idents=idents[-c(8,20)]) 
names(table(Idents(skins)))

rl <- read.delim('Homo Sapiens.txt')
rl$unique <- paste0(rl$Ligand_id, rl$Receptor_id)
rl<-rl %>% distinct(unique, .keep_all = TRUE)

idents<-names(table(Idents(skins)))
rl_table <- data.frame(matrix(NA, nrow=length(idents), ncol = length(idents)))
colnames(rl_table) <- idents
rownames(rl_table) <- idents
rl_table <- rl_table %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

colnames(rl_table) <- c("Ligand","Receptor","Interactions")

count=0
df<-data.frame(Index='empty', Ligand='empty', Receptor='empty', Ligand.Clust = 'empty', Receptor.Clust = 'empty', 
               r.avgLog2FC = 'empty', r.receptor.pct = 'empty', r.ligand.pct='empty', r.pvaladj = 'empty',
               l.avgLog2FC = 'empty', l.receptor.pct = 'empty', l.ligand.pct='empty', l.pvaladj = 'empty', Type = 'empty')

for (n in 1:nrow(rl_table)){
  ident1 <- rl_table$Ligand[n]
  ident2 <- rl_table$Receptor[n]
  count=0
  print(paste0(ident1," & ", ident2, " - ", n,"/",nrow(rl_table)))
  if(ident1 != ident2) {
    #find differential markers + sort for significant ones
    x<-FindMarkers(skins, ident.1=ident1, ident.2=ident2, logfc.threshold = 0.4, verbose = T)
    x<-x[x$p_val_adj<=0.05,]
    ident1_gene <- x[x$avg_log2FC>0,]
    ident2_gene <- x[x$avg_log2FC<0,]
    ident1_gene$gene <- rownames(ident1_gene)
    ident2_gene$gene <- rownames(ident2_gene)
    
    for(i in 1:nrow(rl)){
      if(rl$Receptor_id[i] %in% rownames(ident1_gene)){
        if(rl$Ligand_id[i] %in% rownames(ident2_gene)){
          count=count+1
          df<-rbind(df, data.frame(Index=i, Ligand=rl$Ligand_id[i], Receptor=rl$Receptor_id[i], Ligand.Clust = ident1, Receptor.Clust = ident2, 
                                   l.avgLog2FC = as.numeric(ident1_gene$avg_log2FC[ident1_gene$gene==rl$Ligand_id[i]]), 
                                   l.receptor.pct = as.numeric(ident1_gene$pct.1[ident1_gene$gene==rl$Ligand_id[i]]), 
                                   l.ligand.pct = as.numeric(ident1_gene$pct.2[ident1_gene$gene==rl$Ligand_id[i]]),
                                   l.pvaladj = as.numeric(ident1_gene$p_val_adj[ident1_gene$gene==rl$Ligand_id[i]]),
                                   r.avgLog2FC = as.numeric(ident2_gene$avg_log2FC[ident2_gene$gene==rl$Receptor_id[i]]), 
                                   r.receptor.pct = as.numeric(ident2_gene$pct.1[ident2_gene$gene==rl$Receptor_id[i]]), 
                                   r.ligand.pct=as.numeric(ident2_gene$pct.2[ident2_gene$gene==rl$Receptor_id[i]]), 
                                   r.pvaladj = as.numeric(ident2_gene$p_val_adj[ident2_gene$gene==rl$Receptor_id[i]]), type = rl$Type[i]))
        }
      }
    }
  }
  rl_table$`Interactions`[n] <- count
}

ggplot(rl_table, aes(x = Receptor, y = Ligand, fill = Interactions)) +
  geom_tile()+ scale_fill_gradientn(colors = c("dark blue", "yellow"), limits = c(0,27)) +RotatedAxis()

# separate by adhesion 
cell.ad <- df[df$Type=="Cell adhesion",]
rl_table2 <- rl_table
rl_table2$Interactions <- NA

for (n in 1:nrow(rl_table2)) {
  rl_table2$Interactions[n] <- nrow(cell.ad[cell.ad$Ligand.Clust==rl_table2$Ligand[n] & cell.ad$Receptor.Clust==rl_table2$Receptor[n],])
}

ggplot(rl_table2, aes(x = Receptor, y = Ligand, fill = Interactions)) +
  geom_tile()+ scale_fill_gradientn(colors = c("dark blue", "yellow"), limits = c(0,15)) +RotatedAxis()


# separate by adhesion 
secrete <- df[df$Type=="Secreted protein to receptor interaction" | df$Type=="Cytokine-cytokine receptor interaction" ,]
rl_table3 <- rl_table
rl_table3$Interactions <- NA

for (n in 1:nrow(rl_table3)) {
  rl_table3$Interactions[n] <- nrow(secrete[secrete$Ligand.Clust==rl_table3$Ligand[n] & secrete$Receptor.Clust==rl_table3$Receptor[n],])
}

ggplot(rl_table3, aes(x = Receptor, y = Ligand, fill = Interactions)) +
  geom_tile()+ scale_fill_gradientn(colors = c("dark blue", "yellow"), limits = c(0,15)) +RotatedAxis()

