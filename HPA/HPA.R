# ## set up 
# prot <-read.delim("/Users/thomasoneil/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/zOrganisedExperiments/normal_tissue.tsv",row.names = NULL)
# prot$Exp <- ifelse(prot$Level == "High", 3, ifelse(prot$Level=="Medium", 2, ifelse(prot$Level == "Low", 1, ifelse(prot$Level == 'Not detected', 0, 4))))
# prot = prot[prot$Exp<4,]
# 
# prot2 <- prot[prot$Tissue %in% c("cervix", "colon","duodenum", 'endometrium',     'endometrium 1',     'endometrium 2', 'lymph node', 'skin' ,          'skin 1'      ,      'skin 2', "rectum",
#                                  'small intestine'   ,  'smooth muscle' , 'vagina'),]
# 
# 
# 
# prot2$Tissue[prot2$Tissue=="skin 1" | prot2$Tissue=="skin 2"] = "skin"
# prot2$Tissue[prot2$Tissue=="endometrium 1" | prot2$Tissue=="endometrium 2"] = "endometrium"
# 
# prot2$Tissue <- factor(prot2$Tissue, levels=c("colon", "rectum", "small intestine", "duodenum", 'vagina', 'cervix', "endometrium", "skin", "lymph node", "smooth muscle"))
# 
# prot3 <- prot2[prot2$Cell.type %in% c("Langerhans","langerhans cells", "macrophages","lymphocytes","mucosal lymphoid cells",
# 
#                                       "keratinocytes" ,"epidermal cells","cells in basal layer","cells in granular layer","squamous epithelial cells","ciliated epithelial cells","nonciliated glandular epithelial cells","nonciliated luminal epithelial cells","smooth muscle cells", "fibroblasts","stromal fibroblasts", "fibrohistiocytic cells","endothelial cells","vascular mural cells",
#                                       "enterocytes", "enterocytes - Microvilli" ,"goblet cells","hair follicles", "glands of Brunner", "sebaceous glands","sebaceous cells", "secretory cells","sweat ducts",
#                                       "melanocytes","glandular cells", "peripheral nerve/ganglion","germinal center cells","non-germinal center cells","extracellular matrix"),]
# prot3$Cell.type <-factor(prot3$Cell.type, levels=c("Langerhans","langerhans cells", "macrophages","lymphocytes","mucosal lymphoid cells",
# 
#                                                    "keratinocytes" ,"epidermal cells","cells in basal layer","cells in granular layer","squamous epithelial cells","ciliated epithelial cells","nonciliated glandular epithelial cells","nonciliated luminal epithelial cells","smooth muscle cells", "fibroblasts","stromal fibroblasts", "fibrohistiocytic cells","endothelial cells","vascular mural cells",
#                                                    "enterocytes", "enterocytes - Microvilli" ,"goblet cells","hair follicles", "glands of Brunner", "sebaceous glands","sebaceous cells", "secretory cells","sweat ducts",
#                                                    "melanocytes","glandular cells", "peripheral nerve/ganglion","germinal center cells","non-germinal center cells","extracellular matrix"))
# 
# rna  <-  read.delim("/Users/thomasoneil/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/zOrganisedExperiments/rna_single_cell_type.tsv",row.names = NULL)
# 
# useful <- rna[rna$Cell.type %in% c("B-cells", "Basal keratinocytes", "Basal squamous epithelial cells", "dendritic cells", "Endothelial cells", "granulocytes", "Fibroblasts",
#                                    "Intestinal goblet cells", "Macrophages", "Lymphatic endothelial cells", "Langerhans cells", "Melanocytes", "monocytes", "NK-cells", "Paneth cells",
#                                    "Plasma cells", "Schwann cells", "Smooth muscle cells","Squamous epithelial cells", 'Suprabasal keratinocytes', "T-cells"),]
# 
# useful$Cell.type <- factor(useful$Cell.type, levels = c("Langerhans cells","dendritic cells" ,"Macrophages"  ,"monocytes","T-cells","NK-cells","B-cells","Plasma cells","granulocytes",
#                                                         "Basal keratinocytes","Suprabasal keratinocytes","Basal squamous epithelial cells","Squamous epithelial cells","Smooth muscle cells" ,"Fibroblasts","Endothelial cells", "Lymphatic endothelial cells",
#                                                         "Melanocytes", "Intestinal goblet cells", "Paneth cells","Schwann cells" ))
# 
# rna1  <-  read.delim("/Users/thomasoneil/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/zOrganisedExperiments/rna_single_cell_type_tissue.tsv",row.names = NULL)
# 
# gut <- rna1[rna1$Tissue %in% c("colon", "small intestine", "skin", "rectum", "lymph node"),]
# gut$Cell.type <- factor(gut$Cell.type, levels = c("langerhans cells", "macrophages", 't-cells',"b-cells",  "granulocytes", "mixed immune cells",
#                                                   'basal keratinocytes', "suprabasal keratinocytes","smooth muscle cells",  "fibroblasts", 'endothelial cells','distal enterocytes', 'enteroendocrine cells',
#                                                   "intestinal goblet cells","melanocytes", "paneth cells","proximal enterocytes", "mixed cell types","undifferentiated cells"))
# 
# write.csv(prot3, '/Users/thomasoneil/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/Scripts/ProteinAtlas/prot.tsv')
# write.csv(useful, '/Users/thomasoneil/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/Scripts/ProteinAtlas/celltype.csv')
# write.csv(gut, '/Users/thomasoneil/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/Scripts/ProteinAtlas/celltypetissue.tsv')
# HPA <- list(Protein = prot3, 
#             CellType=useful, 
#             CellTypeTissue=gut)

library(ggplot2);library(Seurat);library(cowplot)
HPA <- readRDS(url("https://raw.githubusercontent.com/tomoneil58/LabCode/main/HPA/HPA.rds"))
hpa <- function(gene="CD14", log=F) {
  if(log==T) {
    p2 <- ggplot(HPA[[1]][HPA[[1]]$Gene.name==gene,], aes(Cell.type, log(Exp), fill=Tissue))+geom_col()+theme_classic()+RotatedAxis()+NoLegend()+ylab("Protein Exp")
    p1 <- ggplot(HPA[[1]][HPA[[1]]$Gene.name==gene,], aes(Tissue, log(Exp), fill=Tissue))+geom_col()+theme_classic()+RotatedAxis()+NoLegend()+geom_vline(xintercept = c(4.5,7.5))+ggtitle(gene)+ylab("Protein Exp")
    p3 <-ggplot(HPA[[2]][HPA[[2]]$Gene.name==gene,], aes(Cell.type, log(nTPM+1), fill=Cell.type))+geom_col()+theme_classic()+RotatedAxis()+NoLegend()+geom_vline(xintercept = c(4.5,6.5,8.5,13.5,15.5,17.5))
    p4 <-ggplot(HPA[[3]][HPA[[3]]$Gene.name==gene,], aes(Cell.type, log(nTPM+1), fill=Tissue))+geom_col()+theme_classic()+RotatedAxis()+theme(legend.position = 'top', legend.direction= "horizontal")+geom_vline(xintercept = c(2.5,4.5,6.5,8.5,10.5))  
  }
  else{
    p2 <- ggplot(HPA[[1]][HPA[[1]]$Gene.name==gene,], aes(Cell.type, Exp, fill=Tissue))+geom_col()+theme_classic()+RotatedAxis()+NoLegend()+ylab("Protein Exp")
    p1 <- ggplot(HPA[[1]][HPA[[1]]$Gene.name==gene,], aes(Tissue, Exp, fill=Tissue))+geom_col()+theme_classic()+RotatedAxis()+NoLegend()+geom_vline(xintercept = c(4.5,7.5))+ggtitle(gene)+ylab("Protein Exp")
    p3 <-ggplot(HPA[[2]][HPA[[2]]$Gene.name==gene,], aes(Cell.type, nTPM, fill=Cell.type))+geom_col()+theme_classic()+RotatedAxis()+NoLegend()+geom_vline(xintercept = c(4.5,6.5,8.5,13.5,15.5,17.5))
    p4 <-ggplot(HPA[[3]][HPA[[3]]$Gene.name==gene,], aes(Cell.type, nTPM, fill=Tissue))+geom_col()+theme_classic()+RotatedAxis()+theme(legend.position = 'top', legend.direction= "horizontal")+geom_vline(xintercept = c(2.5,4.5,6.5,8.5,10.5))
  }
  print(plot_grid(p1,p2,p3,p4))
}
