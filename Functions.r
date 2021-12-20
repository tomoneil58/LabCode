## get CSV
# specific data object, genes and initals
# e.g. getCSV(cpm, c("CD207", "CD209"), "TO")
# e.g. getCSV(cpm, top100, "TO") 
# Will output 



## cpmBox
# takes a gene and heading
# requires a meta-file called 
cpmBox <- function(gene, heading = "Bulk RNA Expression (cpm)") {
  temp <- meta
  temp$gene <- cpm[gene,]
  ggplot(temp, 
  aes(x=group, #change this
      y=gene, #change this
      colour=tcell) #change this
   ) + ggtitle(heading)+geom_boxplot()+geom_point()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()+facet_wrap(~tissue) #change this
}

## get CSV
# specific data object, genes and initals
# e.g. getCSV(cpm, c("CD207", "CD209"), "TO")
# e.g. getCSV(cpm, top100, "TO") 
# Will output cpmBox boxplots as a pdf and label the csv with the date and your initals

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

DiffGenes <- function(data, subset=c(), split.by=c(), n.genes=10) {
  temp <- subset(data#
                 , idents=c(subset))#
  Idents(temp) <- split.by
  
  avg.cells <- as.data.frame(log1p(AverageExpression(temp, verbose = FALSE)$RNA))
  avg.cells$gene <- rownames(avg.cells)
  
  marks <- FindAllMarkers(temp, assay = 'RNA', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)
  topm <- marks %>% group_by(cluster) %>% top_n(n = n.genes, wt = avg_log2FC)
  
  nam <- colnames(avg.cells)
  colnames(avg.cells) <- c("ident.1", "ident.2")
  
  p1 <- ggplot(avg.cells, aes(ident.1, ident.2)) + geom_point() + ggtitle(paste("Differential expressions of", subset))+xlab(nam[1])+ylab(nam[2])
  LabelPoints(plot = p1, points = topm$gene, repel = TRUE)
}
