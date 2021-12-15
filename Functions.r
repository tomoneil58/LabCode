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
