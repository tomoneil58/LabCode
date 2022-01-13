#Jacinta's code


setwd("~/Desktop/Jacinta/")
oldList<- readRDS("ListForHarmanBulk.rds")
newList <- readRDS("ListofObject_BulkNewDCData.rds")

oldmeta <- oldList[[6]]
newmeta <- newList[[6]]

oldcpm <- oldList[[4]]
newcpm <- newList[[4]]

oldlcpm <- oldList[[5]]
newlcpm <- newList[[5]]

# getting just dermal DCs
#from old:
# cDC1, cDC2, Lang+cDC2

ref <- grep("CDC1", oldmeta$group)
ref <- c(ref, grep("CDC2", oldmeta$group))
ref <- c(ref, grep("Lang", oldmeta$group))

oldmetacut <- oldmeta[ref,]
oldcpm <- oldcpm[,ref]
oldlcpm <- oldlcpm[,ref]
oldmetacut <- oldmetacut[,-c(1,2)]
oldmetacut <- oldmetacut[,c(2,1,3)]
oldmetacut$data <- 'old'

oldmetacut$group <- ifelse(oldmetacut$group == 'CDC1', "cDC1", #to match the new meta format
                           ifelse(oldmetacut$group == 'CDC2', "cDC2LangN", "cDC2LangP")) 
colnames(oldmetacut) <- c("ident", 'group', 'tissue', 'data')
#from new:
# cDC1, cDC2LangN, cDC2LangP
ref <- grep("cDC1", newmeta$group)
ref <- c(ref, grep("cDC2LangN", newmeta$group))
ref <- c(ref, grep("cDC2LangP", newmeta$group))

newmetacut <- newmeta[ref,]
newcpm <- newcpm[,ref]
newlcpm <- newlcpm[,ref]
newmetacut$data <- 'new'


#function
cpmBox <- function(gene, heading = "Bulk RNA Expression (cpm) - N=4") {
  tempold <- oldmetacut
  tempold$gene <- oldcpm[gene,]
  tempnew <- newmetacut
  tempnew$gene <- newcpm[gene,]
  temp <- rbind(tempold, tempnew)
  
  ggplot(temp, aes(x=group, y=gene, colour=data)) + ggtitle(heading)+geom_boxplot()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()
}
lcpmBox <- function(gene, heading = "Bulk RNA Expression (lcpm) - N=4") {
  tempold <- oldmetacut
  tempold$gene <- oldlcpm[gene,]
  tempnew <- newmetacut
  tempnew$gene <- newlcpm[gene,]
  temp <- rbind(tempold, tempnew)
  
  ggplot(temp, aes(x=group, y=gene, colour=data)) + ggtitle(heading)+geom_boxplot()+ylab(paste(gene, 'Expression'))+ xlab("") +theme_classic()+RotatedAxis()
}

