
setwd("~/Desktop/CSV/")
time1 <- Sys.time()
temp <- read.csv(list.files()[1])
label=nchar(temp$Label[6])-23
df <- data.frame(
  Sample = substr(temp$Label[6], 1,label), 
  TotalArea = temp$Area[2], 
  TotalDABMean = temp$Mean[4], 
  TotalNuclearMean = temp$Mean[2],
  LenEpi = temp$Area[5], 
  LenSub = temp$Area[6],
  
  DABArea = temp$Area[3],
  DABmean = temp$Mean[3],
  DABInt = temp$IntDen[3], 
  DABpc = 100*(temp$Area[3]/temp$Area[2]), 
  DABMeanToTotal = temp$Mean[3]/temp$Mean[4],
  
  NuclearArea = temp$Area[1], 
  NuclearMean = temp$Mean[1],
  NuclearInt = temp$IntDen[1], 
  NuclearPc = 100*(temp$Area[1]/temp$Area[2]),
  NuclearMeanToTotal = temp$Mean[1]/temp$Mean[2],
  
  DABtoNuclear = 100*(temp$Area[3]/temp$Area[1]),
  NuclearAreaToEpi = temp$Area[1]/temp$Area[5],
  NuclearAreaToSub = temp$Area[1]/temp$Area[6],
  DABAreaToEpi = temp$Area[3]/temp$Area[5],
  DABAreaToSub = temp$Area[3]/temp$Area[6],
  
  TotalAreaToEpi = temp$Area[2]/temp$Area[5],
  TotalAreaToSub = temp$Area[2]/temp$Area[6],
  
  Error = ifelse(temp$Area[1] >= temp$Area[4], "Error", "")
)


for(i in 2:length(list.files())){
  temp <- read.csv(list.files()[i])
  label=nchar(temp$Label[6])-23
  data <- data.frame(
    Sample = substr(temp$Label[6], 1,label), 
    TotalArea = temp$Area[2], 
    TotalDABMean = temp$Mean[4], 
    TotalNuclearMean = temp$Mean[2],
    LenEpi = temp$Area[5], 
    LenSub = temp$Area[6],
    
    DABArea = temp$Area[3],
    DABmean = temp$Mean[3],
    DABInt = temp$IntDen[3], 
    DABpc = 100*(temp$Area[3]/temp$Area[2]), 
    DABMeanToTotal = temp$Mean[3]/temp$Mean[4],
    
    NuclearArea = temp$Area[1], 
    NuclearMean = temp$Mean[1],
    NuclearInt = temp$IntDen[1], 
    NuclearPc = 100*(temp$Area[1]/temp$Area[2]),
    NuclearMeanToTotal = temp$Mean[1]/temp$Mean[2],
    
    DABtoNuclear = 100*(temp$Area[3]/temp$Area[1]),
   
    NuclearAreaToEpi = temp$Area[1]/temp$Area[5],
    NuclearAreaToSub = temp$Area[1]/temp$Area[6],
    DABAreaToEpi = temp$Area[3]/temp$Area[5],
    DABAreaToSub = temp$Area[3]/temp$Area[6],
    
    TotalAreaToEpi = temp$Area[2]/temp$Area[5],
    TotalAreaToSub = temp$Area[2]/temp$Area[6],
    
    Error = ifelse(temp$Area[1] >= temp$Area[4], "Error", "")
    
  )
  df <- rbind(df, data)
}

meta <-read.csv("meta.csv")
meta$Image <- paste0(meta$Image, ".tiff")
for( i in 1:nrow(df)) {
  if (df$Sample[i] %in% meta$Image) {
    df$Cage[i] <- meta$Cage[meta$Image == df$Sample[i]]
    df$Group[i] <-  meta$Group[meta$Image == df$Sample[i]]
    df$Mice[i] <- meta$Mice[meta$Image == df$Sample[i]]
    df$Order[i] <- meta$Order[meta$Image == df$Sample[i]]
  }else{
    df$Error[i] <- 'Error'
  }
}
message("There are ",nrow(df)-table(df$Error[df$Error == '']), " Errors")
time2 <- Sys.time()
time2-time1

ggplot(df, aes(Group, DABtoNuclear, color=Group))+geom_boxplot()+geom_boxplot()+geom_point(aes(shape=Mice))
ggplot(df, aes(Mice, DABmean))+geom_boxplot()+geom_point()+facet_wrap(~Group)
ggplot(df, aes(Group,DABtoNuclear))+geom_boxplot()+geom_point()
ggplot(df, aes(Group,DABpc))+geom_boxplot()+geom_point()
ggplot(df, aes(Mice, TotalAreaToSub))+geom_boxplot()+geom_point()+facet_wrap(~Group)
ggplot(df, aes(Mice, NuclearPc))+geom_boxplot()+geom_point()+facet_wrap(~Group)
ggplot(df, aes(Mice, NuclearAreaToSub))+geom_boxplot()+geom_point()+facet_wrap(~Group)
ggplot(df, aes(as.factor(Order), DABAreaToSub))+geom_boxplot()+geom_point()+facet_wrap(~Group)
ggplot(df, aes(Mice, DABInt))+geom_boxplot()+geom_point()+facet_wrap(~Group)

ggplot(df, aes(as.factor(Order), DABInt/LenSub))+geom_boxplot()+geom_line(aes(group=Mice))+facet_wrap(~Group)


ggplot(df, aes(Mice,NuclearArea))+geom_boxplot()+geom_point()+facet_wrap(~Group)













