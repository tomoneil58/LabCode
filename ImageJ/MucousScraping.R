
setwd("~/Desktop/CSV/")
time1 <- Sys.time()
temp <- read.csv(list.files()[1])
label=nchar(temp$Label)-19
df <- data.frame(
  Sample = substr(temp$Label[1], 1,label), 
  Area = temp$Area[1],
  BlueIntensityTotal = temp$Mean[1],
  BlueIntensityThreshold = temp$Mean[2],
  BlueIntensityRatio = temp$Mean[2]/temp$Mean[1],
  BlueArea = temp$Area[2],
  BlueAreaPercent = 100*(temp$Area[2]/temp$Area[1]),
  PurpleIntensityTotal = temp$Mean[3],
  PurpleIntensityThreshold = temp$Mean[4],
  PurpleIntensityRatio = temp$Mean[4]/temp$Mean[3],
  PurpleArea = temp$Area[4],
  PurpleAreaPercent = 100*(temp$Area[4]/temp$Area[3]), 
  BlueAreaOverPurpleArea = temp$Area[2]/temp$Area[4], 
  Error = ifelse(temp$Area[4] >= temp$Area[1], "Error", "")
)

for(i in 2:length(list.files())){
  temp <- read.csv(list.files()[i])
  label=nchar(temp$Label)-19
  data <- data.frame(
    Sample = substr(temp$Label[1], 1,label), 
    Area = temp$Area[1],
    BlueIntensityTotal = temp$Mean[1],
    BlueIntensityThreshold = temp$Mean[2],
    BlueIntensityRatio = temp$Mean[2]/temp$Mean[1],
    BlueArea = temp$Area[2],
    BlueAreaPercent = 100*(temp$Area[2]/temp$Area[1]),
    PurpleIntensityTotal = temp$Mean[3],
    PurpleIntensityThreshold = temp$Mean[4],
    PurpleIntensityRatio = temp$Mean[4]/temp$Mean[3],
    PurpleArea = temp$Area[4],
    PurpleAreaPercent = 100*(temp$Area[4]/temp$Area[3]),
    BlueAreaOverPurpleArea = temp$Area[2]/temp$Area[4],
    Error = ifelse(temp$Area[4] >= temp$Area[1], "Error", "")
    
  )
  df <- rbind(df, data)
}
meta <-read.csv("Meta1.csv")
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

time2 <- Sys.time()
time2-time1

ggplot(df, aes(x=Group, y=BlueAreaPercent))+geom_boxplot()
ggplot(df, aes(x=BlueIntensityThreshold, y=BlueAreaPercent, shape=Mice))+
  geom_point(alpha=0.5)+geom_smooth(method='lm', formula= y~x, aes(color= Mice))+theme_bw()

ggplot(df, aes(x=BlueIntensityThreshold, y=BlueAreaPercent, color=Group))+
  geom_point(alpha=0.9)+geom_smooth(method='lm', formula= y~x)+theme_bw()

ggplot(df, aes(Order, BlueAreaPercent)) + geom_point() + facet_wrap(~Group)

#averaging

df2 <- data.frame(Group = NA, Mice = NA, Intensity = NA, AreaPc=NA)

df$uniqueGroup <- paste0(df$Group, df$Mice)

unique(df$uniqueGroup)

for (i in 1:length(unique(df$uniqueGroup))) {
  temp <-df[df$uniqueGroup == unique(df$uniqueGroup)[i],]
  tempdf <- data.frame(Group = temp$Group[1],
    Mice = temp$Mice[1],
    Intensity = mean(temp$BlueIntensityThreshold),
    AreaPc = mean(temp$BlueAreaPercent)
  )
  df2 <-rbind(df2, tempdf)
}

df2 <- df2[-1,]
ggplot(df2, aes(Group, AreaPc))+geom_point(aes(color=Mice))
ggplot(df2, aes(Intensity, AreaPc, color=Group)) + geom_point()#+geom_smooth(method='lm', formula= y~x)

write.csv(df,"results.csv");write.csv(df2, "averageresults.csv")


