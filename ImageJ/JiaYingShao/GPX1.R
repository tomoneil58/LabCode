
setwd("~/Desktop/CSV/")
time1 <- Sys.time()
temp <- read.csv(list.files()[1])
label=nchar(temp$Label)-20

df <- data.frame(
  Sample = substr(temp$Label[1], 1,label), 
  EpiArea = temp$Area[3],
  EpiBrown = temp$Mean[3], 
  EpiBrownSum = temp$IntDen[3]
)

for(i in 2:length(list.files())){
  temp <- read.csv(list.files()[i])
  label=nchar(temp$Label)-20
  data <- data.frame(
    Sample = substr(temp$Label[1], 1,label), 
    EpiArea = temp$Area[3],
    EpiBrown = temp$Mean[3], 
    EpiBrownSum = temp$IntDen[3]
    
  )
  df <- rbind(df, data)
}
grep("Cont", df$Sample)
grep("Met", df$Sample)
grep("Nano", df$Sample)
grep("NaSe", df$Sample)

group <- c(rep("Control", each=14), rep("Met", each=6), rep("Nano", each=15), rep("NaSe", each=6))
df$Group <- group

mice <- c()
for(i in 1:41) {
  x <- df$Sample
  label=nchar(df$Sample[i]) - 5
  df$Mice[i]<-substr(df$Sample[i], 1,label)
}

df$OD <- log(255/(255-df$EpiBrown))
ggplot(df, aes(Group, OD))+geom_boxplot()+geom_point()


write.csv(df, "GPX1.csv")







