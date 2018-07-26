library(plyr)

library(ggplot2)

ddd <- read.csv("DDD-16nGroups.mem", header=F)
sss <- read.csv("SSS-32nGroups.mem", header=F)

names(ddd) <- "memory"
names(sss) <- "memory"
#u <- max(nrow(df1), nrow(df2))

ddd$Variant <- "double"
sss$Variant <- "single"
ddd$time <- seq(1, nrow(ddd))
sss$time <- seq(1, nrow(sss))

memdf <- rbind(ddd[1:100,], sss[1:100,])

head(sss)
head(ddd)
#y <- ggplot(data=ddd, aes(x=time, y=memory))+geom_point()

memplot <- ggplot(data=memdf,aes(x=time, y = memory, group=Variant, color=Variant))+geom_point()+
#geom_line(data=df1, aes(y=dbl, color="double"))+
#geom_line(data=df2, aes(y=single, color="single"))+
ggtitle("Memory use over time for 16 nGroups")

ggsave("memplot.pdf", memplot)




