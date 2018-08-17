library(ggplot2)

#df <- read.csv("Data_for_1node_and_8_node_runs/averaged_values.csv", header=TRUE, sep=",")
df <- read.csv("../tycho2_plotter.csv", header=TRUE, sep=",")

#ggplot(data =df, aes(x = nodes, y=L2RelErr, group=1))+geom_line()+geom_point()
#ggplot(data =df, aes(nodes)+ geom_line(aes(y=L2RelErr, colour="red")) + geom_line(aes(y=L2RelErr, colour="blue"))
df$MPI <- sapply(df$MPI, as.factor)
#df$nodes- sapply(df$nodes, as.numeric)
l2plot <- ggplot(data =df, aes(x = nodes, y=L2RelErr, colour=precision))+geom_point(size=4)+ggtitle("Number of nodes versus L2 Error to Answer")+
scale_x_log10(breaks=df$nodes, labels=df$nodes)+
geom_line(aes(linetype=precision))
              
timeplot <- ggplot(data =df, aes(x = nodes, y=TotSourceIterTime, colour=precision))+geom_point(size=4)+ggtitle("Number of nodes versus total source iteration time to answer")+
scale_x_log10(breaks=df$nodes, labels=df$nodes)+
geom_line(aes(linetype=precision))+
ylab("time(sec)")


ggsave("l2errors-1e5.pdf", l2plot)
ggsave("runtime-1e5.pdf", timeplot)

dfs <- split(df, df$precision)
z<- sapply(dfs, function(x) tapply(x$TotSourceIterTime, x$nodes, max))

print(z)
