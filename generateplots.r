library(ggplot2)

#df <- read.csv("Data_for_1Node_and_8_Node_runs/averaged_values.csv", header=TRUE, sep=",")
df <- read.csv("../tycho2results_plotter.csv", header=TRUE, sep=",")

df$MPI <- sapply(df$MPI, as.factor)
#df$Nodes- sapply(df$Nodes, as.numeric)
l2plot <- ggplot(data =df, aes(x = Nodes, y=L2RelErr, shape= MPI, colour=Variant, group=interaction(MPI, Variant)))+geom_line(aes(linetype=Variant))+geom_point()+ggtitle("Number of Nodes versus L2 Error to Answer")+
scale_x_log10(breaks=df$Nodes, labels=df$Nodes)
              
timeplot <- ggplot(data =df, aes(x = Nodes, y=TotSourceIterTime, shape= MPI, colour=Variant, group=interaction(MPI, Variant)))+geom_line(aes(linetype=Variant))+geom_point()+ggtitle("Number of Nodes versus Total Source Iteration Time to Answer")+
scale_x_log10(breaks=df$Nodes, labels=df$Nodes)

ggsave("l2errors.pdf", l2plot)
ggsave("runtime.pdf", timeplot)
