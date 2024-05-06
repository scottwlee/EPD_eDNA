############################################################
############################################################
### Package: NetSwan
.libPaths()

# ##### Install packages
# install.packages("igraph")
# install.packages("NetSwan")

library(igraph)
library(NetSwan)
library(dplyr)

getwd()
setwd("D:/")

############################################################
############################################################
### Southern region Edge-table
elec2 <- read.csv("Southern_EdgeTable2.csv") %>% 
  select(NF,NT) %>% 
  as.matrix()

head(elec2)
class(elec2)
View(elec2)

gra<-graph.edgelist(elec2, directed=FALSE)
f<-swan_efficiency(gra)
vertex_attr(gra, "efficiency_loss", index = V(gra))<-f
vertex_attr(gra)
f2<-swan_closeness(gra)
bet<-betweenness(gra)
reg<-lm(bet~f2)
summary(reg)
f3<-swan_connectivity(gra)
f4<-swan_combinatory(gra,10)

### Eastern region Edge-table
elec3 <- read.csv("Eastern_EdgeTable2.csv") %>% 
  select(NF,NT) %>% 
  as.matrix()

head(elec3)
class(elec3)
View(elec3)

gra2<-graph.edgelist(elec3, directed=FALSE)
f5<-swan_efficiency(gra2)
vertex_attr(gra2, "efficiency_loss", index = V(gra2))<-f5
vertex_attr(gra2)
f6<-swan_closeness(gra2)
bet2<-betweenness(gra2)
reg2<-lm(bet2~f6)
summary(reg2)
f7<-swan_connectivity(gra2)
f8<-swan_combinatory(gra2,10)

############################################################
############################################################
##### Plot a figure
robustness<-plot(f4[,1],f4[,4], type='o', col='dark green',xlab="Fraction of nodes removed",
                 ylab="Connectivity loss", cex=0.5)
lines(f4[,1],f4[,2], type='o', col='green', cex=0.5)
lines(f8[,1],f8[,4], type='o', col='brown', cex=0.5)
lines(f8[,1],f8[,2], type='o', col='orange', cex=0.5)
legend('bottomright',c("Cascading Light-yellow (Eastern)", "Betweenness Light-yellow (Eastern)", "Cascading Dark-green (Southern)", "Betweenness Dark-green (Southern)"),
       lty=c(1,1,1,1), pch=c(1,1,1,1),
       col=c("brown","orange","dark green", "green"))

############################################################
############################################################
# ##### The following command: Extract image in "tiff" format
# tiff("Robustness_final.tiff", width = 6, height = 6, units = "in", res = 600)
# 
# robustness <- plot(f4[,1],f4[,4], type='o', col='dark green',xlab="Fraction of nodes removed",
#                    ylab="Connectivity loss", cex=0.5)
# lines(f4[,1],f4[,2], type='o', col='green', cex=0.5)
# lines(f8[,1],f8[,4], type='o', col='brown', cex=0.5)
# lines(f8[,1],f8[,2], type='o', col='orange', cex=0.5)
# legend('bottomright',c("Cascading Light-yellow (Eastern)", "Betweenness Light-yellow (Eastern)", "Cascading Dark-green (Southern)", "Betweenness Dark-green (Southern)"),
#        lty=c(1,1,1,1), pch=c(1,1,1,1),
#        col=c("brown","orange","dark green", "green"))
# dev.off()