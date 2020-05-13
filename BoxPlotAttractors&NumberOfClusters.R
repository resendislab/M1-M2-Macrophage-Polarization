BasinMacrophages<-read.table(file.choose(), header=T)
library(RColorBrewer)
library(ggplot2)
colourCount.stability = length(unique(BasinMacrophages$Phenotypes))
getPalette2 = colorRampPalette(brewer.pal(9, "Set1"))
p <- ggplot(BasinMacrophages, aes(x=Phenotypes, y=Value, fill=Phenotypes )) + geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired")  +
  theme(legend.position = "none") +
  labs( y = "Basin of Attraction")
  q <- p + coord_flip() 
  q
  svg(file="StabilityAttractorsAll.pdf")
  q
  dev.off()
  
  
##############################################
#            Basin of attraction with        # 
#           HIF1A =0 and NFKB=1              #
#                                            #
#                                            # 
##############################################
BasinMacrophagesNoHIF1AandNFKB<-read.table(file.choose(), header=T)
library(ggplot2)
t <- ggplot(BasinMacrophagesNoHIF1AandNFKB, aes(x=Phenotypes, y=Value, fill=Phenotypes )) + geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired")  +
  theme(legend.position = "none") + theme(text=element_text(family="Arial", size=10))+
  labs( y = "Basin of Attraction")
  u <- t + coord_flip()
  ggsave("Fig6A.png", u, width=4, height=4)
  
##############################################
#            Number of Clusters              # 
#                                            #
#                                            #
#                                            # 
##############################################
NumberOfClustersMacrophages<-read.table(file.choose(), header=T)
NumberOfClustersMacrophages <- data.frame(
  Clusters=c(10,11,12,13,23,25,28,30),
  Frecuency=c(6,2,4,2,4,1,1,3)
)
Clusters.Colors <- brewer.pal(8, "Set2")
barplot(height=NumberOfClustersMacrophages$Frecuency, names=NumberOfClustersMacrophages$Clusters, 
        col=Clusters.Colors,
        xlab="Number of Clusters", 
        ylab="Frecuency among all indices",
        names.arg = c("10", "11", "12", "13", "23", "25", "28", "30") ,
        main="Best Partition Analysis with K-Means", 
        ylim=c(0,6)
        )
