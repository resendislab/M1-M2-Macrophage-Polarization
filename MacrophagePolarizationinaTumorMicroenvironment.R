#This code will let you simulate the polarization of macrophages in a tumor microenvironment. 
#Remember to download all txt files and csv files from github and save all of them in the same folder.
#Install packages needed in this work 
rm(list = ls())
install.packages(c("devtools","BoolNet", "BoolNetPerturb", "ggplot2", "dplyr", "magrittr", "ggrepel", "Rtsne", "mclust", "fpc"))
install_github("mar-esther23/boolnet-perturb")
##Once installed, invoked them
#####################################################################################
library(BoolNet)
library(BoolNetPerturb)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggrepel)
library(Rtsne)
library(mclust)
library(fpc)
####################################################################################
#Lets invoke the set of logical rules of figure 1 in the article. 
macrophage <- read.table(file = "PathWhereYouDownloadAllTheCode/macrofago2.txt", header = TRUE)
#Upload a data frame where the phenotypes of the macrophages are saved. 
Macrophage.Phenotypes<-data.frame(labels=c('M0', 'M1', 'M2a','M2b','M2c', 'M2d'), rules=c('!(NFKB | STAT1 |STAT3 |STAT6|HIF1A|ERK|Fra1|AP1)', 'NFKB|STAT1|TNFA & AP1| TNFAe & AP1', 'STAT6', '  (IL1B & AP1)| ERK','STAT3', '(TLR4 & A2a) | (Fra1 & AP1) | HIF1A| Fra1'))

#Once the data and the labels are in R, we call it macrophage. This is important for BoolNet
macrophage<-loadNetwork("macrofago2.txt")
#We obtain the attractors of the network. This part will take a while due to the method chosen. So grab a coffee and read our article.
attractors <- getAttractors(macrofago, method = "exhaustive", type="synchronous")
#Once we have the attractors, we want to put labels based on experimental evidence
Macrophages.labels <- labelAttractors(attractors, Macrophage.Phenotypes)
#We save the attractors in a data frame. 
Macrophages.Phenotypes<-attractorToDataframe(Macrophages.labels, Boolean = TRUE)
##We save it as a csv file in the same where all file are saved. 
write.table(Macrophage.Phenotypes, "atractores9deNov2020Boolean.txt")


###################################################################################
#We used the tsne to try to visualize the phenotypes in a two dimension space 
#Check where the script is kept, redirect the path to the script
source("PathWhereYouDownloadAllTheCode/Funn.R")
#Same here, put the path where the archive is being kept
Data <- read.table(file = "PathWhereYouDownloadAllTheCode/atractores9deNov2020Boolean.txt",
                   header = TRUE)
rnames <- as.vector(Data[,1])
Data <- Data[,-1]
Pheno.names <- unique(rnames)
n.Pheno <- length(Pheno.names)


Pheno.asig <- lapply(1:n.Pheno, FUN = function(x) which(rnames==Pheno.names[x]))
Pheno.colors <- sample(367:657,n.Pheno,replace = F)
Col <- vector(mode = "integer", length = length(rnames))
for (i in 1:n.Pheno){
  Col[Pheno.asig[[i]]] <- Pheno.colors[i] 
}
library(Rtsne)
set.seed(1)

tSNE.Macro <- Rtsne(Data, perplexity=30,max_iter=5000,verbose=0)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(tSNE.Macro$Y, pch=21,bg=Col)
legend("bottomleft",inset=c(1,.25), title = "Phenotypes",
       Pheno.names,
       pch = rep(21,n.Pheno),
       pt.bg = Pheno.colors,
       bty = "n",
       pt.cex = 1,
       cex = 1
)
##### kmeans utilizando consenso para determinar el número de clusters se 
##### puede cargar el análsisi para no volverloa  correr

# library(NbClust)
# Res <- NbClust(data = tSNE.Macro$Y, distance = "euclidean",diss = NULL, 
#                min.nc = 10, max.nc = 30, 
#                method = "kmeans", index = "all")

### Número de clusters
n <- max(Res$Best.partition)

Cluster.asig <- lapply(1:n, FUN = function(x) which(x==Res$Best.partition))
### Asignación de color por clusters
Cluster.colors<-c(brewer.pal(name="Set1", n=5), brewer.pal(name="Set3", n=5))
Col <- vector(mode = "integer", length = length(rnames))

for (i in 1:n){
  Col[Cluster.asig[[i]]] <- Cluster.colors[i] 
}

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(tSNE.Macro$Y, pch=21,bg=Col)
legend("bottomleft",inset=c(1,.25), title = "Clusters",
       LETTERS[1:n],
       pch = rep(21,n),
       pt.bg =Cluster.colors,
       bty = "n",
       pt.cex = 1,
       cex = 1
)

##### Cálculo de proporciones
Prop <- Proportion(Res$Best.partition,
                   n,rnames,n.Pheno,Pheno.names,tot)
##### By Clusters
MyplotBar(Prop[[1]],"Kmeans_Proportion_Macrophages_by_Clusters",
          "Groups percentage","Phenotypes",Pheno.colors)
MyplotBar(Prop[[2]],"Kmeans_Proportion_Macrophages_by_clusters",
          "Sample percentage","Phenotypes",Pheno.colors)

##### By Phenotypes

MyplotBar(t(Prop[[1]]),"Kmeans_Proportion_Macrophages_by_Phenotypes",
          "Groups percentage","Clusters",Cluster.colors)
MyplotBar(t(Prop[[2]]),"Kmeans_Proportion_Macrophages_by_Phenotypes",
          "Sample percentage","Clusters",Cluster.colors)
          
#### We evaluate the knockout scenario of the nodes or transcriptional factors
#Now we simulate a KNOCK-OUT of transcriptional Factors 
KOAP1.macrofago<-fixGenes(macrophage,"AP1", 0)
KOAP1.attr<-getAttractors(KOAP1.macrofago)
KOAP1.labels<-labelAttractors(KOAP1.attr,Macrophage.Phenotypes)
table(KOAP1.labels)  #We write the number of phenotypes
#########################################################################
KOErk.macrofago<-fixGenes(macrophage,"ERK", 0)
KOErk.attr<-getAttractors(KOErk.macrofago)
KOErk.labels<-labelAttractors(KOErk.Macrophage.Phenotypes)
table(KOErk.labels) #We write the number of phenotypes
###########################################################################
KOFra1.macrofago<-fixGenes(macrofago,"Fra1", 0)
KOFra1.attr<-getAttractors(KOFra0.macrofago)
KOFra1.labels<-labelAttractors(KOFra0.attr,Macrophage.Phenotypes)
table(KOFra1.labels) #We write the number of phenotypes
########################################################################################
KOHIF1A.macrofago<-fixGenes(macrofago,"HIF1A", 0)
KOHIF1A.attr<-getAttractors(KOHIF0A.macrofago)
KOHIF1A.labels<-labelAttractors(KOHIF0A.attr, Macrophage.Phenotypes)
table(KOHIF1A.labels)   #We write the number of phenotypes
###############################################################################################
KOIL10.macrofago<-fixGenes(macrofago,"IL10", 0)
KOIL10.attr<-getAttractors(KOIL10.macrofago)
KOIL10.labels<-labelAttractors(KOIL10.attr,Macrophage.Phenotypes)
table(IL10.labels)
##########################################################################################################
KONFKB.macrofago<-fixGenes(macrofago,"NFKB", 0)
KONFKB.attr<-getAttractors(KONFKB.macrofago)
KONFKB.labels<-labelAttractors(KONFKB.attr, Macrophage.Phenotypes)
table(KONFKB.labels) #We write the number of phenotypes
############################################################################################################
KOSMAD23.macrofago<-fixGenes(macrophage,"SMAD23", 0)
KOSMAD23.attr<-getAttractors(KOSMAD23.macrofago)
KOSMAD23.labels<-labelAttractors(KOSMAD23.attr, Macrophage.Phenotypes)
table(KOSMAD23.labels)
#################################################################################################################
KOSOCS3.macrofago<-fixGenes(macrophage,"SOCS3", 0)
KOSOCS3.attr<-getAttractors(KOSOCS0.macrofago, method = "exhaustive", type="synchronous")
KOSOCS3.labels<-labelAttractors(KOSOCS0.attr, Macrophage.Phenotypes)
table(KOSOCS3.labels)
#################################################################################################################
KOSOCS1.macrofago<-fixGenes(macrofago,"SOCS1", 0)
KOSOCS1.attr<-getAttractors(KOSOCS1.macrofago, method = "exhaustive", type="synchronous")
KOSOCS1.labels<-labelAttractors(KOSOCS1.attr, Macrophage.Phenotypes)
table(KOSOCS1.labels)
#######################################################################################################################
KOSTAT1.macrofago<-fixGenes(macrophage,"STAT1", 0)
KOSTAT1.attr<-getAttractors(KOSTAT1.macrofago, method = "exhaustive", type="synchronous")
KOSTAT1.labels<-labelAttractors(KOSTAT1.attr, Macrophage.Phenotypes)
table(KOSTAT1.labels)
##############################################################################################################################
KOSTAT3.macrofago<-fixGenes(macrophage,"STAT3", 0)
KOSTAT3.attr<-getAttractors(KOSTAT3.macrofago, method = "exhaustive", type="synchronous")
KOSTAT3.labels<-labelAttractors(KOSTAT3.attr, Macrophage.Phenotypes)
table(KOSTAT3.labels)
################################################################################################################################
KOSTAT6.macrofago<-fixGenes(macrophage,"STAT6", 0)
KOSTAT6.attr<-getAttractors(KOSTAT6.macrofago, method = "exhaustive", type="synchronous")
KOSTAT6.labels<-labelAttractors(KOSTAT6.attr, Macrophage.Phenotypes)
table(KOSTAT6.labels)
#################################################################################################################################
KOTGFB.macrofago<-fixGenes(macrophage,"TGFB", 0)
KOTGFB.attr<-getAttractors(KOTGFB.macrofago, method = "exhaustive", type="synchronous")
KOTGFB.labels<-labelAttractors(KOTGFB.attr, Macrophage.Phenotypes)
table(KOTGFB.labels)
##################################################################################################################################
KOTNFA.macrofago<-fixGenes(macrophage,"TNFA", 0)
KOTNFA.attr<-getAttractors(KOTNFA.macrofago, method = "exhaustive", type="synchronous")
KOTNFA.labels<-labelAttractors(KOTNFA.attr, Macrophage.Phenotypes)
table(KOTNFA.labels)


#### We simulate an over expression of the transcriptional factors in our model
##
OverAP1.macrofago<-fixGenes(macrophage,"AP1", 1)
OverAP1.attr<-getAttractors(OverAP1.macrofago, method = "exhaustive", type="synchronous")
OverAP1.labels<-labelAttractors(OverAP1.attr, Macrophage.Phenotypes)
table(OverAP1.labels)
#########################################################################
OverErk.macrofago<-fixGenes(macrophage,"ERK", 1)
OverErk.attr<-getAttractors(OverErk.macrofago, method = "exhaustive", type="synchronous")
OverErk.labels<-labelAttractors(OverErk.attr, Macrophage.Phenotypes)
table(OverErk.labels)
###########################################################################
OverFra1.macrofago<-fixGenes(macrophage,"Fra1", 1)
OverFra1.attr<-getAttractors(OverFra1.macrofago, method = "exhaustive", type="synchronous")
OverFra1.labels<-labelAttractors(OverFra1.attr, Macrophage.Phenotypes)
table(OverFra1.labels)
########################################################################################
OverHIF1A.macrofago<-fixGenes(macrophage,"HIF1A", 1)
OverHIF1A.attr<-getAttractors(OverHIF1A.macrofago, method = "exhaustive", type="synchronous")
OverHIF1A.labels<-labelAttractors(OverHIF1A.attr, Macrophage.Phenotypes)
table(OverHIF1A.labels)
###############################################################################################
OverIL10.macrofago<-fixGenes(macrophage,"IL10", 1)
OverIL10.attr<-getAttractors(OverIL10.macrofago, method = "exhaustive", type="synchronous")
OverIL10.labels<-labelAttractors(OverIL10.attr, Macrophage.Phenotypes)
table(OverIL10.labels)
##########################################################################################################
OverNFKB.macrofago<-fixGenes(macrophage,"NFKB", 1)
OverNFKB.attr<-getAttractors(OverNFKB.macrofago, method = "exhaustive", type="synchronous")
OverNFKB.labels<-labelAttractors(OverNFKB.attr, Macrophage.Phenotypes)
table(OverNFKB.labels)
############################################################################################################
OverSMAD23.macrofago<-fixGenes(macrophage,"SMAD23", 1)
OverSMAD23.attr<-getAttractors(OverSMAD23.macrofago, method = "exhaustive", type="synchronous")
OverSMAD23.labels<-labelAttractors(OverSMAD23.attr, Macrophage.Phenotypes)
table(OverSMAD23.labels)
#################################################################################################################
OverSOCS1.macrofago<-fixGenes(macrofago,"SOCS1", 1)
OverSOCS1.attr<-getAttractors(OverSOCS1.macrofago, method = "exhaustive", type="synchronous")
OverSOCS1.labels<-labelAttractors(OverSOCS1.attr, Macrophage.Phenotypes)
table(OverSOCS1.labels)
#################################################################################################################
OverSOCS3.macrofago<-fixGenes(macrophage,"SOCS3", 1)
OverSOCS3.attr<-getAttractors(OverSOCS3.macrofago, method = "exhaustive", type="synchronous")
OverSOCS3.labels<-labelAttractors(OverSOCS3.attr, Macrophage.Phenotypes)
table(OverSOCS3.labels)
#######################################################################################################################
OverSTAT1.macrofago<-fixGenes(macrophage,"STAT1", 1)
OverSTAT1.attr<-getAttractors(OverSTAT1.macrofago, method = "exhaustive", type="synchronous")
OverSTAT1.labels<-labelAttractors(OverSTAT1.attr, Macrophage.Phenotype)
table(OverSTAT1.labels)
##############################################################################################################################
OverSTAT3.macrofago<-fixGenes(macrophage,"STAT3", 1)
OverSTAT3.attr<-getAttractors(OverSTAT3.macrofago, method = "exhaustive", type="synchronous")
OverSTAT3.labels<-labelAttractors(OverSTAT3.attr, Macrophage.Phenotypes)
table(OverSTAT3.labels)
################################################################################################################################
OverSTAT6.macrofago<-fixGenes(macropage,"STAT6", 1)
OverSTAT6.attr<-getAttractors(OverSTAT6.macrofago, method = "exhaustive", type="synchronous")
OverSTAT6.labels<-labelAttractors(OverSTAT6.attr, Macrophage.Phenotypes)
table(OverSTAT6.labels)
#################################################################################################################################
OverTGFB.macrofago<-fixGenes(macrophage,"TGFB", 1)
OverTGFB.attr<-getAttractors(OverTGFB.macrofago, method = "exhaustive", type="synchronous")
OverTGFB.labels<-labelAttractors(OverTGFB.attr, Macrophage.Phenotypes)
table(OverTGFB.labels)
##################################################################################################################################
OverTNFA.macrofago<-fixGenes(macrophage,"TNFA", 1)
OverTNFA.attr<-getAttractors(OverTNFA.macrofago, method = "exhaustive", type="synchronous")
OverTNFA.labels<-labelAttractors(OverTNFA.attr, Macrophage.Phenotypes)
table(OverTNFA.labels)

###Once we have all the information from the knockouts or overexpression in a csv file, we uploaded and proceed to the heatmap

#####Código para obtener figuras del artículo. 
#Figura 3 
###############################################################
#                                                             #
#              HeatMap of the Knock-Outs of                   #
#                 Transcriptional Factors                     #
###############################################################
library(gplots)
library(RColorBrewer)
HeatMapKOTranscriptional<- read.csv("PathWhereYouDownloadAllTheCode/HeatMapKnockOutTranscriptionalFactors.csv", comment.char="#")
#Obtaining the phenotypes from the csv format
Phenotypes <- HeatMapKOTranscriptional[,1] 
#Transforming the matrix into a data fram and naming rownames with the phenotypes obtained earlier                           
mat_HeatMapKOTranscriptional<-data.frame(HeatMapKOTranscriptional[,2:ncol(HeatMapKOTranscriptional)])
rownames(mat_HeatMapKOTranscriptional)<- Phenotypes
#Customizing our scale for the heatmap
MyScale<-function(a,...){
a.min<-min(a,...)
a.max<-max(a,...)
d<-a.max-a.min
A<-(a)/(log10(d+1.0001))
return(A)
}
#Scaling the data frame
scaleKOTransc = lapply(mat_HeatMapKOTranscriptional, MyScale, na.rm = T)
#Transforming the list obtain earlier into a data frame
mat_HeatMapKOTranscriptional<-data.frame(AP1=unlist(scaleKOTransc[1]),ERK=unlist(scaleKOTransc[2]), Fra1=unlist(scaleKOTransc[3]), HIF1A=unlist(scaleKOTransc[4]),NFKB=unlist(scaleKOTransc[5]), IL10=unlist(scaleKOTransc[6]), SMAD23=unlist(scaleKOTransc[7]), SOCS3=unlist(scaleKOTransc[8]), SOCS1=unlist(scaleKOTransc[9]), STAT1=unlist(scaleKOTransc[10]), STAT3=unlist(scaleKOTransc[11]),STAT6=unlist(scaleKOTransc[12]), TGFB=unlist(scaleKOTransc[13]), TNFA=unlist(scaleKOTransc[14]))
rownames(mat_HeatMapKOTranscriptional)<- Phenotypes
#Verifing the obtain data frame with row names and new value from the scale
mat_HeatMapKOTranscriptional
#Transforming the data frame into a matrix for we could draw the heatmap
mat_HeatMapKOTranscriptional<-data.matrix(mat_HeatMapKOTranscriptional)
#Creating a new window where we are going to plot our heatmap of the OV of Extracellular factors
library("extrafont")
loadfonts()
pdf("Fig3B.pdf",family="Times New Roman", width=6.83, height=6)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
HeatMapKOextracel<- heatmap.2(mat_HeatMapKOTranscriptional,
   main = "Gene Deletion Intracellular", 
  notecol="black",      
  density.info="none",  
  trace="none",     
  col=redgreen(75),
  margins=c(12,8),
  na.color="gray",  
   cexRow=1, cexCol=1,    
  hclustfun=function(x)hclust(x,method="ward.D2"))
 par(op)               
  dev.off()
  HeatMapKOextracel

###############################################################
#                                                             #
#              HeatMap of the Over-expression of              #
#              Transcriptional Factors                        #
###############################################################
#Reading the data and tranforming it to a matrix 
HeatMapOverTransc<- read.csv("PathWhereYouDownloadAllTheCode/OverExpressionTranscriptionalFactors.csv", comment.char="#")
#Obtaining the phenotypes from the csv format
Phenotypes <- HeatMapOverTransc[,1] 
#Transforming the matrix into a data fram and naming rownames with the phenotypes obtained earlier                           
mat_HeatMapOverTranscriptional<-data.frame(HeatMapOverTransc[,2:ncol(HeatMapOverTransc)])
rownames(mat_HeatMapOverTranscriptional)<- Phenotypes
#Verifing the obtain data frame with row names and new value from the scale
mat_HeatMapOverTranscriptional
#Transforming the data frame into a matrix for we could draw the heatmap
mat_HeatMapOverTranscriptional<-data.matrix(mat_HeatMapOverTranscriptional)
dist_no_na <- function(mat) {
    edist <- dist(mat)
    edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
    return(edist)
}
#Creating a new window where we are going to plot our heatmap of the OV of Extracellular factors
library("extrafont")
loadfonts()
pdf("Fig3A.pdf",family="Times New Roman", width=6.83, height=6)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
HeatMapKOextracel<- heatmap.2(mat_HeatMapOverTranscriptional,
   main = "Gene Activation Intracellular", 
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins =c(12, 7), 
  cexRow=0.95, cexCol=0.95 ,    
  col=redgreen(75), 
  na.color="gray", 
  distfun=dist_no_na) 
  par(op) 
dev.off()

##########################################################################################
#Hheatmap comparing the basin of the microenvironment with the WT
#Reading the data and tranforming it to a matrix (Fig4)
#########################################################################################
library(gplots)
library(RColorBrewer)
HeatMapMircoE<- read.csv("PathWhereYouDownloadAllTheCode/Microambientes2Final.csv", comment.char="#")
Phenotypes <- HeatMapMircoE[,1]                            
mat_HeatMapMircoE <- data.matrix(HeatMapMircoE[,2:ncol(HeatMapMircoE)])  
rownames(mat_HeatMapMircoE)<- Phenotypes
 #We create a png image
library("extrafont")
loadfonts()
pdf("PathWhereYouDownloadAllTheCod/Fig4.pdf",family="Arial", width=6.83, height=6)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1) 
HeatMapMicroE<- heatmap.2(mat_HeatMapMircoE,
   main = "Microenvironments", 
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins =c(12, 7),     
  col=redgreen(75),  
  na.color="gray", 
   cexRow=.95, cexCol=.95,    
   distfun=dist_no_na)
par(op)
dev.off()

##########################################
#         Cell Fate Map of               #
#         Macrophage Polarization        #
##########################################
##We took the cell fate map function from BoolNetPerturb and ajusted it to our analysis of macrophage polarization. The only modification we
## did was on the method obtaining the attractors
CellFateMap <- function(net,states,genes,time=1, label.rules, ...) {    
    # validate and select random if default
    if (!is(net, "BooleanNetwork")) { stop("Error: non-valid network") }
    if (missing(states)) { #run for all attractors
        states <- getAttractors(net, method="exhaustive")
        states <- attractorToDataframe(states, Boolean=T)
        states["attractor"] <-NULL
        states["state"] <- NULL
    }
    if (missing(genes)) {genes = net$genes}
    res <- apply(states, 1, function(s) {
        r <- sapply(genes, function(g) {
            perturbState(net,s,g,time,all.data=T,int=T) 
        })
        t(data.frame(r))
    })
    res <- do.call("rbind", res)
    row.names(res) <- NULL
    res <- data.frame(res, row.names=NULL)
    
    if (!missing(label.rules)) { #label
        res["initial"] <- suppressWarnings( 
            labelList(res["initial"], label.rules,net$genes) )
        res["final"] <- suppressWarnings( 
            labelList(res["final"], label.rules,net$genes) )
    }       
    
    return(res)
}  
#This part will take a while to finish as well. 
CellFateMapMacrophagesPlasticity<-CellFateMap(macrophage, label.rules = Macrophage.Phenotypes, time=NULL) 

#################################################################
#                                                               #
# Theoretical pharmeceutical approach                           #
#################################################################
NfkbNoHif1a<-fixGenes(macrofago, c('NFKB', 'HIF1A'), c(1,0))
NfkbNoHif1a.attr<-getAttractors(NfkbNoHif1a, method = "exhaustive", type="synchronous")
NfkbNoHif1a.labels<-labelAttractors(NfkbNoHif1a.attr,labels.rules)
table(NfkbNoHif1a.labels)
## Cell fate mao of our theoretical pharmaceutical approach
CellFateMapTGEM<-CellFateMap(NfkbNoHif1a, label.rules = Macrophage.Phenotypes, time=NULL) 
##############################################
#            Basin of attraction with        # 
#           HIF1A =0 and NFKB=1 Fig 6A       #
#                                            #
#                                            # 
##############################################
BasinMacrophagesNoHIF1AandNFKB<- read.csv("PathWhereYouDownloadAllTheCode/BasinOfAttractionNfbkNoHif1a.csv", comment.char="#")
library(ggplot2)
library("extrafont")
loadfonts()
pdf("PathWhereYouDownloadAllTheCode/Fig6A.pdf",family="Arial", width=6.83, height=6)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggplot(BasinMacrophagesNoHIF1AandNFKB, aes(x=Phenotypes, y=Value, fill=Phenotypes )) + geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired")  +
  theme(legend.position = "none") + theme(text=element_text(family="Arial", size=10))+
  labs( y = "Basin of Attraction") + coord_flip()
par(op)
dev.off()
#################################################################
#                                                               #
# Theoretial approach in breast cancer microenvironment         #
#################################################################
#First let´s test this approach by applying a possible breast cancer microenvironment
# We are going to fixed th behaviour of IFNg and INFb as cero for all the simulations 
IL10TGFB<-fixGenes(NfkbNoHif1a, c('IFNG', 'IFNB', 'IL10e', 'TGFBe'), c(0,0,1,1))
IL10TGFB.attr<-getAttractors(IL10TGFB, method = "exhaustive", type="synchronous")
IL10TGFB.labels<-labelAttractors(IL10TGFB.attr,labels.rules)
table(IL10TGFB.labels)
######################################################################################
# Now let´s see the importance of IgG and adenosines in the behaviour of our model
IgGA2a<-fixGenes(NfkbNoHif1a, c('IFNG', 'IFNB', 'IgG', 'A2a'), c(0,0,1,1))
IgGA2a.attr<-getAttractors(IgGA2a, method = "exhaustive", type="synchronous")
IgGA2a.labels<-labelAttractors(IgGA2a.attr,labels.rules)
table(IgGA2a.labels)
#######################################################################################
HipGCGCR<-fixGenes(NfkbNoHif1a, c('IFNG', 'IFNB', 'Hipoxia', 'GCGCR'), c(0,0,1,1))
HipGCGCR.attr<-getAttractors(HipGCGCR, method = "exhaustive", type="synchronous")
HipGCGCR.labels<-labelAttractors(HipGCGCR.attr,labels.rules)
table(HipGCGCR.labels)
#######################################################################################
IL1BIL6<-fixGenes(NfkbNoHif1a, c('IFNG', 'IFNB', 'IL1B', 'IL6e'), c(0,0,1,1))
IL1BIL6.attr<-getAttractors(IL1BIL6, method = "exhaustive", type="synchronous")
IL1BIL6.labels<-labelAttractors(IL1BIL6.attr,labels.rules)
table(IL1BIL6.labels)
############################################################################################################################
#                                                                                                                          #
#                                          HeatMaps   fig7A. Breast Cancer Microenvironments                               #
############################################################################################################################
#Now let´s make the heatmap of the results of the microenvironment
library(gplots)
library(RColorBrewer)
HeatMapBreastCancerMicro<- read.csv("PathWhereYouDownloadAllTheCode/BreastCancerMacrophage.csv", comment.char="#")
Phenotypes <- HeatMapBreastCancerMicro[,1]                            
mat_HeatMapBreastCancerMicro <- data.matrix(HeatMapBreastCancerMicro[,2:ncol(HeatMapBreastCancerMicro)])  
rownames(mat_HeatMapBreastCancerMicro)<- Phenotypes
library("extrafont")
loadfonts()
pdf("Fig7A.pdf",family="Arial", width=6.83, height=6)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
HeatMapBreastCancer<- heatmap.2(mat_HeatMapBreastCancerMicro,
   main = " Breast Cancer Microenvironments", 
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins =c(11, 6),     
  col=redgreen(75),  
  cexRow=1.1, cexCol=1.1,    
  na.color="gray", 
  hclustfun=function(x)hclust(x,method="ward.D2")) 
par(op)
dev.off()
#####################################################################################
#                                                                                   #
#      Cell fate map of said microenvironments                                      #
#####################################################################################
CellFateMapTGEM.1<-CellFateMap(IL10TGFB, label.rules = Macrophage.Phenotypes, time=NULL) 
CellFateMapTGEM.2<-CellFateMap(IgGA2a, label.rules = Macrophage.Phenotypes, time=NULL) 
CellFateMapTGEM.3<-CellFateMap(HipGCGCR, label.rules = Macrophage.Phenotypes, time=NULL) 
CellFateMapTGEM.4<-CellFateMap(IL1BIL6, label.rules = Macrophage.Phenotypes, time=NULL) 

#################################################################
#                                                               #
#                Derrida Curves                                 #
#################################################################
##We took the derrida function from BoolNetPerturb and ajusted it to our analysis of macrophage polarization. The only modification we
## did was on the method obtaining the attractors
derrida <- function(net, repetitions=10000) {
  states <- getAttractors(net, method="exhaustive")
  repetitions = round(repetitions/length(net$genes))
  res = 1:length(net$genes)
  for (i in res) { 
    hamming <- sapply(1:repetitions, function(n) {
      state <- sample.int(2^length(net$genes),1) #random initial state
      state <- validateState(state, net$genes)
      genes <- sample(net$genes,i) #random nodes to perturb
      ori <- stateTransition(net,state)
      per <- perturbState(net, state, genes, result = "nextState")
      sum(ori!=per) #hamming between two states
    })
    res[[i]] <- mean(hamming)
  }
  res <- c(0, res)
  names(res) <- 0:length(net$genes)
  res
}
macrophage.derrida<-derrida(macrophage)
macrophage.derrida<- read.csv("PathWhereYouDownloadAllTheCode/Macrophage.derrida.csv", comment.char="#") 
library("extrafont")
loadfonts()
pdf("Fig8A.pdf",family="Arial", width=8.83, height=8)
op <- par(mar = c(8, 7, 0.5, 0.5) + 0.1)
plot(macrophage.derrida$x, macrophage.derrida$y,
   ylim=c(0, length(macrophage.derrida$y)), 
     pch=18, 
     cex=2, 
     col="#69b3a2",
     xlab="h_t", ylab="h_t+1")
abline(0,1)
par(op)
dev.off()
#################################################################
#                                                               #
#                Sensitiviy Analysis                            #
#################################################################
# A function to analize the sensitivity values of the nodes or variables in my mathematical model. 
SensitivityMacrophage <- function(macrophage, nSamples = 10000){
  nodes <- macrophage$nodes
  SensitivitiesValues <- c()
  for(n in nodes){
    sensi <- perturbTrajectories(macrophage, measure = "sensitivity", numSamples = nSamples, flipBits = 1, gene = n)
    sensitivities <- append(SensitivitiesValue,sensi$value)}
  names(SensitivitiesValues) <- nodes
  return(SensitivitiesValues)
} 
Sensibilidad<-SensitivityMacrophage(macrophage, nSamples = 100000)  
##Sensitivity Analysis for Our Macrophage Polarization in a Tumor Microenvironment 
Sensibilidad<-read.csv("PathWhereYouDownloadAllTheCode/AnalisisDeSensibilidad.csv", header=T)
library(ggplot2)
library("extrafont")
loadfonts()
pdf("Fig8B.pdf",family="Arial", width=6.83, height=6) 
ggplot(Sensibilidad, aes(x=Node, y=Sensitivity)) +
  geom_segment( aes(x=Node, xend=Node, y=0, yend=Sensitivity)) +
  geom_point( size=5, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=2)+
  coord_flip() + geom_hline(yintercept = 0.028, linetype="solid", 
                            color = "green", size=0.5)
dev.off()                        
                                                      

