Proportion <- function(clust,n,col,n.Pheno,Pheno,tot){
  nam= Pheno
  prop <- vector(mode = "double", length = n*n.Pheno)
  prop1 <- vector(mode = "double", length = n*n.Pheno)
  id.fl <- 1
  id1 <- 0
  #### id <- Cluster
  #### id1 <- Phenotype
  for(i in 0:(n*n.Pheno-1)){ #n*n.Pheno
    id <- i%/%n.Pheno+1
    if (id.fl == id){
      id1 = id1+1
    }else{ 
      id1=1
      id.fl= id}
#    print(c(id, id.fl,id1))
    b=which(clust %in% id)
    prop[i+1]<-length(grep(paste("^",nam[id1],"$",sep = ""),
                           col[b],value = FALSE))/tot[id1]*100
    if (id1 == n.Pheno){
      prop1[(i-n.Pheno+2):(i+1)]=100*prop[(i-n.Pheno+2):(i+1)]/sum(prop[(i-n.Pheno+2):(i+1)])
    } 
  }
  prop <-array(prop,dim=c(n.Pheno,n))
  prop1 <-array(prop1,dim=c(n.Pheno,n))
  prop <- as.table(prop)
  prop1 <- as.table(prop1)
  rownames(prop) <- Pheno
  #  colnames(prop) <- c(toString(1:n))
  
  rownames(prop1) <- Pheno
  #  colnames(prop1) <- c(toString(1:n))
  return(list(prop,prop1))
}

MyplotBar <- function(prop,tit,comp,lg,col){
  png(file=paste('~/Documents/UGO/Proportions_',tit,'_',
                 comp,'.png',
                 sep = ""),width = 2500, height = 700,
      units = "px", pointsize = 24)
  par(mar= c(5.1, 2.0, 4.2,4.5),xpd=TRUE)
#      oma=c(0, 2, 0, 0))
#  barplot(prop, ylim = c(0,100), col=1:col)#,
#          legend = rownames(prop[[1]]),main = 'Sample percentage')
  barplot(prop, col=col,
          main = comp, cex.names = 0.9)
  legend("bottomleft",inset=c(0.97,0), title = lg,
         legend = rownames(prop),
         pch = rep(21,dim(prop)[1]),
         pt.bg = col,
         bty = "n",
         pt.cex = 1,
         cex = 0.8)
  dev.off()
}
