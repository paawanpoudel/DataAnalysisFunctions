hypergeometric_test <-
function(b){
  
  # perform the hypergeometric test
  pv<-matrix(NA,nrow(b),ncol(b))
  for(i in 1:ncol(b)) {
    for(j in 1:nrow(b)) {
      pv[j,i]<-phyper((b[j,i]-1),sum(b[,i]),sum(b)-sum(b[,i]),sum(b[j,]),lower.tail=FALSE)
    }
  }
  
  
  apval<-matrix(p.adjust(c(pv),"fdr"),nrow(pv),ncol(pv))
  rownames(apval)<-rownames(pv)<-rownames(b)
  colnames(apval)<-colnames(pv)<-colnames(b)
  
  return(apval)
}
