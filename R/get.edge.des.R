#' @export
get.edge.des<-function(node,tree){
  if(length(node)>1){
    stop('get.edge.des only does one node at a time!')
  }
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  out<-new.edges<-which(tree$edge[,1]==node)
  if(length(out)==0){
    stop('given node is a tip: it has no descendants!')
  }
  while(length(new.edges)>0){
    new.edges<-unlist(lapply(tree$edge[new.edges,2],function(ii) which(tree$edge[,1]==ii)))
    out<-c(out,new.edges)
  }
  out
}

tree<-pbtree(n=50)
sim<-gen.corateBM(tree,Rsig2=0,Rmu=-0.5)
fit<-fit.corateBM(tree,sim$X,constrain.Rsig2=T,chains=1,trend=T,iter=2000)


