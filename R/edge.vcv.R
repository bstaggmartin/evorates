#get the edge-wise expected variance covariance evolving under BM with branch lengths equal to evolutionary rate
#' @export
old.edge.vcv<-function(tree){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  ancs<-vector(mode='list',length=nrow(tree$edge))
  for(i in 1:nrow(tree$edge)){
    adj.edge<-which(tree$edge[,2]==tree$edge[i,1])
    if(length(adj.edge)>0){
      ancs[[i]]<-c(ancs[[adj.edge]],i)
    }else{
      ancs[[i]]<-i
    }
  }
  edge.hgts<-node.depth.edgelength(tree)[tree$edge[,2]]
  mat<-matrix(NA,nrow=nrow(tree$edge),ncol=nrow(tree$edge))
  get.mrca<-function(other.edge,ref.edge,ancs){
    ind<-which(ancs[[ref.edge]]%in%ancs[[other.edge]])
    if(length(ind)==0){
      NA
    }else{
      ancs[[ref.edge]][max(ind)]
    }
  }
  for(i in 1:nrow(tree$edge)){
    tmp.edge<-unlist(lapply((1:i)[-i],get.mrca,ref.edge=i,ancs=ancs))
    if(all(is.na(tmp.edge))){
      tmp.hgts<-rep(0,length(tmp.edge))
    }else{
      anc.des<-as.numeric(tmp.edge==1:length(tmp.edge))
      tmp.hgts<-edge.hgts[tmp.edge]-anc.des*tree$edge.length[tmp.edge]/2
      tmp.hgts[is.na(tmp.hgts)]<-0
    }
    mat[i,1:i]<-c(tmp.hgts,edge.hgts[i]-2*tree$edge.length[i]/3)
  }
  mat[upper.tri(mat)]<-t(mat)[upper.tri(mat)]
  mat
}
