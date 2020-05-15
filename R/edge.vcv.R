#get the edge-wise expected variance covariance evolving under BM with branch lengths equal to evolutionary rate
#' @export
edge.vcv<-function(tree){
  ancs<-vector(mode='list',length=nrow(tree$edge))
  for(i in 1:nrow(tree$edge)){
    adj.edge<-which(tree$edge[,2]==tree$edge[i,1])
    if(length(adj.edge)>0){
      ancs[[i]]<-c(ancs[[adj.edge]],i)
    }else{
      ancs[[i]]<-i
    }
  }
  for(i in 1:length(ancs)){
    ancs[[i]]<-ancs[[i]][-which(ancs[[i]]==i)]
  }
  edge.hgts<-ape::node.depth.edgelength(tree)[tree$edge[,2]]
  mrcas<-ape::mrca(tree,full=T)
  mat<-matrix(NA,nrow=nrow(tree$edge),ncol=nrow(tree$edge))
  mat[1,1]<-edge.hgts[1]-2*tree$edge.length[1]/3
  for(i in 2:nrow(tree$edge)){
    tmp.edge<-sapply(mrcas[tree$edge[i,2],tree$edge[1:(i-1),2]],function(ii) which(tree$edge[,2]==ii))
    tmp.edge[lengths(tmp.edge)==0]<-NA;tmp.edge<-unlist(tmp.edge)
    anc.des<-sapply(1:(i-1),function(ii) if(ii%in%ancs[[i]]) 1 else 0)
    mat[i,1:i]<-c(ifelse(is.na(tmp.edge),0,edge.hgts[tmp.edge]-anc.des*tree$edge.length[tmp.edge]/2),
                  edge.hgts[i]-2*tree$edge.length[i]/3)
  }
  mat[upper.tri(mat)]<-t(mat)[upper.tri(mat)]
  mat
}