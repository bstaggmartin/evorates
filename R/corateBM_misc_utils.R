#' @export
prep.tree<-function(tree,tol=0){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladwise')
  tree<-di2multi(tree,tol=tol)
}

#' @export
coerce.to.cov.mat<-function(in.mat,vars){
  dims<-dim(in.mat)
  if(length(dims)>0){
    new.diag<-rep(diag(in.mat),length.out=vars)
    mat<-do.call(rbind,c(list(in.mat),rep(0,vars-nrow(in.mat))))
    if(ncol(mat)>vars){
      mat<-mat[,1:vars]
    }else{
      mat<-do.call(cbind,c(list(mat),rep(0,vars-ncol(mat))))
    }
    diag(mat)<-new.diag
    if(!isSymmetric(mat)){
      warning(paste(deparse(substitute(in.mat)),'is not symmetric: reflected lower triangle into upper triangle'))
      mat[upper.tri(mat)]<-t(mat)[upper.tri(mat)]
    }
    if(any(eigen(mat)$values<=0)){
      stop(paste('failed to create properly-formed covariance matrix from ',deparse(substitute(in.mat)),': make sure variance-covariance structure makes sense',sep=''))
    }
    mat
  }else{
    mat<-matrix(0,vars,vars)
    diag(mat)<-rep(in.mat,length.out=vars)
    mat
  }
}

#' @export
quick.recon<-function(X,tree){
  if(!is.null(names(X))){
    X<-X[tree$tip.label]
  }
  XX<-rep(NA,nrow(tree$edge)+1)
  names(XX)<-c(tree$tip.label,length(tree$tip.label)+1:tree$Nnode)
  PP<-rep(NA,nrow(tree$edge)+1)
  for(e in nrow(tree$edge):0){
    if(e==0){
      n<-length(tree$tip.label)+1
      d<-tree$edge[which(tree$edge[,1]==n),2]
      sum.P<-sum(PP[d])
      XX[n]<-sum(XX[d]*PP[d]/sum.P)
      PP[n]<-sum(sum.P)
      break
    }
    n<-tree$edge[e,2]
    if(length(which(tree$edge[,1]==n))==0){
      XX[n]<-X[n]
      PP[n]<-1/tree$edge.length[e]
    }
    if(length(which(tree$edge[,1]==n))>0){
      d<-tree$edge[which(tree$edge[,1]==n),2]
      sum.P<-sum(PP[d])
      XX[n]<-sum(XX[d]*PP[d]/sum.P)
      PP[n]<-sum.P/(1+tree$edge.length[e]*sum.P)
    }
  }
  for(e in 1:nrow(tree$edge)){
    n<-tree$edge[e,2]
    if(length(which(tree$edge[,1]==n))!=0){
      a<-tree$edge[e,1]
      XX[n]<-XX[n]*PP[n]*tree$edge.length[e]+XX[a]-XX[a]*PP[n]*tree$edge.length[e]
      #no need for variance if only goal is ancestral states
    }
  }
  XX[-(1:length(tree$tip.label))]
}

#helpful in fit.corateBM function -- excludes rows from a 3D array without any simplification or destroying
#names
#' @export
reduce.array<-function(arr,excl.inds){
  new.dims<-dim(arr)
  new.dims[1]<-new.dims[1]-length(excl.inds)
  array(arr[-excl.inds,,],new.dims,dimnames(arr))
}
