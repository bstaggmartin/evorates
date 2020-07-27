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

#newest version of check.n.proc and coerce.to.array--expands any given element to a 3D array (if chains,
#quantiles, sampler, or parameter diagnostics) or a 2D matrix (if means or MAPs)
#finally got it to a point where it should work with just about anything...
#' @export
expand.element<-function(arr){
  if(length(dim(arr))==0){
    new.arr<-as.matrix(arr)
    new.dimnames<-dimnames(new.arr)
    if(is.null(dimnames(new.arr)[[1]])){
      names(new.dimnames)[1]<-'iterations'
    }else if(sum(grepl('%',dimnames(new.arr)[[1]]))!=0){
      names(new.dimnames)[1]<-'quantiles'
    }else if(sum(grepl('_',dimnames(new.arr)[[1]]))!=0){
      if(all(dimnames(new.arr)[[1]]==c('inits','bulk_ess','tail_ess','Rhat'))){
        names(new.dimnames)<-''
      }else{
        names(new.dimnames)[1]<-'parameters'
      }
    }else if(sum(grepl('chain',dimnames(new.arr)[[1]]))!=0){
      names(new.dimnames)[1]<-'chains'
    }
    dimnames(new.arr)<-new.dimnames
    for(i in c('quantiles','parameters','chains')){
      if(!is.null(attr(arr,i))){
        attr(new.arr,i)<-attr(arr,i)
      }
    }
    arr<-new.arr
  }
  if(length(dim(arr))==2){
    if('quantiles'%in%names(dimnames(arr))|'quantiles'%in%names(attributes(arr))){
      first.dim<-'quantiles'
    }else if(''%in%names(dimnames(arr))){
      first.dim<-''
    }else{
      first.dim<-'iterations'
    }
    dim.map<-lapply(c(first.dim,'parameters','chains'),
                    function(ii) which(names(dimnames(arr))==ii))
    out.dim<-rep(NA,3)
    out.dimnames<-setNames(vector('list',3),c(first.dim,'parameters','chains'))
    for(i in 1:3){
      if(length(dim.map[[i]])==0){
        out.dim[i]<-1
        out.dimnames[i]<-list(attr(arr,names(out.dimnames)[i]))
      }else{
        out.dim[i]<-dim(arr)[dim.map[[i]]]
        out.dimnames[i]<-list(dimnames(arr)[[dim.map[[i]]]])
      }
    }
    dim.map[lengths(dim.map)==0]<-NA
    perm<-unique(unlist(dim.map))
    if(length(perm)==2){
      perm[is.na(perm)]<-if(perm[!is.na(perm)]==1) 2 else 1
    }else{
      perm<-perm[!is.na(perm)]
    }
    arr<-array(aperm(arr,perm),out.dim,out.dimnames)
  }
  # if(dim(arr)[1]==1&names(dimnames(arr)[1])!='quantiles'){
  #   arr<-array(arr,dim(arr)[-1],dimnames(arr)[-1])
  # }
  arr
}

#opposite of expand.element--collapses dimensions of length 1 and adds them to parameters. In cases of
#length 1 output, prioritizes parameter names over everything else
##NEED TO UPDATE TO WORK WITH 2D ELEMENTS (like means and MAPs)
#' @export
simplify.element<-function(arr){
  o.dimnames<-dimnames(arr)
  o.dims<-dim(arr)
  out<-arr[1:o.dims[1],1:o.dims[2],1:o.dims[3]]
  for(i in 1:3){
    if(o.dims[i]==1){
      if(!is.null(o.dimnames[[i]])){
        attr(out,names(o.dimnames)[i])<-o.dimnames[[i]]
      }
    }
  }
  if(length(out)==1){
    names(out)<-attr(out,names(o.dimnames)[2])
    attr(out,names(o.dimnames)[2])<-NULL
  }
  out
}

#newest, generalized version of reduce.array--selects indices from arbitrary arrays without
#simplification or destroying names
#' @export
index.array<-function(arr,inds,dims,invert=F){
  inds.list<-vector('list',length(dim(arr)))
  if(is.numeric(inds)){
    inds<-list(inds)
  }
  inds.list[dims]<-inds
  if(invert){
    inds.list<-lapply(inds.list,function(ii) -1*ii)
  }else{
    inds.list<-lapply(1:length(inds.list),
                      function(ii) if(is.null(inds.list[[ii]])) 1:dim(arr)[ii] else inds.list[[ii]])
  }
  new.dims<-dim(arr)
  if(invert){
    new.dims<-new.dims-lengths(inds.list)
    inds.list[lengths(inds.list)==0]<- -dim(arr)[lengths(inds.list)==0]-1
  }else{
    new.dims<-lengths(inds.list)
  }
  inds.list<-lapply(inds.list,sort)
  new.dimnames<-dimnames(arr)
  new.dimnames<-lapply(1:length(inds.list),
                       function(ii) new.dimnames[[ii]][inds.list[[ii]]])
  names(new.dimnames)<-names(dimnames(arr))
  array(arr[inds.list[[1]],inds.list[[2]],inds.list[[3]]],new.dims,new.dimnames)
}



#helpful in fit.corateBM function -- excludes rows from a 3D array without any simplification or destroying
#names
#' @export
reduce.array<-function(arr,excl.inds){
  new.dims<-dim(arr)
  new.dims[1]<-new.dims[1]-length(excl.inds)
  array(arr[-excl.inds,,],new.dims,dimnames(arr))
}