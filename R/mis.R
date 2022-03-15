.coerce.to.cov.mat<-function(in.mat,vars){
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

.anc.recon<-function(tree,XX,LL,PP,stochastic){
  ntips<-length(tree$tip.label)
  nedges<-nrow(tree$edge)
  for(e in nedges:0){
    if(e==0){
      n<-ntips+1
    }else{
      n<-tree$edge[e,2]
    }
    d<-tree$edge[which(tree$edge[,1]==n),2]
    if(length(d)==0){
      if(is.infinite(PP[1,n,1])){
        PP[,n,]<-1/LL[,n,]
      }else{
        PP[,n,]<-PP[,n,]/(1+LL[,n,]*PP[,n,])
      }
    }else{
      des_P<-PP[,d,,drop=F]
      des_X<-XX[,d,,drop=F]
      infs<-is.infinite(des_P[1,,1])
      inf.code<-sum(infs)
      if(inf.code>1){
        stop('Two tips are occupying the exact same phylogenetic position: something went very wrong here.')
      }else if(inf.code==1){
        XX[,n,]<-des_X[,infs,]
        PP[,n,]<-1/LL[,n,]
      }else{
        sum_P<-apply(des_P,c(1,3),sum)
        #this could be improved...assumes the same thing the infinite does earlier--that iterations won't have mixes of inf and
        #non-infs at the same node
        if(sum_P[1,1]==0){
          XX[,n,]<-0
        }else{
          XX[,n,]<-apply(XX[,d,,drop=F]*PP[,d,,drop=F],c(1,3),sum)/sum_P
        }
        PP[,n,]<-sum_P/(1+LL[,n,]*sum_P)
      }
    }
  }
  if(stochastic){
    niter<-dim(XX)[1]
    nchains<-dim(XX)[3]
    XX[,ntips+1,]<-rnorm(niter*nchains,XX[,ntips+1,],sqrt(1/PP[,ntips+1,]))
    sds<-PP
    sds[,,]<-1/sqrt(PP/(1-LL*PP)+1/LL)
    for(e in 1:nedges){
      n<-tree$edge[e,2]
      a<-tree$edge[e,1]
      if(is.infinite(PP[1,n,1])){
        XX[,n,]<-XX[,n,]
      }else{
        PL<-PP[,n,]*LL[,n,]
        XX[,n,]<-(XX[,n,]-XX[,a,])*PL+XX[,a,]
        XX[,n,]<-rnorm(niter*nchains,XX[,n,],sds[,n,])
      }
    }
    XX
  }else{
    for(e in 1:nedges){
      n<-tree$edge[e,2]
      a<-tree$edge[e,1]
      if(is.infinite(PP[1,n,1])){
        XX[,n,]<-XX[,n,]
      }else{
        PL<-round(PP[,n,]*LL[,n,],15)
        if(PL[1]==1){
          PP[,n,]<-Inf
        }else{
          P_diff<-PP[,a,]-PP[,n,]
          XX[,n,]<-(XX[,n,]-XX[,a,])*PL+XX[,a,]
          if(is.infinite(P_diff[1])){
            PP[,n,]<-PP[,n,]/(1-PL)+1/LL[,n,]
          }else{
            PP[,n,]<-PP[,n,]/(1-PL)+P_diff/(1+LL[,n,]*P_diff)
          }
        }
      }
    }
    list(XX,1/PP)
  }
}

.lin.interp<-function(x,length.out){
  xx<-seq(1,length(x),length.out=length.out)
  in.inds<-xx%in%(1:length(xx))
  out<-rep(NA,length.out)
  out[in.inds]<-x[xx[in.inds]]
  xx<-xx[!in.inds]
  out[!in.inds]<-(x[ceiling(xx)]-x[floor(xx)])/(ceiling(xx)-floor(xx))*(xx-floor(xx))+x[floor(xx)]
  out
}