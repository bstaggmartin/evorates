#fit an autocorrelated Brownian motion model to a tree and trait data
#' @export
fit.corateBM<-function(tree,X,niter=2000,report=niter/10,nchain=1,mcmc.pars=list(adapt_delta=0.8,max_treedepth=10)){
  X<-X[tree$tip.label]
  n<-length(tree$tip.label)
  n_e<-nrow(tree$edge)
  eV<-edge.vcv(tree)
  o.hgt<-max(eV)
  tree$edge.length<-tree$edge.length/o.hgt
  eV<-eV/o.hgt
  o.Xsig2<-sum(pic(X,multi2di(tree))^2)/n
  xx<-c(X,tree$Node)
  for(e in nrow(tree$edge):0){
    if(e==0){
      des<-which(tree$edge[,1]==n+1)
      xx[n+1]<-sum((1/tree$edge.length[des])/sum(1/tree$edge.length[des])*xx[tree$edge[des,2]])
      break
    }
    nn<-tree$edge[e,2]
    if(nn<=n){
      next
    }
    des<-which(tree$edge[,1]==nn)
    xx[nn]<-sum((1/tree$edge.length[des])/sum(1/tree$edge.length[des])*xx[tree$edge[des,2]])
  }
  o.X0<-xx[n+1]
  X<-(X-o.X0)/sqrt(o.Xsig2)
  
  if(!is.binary(tree)){
    poly.nodes<-which(sapply(1:max(tree$edge),function(ii) length(which(ii==tree$edge[,1])))>2)
    d_poly<-lapply(poly.nodes,function(ii) which(ii==tree$edge[,1]))
    tmp<-sort(unlist(lapply(d_poly,function(ii) ii[seq(2,length(ii)-1)])))
    tmp<-tmp+0:(length(tmp)-1)
    real_e<-(1:(2*n-2))[-tmp]+1
    tree<-multi2di(tree,random=F)
  }else{
    real_e<-1:n_e+1
  }
  des_e<-sapply(1:nrow(tree$edge),function(ii) which(tree$edge[,1]==tree$edge[ii,2])+1)
  tip_e<-which(lengths(des_e)==0)
  tip_e<-tip_e[order(tree$edge[tip_e,2])]
  des_e[tip_e]<-list(rep(-1,2))
  des_e<-matrix(unlist(des_e),ncol=2,byrow=T)
  root.edges<-which(tree$edge[,1]==n+1)+1
  des_e<-rbind(root.edges,des_e)
  prune_T<-c(0,tree$edge.length)
  prune_seq<-((2*n-2):1)[-((2*n-2)-tip_e)]
  tip_e<-tip_e+1
  
  dat<-list('n'=n,'n_e'=n_e,'X'=X,'eV'=eV,'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,'prune_seq'=prune_seq)
  
  # if(optimize){
  #   ret<-rstan::optimizing(object=stanmodels$univar_corateBM,data=dat)
  #   
  #   out<-ret$par[-((1:n_e)+2)]
  #   
  #   out[1]<-out[1]+log(o.Xsig2)-log(o.hgt)
  #   out[2]<-out[2]/o.hgt
  #   out[3]<-out[3]*sqrt(o.Xsig2)+o.X0
  #   out[-(1:3)]<-out[-(1:3)]+log(o.Xsig2)-log(o.hgt)
  #   
  #   return(out)
  # }else{
  ret<-rstan:::sampling(object=stanmodels$univar_corateBM,data=dat,iter=niter,chains=nchain,refresh=report,control=mcmc.pars)
  
  R0<-rstan::extract(ret,"R0",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
  R<-rstan::extract(ret,"R",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
  Rsig2<-rstan::extract(ret,"Rsig2",permute=FALSE,inc_warmup=FALSE)/o.hgt
  X0<-rstan::extract(ret,"X0",permute=FALSE,inc_warmup=FALSE)*sqrt(o.Xsig2)+o.X0
  
  return(list('R'=R,'R0'=R0,'Rsig2'=Rsig2,'X0'=X0))
  # }
}