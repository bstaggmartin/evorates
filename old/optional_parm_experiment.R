#' @export

#5/14: add trace and profile plotting functions
fit.corateBM.exp<-function(tree,X,R0.prior=10,Rsig2.prior=20,X0.prior=100,
                       intra.var=F,X.prior=200,Xsig2.prior=50,
                       trend=F,Rmu.prior=Rsig2.prior,
                       report.quantiles=c(0.025,0.5,0.975),report.means=T,report.devs=T,
                       constrain.Rsig2=F,constrain.Rmu=F,...){
  
  
  if(hasArg(chains)){
    nchain<-list(...)$chains
  }else{
    nchain<-4
  }
  
  
  if(is.null(names(X))){
    stop('trait data (X) is unlabelled--please name each element with its corresponding tip label')
  }
  X<-X[unlist(sapply(tree$tip.label,function(ii) which(names(X)==ii)))]
  n_obs<-sapply(tree$tip.label,function(ii) sum(names(X)==ii))
  n<-length(tree$tip.label)
  n_e<-nrow(tree$edge)
  eV<-edge.vcv(tree)
  o.hgt<-max(eV)
  tree$edge.length<-tree$edge.length/o.hgt
  eV<-eV/o.hgt
  o.Xsig2<-sum(ape::pic(tapply(X,names(X),mean)[tree$tip.label],ape::multi2di(tree))^2)/n
  xx<-c(tapply(X,names(X),mean)[tree$tip.label],tree$Nnode)
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
  
  
  dat<-list('n'=n,'n_e'=n_e,'eV'=eV,
            'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,'prune_seq'=prune_seq,
            'R0_prior'=R0.prior,'Rsig2_prior'=Rsig2.prior,'X0_prior'=X0.prior,'Rmu_prior'=Rmu.prior,
            'constr_Rsig2'=as.numeric(constrain.Rsig2),'constr_Rmu'=as.numeric(constrain.Rmu))
  
  ret<-rstan::sampling(object=stanmodels$optional_parm_experiment,data=dat,...)
  R0<-rstan::extract(ret,"R0",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
  X0<-rstan::extract(ret,"X0",permute=FALSE,inc_warmup=FALSE)*sqrt(o.Xsig2)+o.X0
  out<-list('R0'=R0,'X0'=X0)
  if(!constrain.Rmu){
    Rmu<-rstan::extract(ret,"Rmu",permute=FALSE,inc_warmup=FALSE)/o.hgt
    out$Rmu<-Rmu
  }
  if(!constrain.Rsig2){
    Rsig2<-rstan::extract(ret,"Rsig2",permute=FALSE,inc_warmup=FALSE)/o.hgt
    R<-rstan::extract(ret,"R",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
    out$Rsig2<-Rsig2
    out$R<-R
  }
  out
}