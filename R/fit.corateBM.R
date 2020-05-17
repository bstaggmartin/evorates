#fit an autocorrelated Brownian motion model to a tree and trait data


#' Fit an autocorrelated rate Brownian motion model
#'
#' This 
#'
#' @param tree A phylogenetic tree of class 'phylo'; function does not yet support multiple phylogenies
#' but will handle polytomous trees just fine.
#' @param X A vector of trait data, with each element labelled according to its corresponding tip label.
#' @param R0.prior The standard deviation of the cauchy prior for the \code{R0} parameter
#' @param Rsig2.prior The standard deviation of the cauchy prior for the \code{Rsig2} parameter.
#' @param X0.prior The standard deviation of the normal prior for the \code{X0} parameter.
#' @param report.quantiles A vector of values between 0 and 1, telling the function which quantiles of
#' of the marginal posterior distributions for each parameter to calculate. Use this to extract
#' 'credible intervals' for parameter values. Set to NULL if no quantiles are desired.
#' @param report.MAPs TRUE or FALSE value: should the function calculate mean a posteriori parameter
#' estimates?
#' @param report.devs TRUE or FALSE value: should the function calculate differences between branch-wise
#' rates and the 'background rate' (see below)?
#' @param ... Extra arguments to pass to the Stan's mcmc sampler. Use \code{chains} to set the number
#' of chains to run (4 by default), \code{iter} for the number of iterations (2000 by default),
#' \code{warmup} for the number of iterations devoted to the warmup period (\code{iter/2} by default),
#' \code{thin} to set the thinning rate (1 by default), and \code{refresh} to set the frequency of
#' reporting in the R console (\code{iter/10} by default). See \code{rstan::sampling} for more options.
#' @return 
#' @export

#5/14: add trace and profile plotting functions
fit.corateBM<-function(tree,X,R0.prior=10,Rsig2.prior=20,X0.prior=100,
                       report.quantiles=c(0.025,0.5,0.975),report.MAPs=T,report.devs=T,...,
                       optimize=F){
  if(hasArg(chains)){
    nchain<-list(...)$chains
  }else{
    nchain<-4
  }
  X<-X[tree$tip.label]
  n<-length(tree$tip.label)
  n_e<-nrow(tree$edge)
  eV<-edge.vcv(tree)
  o.hgt<-max(eV)
  tree$edge.length<-tree$edge.length/o.hgt
  eV<-eV/o.hgt
  o.Xsig2<-sum(ape::pic(X,ape::multi2di(tree))^2)/n
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
  dat<-list('n'=n,'n_e'=n_e,'X'=X,'eV'=eV,
            'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,'prune_seq'=prune_seq,
            'R0_prior'=R0.prior,'Rsig2_prior'=Rsig2.prior,'X0_prior'=X0.prior)
  if(optimize){
    ret<-rstan::optimizing(object=stanmodels$univar_corateBM,data=dat,...)
    
    out<-ret$par[-((1:n_e)+2)]
    
    out[1]<-out[1]+log(o.Xsig2)-log(o.hgt)
    out[2]<-out[2]/o.hgt
    out[3]<-out[3]*sqrt(o.Xsig2)+o.X0
    out[-(1:3)]<-out[-(1:3)]+log(o.Xsig2)-log(o.hgt)
    
    return(out)
  }
  ret<-rstan::sampling(object=stanmodels$univar_corateBM,data=dat,...)
  R0<-rstan::extract(ret,"R0",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
  R<-rstan::extract(ret,"R",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
  Rsig2<-rstan::extract(ret,"Rsig2",permute=FALSE,inc_warmup=FALSE)/o.hgt
  X0<-rstan::extract(ret,"X0",permute=FALSE,inc_warmup=FALSE)*sqrt(o.Xsig2)+o.X0
  wgts<-tree$edge.length/sum(tree$edge.length)
  if(nchain==1){
    bg.rate=apply(R,1,function(ii) sum(ii*wgts))
  }else{
    bg.rate=apply(R,c(1,2),function(ii) sum(ii*wgts))
  }
  out<-list(chains=array(dim=c(dim(R)[1],4+dim(R)[3],nchain)))
  dimnames(out$chains)<-list(iterations=NULL,
                             parameters=c('R0','Rsig2','X0','bg.rate',
                                          paste('R[',1:dim(R)[3],']',sep='')),
                             chains=paste('chain',1:nchain))
  out$chains[,'R0',]<-R0
  out$chains[,'Rsig2',]<-Rsig2
  out$chains[,'X0',]<-X0
  out$chains[,'bg.rate',]<-bg.rate
  out$chains[,5:dim(out$chains)[2],]<-aperm(R,c(1,3,2))
  class(out)<-'corateBM_fit'
  if(report.devs){
    if(nchain==1){
      rate.devs<-apply(R[,1,],2,function(ii) ii-bg.rate)
    }else{
      rate.devs<-R
      for(i in 1:dim(R)[2]){
        rate.devs[,i,]<-apply(R[,i,],2,function(ii) ii-bg.rate[,i])
      }
    }
    tmp<-dimnames(out$chains)
    tmp[[2]]<-c(tmp[[2]],paste('R[',1:dim(R)[3],'] dev',sep=''))
    out$chains<-aperm(array(c(aperm(out$chains,c(1,3,2)),rate.devs),
                            dim=c(dim(R)[1],nchain,4+dim(R)[3]*2),dimnames=tmp[c(1,3,2)]),
                      c(1,3,2))
  }
  if(!is.null(report.quantiles)){
    out$quantiles<-apply(out$chains,c(2,3),quantile,probs=report.quantiles)
  }
  if(report.MAPs){
    out$MAPs<-apply(out$chains,c(2,3),mean)
  }
  if(report.devs){
    if(nchain==1){
      out$post.probs<-apply(out%chains%'dev',2,function(ii) sum(ii>0)/length(ii))
    }else{
      out$post.probs<-apply(out%chains%'dev',c(2,3),function(ii) sum(ii>0)/length(ii))
    }
  }
  out
}