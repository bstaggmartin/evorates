#' Fit an autocorrelated rates Brownian motion model to a tree and trait data
#'
#' This function processes tree and trait data and fits an autocorrelated rates Brownian motion model to
#' it using a Stan-based mcmc sampler. Under an autocorrelated rates Brownian motion model, the rates at
#' trait evolution themselves evolve according to a Brownian motion process, allowing rates of
#' evolution to vary across the tree in a phylogenetically autocorrelated manner. The results are
#' formatted into a (relatively) user-friendly list of arrays for further analysis and diagnostics.
#'
#' @param tree A phylogenetic tree of class 'phylo'; function does not yet support multiple phylogenies
#' but will handle polytomous trees just fine.
#' @param X A vector of trait data, with each element labelled according to its corresponding tip label.
#' @param R0.prior The standard deviation of the cauchy prior for the \code{R0} parameter determining
#' the rate of trait evolution at the root.
#' @param Rsig2.prior The standard deviation of the cauchy prior for the \code{Rsig2} parameter
#' determining the rate at which rate variance accumulates over time (inversely related to phylogenetic
#' autocorrelation of rates).
#' @param X0.prior The standard deviation of the cauchy prior for the \code{X0} parameter determining
#' the trait value at the root.
#' @param intra.var TRUE or FALSE value: should tip means be estimated with an error term or fixed
#' according to the data?
#' @param X.prior The standard deviation of the cauchy prior for the \code{X} parameter determining the
#' mean trait values at the tips. Only matters when \code{intra.var} is set to TRUE.
#' @param Xsig2.prior The standard deviation of the cauchy prior for the \code{Xsig2} parameter
#' determining the variance of observed trait values due to both intraspecific variation and/or
#' measurement error. Only matters when \code{intra.var} is set to TRUE.
#' @param trend TRUE or FALSE value: should the model fit an overall temporal trend to rates, as in an
#' early burst or accelarating rates model?
#' @param Rmu.prior The standard deviation of the cauchy prior for the \code{Rmu} parameter determining
#' the rate at which rates of trait evolution increase or decrease over time. Only matters when
#' \code{trend} is set to TRUE.
#' @param constrain.Rsig2 TRUE or FALSE value: should the model constrain Rsig2 to be 0, reducing the
#' model to either a plain single-rate Brownian motion model (if \code{trend} is set to FALSE) or a
#' early burst or accelarating rates model (if \code{trend} is set to TRUE)?
#' @param report.quantiles A vector of values between 0 and 1, telling the function which quantiles of
#' of the marginal posterior distributions for each parameter to calculate. Use this to extract
#' 'credible intervals' for parameter values. Set to NULL if no quantiles are desired.
#' @param report.means TRUE or FALSE value: should the function calculate mean posterior parameter
#' estimates?
#' @param report.devs TRUE or FALSE value: should the function calculate differences between average
#' branch-wise rates and the 'background rate' (see below)?
#' @param ... Extra arguments to pass to the Stan's mcmc sampler. Use \code{chains} to set the number
#' of chains to run (4 by default), \code{iter} for the number of iterations (2000 by default),
#' \code{warmup} for the number of iterations devoted to the warmup period (\code{iter/2} by default),
#' \code{thin} to set the thinning rate (1 by default), and \code{refresh} to set the frequency of
#' reporting in the R console (\code{iter/10} by default). See \code{rstan::sampling} for more options.
#' 
#' @return A list of class 'corateBM_fit', containing an array 'chains', with rows corresponding to
#' iterations in the mcmc sampler and columns to parameters. The array is 3D, with each
#' matrix slice corresponding to a particular mcmc chain. Parameters may include \code{R0}, \code{Rsig2},
#' \code{X0}, \code{Xsig2}, and \code{Rmu}, as defined above. Average branch-wise rates are also included
#' and denoted
#' \code{R[i]} where \code{i} is the index of the edge as given in \code{tree}. If \code{intra.var} is
#' set to TRUE, the array will also include the tip mean parameters, denoted by their tip labels as given
#' in \code{tree}. Additionally, the function calculates the posterior distribution of 'background rates'
#' by taking the average of all average branch-wise rate parameters, weighted by the relative lengths of 
#' their respective branches, denoted \code{bg.rate}.
#' 
#' Optionally, the object may contain additional arrays including 'quantiles', 'means', and 'post.probs'.
#' 'quantiles' is an array formatted similarly to 'chains' with rows corresponding to quantiles of the
#' marginal posterior distributions of parameters. 'means' is a matrix of mean posterior parameter
#' estimates with rows corresponding to parameters and columns to separate mcmc chains. If 
#' \code{report.devs} is set to TRUE, the function additionally calculates the posterior distributions of
#' differences between each average branch-wise rate parameter and the background rate ('rate
#' deviations', denoted \code{R[i] dev}), and creates a matrix, 'post.probs', of posterior probabilities
#' that each branch-wise rate deviation is greater than 0. These can be interpreted as the probability
#' that a given branch in the phylogeny exhibits an anomalously high rate of trait evolution (note that
#' low probabilities indicate anomalously low rates). Quantiles and mean posterior estimates will also be
#' calculated for rate deviation parameters, if applicable.
#' 
#' @export

#5/14: add trace and profile plotting functions
fit.corateBM<-function(tree,X,R0.prior=10,Rsig2.prior=20,X0.prior=100,
                       intra.var=F,X.prior=200,Xsig2.prior=50,
                       trend=F,Rmu.prior=Rsig2.prior,
                       constrain.Rsig2=F,
                       report.quantiles=c(0.025,0.5,0.975),report.means=T,report.devs=T,
                       return.stanfit=F,...){
  #initial checks
  if(is.null(names(X))){
    stop('trait data (X) is unlabelled--please name each element with its corresponding tip label')
  }
  if(constrain.Rsig2&!trend&report.devs){
    report.devs<-F
    warning('report.devs was set to FALSE: rate deviations are meaningless with no rate heterogeneity')
  }
  if(hasArg(chains)){
    nchain<-list(...)$chains
  }else{
    nchain<-4
  }
  
  #process/format data
  X<-X[unlist(sapply(tree$tip.label,function(ii) which(names(X)==ii)))]
  n_obs<-sapply(tree$tip.label,function(ii) sum(names(X)==ii))
  n<-length(tree$tip.label)
  n_e<-nrow(tree$edge)
  eV<-edge.vcv(tree)
  o.hgt<-max(eV)
  tree$edge.length<-tree$edge.length/o.hgt
  eV<-eV/o.hgt
  if(any(n_obs==0)){
    missing.dat<-tree$tip.label[n_obs==0]
    tmp.tree<-ape::drop.tip(tree,missing.dat)
  }else{
    tmp.tree<-tree
  }
  o.Xsig2<-sum(ape::pic(tapply(X,names(X),mean)[tmp.tree$tip.label],ape::multi2di(tmp.tree))^2)/n
  xx<-c(tapply(X,names(X),mean)[tmp.tree$tip.label],rep(NA,tmp.tree$Nnode))
  for(e in nrow(tmp.tree$edge):0){
    if(e==0){
      des<-which(tmp.tree$edge[,1]==n+1)
      xx[n+1]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
      break
    }
    nn<-tmp.tree$edge[e,2]
    if(nn<=n){
      next
    }
    des<-which(tmp.tree$edge[,1]==nn)
    xx[nn]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
  }
  o.X0<-xx[n+1]
  X<-(X-o.X0)/sqrt(o.Xsig2)
  if(!ape::is.binary(tree)){
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
  dat<-list('obs'=sum(n_obs),'n'=n,'n_e'=n_e,'X_obs'=X,'n_obs'=n_obs,'eV'=eV,
            'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,'prune_seq'=prune_seq,
            'R0_prior'=R0.prior,'Rsig2_prior'=Rsig2.prior,'X0_prior'=X0.prior,'Rmu_prior'=Rmu.prior,
            'constr_Rsig2'=as.numeric(constrain.Rsig2),'constr_Rmu'=as.numeric(!trend))
  
  #run mcmc
  if(intra.var){
    dat$X_prior=X.prior
    dat$Xsig2_prior=Xsig2.prior
    ret<-rstan::sampling(object=stanmodels$hier_intravar_univar_corateBM,data=dat,...)
  }else{
    if(any(n_obs>1)){
      warning('multiple observations per tip: observations were averaged, but we strongly recommend running function with intra.var set to TRUE')
      X<-tapply(X,names(X),mean)[tree$tip.label]
    }
    dat$X_obs<-X
    ret<-rstan::sampling(object=stanmodels$univar_corateBM,data=dat,...)
  }

  #process/format output
  R0<-rstan::extract(ret,"R0",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
  X0<-rstan::extract(ret,"X0",permute=FALSE,inc_warmup=FALSE)*sqrt(o.Xsig2)+o.X0
  out.X<-rstan::extract(ret,"X",permute=FALSE,inc_warmup=FALSE)*sqrt(o.Xsig2)+o.X0
  out<-list(chains=array(dim=c(dim(R0)[1],6+n_e+n,nchain)))
  param.names<-c('R0','X0',paste('R[',1:n_e,']',sep=''),'Rsig2','Rmu','bg.rate','Xsig2',tree$tip.label)
  dimnames(out$chains)<-list(iterations=NULL,
                             parameters=param.names,
                             chains=paste('chain',1:nchain))
  out$chains[,'R0',]<-R0
  out$chains[,'X0',]<-X0
  out$chains[,tree$tip.label,]<-aperm(out.X,c(1,3,2))
  if(intra.var){
    Xsig2<-rstan::extract(ret,"Xsig2",permute=FALSE,inc_warmup=FALSE)*o.Xsig2
    out$chains[,'Xsig2',]<-Xsig2
  }
  if(!constrain.Rsig2){
    Rsig2<-rstan::extract(ret,"Rsig2",permute=FALSE,inc_warmup=FALSE)/o.hgt
    out$chains[,'Rsig2',]<-Rsig2
  }
  if(!constrain.Rsig2|trend){
    R<-rstan::extract(ret,"R",permute=FALSE,inc_warmup=FALSE)+log(o.Xsig2)-log(o.hgt)
    wgts<-tree$edge.length/sum(tree$edge.length)
    if(nchain==1){
      bg.rate=apply(R,1,function(ii) sum(ii*wgts))
    }else{
      bg.rate=apply(R,c(1,2),function(ii) sum(ii*wgts))
    }
    out$chains[,'bg.rate',]<-bg.rate
    out$chains[,2+1:n_e,]<-aperm(R,c(1,3,2))
  }
  if(trend){
    Rmu<-rstan::extract(ret,"Rmu",permute=FALSE,inc_warmup=FALSE)/o.hgt
    out$chains[,'Rmu',]<-Rmu
  }
  incl.inds<-which(apply(out$chains[,,1],2,function(ii) !all(is.na(ii))))
  new.dim<-dim(out$chains);new.dim[2]<-length(incl.inds)
  new.dimnames<-dimnames(out$chains);new.dimnames[[2]]<-dimnames(out$chains)[[2]][incl.inds]
  out$chains<-array(out$chains[,incl.inds,],dim=new.dim,dimnames=new.dimnames)
  class(out)<-'corateBM_fit'

  #add rate deviation chains
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
                            dim=c(dim(out$chains)[1],nchain,dim(out$chains)[2]+n_e),
                            dimnames=tmp[c(1,3,2)]),
                      c(1,3,2))
  }

  #add quantiles
  if(!is.null(report.quantiles)){
    out$quantiles<-apply(out$chains,c(2,3),quantile,probs=report.quantiles)
    if(is.matrix(out$quantiles)){
      out$quantiles<-array(out$quantiles,dim=c(1,dim(out$quantiles)))
      dimnames(out$quantiles)<-c('quantiles'=paste(report.quantiles*100,'%',sep=''),dimnames(out$chains)[-1])
    }
  }

  #add mean posterior estimates
  if(report.means){
    out$means<-apply(out$chains,c(2,3),mean)
  }

  #add rate deviation posterior probabilities
  if(report.devs){
    if(nchain==1){
      out$post.probs<-apply(out%chains%'dev',2,function(ii) sum(ii>0)/length(ii))
      out$post.probs<-as.matrix(out$post.probs)
      dimnames(out$post.probs)<-list(parameters=rownames(out$post.probs),chains='chain 1')
    }else{
      out$post.probs<-apply(out%chains%'dev',c(2,3),function(ii) sum(ii>0)/length(ii))
    }
  }
  
  if(return.stanfit){
    out$stanfit<-ret
  }

  out
}