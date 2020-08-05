#' Fit a correlated rates Brownian motion model to a tree and trait data
#'
#' This function processes tree and trait data and fits an correlated rates Brownian motion model to
#' it using a Stan-based mcmc sampler. Under an correlated rates Brownian motion model, the rates at
#' trait evolution themselves evolve according to a Brownian motion process, allowing rates of
#' evolution to vary across the tree in a phylogenetically autocorrelated manner. The results are
#' formatted into a (relatively) user-friendly list of arrays for further analysis and diagnostics.
#'
#' @param tree A phylogenetic tree of class 'phylo'; function does not yet support multiple phylogenies
#' but will handle polytomous trees just fine. Tree must be in cladewise order with all edges of length
#' of length 0 collapsed into polytomies; running your tree through the function \code{prep.tree} first
#' takes care of this.
#' @param trait.data A numeric vector/matrix or data.frame of trait data. Rows of data (or elements of vector) must
#' correspond to observations and be appropriately labelled with their corresponding tip label as given in 
#' \code{tree}. Columns must correspond to different traits. Missing data must be coded as \code{NA}.
#' If a vector, labels are taken from the vector's \code{names} attribute. If a
#' matrix, labels are taken from the matrix's \code{rownames} attribute. If a data.frame, it is
#' best practice to have a column corresponding to tip labels with \code{colname} 'tip.label', as R doesn't
#' allow non-unique rownames in data.frame objects (thus one can't have multiple observations per tip).
#' If no 'tip.label' column is found, labels are taken from the data.frame's \code{rownames} attribute.
#' @param R0.prior The standard deviation of the cauchy prior for the \code{R_0} parameter determining
#' the rate of trait evolution at the root.
#' @param Rsig2.prior The standard deviation of the cauchy prior for the \code{R_sig2} parameter
#' determining the rate at which rate variance accumulates over time (inversely related to the apparent
#' 'phylogenetic correlation' of rates).
#' @param X0.prior The standard deviation of the cauchy prior(s) for the parameters determining
#' the trait value(s) at the root. In the case of multiple traits, can be a vector to specify different
#' standard deviations for each trait, in the same order as the corresponding columns in \code{trait.data}.
#' Will be recycled when length of vector and number of traits don't match.
#' 
#' 
#' @param intra.var TRUE or FALSE value: should tip means be estimated with an error term or fixed
#' according to the data?
#' @param intrasig2.prior The standard deviation of the cauchy prior for the \code{Ysig2} parameter
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
#' reporting in the R console (\code{iter/10} by default). See \code{sampling} for more options.
#' 
#' @return A list of class 'corateBM_fit', containing an array 'chains', with rows corresponding to
#' iterations in the mcmc sampler and columns to parameters. The array is 3D, with each
#' matrix slice corresponding to a particular mcmc chain. Parameters may include \code{R0}, \code{Rsig2},
#' \code{X0}, \code{Ysig2}, and \code{Rmu}, as defined above. Average branch-wise rates are also included
#' and denoted
#' \code{R_i} where \code{i} is the index of the edge as given in \code{tree}. If \code{intra.var} is
#' set to TRUE, the array will also include the tip mean parameters, denoted by their tip labels as given
#' in \code{tree}. Additionally, the function calculates the posterior distribution of 'background rates'
#' by taking the average of all average branch-wise rate parameters, weighted by the relative lengths of 
#' their respective branches, denoted \code{bg_rate}.
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
#' 
#' @export
#5/14: add trace and profile plotting functions
#rename Y to dat, Ysig2.prior to intrasig2.prior, Xsig2.prior to evosig2.prior,
#Xcor.prior to evocor.prior, and Ycor.prior to intracor.prior... --done!
#additional out components:
##diagnostics--Rhats, effective sample sizes, stepsizes, energy, lp__,... --done!
##MAPs? --done!
##data used for MCMC --done!
#report warnings BEFORE starting stan sampler? --done!
fit.corateBM<-function(tree,trait.data,R0.prior=10,Rsig2.prior=20,X0.prior=100,
                       intra.var=F,intrasig2.prior=50,
                       trend=F,Rmu.prior=Rsig2.prior,
                       evosig2.prior=1,evocor.prior=1,intracor.prior=1,
                       constrain.Rsig2=F,
                       report.quantiles=c(0.025,0.5,0.975),report.means=TRUE,report.MAPs=TRUE,report.devs=TRUE,
                       return.stanfit=F,
                       include.warmup=F,
                       ...){
  #initial checks
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  Y<-trait.data
  Ysig2.prior<-intrasig2.prior
  Xsig2.prior<-evosig2.prior
  Xcor.prior<-evocor.prior
  Ycor.prior<-intracor.prior
  if(is.data.frame(Y)){
    if(!is.null(Y$tip.label)){
      labels<-Y$tip.label
      labels.col<-which(colnames(Y)=='tip.label')
    }else if(!is.null(rownames(Y))){
      labels<-rownames(Y)
      warning("trait data is a data.frame with no 'tip.label' column: used rownames as labels instead",
              immediate.=T)
    }else{
      stop('trait data is unlabelled: please name each row with its corresponding tip label or add tip.label column')
    }
    non.num.cols<-sapply(1:ncol(Y),function(ii) !is.numeric(Y[,ii]))
    if(any(non.num.cols)){
      if(all(non.num.cols)){
        stop('no numeric trait data provided')
      }
      non.num.cols.inds<-which(non.num.cols)
      if(exists('labels.col')){
        non.num.cols.inds<-non.num.cols.inds[non.num.cols.inds!=labels.col]
        tmp.non.num.cols.inds<-c(non.num.cols.inds[non.num.cols.inds<labels.col],
                                 non.num.cols.inds[non.num.cols.inds>labels.col]-1)
        if(length(X0.prior)==ncol(Y)-1){
          X0.prior<-X0.prior[-tmp.non.num.cols.inds]
        }
        if(length(Xsig2.prior)==ncol(Y)-1){
          Xsig2.prior<-Xsig2.prior[-tmp.non.num.cols.inds]
        }
        if(length(Ysig2.prior)==ncol(Y)-1){
          Ysig2.prior<-Ysig2.prior[-tmp.non.num.cols.inds]
        }
        Y<-Y[,-c(labels.col,non.num.cols.inds)]
        if(length(non.num.cols.inds)>0){
          warning('trait data in column(s) ',
                  paste(non.num.cols.inds,collapse=', '),' is non-numeric: these data were removed',
                  immediate.=T)
        }
      }else{
        if(length(X0.prior)==ncol(Y)-1){
          X0.prior<-X0.prior[-non.num.cols.inds]
        }
        if(length(Xsig2.prior)==ncol(Y)-1){
          Xsig2.prior<-Xsig2.prior[-non.num.cols.inds]
        }
        if(length(Ysig2.prior)==ncol(Y)-1){
          Ysig2.prior<-Ysig2.prior[-non.num.cols.inds]
        }
        Y<-Y[,-non.num.cols.inds]
        warning('trait data in column(s) ',
                paste(non.num.cols.inds,collapse=', '),'is non-numeric: these data were removed',
                immediate.=T)
      }
    }
    Y<-as.matrix(Y)
    rownames(Y)<-labels
  }else if(is.numeric(Y)){
    is.mat<-is.matrix(Y)
    Y<-as.matrix(Y)
    if(is.null(rownames(Y))){
      stop('trait data is unlabelled: please name each ',
           if(is.mat) 'row' else 'element',' with its corresponding tip label')
    }
  }else{
    stop("trait data format isn't recognized: please format it as a numeric vector/matrix or data.frame: see help file for more details")
  }
  k<-ncol(Y)
  colname.lens<-nchar(colnames(Y))
  if(is.null(colnames(Y))){
    colnames(Y)<-paste('X',1:k,sep='')
  }else if(any(colname.lens==0)){
    colnames(Y)<-ifelse(colname.lens==0,paste('X',1:k,sep=''),colnames(Y))
  }
  unmatched.names<-!(rownames(Y)%in%tree$tip.label)
  if(sum(unmatched.names)>0){
    if(all(unmatched.names)){
      stop('no matches between tree tip and trait data labels')
    }
    ind.unmatched.names<-which(unmatched.names)
    warning('could not find tree tip labels matching with ',
            paste(unique(rownames(Y)[unmatched.names]),collapse=', '),': these data were removed',
            immediate.=T)
    Y<-as.matrix(Y[-ind.unmatched.names,])
  }
  if(constrain.Rsig2&!trend&report.devs){
    report.devs<-F
    warning('report.devs was set to FALSE: rate deviations are meaningless with no rate heterogeneity',
            immediate.=T)
  }
  
  #process/format data
  n<-length(tree$tip.label)
  e<-nrow(tree$edge)
  eV<-edge.vcv(tree)
  if(intra.var){
    is.obs<-!as.matrix(apply(Y,1,is.na))
    if(k>1){
      is.obs<-t(is.obs)
    }
    codes<-apply(is.obs,1,function(ii) sum(ii*10^((length(ii)-1):0)))
    if(any(codes==0)){
      null.obs<-which(codes==0)
      tmp.codes<-codes
      Y<-Y[-null.obs,]
      is.obs<-is.obs[-null.obs,]
      codes<-codes[-null.obs]
      for(i in sort(ind.unmatched.names)){
        tmp.codes<-append(tmp.codes,-1,i)
      }
      report.null.obs<-which(tmp.codes==0)
      warning('found no trait data for observation(s) ',
              paste(report.null.obs,collapse=', '),': these data were removed')
    }
    X_id<-match(rownames(Y),tree$tip.label)
    code_sizes<-tapply(codes,codes,length)
    obs_code<-unlist(lapply(names(code_sizes),function(ii) which(codes==as.numeric(ii))))
    names(code_sizes)<-paste(lapply(k-nchar(names(code_sizes)),function(ii) paste(rep(0,ii),collapse='')),
                             names(code_sizes),sep='')
    n_code<-length(code_sizes)
    parsed.codes<-gregexpr('1',names(code_sizes))
    code_ks<-lengths(parsed.codes)
    code_key<-unlist(parsed.codes)
  }else{
    col.names<-colnames(Y)
    n.obs<-sapply(rownames(Y),function(ii) sum(rownames(Y)==ii))
    if(any(n.obs>1)){
      warning("multiple observations per tip found: observations were averaged, but it's strongly recommend to run function with intra.var set to TRUE",
              immediate.=T)
      Y<-as.matrix(sapply(1:k,function(ii) tapply(Y[,ii],rownames(Y),mean,na.rm=T)))
      Y[is.nan(Y)]<-NA
    }
    missing.Y<-names(which(sapply(tree$tip.label,function(ii) sum(rownames(Y)==ii))==0))
    Y<-do.call(rbind,c(list(Y),setNames(rep(list(NA),length(missing.Y)),missing.Y)))
    Y<-as.matrix(Y[tree$tip.label,])
    colnames(Y)<-col.names
    is.obs<-!as.matrix(apply(Y,1,is.na))
    if(k>1){
      is.obs<-t(is.obs)
    }
    which_mis<-lapply(asplit(Y,2),function(ii) which(is.na(ii)))
    if(length(which_mis)==0){
      k_mis<-rep(0,k)
    }else{
      k_mis<-lengths(which_mis)
      which_mis<-unlist(which_mis)
    }
  }
  untrans.tree<-tree
  untrans.Y<-Y
  o.hgt<-max(eV)
  tree$edge.length<-tree$edge.length/o.hgt
  eV<-eV/o.hgt
  o.Xsig2<-rep(NA,k)
  o.X0<-rep(NA,k)
  for(i in 1:k){
    tmp.Y<-Y[,i]
    missing.dat<-names(which(tapply(tmp.Y,names(tmp.Y),function(ii) sum(!is.na(ii)))==0))
    missing.dat<-c(missing.dat,tree$tip.label[!(tree$tip.label%in%names(tmp.Y))])
    if(length(missing.dat)>0){
      tmp.tree<-drop.tip(tree,missing.dat)
    }else{
      tmp.tree<-tree
    }
    tmp.X<-tapply(tmp.Y,names(tmp.Y),mean,na.rm=T)[tmp.tree$tip.label]
    o.Xsig2[i]<-sum(pic(tmp.X,multi2di(tmp.tree))^2)/n
    xx<-c(tmp.X,rep(NA,tmp.tree$Nnode))
    for(ee in nrow(tmp.tree$edge):0){
      if(ee==0){
        des<-which(tmp.tree$edge[,1]==n+1)
        xx[n+1]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
        break
      }
      nn<-tmp.tree$edge[ee,2]
      if(nn<=n){
        next
      }
      des<-which(tmp.tree$edge[,1]==nn)
      xx[nn]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
    }
    o.X0[i]<-xx[n+1]
    Y[,i]<-(Y[,i]-o.X0[i])/sqrt(o.Xsig2[i])
  }
  if(!is.binary(tree)){
    poly.nodes<-which(sapply(1:max(tree$edge),function(ii) length(which(ii==tree$edge[,1])))>2)
    d_poly<-lapply(poly.nodes,function(ii) which(ii==tree$edge[,1]))
    tmp<-sort(unlist(lapply(d_poly,function(ii) ii[seq(2,length(ii)-1)])))
    tmp<-tmp+0:(length(tmp)-1)
    real_e<-(1:(2*n-2))[-tmp]+1
    tree<-multi2di(tree,random=F)
  }else{
    real_e<-1:e+1
  }
  des_e<-sapply(1:nrow(tree$edge),function(ii) which(tree$edge[,1]==tree$edge[ii,2])+1)
  tip_e<-which(lengths(des_e)==0)
  tip_e<-tip_e[order(tree$edge[tip_e,2])]
  des_e[tip_e]<-list(rep(-1,2))
  des_e<-matrix(unlist(des_e),ncol=2,byrow=T)
  root.edges<-which(tree$edge[,1]==n+1)+1
  des_e<-rbind(root.edges,des_e)
  prune_T<-c(0,tree$edge.length)
  postorder<-((2*n-2):1)[-((2*n-2)-tip_e)]
  tip_e<-tip_e+1
  Y[!is.obs]<-0
  dat<-list('n'=n,'e'=e,'Y'=t(Y),'eV'=eV,
            'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,
            'R0_prior'=R0.prior,'Rsig2_prior'=Rsig2.prior,'Rmu_prior'=Rmu.prior,
            'X0_prior'=rep(X0.prior,length.out=k),
            'constr_Rsig2'=as.numeric(constrain.Rsig2),'constr_Rmu'=as.numeric(!trend))
  
  #run mcmc
  if(intra.var){
    dat$obs<-nrow(Y)
    datX_id<-X_id
    dat$Ysig2_prior<-rep(Ysig2.prior,length.out=k)
    if(k>1){
      dat$k<-k
      dat$obs_code<-obs_code
      dat$n_code<-n_code
      dat$code_sizes<-code_sizes
      dat$code_ks<-as.array(code_ks)
      dat$code_key<-code_key
      dat$Xsig2_prior=rep(Xsig2.prior,length.out=k)
      dat$Xcor_prior=Xcor.prior
      dat$Ycor_prior=Ycor.prior
      ret<-sampling(object=stanmodels$intravar_multivar_corateBM,data=dat,...)
    }else{
      dat$Y<-as.vector(Y)
      ret<-sampling(object=stanmodels$intravar_univar_corateBM,data=dat,...)
    }
  }else{
    dat$which_mis<-which_mis
    dat$postorder<-postorder
    if(k>1){
      dat$k<-k
      dat$k_mis<-k_mis
      dat$Xsig2_prior=rep(Xsig2.prior,length.out=k)
      dat$Xcor_prior=Xcor.prior
      ret<-sampling(object=stanmodels$multivar_corateBM,data=dat,...)
    }else{
      dat$Y<-as.vector(Y)
      dat$mis<-length(which_mis)
      ret<-sampling(object=stanmodels$univar_corateBM,data=dat,...)
    }
  }

  #get sampler arguments/parameters for diagnostics
  sampler.args<-attr(ret@sim[[1]][[1]],'args')
  sampler.args$chains<-length(ret@sim[[1]])
  nchain<-sampler.args$call.chains<-sampler.args$chains
  sampler.args$call.iter<-sampler.args$iter
  sampler.args$call.warmup<-sampler.args$warmup
  sampler.args$iter<-floor(sampler.args$iter/sampler.args$thin)
  sampler.args$warmup<-floor(sampler.args$warmup/sampler.args$thin)
  sampler.args$call.thin<-sampler.args$thin
  if(include.warmup){
    niter<-sampler.args$iter
    excl.iter<-numeric(0)
  }else{
    niter<-sampler.args$iter-sampler.args$warmup+1
    excl.iter<-2:sampler.args$warmup
  }
  lp<-extract(ret,"lp__",permute=F,inc_warmup=T)
  sampler.params<-array(NA,c(sampler.args$iter,7,nchain))
  dimnames(sampler.params)<-list(iterations=NULL,
                                parameters=c(names(attr(ret@sim$samples[[1]],'sampler_params')),'lp__'),
                                chains=paste('chain',1:nchain))
  for(i in 1:nchain){
    sampler.params[,,i]<-as.matrix(do.call(cbind,
                                           c(attr(ret@sim$samples[[i]],'sampler_params'),list(lp[,i,]))))
  }
  out<-list(sampler.control=sampler.args[-2],sampler.params=sampler.params)

  #process/format output
  out$chains<-array(dim=c(niter,4+e+k*(n+3+(k-1)),nchain))
  X.param.names<-paste(rep(colnames(Y),n),'_',rep(tree$tip.label,each=k),sep='')
  Xcov.param.names<-paste(rep(colnames(Y),each=k),',',
                          rep(colnames(Y),k),'_evocov',
                          sep='')[lower.tri(matrix(NA,k,k),diag=T)]
  Ycov.param.names<-paste(rep(colnames(Y),each=k),',',
                          rep(colnames(Y),k),'_intracov',
                          sep='')[lower.tri(matrix(NA,k,k),diag=T)]
  param.names<-c('R_0','R_sig2','R_mu','bg_rate',
                 paste('R_',1:e,sep=''),
                 paste(colnames(Y),'_0',sep=''),Xcov.param.names,Ycov.param.names,X.param.names)
  dimnames(out$chains)<-list(iterations=NULL,
                             parameters=param.names,
                             chains=paste('chain',1:nchain))
  X0<-.index.element(extract(ret,"X0",permute=F,inc_warmup=T),excl.iter,1,T)*
    rep(sqrt(o.Xsig2),each=niter*nchain)+rep(o.X0,each=niter*nchain)
  out$chains[,paste(colnames(Y),'_0',sep=''),]<-aperm(X0,c(1,3,2))
  out$call$X0.prior<-dat$X0_prior
  R0<-.index.element(extract(ret,"R0",permute=F,inc_warmup=T),excl.iter,1,T)-
    log(o.hgt)
  if(k>1){
    chol_Xcov<-.index.element(extract(ret,"chol_Xcov",permute=F,inc_warmup=T),excl.iter,1,T)
    Xcov<-array(0,c(k,k,niter*nchain))
    for(i in 1:k){
      for(j in 1:i){
        Xcov[i,j,]<-as.vector(chol_Xcov[,,paste('chol_Xcov[',i,',',j,']',sep='')])
      }
    }
    Xcov<-lapply(asplit(Xcov,3),function(ii) ii%*%t(ii))
    Xcov<-lapply(Xcov,function(ii) diag(sqrt(o.Xsig2))%*%ii%*%diag(sqrt(o.Xsig2)))
    new.fac<-sapply(1:length(Xcov),function(ii) 1/mean(diag(Xcov[[ii]])))
    Xcov<-lapply(1:length(Xcov),function(ii) Xcov[[ii]]*new.fac[ii])
    Xcov<-array(unlist(Xcov),dim=c(k,k,length(Xcov)))
    tmp<-which(lower.tri(matrix(NA,k,k),diag=T),arr.ind=T)
    for(i in 1:length(Xcov.param.names)){
      out$chains[,Xcov.param.names[i],]<-Xcov[tmp[i,1],tmp[i,2],]
    }
    out$call$evosig2.prior<-dat$Xsig2_prior
    out$call$evocor.prior<-dat$Xcor_prior
    R0<-R0-log(new.fac)
    if(intra.var){
      chol_Ycov<-.index.element(extract(ret,"chol_Ycov",permute=F,inc_warmup=T),excl.iter,1,T)
      Ycov<-array(0,c(k,k,niter*nchain))
      for(i in 1:k){
        for(j in 1:i){
          Ycov[i,j,]<-as.vector(chol_Ycov[,,paste('chol_Ycov[',i,',',j,']',sep='')])
        }
      }
      Ycov<-lapply(asplit(Ycov,3),function(ii) ii%*%t(ii))
      Ycov<-lapply(Ycov,function(ii) diag(sqrt(o.Xsig2))%*%ii%*%diag(sqrt(o.Xsig2)))
      Ycov<-array(unlist(Ycov),dim=c(k,k,length(Ycov)))
      for(i in 1:length(Ycov.param.names)){
        out$chains[,Ycov.param.names[i],]<-Ycov[tmp[i,1],tmp[i,2],]
      }
      out$call$intrasig2.prior<-dat$Ysig2_prior
      out$call$intracor.prior<-dat$Ycor_prior
    }
  }else{
    if(intra.var){
      Ysig2<-.index.element(extract(ret,"Ysig2",permute=F,inc_warmup=T),excl.iter,1,T)*
        sqrt(o.Xsig2)
      out$chains[,Ycov.param.names,]<-Ysig2
      out$call$intrasig2.prior<-dat$Ysig2_prior
    }
    R0<-R0+log(o.Xsig2)
  }
  out$chains[,'R_0',]<-R0
  out$call$R0.prior<-dat$R0_prior
  if(intra.var){
    out.X<-.index.element(extract(ret,"X",permute=F,inc_warmup=T),excl.iter,1,T)*
      rep(sqrt(o.Xsig2),each=niter*nchain)+rep(o.X0,each=niter*nchain)
    out$chains[,X.param.names,]<-aperm(out.X,c(1,3,2))
  }else if(length(which_mis)>0){
    mis.Y<-.index.element(extract(ret,"mis_Y",permute=F,inc_warmup=T),excl.iter,1,T)*
      rep(sqrt(o.Xsig2),niter*nchain*k_mis)+rep(o.X0,niter*nchain*k_mis)
    tmp<-(which_mis-1)*k+rep(1:k,k_mis)
    out$chains[,X.param.names[tmp],]<-aperm(mis.Y,c(1,3,2))
  }
  if(!constrain.Rsig2){
    Rsig2<-.index.element(extract(ret,"Rsig2",permute=F,inc_warmup=T),excl.iter,1,T)/
      o.hgt
    out$chains[,'R_sig2',]<-Rsig2
    out$call$Rsig2.prior<-dat$Rsig2_prior
  }
  if(trend){
    Rmu<-.index.element(extract(ret,"Rmu",permute=F,inc_warmup=T),excl.iter,1,T)/
      o.hgt
    out$chains[,'R_mu',]<-Rmu
    out$call$Rmu.prior<-dat$Rmu_prior
  }
  if(!constrain.Rsig2|trend){
    R<-.index.element(extract(ret,"R",permute=F,inc_warmup=T),excl.iter,1,T)-
      log(o.hgt)
    if(k>1){
      R<-R-log(new.fac)
    }else{
      R<-R+log(o.Xsig2)
    }
    wgts<-tree$edge.length/sum(tree$edge.length)
    if(nchain==1){
      bg.rate<-apply(R,1,function(ii) sum(ii*wgts))
    }else{
      bg.rate<-apply(R,c(1,2),function(ii) sum(ii*wgts))
    }
    out$chains[,'bg_rate',]<-bg.rate
    out$chains[,paste('R_',1:e,sep=''),]<-aperm(R,c(1,3,2))
  }
  incl.inds<-which(apply(out$chains[,,1],2,function(ii) !all(is.na(ii))))
  out$chains<-.index.element(out$chains,incl.inds,2)
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
    tmp[[2]]<-c(tmp[[2]],paste('R_',1:e,'_dev',sep=''))
    out$chains<-aperm(array(c(aperm(out$chains,c(1,3,2)),rate.devs),
                            dim=c(dim(out$chains)[1],nchain,dim(out$chains)[2]+e),
                            dimnames=tmp[c(1,3,2)]),
                      c(1,3,2))
  }
  
  #create parameter diagnostics table
  out$param.diags<-.index.element(out$chains,1:4,1)
  dimnames(out$param.diags)<-c(diagnostics=list(c('inits','bulk_ess','tail_ess','Rhat')),
                               dimnames(out$chains)[-1])
  if(!include.warmup){
    out$chains<-.index.element(out$chains,1,1,T)
  }
  out$param.diags[2,,]<-apply(out$chains,c(2,3),rstan::ess_bulk)
  out$param.diags[3,,]<-apply(out$chains,c(2,3),rstan::ess_tail)
  out$param.diags[4,,]<-apply(out$chains,c(2,3),rstan::Rhat)

  #add quantiles
  if(!is.null(report.quantiles)){
    out$quantiles<-.int.quantiles(out,c('.|dev'))
  }

  #add means
  if(report.means){
    out$means<-.int.means(out,c('.|dev'))
    out$means<-array(out$means,dim(out$means)[-1],dimnames(out$means)[-1])
  }
  
  #add MAPs
  if(report.MAPs){
    out$MAPs<-.int.MAPs(out,c('.|dev'))
    out$MAPs<-array(out$MAPs,dim(out$MAPs)[-1],dimnames(out$MAPs)[-1])
  }

  #add rate deviation posterior probabilities
  if(report.devs){
    out$post.probs<-apply(.int.chains(out,'R_\\d+_dev'),c(2,3),function(ii) sum(ii>0)/length(ii))
  }

  if(return.stanfit){
    out$stanfit<-ret
  }

  out$call$trait.data<-untrans.Y
  out$call$tree<-untrans.tree

  out
}

#returns error when thin!=1, probs need to distinguish between 'true' iter and call.iter...same for warmup