#need to have some handlers to make sure lik_power greater than 0 and T/F values are T/F
#need to add support for forcing proper priors, i.e., making sure any missing data has an associated prior, coded via an
#X.prior.mu and X.prior.sig, handled in much the same way as X0.prior.mu and X0.prior.sig
#' @export
form.input<-function(tree,trait.data,intra.var=F,ensure.prop.prior=T,
                     constrain.Rsig2=F,trend=F,lik.power=1,
                     X0.prior.mu=NULL,X0.prior.sig=NULL,
                     X.prior.mu=NULL,X.prior.sig=NULL,scale.X.prior=T,
                     evosig2.prior=NULL,evocor.prior=NULL,
                     R0.prior.mu=NULL,R0.prior.sig=NULL,Rsig2.prior=NULL,
                     intrasig2.prior=NULL,intracor.prior=NULL,
                     Rmu.prior.mu=NULL,Rmu.prior.sig=NULL){
  
  ##INITIAL DATA/PRIOR CLEAN-UP##
  #make sure any prior scale parameters are positive
  for(i in c('X0.prior.sig',
             'X.prior.sig',
             'evosig2.prior','evocor.prior',
             'R0.prior.sig','Rsig2.prior',
             'intrasig2.prior','intracor.prior',
             'Rmu.prior.sig')){
    tmp<-get(i)
    if(any(tmp<0)){
      warning(i,' must consist of postive numbers: entries that were negative or 0 set to defaults',
              immediate.=T)
      if(length(tmp)==1){
        tmp<-NULL
      }else{
        tmp[tmp<=0]<-NA
      }
      assign(i,tmp)
    }
  }
  if(lik.power<0){
    warning('likelihood powers must be 0 or positive: likelihood power set to default of 1 (i.e., sampling from normal posterior)',
            immediate.=T)
    lik.power<-1
  }
  #coerce trait.data into numeric matrix, make priors match it
  tmp<-.format.trait.data(trait.data,
                          X0.prior.mu,X0.prior.sig,
                          X.prior.mu,X.prior.sig,scale.X.prior,
                          evosig2.prior,
                          intrasig2.prior,tree)
  for(i in names(tmp)){
    assign(i,tmp[[i]])
  }
  #coerce tree to be compatible (collapse internal 0-length edges to polytomies, identical tips to single tips)
  tmp<-.format.tree(tree,Y)
  for(i in names(tmp)){
    assign(i,tmp[[i]])
  }
  if(!is.null(scale.X.prior)){
    scale.X.prior[is.na(scale.X.prior)]<-T
  }
  #split tip priors out from trait.data and clean everything up (tip prior mus must have sigs, if intra.var is
  #FALSE a tip cannot be assigned data and a tip prior)
  tmp<-.split.tp(Y,tree,intra.var)
  for(i in names(tmp)){
    assign(i,tmp[[i]])
  }
  if(length(Y)==0&length(tp.mu)==0){
    stop('no valid trait data was provided')
  }
  call<-list(tree=tree,trait.data=Y,tip.prior.means=tp.mu,tip.prior.sd=tp.sig)
  
  ##DATA/PRIOR TRANSFORMATIONS##
  tmp<-.get.trans.const(tree,Y,tp.mu,tp.sig,
                        X0.prior.mu,X0.prior.sig,
                        X.prior.mu,X.prior.sig,scale.X.prior,
                        R0.prior.mu,Rsig2.prior,
                        Ysig2.prior,
                        Rmu.prior.mu,Rmu.prior.sig)
  for(i in names(tmp)){
    assign(i,tmp[[i]])
  }
  
  ##FORMAT DATA FOR STAN##
  #basics
  n<-length(tree$tip.label)
  e<-nrow(tree$edge)
  k<-ncol(Y)
  if(ensure.prop.prior){
    if(is.null(X.prior.mu)){
      X.prior.mu<-rep(NA,k)
    }
    X.prior.mu<-ifelse(is.na(X.prior.mu),0,X.prior.mu)
    if(is.null(X.prior.sig)){
      X.prior.sig<-rep(NA,k)
    }
    X.prior.sig<-ifelse(is.na(X.prior.sig),100,X.prior.sig)
    if(!intra.var){
      missing.tp<-tree$tip.label[!(tree$tip.label%in%rownames(Y))&!(tree$tip.label%in%rownames(tp.mu))]
      if(length(missing.tp)>0){
        tmp.mu<-matrix(rep(X.prior.mu,each=length(missing.tp)),length(missing.tp),k)
        tmp.sig<-matrix(rep(X.prior.sig,each=length(missing.tp)),length(missing.tp),k)
        rownames(tmp.mu)<-rownames(tmp.sig)<-missing.tp
        tp.mu<-rbind(tp.mu,tmp.mu)
        tp.sig<-rbind(tp.sig,tmp.sig)
      }
    }else{
      tmp.mu<-matrix(rep(X.prior.mu,each=n),n,k)
      tmp.sig<-matrix(rep(X.prior.sig,each=n),n,k)
      rownames(tmp.mu)<-rownames(tmp.sig)<-tree$tip.label
      tp.mu<-rbind(tp.mu,tmp.mu)
      tp.sig<-rbind(tp.sig,tmp.sig)
    }
    call.tp.mu<-tp.mu
    call.tp.sig<-tp.sig
    for(i in 1:k){
      tp.mu[is.na(tp.mu[,i]),i]<-X.prior.mu[i]
      tp.mu[is.na(tp.sig[,i]),i]<-X.prior.sig[i]
      call.tp.mu[,i]<-tp.mu[,i]*sqrt(trans.const$X_sig2[i])+trans.const$X_0[i]
      call.tp.sig[,i]<-tp.sig[,i]*sqrt(trans.const$X_sig2[i])
    }
    call$tip.prior.means<-call.tp.mu
    call$tip.prior.sd<-call.tp.sig
  }
  eV<-edge.vcv(tree)
  #missing data
  if(intra.var){
    if(length(Y)==0){
      obs<-n_code<-0
      X_id<-code_sizes<-obs_code<-code_ks<-code_key<-integer(0)
    }else{
      obs<-nrow(Y)
      is.obs<-!as.matrix(apply(Y,1,is.na))
      if(k>1){
        is.obs<-t(is.obs)
      }
      if(obs==1){
        rownames(is.obs)<-rownames(Y)
      }
      codes<-apply(is.obs,1,function(ii) sum(ii*10^((length(ii)-1):0)))
      X_id<-match(rownames(Y),tree$tip.label)
      code_sizes<-tapply(codes,codes,length)
      obs_code<-unlist(lapply(names(code_sizes),function(ii) which(codes==as.numeric(ii))))
      names(code_sizes)<-paste(lapply(k-nchar(names(code_sizes)),function(ii) paste(rep(0,ii),collapse='')),
                               names(code_sizes),sep='')
      n_code<-length(code_sizes)
      parsed.codes<-gregexpr('1',names(code_sizes))
      code_ks<-lengths(parsed.codes)
      code_key<-unlist(parsed.codes)
      Y[!is.obs]<-0
    }
  }else{
    missing.Y<-names(which(sapply(tree$tip.label,function(ii) sum(rownames(Y)==ii))==0))
    Y<-do.call(rbind,c(list(Y),setNames(rep(list(NA),length(missing.Y)),missing.Y)))
    Y<-Y[tree$tip.label,,drop=F]
    is.obs<-!as.matrix(apply(Y,1,is.na))
    if(k>1){
      is.obs<-t(is.obs)
    }
    if(nrow(Y)==1){
      rownames(is.obs)<-rownames(Y)
    }
    which_mis<-lapply(asplit(Y,2),function(ii) which(is.na(ii)))
    if(length(which_mis)==0){
      k_mis<-rep(0,k)
    }else{
      k_mis<-lengths(which_mis)
      which_mis<-unlist(which_mis)
    }
    Y[!is.obs]<-0
  }
  #tip priors
  if(length(tp.mu)>0){
    if(intra.var){
      n_tp<-rep(NA,k)
      which_tp<-vector('list',k)
      for(i in 1:k){
        tmp<-!is.na(tp.mu[,i])
        n_tp[i]<-sum(tmp)
        which_tp[[i]]<-match(rownames(tp.mu)[tmp],tree$tip.label)
      }
      which_tp<-unlist(which_tp)
      tp_mu<-as.vector(tp.mu[!is.na(tp.mu)])
      tp_sig<-as.vector(tp.sig[!is.na(tp.mu)])
    }else{
      inds<-which(!is.na(tp.mu),arr.ind=T)
      inds<-paste(colnames(tp.mu)[inds[,2]],rownames(tp.mu)[inds[,1]],sep='.')
      n_tp<-length(inds)
      which_tp<-match(inds,names(which_mis))
      tp_mu<-as.vector(tp.mu[!is.na(tp.mu)])
      tp_sig<-as.vector(tp.sig[!is.na(tp.mu)])
    }
  }else{
    if(intra.var){
      n_tp<-rep(0,k)
    }else{
      n_tp<-0
    }
    which_tp<-integer(0)
    tp_mu<-numeric(0)
    tp_sig<-numeric(0)
  }
  #tree
  #ah, figured out the source of the issue, the below function does not work when some of the polytomies
  #lead to tips...
  #this actually seems to work, not sure what was going on earlier?
  #figured it out! Needs better way to handle the root edge if it's a polytomy...
  #is.binary is FALSE for a tree with three edges coming off root? Should work now...
  poly.nodes<-which(sapply(1:max(tree$edge),function(ii) length(which(ii==tree$edge[,1])))>2)
  if(length(poly.nodes)>0){
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
  #misc
  constr_Rsig2<-as.numeric(constrain.Rsig2)
  constr_Rmu<-as.numeric(!trend)
  lik_power<-lik.power
  #final data list
  dat<-list(n=n,e=e,k=k,Y=t(Y),eV=eV,
            n_tp=n_tp,which_tp=which_tp,tp_mu=tp_mu,tp_sig=tp_sig,
            prune_T=prune_T,des_e=des_e,tip_e=tip_e,real_e=real_e,postorder=postorder,
            X0_prior_mu=X0.prior.mu,X0_prior_sig=X0.prior.sig,Xsig2_prior=Xsig2.prior,Xcor_prior=evocor.prior,
            R0_prior_mu=R0.prior.mu,R0_prior_sig=R0.prior.sig,Rsig2_prior=Rsig2.prior,
            Ysig2_prior=Ysig2.prior,Ycor_prior=intracor.prior,
            Rmu_prior_mu=Rmu.prior.mu,Rmu_prior_sig=Rmu.prior.sig,
            constr_Rsig2=constr_Rsig2,constr_Rmu=constr_Rmu,
            lik_power=lik_power)
  if(intra.var){
    extra.dat<-list(obs=obs,X_id=X_id)
    if(k>1){
      extra.dat<-c(extra.dat,
                   list(n_code=n_code,code_sizes=code_sizes,obs_code=obs_code,code_ks=array(code_ks),code_key=code_key))
    }
  }else{
    extra.dat<-list(k_mis=k_mis,which_mis=which_mis)
  }
  dat<-c(dat,extra.dat)
  
  ##OUTPUT##
  list(call=call,trans.const=trans.const,dat=dat)
}

#kinda a shitty output when sampling prior --> maybe change some priors to normal dists or tighten them?
#Ysig2 --> 5/10, X0 --> 10/20, R0 --> 7, Rsig2 --> Rmu --> 10
#maybe keep things as they are but drop the cauchy...
#might want to add better handling of extra stan arguments (don't allow exclude_pars, include, or chain_id 
#[I think this will break later functions?])
#So new behavior: if return.as.obj is FALSE, results are always saved as files. If out.file isn't NULL,
#results are always saved as files. To skip exporting files, return.as.obj must be TRUE and out.file must
#be NULL
#better way to explain it: files will not be saved unless out.file isn't NULL; however, if return.as.obj is
#FALSE, a name for out.file will automatically be picked based on time/date
#' @export
run.corateBM<-function(corateBM.form.dat,return.as.obj=T,out.file=NULL,check.overwrite=T,...,
                       constrain.Rsig2=NULL,trend=NULL,lik.power=NULL,
                       X0.prior.mu=NULL,X0.prior.sig=NULL,
                       evosig2.prior=NULL,evocor.prior=NULL,
                       R0.prior.mu=NULL,R0.prior.sig=NULL,Rsig2.prior=NULL,
                       intrasig2.prior=NULL,intracor.prior=NULL,
                       Rmu.prior.mu=NULL,Rmu.prior.sig=NULL,
                       centered=F){
  if(hasArg(chains)){
    nchain<-list(...)$chains
  }else{
    nchain<-4
  }
  prep<-.prep.run(corateBM.form.dat,return.as.obj,out.file,check.overwrite,nchain=nchain,
                  constrain.Rsig2,trend,lik.power,
                  X0.prior.mu,X0.prior.sig,
                  evosig2.prior,evocor.prior,
                  R0.prior.mu,R0.prior.sig,Rsig2.prior,
                  intrasig2.prior,intracor.prior,
                  Rmu.prior.mu,Rmu.prior.sig)
  if(centered&grepl('intravar',prep$stanobj)){
    prep$stanobj<-paste(prep$stanobj,'_centered',sep='')
  }else{
    warning('centered parameterizations not allowed for models without intraspecific variation')
  }
  if(return.as.obj){
    ret<-sampling(object=stanmodels[[prep$stanobj]],data=prep$dat,
                  pars=prep$exclude.pars,include=F,sample_file=prep$out.file,...)
  }else{
    sampling(object=stanmodels[[prep$stanobj]],data=prep$dat,
             pars=prep$exclude.pars,include=F,sample_file=prep$out.file,...)
  }
  if(!is.null(prep$out.file)){
    saveRDS(prep$corateBM.form.dat,prep$file.strings[length(prep$file.strings)])
    message('MCMC samples were saved as:\n\t',
            paste(basename(prep$file.strings[-length(prep$file.strings)]),'\n',collapse='\t'),
            'and input data/transformation constants were saved as:\n\t',
            basename(prep$file.strings[length(prep$file.strings)]),'\n',
            'in ',normalizePath(prep$directory,winslash='/'))
    #return output file basename
    gsub('_\\d+\\.csv$','',prep$file.strings[1])
  }
  if(return.as.obj){
    c(stanfit=ret,prep$corateBM.form.dat)
  }
}

#corateBM.run can either be a character specifying the name of sampling files or an object (the output from corateBM.run)
#filenames ending in _i.csv, .csv, or _ifno will be truncated
#specifying any of the specific components will overwrite those found in corateBM.run
#check to see if read_stan_csv checks for compatibility between chains...
#' @export
form.output<-function(corateBM.run,stanfit=NULL,call=NULL,trans.const=NULL,dat=NULL,include.warmup=F,
                      report.quantiles=c(0.025,0.5,0.975),report.means=TRUE,report.MAPs=TRUE,report.devs=TRUE){
  ##GATHER UP COMPONENTS##
  if(is.character(corateBM.run)){
    simple.name<-gsub('_\\d+\\.csv$|\\.csv$|_info$','',corateBM.run)
    directory<-dirname(simple.name)
    file.list<-list.files(directory)
    file.strings<-grep(paste0(basename(simple.name),'_\\d+\\.csv'),file.list,value=T)
    if(length(file.strings)>0){
      file.strings<-paste0(simple.name,
                           substr(file.strings,regexpr('_\\d+\\.csv',file.strings),nchar(file.strings)))
    }else{
      file.strings<-NULL
    }
    file.strings<-c(file.strings,paste0(simple.name,'_info'))
    if(length(file.strings>1)){
      tmp.stanfit<-rstan::read_stan_csv(file.strings[-length(file.strings)])
    }else{
      tmp.stanfit<-NULL
    }
    info<-try(readRDS(file.strings[length(file.strings)]))
    if(inherits(info,'try-error')){
      tmp.call<-tmp.trans.const<-tmp.dat<-NULL
    }else{
      for(i in c('call','trans.const','dat')){
        tmp<-info[[i]]
        if(!is.null(tmp)){
          assign(paste0('tmp.',i),tmp)
        }else{
          assign(paste0('tmp.',i),NULL)
        }
      }
    }
  }else if(is.list(corateBM.run)){
    for(i in c('stanfit','call','trans.const','dat')){
      tmp<-corateBM.run[[i]]
      if(!is.null(tmp)){
        assign(paste0('tmp.',i),tmp)
      }else{
        assign(paste0('tmp.',i),NULL)
      }
    }
  }
  for(i in c('stanfit','call','trans.const','dat')){
    cur.tmp<-get(i)
    pot.tmp<-get(paste0('tmp.',i))
    if(is.null(cur.tmp)){
      assign(i,pot.tmp)
    }
  }
  if(is.null(stanfit)){
    stop('no MCMC samples were provided: please recheck provided filenames and/or R objects')
  }
  if(is.null(call)){
    stop('no information on the call (i.e., trait data, tip priors, and tree) were provided: please recheck provided filenames and/or R objects')
  }
  if(is.null(trans.const)){
    trans.const<-.get.trans.const(call$tree,call$trait.data,call$tip.prior.means,call$tip.prior.sd)$trans.const
  }
  if(is.null(dat)){
    stop('no information on the inputted stan data were provided: please recheck provided filenames and/or R objects')
  }
  
  n<-length(call$tree$tip.label)
  e<-nrow(call$tree$edge)
  k<-ncol(call$trait.data)
  trait.names<-colnames(call$trait.data)
  tip.names<-call$tree$tip.label
  checks<-setNames(vector('logical',3),c(if(k>1) 'chol_Ycov' else 'Ysig2','Rsig2','Rmu'))
  for(i in 1:length(checks)){
    tmp.check<-try(extract(stanfit,names(checks)[i]),silent=T)
    if(inherits(tmp.check,'try-error')){
      checks[i]<-F
    }else{
      checks[i]<-T
    }
  }
  checks[2]<-!checks[2]
  #checks --> is there intraspecific variation?, is Rsig2 constrained?, is there a trend?
  if(checks[2]&!checks[3]&report.devs){
    report.devs<-F
    warning('report.devs was set to FALSE: rate deviations are meaningless with no rate heterogeneity',
            immediate.=T)
  }
  
  #consider missing Y (other than ones with tip priors) a form a data imputation--they are ignored when
  #lik_power is 0... (BUT, if this is a hierarchical model, shouldn't the multinormal still be influential as
  #a prior??? Ugh...)
  #alright, no--unfortunately, missing data needs a proper prior to be valid. In practice, this means that each missing
  #observation has to be initialized with a prior--might better hardcode this into the model later, but for now I'll just
  #modify the  form.input function
  if(!checks[1]){
    if(is.null(dat)){
      missing.Y<-names(which(sapply(tip.names,function(ii) sum(rownames(call$trait.data)==ii))==0))
      Y<-do.call(rbind,c(list(call$trait.data),setNames(rep(list(NA),length(missing.Y)),missing.Y)))
      Y<-Y[tip.names,,drop=F]
      is.obs<-!as.matrix(apply(Y,1,is.na))
      if(k>1){
        is.obs<-t(is.obs)
      }
      if(nrow(Y)==1){
        rownames(is.obs)<-rownames(Y)
      }
      which_mis<-lapply(asplit(Y,2),function(ii) which(is.na(ii)))
      if(length(which_mis)==0){
        k_mis<-rep(0,k)
      }else{
        k_mis<-lengths(which_mis)
        which_mis<-unlist(which_mis)
      }
    }else{
      k_mis<-dat$k_mis
      which_mis<-dat$which_mis
    }
  }
  
  #get sampler arguments/parameters for diagnostics
  sampler.args<-stanfit@stan_args[[1]]
  sampler.args<-sampler.args[-which(names(sampler.args)=='chain_id')]
  sampler.args$chains<-length(stanfit@sim[[1]])
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
  sampler.params<-array(NA,c(sampler.args$iter,9,nchain))
  dimnames(sampler.params)<-list(iterations=NULL,
                                 parameters=c(names(attr(stanfit@sim$samples[[1]],'sampler_params')),
                                              'prior','lik','post'),
                                 chains=paste('chain',1:nchain))
  for(i in 1:nchain){
    sampler.params[,1:6,i]<-unlist(attr(stanfit@sim$samples[[i]],'sampler_params'))
  }
  sampler.params[,7,]<-extract(stanfit,"prior",permute=F,inc_warmup=T)
  sampler.params[,8,]<-extract(stanfit,"lik",permute=F,inc_warmup=T)
  sampler.params[,9,]<-sampler.params[,7,]+dat$lik_power*sampler.params[,8,]
  out<-list(sampler.control=sampler.args[-2],sampler.params=sampler.params)
  out$call$intra.var<-unname(checks[1])
  out$call$constrain.Rsig2<-unname(checks[2])
  out$call$trend<-unname(checks[3])
  out$call$tree<-call$tree
  out$call$trait.data<-call$trait.data
  out$call$tip.prior.means<-call$tip.prior.means
  out$call$tip.prior.sd<-call$tip.prior.sd
  #hmmmm...this relies on dat--it may be best just to prevent doing this when lacking dat altogether
  out$call$lik_power<-dat$lik_power
  
  #process/format output
  out$chains<-array(dim=c(niter,4+e+k*(n+3+(k-1)),nchain))
  X.param.names<-paste(rep(trait.names,n),'_',rep(tip.names,each=k),sep='')
  Xcov.param.names<-paste(rep(trait.names,each=k),',',
                          rep(trait.names,k),'_evocov',
                          sep='')[lower.tri(matrix(NA,k,k),diag=T)]
  Ycov.param.names<-paste(rep(trait.names,each=k),',',
                          rep(trait.names,k),'_intracov',
                          sep='')[lower.tri(matrix(NA,k,k),diag=T)]
  param.names<-c('R_0','R_sig2','R_mu','bg_rate',
                 paste('R_',1:e,sep=''),
                 paste(trait.names,'_0',sep=''),Xcov.param.names,Ycov.param.names,X.param.names)
  dimnames(out$chains)<-list(iterations=NULL,
                             parameters=param.names,
                             chains=paste('chain',1:nchain))
  X0<-.index.element(extract(stanfit,"X0",permute=F,inc_warmup=T),excl.iter,1,T)*
    rep(sqrt(trans.const$X_sig2),each=niter*nchain)+rep(trans.const$X_0,each=niter*nchain)
  out$chains[,paste(trait.names,'_0',sep=''),]<-aperm(X0,c(1,3,2))
  out$call$X0_prior_mu<-dat$X0_prior_mu*sqrt(trans.const$X_sig2)+trans.const$X_0
  out$call$X0_prior_sig<-dat$X0_prior_sig*sqrt(trans.const$X_sig2)
  if(!is.null(out$call$X0_prior_mu)){
    names(out$call$X0_prior_mu)<-names(out$call$X0_prior_sig)<-trait.names
  }
  R0<-.index.element(extract(stanfit,"R0",permute=F,inc_warmup=T),excl.iter,1,T)-
    log(trans.const$hgt)
  if(k>1){
    chol_Xcov<-.index.element(extract(stanfit,"chol_Xcov",permute=F,inc_warmup=T),excl.iter,1,T)
    Xcov<-array(0,c(k,k,niter*nchain))
    for(i in 1:k){
      for(j in 1:i){
        Xcov[i,j,]<-as.vector(chol_Xcov[,,paste('chol_Xcov[',i,',',j,']',sep='')])
      }
    }
    Xcov<-lapply(asplit(Xcov,3),function(ii) ii%*%t(ii))
    Xcov<-lapply(Xcov,function(ii) diag(sqrt(trans.const$X_sig2))%*%ii%*%diag(sqrt(trans.const$X_sig2)))
    new.fac<-sapply(1:length(Xcov),function(ii) 1/mean(diag(Xcov[[ii]])))
    Xcov<-lapply(1:length(Xcov),function(ii) Xcov[[ii]]*new.fac[ii])
    Xcov<-array(unlist(Xcov),dim=c(k,k,length(Xcov)))
    tmp<-which(lower.tri(matrix(NA,k,k),diag=T),arr.ind=T)
    for(i in 1:length(Xcov.param.names)){
      out$chains[,Xcov.param.names[i],]<-Xcov[tmp[i,1],tmp[i,2],]
    }
    out$call$evosig2_prior<-dat$Xsig2_prior*sqrt(trans.const$X_sig2)
    if(!is.null(out$call$evosig2_prior)){
      names(out$call$evosig2_prior)<-trait.names
    }
    out$call$evocor_prior<-dat$Xcor_prior
    R0<-R0-log(new.fac)
    if(checks[1]){
      chol_Ycov<-.index.element(extract(stanfit,"chol_Ycov",permute=F,inc_warmup=T),excl.iter,1,T)
      Ycov<-array(0,c(k,k,niter*nchain))
      for(i in 1:k){
        for(j in 1:i){
          Ycov[i,j,]<-as.vector(chol_Ycov[,,paste('chol_Ycov[',i,',',j,']',sep='')])
        }
      }
      Ycov<-lapply(asplit(Ycov,3),function(ii) ii%*%t(ii))
      Ycov<-lapply(Ycov,function(ii) diag(sqrt(trans.const$X_sig2))%*%ii%*%diag(sqrt(trans.const$X_sig2)))
      Ycov<-array(unlist(Ycov),dim=c(k,k,length(Ycov)))
      for(i in 1:length(Ycov.param.names)){
        out$chains[,Ycov.param.names[i],]<-Ycov[tmp[i,1],tmp[i,2],]
      }
      out$call$intrasig2_prior<-dat$Ysig2_prior*sqrt(trans.const$X_sig2)
      if(!is.null(out$call$intrasig2_prior)){
        names(out$call$intrasig2_prior)<-trait.names
      }
      out$call$intracor_prior<-dat$Ycor_prior
    }
  }else{
    if(checks[1]){
      Ysig2<-.index.element(extract(stanfit,"Ysig2",permute=F,inc_warmup=T),excl.iter,1,T)*
        sqrt(trans.const$X_sig2)
      out$chains[,Ycov.param.names,]<-Ysig2
      out$call$intrasig2_prior<-dat$Ysig2_prior*sqrt(trans.const$X_sig2)
      names(out$call$intrasig2_prior)<-trait.names
    }
    R0<-R0+log(trans.const$X_sig2)
  }
  out$chains[,'R_0',]<-R0
  out$call$R0_prior_mu<-dat$R0_prior_mu-log(trans.const$hgt)+log(mean(trans.const$X_sig2))
  out$call$R0_prior_sig<-dat$R0_prior_sig
  if(checks[1]){
    out.X<-.index.element(extract(stanfit,"X",permute=F,inc_warmup=T),excl.iter,1,T)*
      rep(sqrt(trans.const$X_sig2),each=niter*nchain)+rep(trans.const$X_0,each=niter*nchain)
    out$chains[,X.param.names,]<-aperm(out.X,c(1,3,2))
  }else if(length(which_mis)>0){
    mis.Y<-.index.element(extract(stanfit,"mis_Y",permute=F,inc_warmup=T),excl.iter,1,T)*
      rep(sqrt(trans.const$X_sig2),niter*nchain*k_mis)+rep(trans.const$X_0,niter*nchain*k_mis)
    tmp<-(which_mis-1)*k+rep(1:k,k_mis)
    out$chains[,X.param.names[tmp],]<-aperm(mis.Y,c(1,3,2))
  }
  if(!checks[2]){
    Rsig2<-.index.element(extract(stanfit,"Rsig2",permute=F,inc_warmup=T),excl.iter,1,T)/
      trans.const$hgt
    out$chains[,'R_sig2',]<-Rsig2
    out$call$Rsig2_prior<-dat$Rsig2_prior/trans.const$hgt
  }
  if(checks[3]){
    Rmu<-.index.element(extract(stanfit,"Rmu",permute=F,inc_warmup=T),excl.iter,1,T)/
      trans.const$hgt
    out$chains[,'R_mu',]<-Rmu
    out$call$Rmu_prior_mu<-dat$Rmu_prior_mu/trans.const$hgt
    out$call$Rmu_prior_sig<-dat$Rmu_prior_sig/trans.const$hgt
  }
  if(!checks[2]|checks[3]){
    R<-.index.element(extract(stanfit,"R",permute=F,inc_warmup=T),excl.iter,1,T)-
      log(trans.const$hgt)
    if(k>1){
      R<-R-log(new.fac)
    }else{
      R<-R+log(trans.const$X_sig2)
    }
    wgts<-call$tree$edge.length/sum(call$tree$edge.length)
    if(nchain==1){
      bg.rate<-log(apply(R,1,function(ii) sum(exp(ii)*wgts)))
    }else{
      bg.rate<-log(apply(R,c(1,2),function(ii) sum(exp(ii)*wgts)))
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
  
  out
}

#so an abundance of tip priors leads to an erosion of the signal that supports rate heterogenous models,
#collapsing things to simple BM --> this actually makes sense, and renders your model conservative.
#I think that's a good thing?

#will be interesting to test it on fruit syndrome dataset...

#alright, normal priors appear to be fine--probably the way to go (cauchy gets weird)
#>rename unif_par in exclude.pars to std_par
#>change behavior to allow the use of tip priors in datasets with lik_power set to 0? Could always eliminate tip prior
#influence by manually removing them...easier that way, I think
#besides, the way you have it now, you get really odd behavior with sampling tip means...
#might try to eliminate the modeling of missing data in non-intravar models when lik_power set to 0?
###
#Okay, for now, I have decided to make lik_power remove influence of ALL trait data, priors included. This might be changed in
#the future.
#Further, I have decided to make non-intravar models not bother modeling missing trait means when lik_power set to 0, as
#modeling tip means is really not the goal of those models...the only reasonable alternative I can think of is making the
#model model all tips simultaneously (given that they technically have a multinormal distribution...), but that's just
#the intravar model
#Typically, I imagine folks will not be terribly interested in modeling specific tip means--it's more to integrate them out
#as nuisance parameters. If someone is worried about how their inference reacts to the presence of tip priors, it's probably
#better to simply remove tip priors from the call...

#interesting: messages print immediately, and don't tell you what function they're from...
#perhaps replace warnings in form.dat with messages?

#alright, so 2 things: writing to a file is substantially slower than doing it all within R and doesn't select parameters...
#may be best to give folks an option to run within R
#oh, never mind, it was fine --> still might be best to give folks the option, I think
#also, need to figure out a way to add .csv extension when only 1 chain is run
#file.rename should take care of this
#I think I fixed this by making sure a .csv is tacked onto the end of out.file if it doesn't already have it

#' @export
fit.corateBM<-function(tree,trait.data,intra.var=F,ensure.prop.prior=T,
                       constrain.Rsig2=F,trend=F,lik.power=1,
                       return.as.obj=T,out.file=NULL,check.overwrite=T,...,
                       include.warmup=F,report.quantiles=c(0.025,0.5,0.975),report.means=TRUE,report.MAPs=TRUE,report.devs=TRUE,
                       sampling.scale=F,
                       X0.prior.mu=NULL,X0.prior.sig=NULL,
                       X.prior.mu=NULL,X.prior.sig=NULL,scale.X.prior=T,
                       evosig2.prior=NULL,evocor.prior=NULL,
                       R0.prior.mu=NULL,R0.prior.sig=NULL,Rsig2.prior=NULL,
                       intrasig2.prior=NULL,intracor.prior=NULL,
                       Rmu.prior.mu=NULL,Rmu.prior.sig=NULL){
  if(sampling.scale){
    input<-form.input(tree,trait.data,intra.var,ensure.prop.prior,
                      constrain.Rsig2,trend,lik.power,
                      X.prior.mu=X.prior.mu,X.prior.sig=X.prior.sig,scale.X.prior=rep(F,NCOL(trait.data)))
  }else{
    input<-form.input(tree,trait.data,intra.var,ensure.prop.prior,
                      constrain.Rsig2,trend,lik.power,
                      X0.prior.mu,X0.prior.sig,
                      X.prior.mu,X.prior.sig,scale.X.prior,
                      evosig2.prior,evocor.prior,
                      R0.prior.mu,R0.prior.sig,Rsig2.prior,
                      intrasig2.prior,intracor.prior,
                      Rmu.prior.mu,Rmu.prior.sig)
  }
  if(hasArg(chains)){
    nchain<-list(...)$chains
  }else{
    nchain<-4
  }
  prep<-.prep.run(input,return.as.obj,out.file,check.overwrite,nchain,
                  constrain.Rsig2,trend,lik.power,
                  X0.prior.mu,X0.prior.sig,
                  evosig2.prior,evocor.prior,
                  R0.prior.mu,R0.prior.sig,Rsig2.prior,
                  intrasig2.prior,intracor.prior,
                  Rmu.prior.mu,Rmu.prior.sig)
  if(return.as.obj){
    ret<-sampling(object=stanmodels[[prep$stanobj]],data=prep$dat,
                  pars=prep$exclude.pars,include=F,sample_file=prep$out.file,...)
  }else{
    sampling(object=stanmodels[[prep$stanobj]],data=prep$dat,
             pars=prep$exclude.pars,include=F,sample_file=prep$out.file,...)
  }
  if(!is.null(prep$out.file)){
    saveRDS(prep$corateBM.form.dat,prep$file.strings[length(prep$file.strings)])
    message('MCMC samples were saved as:\n\t',
            paste(basename(prep$file.strings[-length(prep$file.strings)]),'\n',collapse='\t'),
            'and input data/transformation constants were saved as:\n\t',
            basename(prep$file.strings[length(prep$file.strings)]),'\n',
            'in ',normalizePath(prep$directory,winslash='/'))
    gsub('_\\d+\\.csv$','',prep$file.strings[1])
  }
  if(return.as.obj){
    run<-c(stanfit=ret,prep$corateBM.form.dat)
    form.output(run,stanfit=NULL,call=NULL,trans.const=NULL,dat=NULL,include.warmup,
                report.quantiles,report.means,report.MAPs,report.devs)
  }
}