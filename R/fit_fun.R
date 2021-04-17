#new plan (only thinking univariate-ly for now):
##trait.data and trait.se are separate objects
##Still must be labelled
##Can be named vectors, rownamed column matrices, or data.frames with a "tip.label" column
##NAs in trait.data are taken to be missing data (coded internally as 0)
##NAs in SE are taken to indicate unfixed SEs to be estimated during model fitting (coded internally as -1)
##Infs in SE are taken to indicate missing data (coded internally as -2)
##After coercing trait.data and trait.se to rownamed matrices...
####First make sure tree and trait.data is alright--lump together tips separated by 0 distance, remove trait.data with no
######corresponding tip
####Priority of conflicts for SE handling:
########tips with intra-tip observations (>1 observations) and non-unfixed SEs return warnings and set SEs to unfixed
########tips with any observations (>0 observations) and infinite SEs return warnings and set SEs to unfixed
########tips with no observations (0 observations) and a non-infinite SE return warnings and set SEs to infinite
####After making these corrections, check for any double SE specification in which there are two valid SEs specified for
######the same tip, and return error if found
##Additionally, might be good to make T/F "switch" that instead estimates fixed SE for each tip using intra-tip observations and
###using these as fixed SEs, rather than estimating them during model fitting
####Will need to have an additional argument specifying how to handle tips with single observations in this case
########Probably just something you can set to NA (for unfixing these) or some number (for default SE)
########Maybe NULL will use max SD found when calculating SEs for other tips?
##Handle priors with ...?
##Technically, standard errors in this context are standard variance --> should probably square given SEs to make this less
###confusing.
#' @export
input.evorates<-function(tree,trait.data,trait.se=NULL,constrain.Rsig2=FALSE,trend=FALSE,lik.power=1,sampling.scale=FALSE,
                         ...){
  ##INITIAL DATA/PRIOR CLEAN-UP##
  #make sure any prior scale parameters are positive
  prior.list<-do.call(.make.prior.list,list(...))
  #coerce trait.data and trait.se to rownamed matrices
  trait.data<-.coerce.trait.mat(trait.data)
  if(!is.null(trait.se)){
    trait.se<-.coerce.trait.mat(trait.se)
  }
  #coerce tree to be compatible (collapse internal 0-length edges to polytomies, identical tips to single tips)
  tmp<-.coerce.tree(tree,trait.data,trait.se)
  tree<-tmp$tree
  trait.data<-tmp$trait.data
  trait.se<-tmp$trait.se
  #check for "nits" (labels in trait info found to be NNNot IIIn TTTree)
  tmp<-sapply(rownames(trait.data),function(ii) ii%in%tree$tip.label)
  conflicts_nits<-rownames(trait.data)[!tmp]
  if(length(conflicts_nits)>0){
    trait.data<-trait.data[tmp,,drop=FALSE]
    warning('No tip in tree matched with label(s) ',
            paste(conflicts_nits,collapse=', '),
            ' in trait.data: values with aforementioned label(s) removed.',
            immediate.=TRUE)
  }
  tmp<-sapply(rownames(trait.se),function(ii) ii%in%tree$tip.label)
  conflicts_nits<-rownames(trait.se)[!tmp]
  if(length(conflicts_nits)>0){
    trait.se<-trait.se[tmp,,drop=FALSE]
    warning('No tip in tree matched with label(s) ',
            paste(conflicts_nits,collapse=', '),
            ' in trait.se: standard errors with aforementioned label(s) removed.',
            immediate.=TRUE)
  }
  
  
  ##TRAIT INFO##
  n<-length(tree$tip.label)
  e<-nrow(tree$edge)
  parsed.trait.data<-split(trait.data,rownames(trait.data))
  #return NA is every trait val is NA, otherwise ignore NA entries
  foo<-function(x){
    if(all(is.na(x))){
      NA
    }else{
      mean(x,na.rm=TRUE)
    }
  }
  X<-as.matrix(sapply(parsed.trait.data,foo)[tree$tip.label])
  n_obs<-lengths(parsed.trait.data)[tree$tip.label]
  X[is.na(X)]<-n_obs[is.na(X)]<-0
  perc_mis<-sum(n_obs==0)/n
  if(perc_mis==1){
    stop('No trait information provided: please double-check input and labelling.')
  }
  #will have to change this a bit for the multivariate case
  if(perc_mis>0.25){
    warning('You are about to run a model with greater than ',
            floor(perc_mis*100),
            '% of the tips missing trait information: this can be done, but it could indicate some problems with input/labelling.',
            immediate.=TRUE)
  }
  has_intra<-n_obs>1
  mis<-n_obs==0
  p_SE<-rep(-1,n)
  p_SE[mis]<- -2
  names(X)<-names(n_obs)<-names(p_SE)<-tree$tip.label
  if(length(trait.se)>0){
    #make sure SE info doesn't have any conflicts with trait info
    p_SE<-.coerce.trait.se(trait.se,p_SE,has_intra,mis)
  }
  p_SE[p_SE>0]<-p_SE[p_SE>0]^2
  inv_n_obs<-1/n_obs
  n_contra<-sum(lengths(parsed.trait.data)-1)
  foo<-function(x){
    l<-length(x)
    if(l==1){
      NULL
    }else{
      tmp.seq<-seq_along(x)
      rbind(x[-1]-cumsum(x[-l])/tmp.seq[-l],tmp.seq[-1]/tmp.seq[-l])
    }
  }
  tmp<-do.call(cbind,lapply(parsed.trait.data,foo))
  if(length(tmp)==0){
    contra<-contra_var<-vector('numeric',0)
  }else{
    contra<-tmp[1,]
    contra_var<-tmp[2,]
  }
  which_mis_SE<-which(p_SE==-1)
  n_mis_SE<-length(which_mis_SE)
  
  
  ##FORM CALL##
  out.trait.se<-as.matrix(p_SE)
  out.trait.se[out.trait.se==-1]<-NA
  out.trait.se[out.trait.se==-2]<-Inf
  call<-list(trait.data=trait.data,trait.se=out.trait.se,tree=tree)
  
  
  ##TRANSFORMATION CONSTANTS##
  if(hasArg(phy.sd)){
    phy.sd<-list(...)$phy.sd
  }else{
    phy.sd<-TRUE
  }
  trans.const<-.get.trans.const(tree,X,mis,phy.sd=phy.sd)
  tree$edge.length<-tree$edge.length/trans.const$hgt
  for(i in 1:ncol(X)){
    X[,i]<-X[,i]/sqrt(trans.const$X_sig2[i])
  }
  #not prepped for multivariate data
  p_SE[p_SE>0]<-p_SE[p_SE>0]/trans.const$X_sig2[1]
  if(length(contra)>1){
    contra<-contra/sqrt(trans.const$X_sig2[1])
  }
  if(!sampling.scale){
    if(!is.null(prior.list$Ysig2_prior_sig)){
      prior.list$Ysig2_prior_sig<-prior.list$Ysig2_prior_sig/trans.const$X_sig2
    }
    if(!is.null(prior.list$R0_prior_mu)){
      prior.list$R0_prior_mu<-prior.list$R0_prior_mu+log(trans.const$hgt)-log(mean(trans.const$X_sig2))
    }
    for(i in c('Rsig2_prior_sig','Rmu_prior_mu','Rmu_prior_sig')){
      if(!is.null(prior.list[[i]])){
        prior.list[[i]]<-prior.list[[i]]*trans.const$hgt
      }
    }
  }
  
  
  ##TREE INFO##
  eV<-edge.vcv(tree)
  poly.nodes<-which(sapply(1:max(tree$edge),function(ii) length(which(ii==tree$edge[,1])))>2)
  if(length(poly.nodes)>0){
    d_poly<-lapply(poly.nodes,function(ii) which(ii==tree$edge[,1]))
    tmp<-sort(unlist(lapply(d_poly,function(ii) ii[seq(2,length(ii)-1)])))
    tmp<-tmp+0:(length(tmp)-1)
    real_e<-(1:(2*n-2))[-tmp]+1
    tree<-multi2di(tree,random=FALSE)
  }else{
    real_e<-1:e+1
  }
  des_e<-sapply(1:nrow(tree$edge),function(ii) which(tree$edge[,1]==tree$edge[ii,2])+1)
  tip_e<-which(lengths(des_e)==0)
  tip_e<-tip_e[order(tree$edge[tip_e,2])]
  des_e[tip_e]<-list(rep(-1,2))
  des_e<-matrix(unlist(des_e),ncol=2,byrow=TRUE)
  root.edges<-which(tree$edge[,1]==n+1)+1
  des_e<-rbind(root.edges,des_e)
  prune_T<-c(0,tree$edge.length)
  which_0tip<-which(prune_T[real_e]==0)
  n_0tip<-length(which_0tip)
  postorder<-((2*n-2):1)[-((2*n-2)-tip_e)]
  tip_e<-tip_e+1
  XX<-vector('logical',2*n-1)
  XX[tip_e]<-ifelse(p_SE==-2,FALSE,TRUE)
  mis_code<-vector('numeric',n-1)
  which_non_mis<-vector('list',n-1)
  counter<-0
  for(i in postorder){
    des_X<-XX[des_e[i,]]
    counter<-counter+1
    mis_code[counter]<-sum(des_X)
    if(mis_code[counter]==0){
      XX[i]<-FALSE
      which_non_mis[[counter]]<-numeric(0)
    }else if(mis_code[counter]==1){
      XX[i]<-TRUE
      which_non_mis[[counter]]<-which(des_X)
    }else if(mis_code[counter]==2){
      XX[i]<-TRUE
      which_non_mis[[counter]]<-c(1,2)
    }
  }
  n_mis_2<-sum(mis_code==2)
  which_non_mis<-unlist(which_non_mis[mis_code==1])
  if(is.null(which_non_mis)){
    which_non_mis<-vector('integer',0)
  }
  n_mis_1<-length(which_non_mis)
  
  
  ##PUTTING IT ALL TOGETHER##
  dat<-list('n'=n,'e'=e,'eV'=eV,'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,'postorder'=postorder,
            'n_0tip'=n_0tip,'which_0tip'=which_0tip,
            'mis_code'=mis_code,'n_mis_2'=n_mis_2,'n_mis_1'=n_mis_1,'which_non_mis'=array(which_non_mis),
            'X'=as.vector(X),'p_SE'=p_SE,'inv_n_obs'=inv_n_obs,'n_contra'=n_contra,'contra'=contra,'contra_var'=contra_var,
            'n_mis_SE'=n_mis_SE,'which_mis_SE'=which_mis_SE)
  dat<-c(dat,prior.list)
  dat$constr_Rsig2<-as.numeric(constrain.Rsig2)
  dat$constr_Rmu<-as.numeric(!trend)
  if(lik.power<0){
    warning('lik.power must be 0 or positive: set to default of 1.',
            immediate.=TRUE)
    lik.power<-1
  }
  dat$lik_power<-lik.power
  list(call=call,trans.const=trans.const,dat=dat)
}

#kinda a shitty output when sampling prior --> maybe change some priors to normal dists or tighten them?
#Ysig2 --> 5/10, X0 --> 10/20, R0 --> 7, Rsig2 --> Rmu --> 10
#maybe keep things as they are but drop the cauchy...
#2/12/21 update: tightening priors isn't a bad idea, but changing from cauchy to normal priors made things
##substantially slower and worse--need the flexibility in priors with such a limited data scenario
#might want to add better handling of extra stan arguments (don't allow exclude_pars, include, or chain_id 
#[I think this will break later functions?])
#So new behavior: if return.as.obj is FALSE, results are always saved as files. If out.file isn't NULL,
#results are always saved as files. To skip exporting files, return.as.obj must be TRUE and out.file must
#be NULL
#better way to explain it: files will not be saved unless out.file isn't NULL; however, if return.as.obj is
#FALSE, a name for out.file will automatically be picked based on time/date
#' @export
run.evorates<-function(input.evorates.obj,return.as.obj=TRUE,out.file=NULL,check.overwrite=TRUE,
                       constrain.Rsig2=FALSE,trend=FALSE,lik.power=1,...){
  stan.args.list<-list(...)
  if(length(stan.args.list)>0){
    prior.args.names<-c('^(R[\\._]?0).*((mu)|(mean))','(R[\\._]?0).*((sig)|(sd))',
                        '^((intra)|(Y))[\\._]?((sig2)|(var))',
                        '^R[\\._]?((sig2)|(var))',
                        '^((R[\\._]?mu)|(trend)).*((mu)|(mean))','^((R[\\._]?mu)|(trend)).*((sig)|(sd))')
    which.prior.args<-apply(do.call(cbind,lapply(prior.args.names,grepl,x=names(stan.args.list))),1,any)
    stan.args.list<-stan.args.list[!which.prior.args]
    stan.args.list<-stan.args.list[!(names(stan.args.list)%in%
                                       c('object','data','pars','include','sample_file','chain_id','save_warmup'))]
  }
  if(!is.null(stan.args.list$chains)){
    nchain<-stan.args.list$chains
  }else{
    nchain<-4
  }
  prior.list<-do.call(.make.prior.list,list(...))
  prep<-.prep.run(input.evorates.obj,return.as.obj,out.file,check.overwrite,
                  constrain.Rsig2,trend,lik.power,prior.list,nchain)
  if(return.as.obj){
    ret<-do.call(sampling,
                 c(object=list(stanmodels[[prep$stanobj]]),
                   data=list(prep$input.evorates.obj$dat),
                   pars=list(prep$exclude.pars),
                   include=list(FALSE),
                   sample_file=list(prep$out.file),
                   stan.args.list))
  }else{
    do.call(sampling,
            c(object=list(stanmodels[[prep$stanobj]]),
              data=list(prep$input.evorates.obj$dat),
              pars=list(prep$exclude.pars),
              include=list(FALSE),
              sample_file=list(prep$out.file),
              stan.args.list))
  }
  if(!is.null(prep$out.file)){
    saveRDS(prep$input.evorates.obj,prep$file.strings[length(prep$file.strings)])
    message('MCMC samples were saved as:\n\t',
            paste(basename(prep$file.strings[-length(prep$file.strings)]),'\n',collapse='\t'),
            'and input data/transformation constants were saved as:\n\t',
            basename(prep$file.strings[length(prep$file.strings)]),'\n',
            'in ',normalizePath(prep$directory,winslash='/'))
    #return output file basename
    gsub('_\\d+\\.csv$','',prep$file.strings[1])
  }
  if(return.as.obj){
    c(stanfit=ret,prep$input.evorates.obj)
  }
}

#corateBM.run can either be a character specifying the name of sampling files or an object (the output from corateBM.run)
#filenames ending in _i.csv, .csv, or _ifno will be truncated
#specifying any of the specific components will overwrite those found in run.evorates
#check to see if read_stan_csv checks for compatibility between chains...
#' @export
output.evorates<-function(run.evorates.obj,stanfit=NULL,call=NULL,trans.const=NULL,dat=NULL,include.warmup=FALSE,
                          report.quantiles=c(0.025,0.5,0.975),report.means=TRUE,report.MAPs=TRUE,report.devs=TRUE){
  ##GATHER UP COMPONENTS##
  list.from.input<-list(stanfit=stanfit,call=call,trans.const=trans.const,dat=dat)
  list.from.run<-list(stanfit=NULL,call=NULL,trans.const=NULL,dat=NULL)
  #get any components specified by run.evorates object
  if(is.character(run.evorates.obj)){
    simple.name<-gsub('_\\d+\\.csv$|\\.csv$|_info$','',run.evorates.obj)
    directory<-dirname(simple.name)
    file.list<-list.files(directory)
    file.strings<-grep(paste0(basename(simple.name),'_\\d+\\.csv'),file.list,value=TRUE)
    if(length(file.strings)>0){
      file.strings<-paste0(simple.name,
                           substr(file.strings,regexpr('_\\d+\\.csv',file.strings),nchar(file.strings)))
    }else{
      file.strings<-NULL
    }
    file.strings<-c(file.strings,paste0(simple.name,'_info'))
    if(length(file.strings>1)){
      list.from.run$stanfit<-rstan::read_stan_csv(file.strings[-length(file.strings)])
    }
    info<-try(readRDS(file.strings[length(file.strings)]))
    if(!inherits(info,'try-error')){
      for(i in names(info)){
        if(!is.null(info[[i]])){
          list.from.run[[i]]<-info[[i]]
        }
      }
    }
  }else if(is.list(run.evorates.obj)){
    for(i in names(run.evorates.obj)){
      if(!is.null(run.evorates.obj[[i]])){
        list.from.run[[i]]<-run.evorates.obj[[i]]
      }
    }
  }
  #merge list from function input with list from evorates run object
  #function input overwrites evorates run object if it isn't NULL
  for(i in names(list.from.input)){
    if(!is.null(list.from.input[[i]])){
      list.from.run[[i]]<-list.from.input[[i]]
    }
  }
  stanfit<-list.from.run$stanfit
  call<-list.from.run$call
  trans.const<-list.from.run$trans.const
  dat<-list.from.run$dat
  if(is.null(stanfit)){
    stop('no MCMC samples were provided: please recheck provided filenames and/or R objects')
  }
  if(is.null(call)){
    stop('no information on the call (i.e., trait data/standard errors and tree) were provided: please recheck provided filenames and/or R objects')
  }
  if(is.null(trans.const)){
    trans.const<-.get.trans.const(call$tree,call$trait.data,call$trait.se)$trans.const
  }
  if(is.null(dat)){
    stop('no information on the inputted stan data were provided: please recheck provided filenames and/or R objects. This can be re-obtained by running input.evorates() with the SAME trait data, tree, and prior/lik.power settings!')
  }
  
  
  ##FORM CALL##
  checks<-setNames(vector('logical',3),c('Ysig2','Rsig2','Rmu'))
  for(i in names(checks)){
    tmp.check<-try(extract(stanfit,i),silent=TRUE)
    if(inherits(tmp.check,'try-error')){
      checks[i]<-FALSE
    }else{
      checks[i]<-TRUE
    }
  }
  has.intra<-checks[1]
  constrain.Rsig2<-!checks[2]
  trend<-checks[3]
  lik.power<-dat$lik_power
  if(constrain.Rsig2&!trend&report.devs){
    report.devs<-FALSE
    warning('report.devs was set to FALSE: rate deviations are meaningless with no rate heterogeneity',
            immediate.=TRUE)
  }
  out<-list()
  class(out)<-'corateBM_fit'
  out$call<-c(call,
              'constrain.Rsig2'=list(constrain.Rsig2),
              'trend'=list(trend),
              'lik.power'=lik.power)
  
  
  ##SAMPLER ARGUMENTS/PARAMETERS##
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
  sampler.params[,7,]<-extract(stanfit,"prior",permute=FALSE,inc_warmup=TRUE)
  sampler.params[,8,]<-extract(stanfit,"lik",permute=FALSE,inc_warmup=TRUE)
  sampler.params[,9,]<-sampler.params[,7,]+dat$lik_power*sampler.params[,8,]
  out$sampler.control<-sampler.args[-2]
  out$sampler.params<-sampler.params
  
  
  ##FORM CHAINS##
  e<-nrow(call$tree$edge)
  #will have to change back to some previous stuff for multivar extension
  out$chains<-array(dim=c(niter,5+e,nchain))
  #might consider changing param names for better intepretability...
  param.names<-c('R_0','Y_sig2','R_sig2','R_mu','bg_rate',paste('R_',1:e,sep=''))
  dimnames(out$chains)<-list(iterations=NULL,
                             parameters=param.names,
                             chains=paste('chain',1:nchain))
  R0<-.index.element(extract(stanfit,"R0",permute=FALSE,inc_warmup=TRUE),excl.iter,1,TRUE)-
    log(trans.const$hgt)+log(mean(trans.const$X_sig2))
  out$chains[,'R_0',]<-R0
  out$call$R0_prior_mu<-dat$R0_prior_mu-log(trans.const$hgt)+log(mean(trans.const$X_sig2))
  out$call$R0_prior_sig<-dat$R0_prior_sig
  if(has.intra){
    #will need to change for multivar
    Ysig2<-.index.element(extract(stanfit,"Ysig2",permute=FALSE,inc_warmup=TRUE),excl.iter,1,TRUE)*
      trans.const$X_sig2
    out$chains[,'Y_sig2',]<-Ysig2
    out$call$Ysig2_prior_sig<-dat$Ysig2_prior_sig*trans.const$X_sig2
  }
  if(!constrain.Rsig2){
    Rsig2<-.index.element(extract(stanfit,"Rsig2",permute=FALSE,inc_warmup=TRUE),excl.iter,1,TRUE)/
      trans.const$hgt
    out$chains[,'R_sig2',]<-Rsig2
    out$call$Rsig2_prior_sig<-dat$Rsig2_prior_sig/trans.const$hgt
  }
  if(trend){
    Rmu<-.index.element(extract(stanfit,"Rmu",permute=FALSE,inc_warmup=TRUE),excl.iter,1,TRUE)/
      trans.const$hgt
    out$chains[,'R_mu',]<-Rmu
    out$call$Rmu_prior_mu<-dat$Rmu_prior_mu/trans.const$hgt
    out$call$Rmu_prior_sig<-dat$Rmu_prior_sig/trans.const$hgt
  }
  if(!constrain.Rsig2|trend){
    R<-.index.element(extract(stanfit,"R",permute=FALSE,inc_warmup=TRUE),excl.iter,1,TRUE)-
      log(trans.const$hgt)+log(mean(trans.const$X_sig2))
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
  
  
  ##PARAMETER DIAGNOSTICS TABLE##
  out$param.diags<-.index.element(out$chains,1:4,1)
  dimnames(out$param.diags)<-c(diagnostics=list(c('inits','bulk_ess','tail_ess','Rhat')),
                               dimnames(out$chains)[-1])
  if(!include.warmup){
    out$chains<-.index.element(out$chains,1,1,TRUE)
  }
  out$param.diags[2,,]<-apply(out$chains,c(2,3),rstan::ess_bulk)
  out$param.diags[3,,]<-apply(out$chains,c(2,3),rstan::ess_tail)
  out$param.diags[4,,]<-apply(out$chains,c(2,3),rstan::Rhat)
  
  
  ##EXTRA POSTERIOR DISTRIBUTION INFO##
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
fit.evorates<-function(tree,trait.data,trait.se=NULL,constrain.Rsig2=FALSE,trend=FALSE,lik.power=1,sampling.scale=FALSE,
                       return.as.obj=TRUE,out.file=NULL,check.overwrite=TRUE,
                       include.warmup=FALSE,
                       report.quantiles=c(0.025,0.5,0.975),report.means=TRUE,report.MAPs=TRUE,report.devs=TRUE,
                       ...){
  input<-input.evorates(tree,trait.data,trait.se,constrain.Rsig2,trend,lik.power,sampling.scale,...)
  stan.args.list<-list(...)
  if(length(stan.args.list)>0){
    prior.args.names<-c('^(R[\\._]?0).*((mu)|(mean))','(R[\\._]?0).*((sig)|(sd))',
                        '^((intra)|(Y))[\\._]?((sig2)|(var))',
                        '^R[\\._]?((sig2)|(var))',
                        '^((R[\\._]?mu)|(trend)).*((mu)|(mean))','^((R[\\._]?mu)|(trend)).*((sig)|(sd))')
    which.prior.args<-apply(do.call(cbind,lapply(prior.args.names,grepl,x=names(stan.args.list))),1,any)
    stan.args.list<-stan.args.list[!which.prior.args]
    stan.args.list<-stan.args.list[!(names(stan.args.list)%in%
                                       c('object','data','pars','include','sample_file','chain_id','save_warmup'))]
  }
  if(!is.null(stan.args.list$chains)){
    nchain<-stan.args.list$chains
  }else{
    nchain<-4
  }
  prep<-.prep.run(input,return.as.obj,out.file,check.overwrite,
                  constrain.Rsig2,trend,lik.power,NULL,nchain)
  if(return.as.obj){
    ret<-do.call(sampling,
                 c(object=list(stanmodels[[prep$stanobj]]),
                   data=list(prep$input.evorates.obj$dat),
                   pars=list(prep$exclude.pars),
                   include=list(FALSE),
                   sample_file=list(prep$out.file),
                   stan.args.list))
  }else{
    do.call(sampling,
            c(object=list(stanmodels[[prep$stanobj]]),
              data=list(prep$input.evorates.obj$dat),
              pars=list(prep$exclude.pars),
              include=list(FALSE),
              sample_file=list(prep$out.file),
              stan.args.list))
  }
  if(!is.null(prep$out.file)){
    saveRDS(prep$input.evorates.obj,prep$file.strings[length(prep$file.strings)])
    message('MCMC samples were saved as:\n\t',
            paste(basename(prep$file.strings[-length(prep$file.strings)]),'\n',collapse='\t'),
            'and input data/transformation constants were saved as:\n\t',
            basename(prep$file.strings[length(prep$file.strings)]),'\n',
            'in ',normalizePath(prep$directory,winslash='/'))
    gsub('_\\d+\\.csv$','',prep$file.strings[1])
  }
  if(return.as.obj){
    run<-c(stanfit=ret,prep$input.evorates.obj)
    output.evorates(run,NULL,NULL,NULL,NULL,include.warmup,
                    report.quantiles,report.means,report.MAPs,report.devs)
  }
}