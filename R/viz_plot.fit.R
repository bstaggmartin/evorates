#later add support for plotting multiple samples from chains
#okay, can use select.extra to specify a particular iteration from chains too
#use chain to specify a specific chain using select.chains()
#if NULL, combine chains (with warning if trying to plot inits)
#god this function, could be SOOOO much nicer and cleaner-->I have to work on it more
#in particular, think about accomodating non-default rate and post.probs-->I think there's a better system!
#' @export
plot.evorates_fit<-function(fit,
                            element=c('means','quantiles','chains','MAPs','diagnostics'),
                            include.post.probs=TRUE,
                            select.extra=NULL,
                            chain=NULL,
                            R=NULL, #use to specify custom rates (deviations, rate residuals, etc.)
                            post.probs=NULL,
                            anc.recon=c('summary','sample'),
                            ...,
                            post.probs.args=NULL,
                            sim=NULL){
  if(is.null(sim)){
    if(!inherits(fit,'evorates_fit')){
      stop("fit must be a fitted evorates model (class 'evorates_fit')")
    }
    in.fit<-fit
    try.element<-try(match.arg(element,c('means','quantiles','chains','MAPs','diagnostics')),silent=TRUE)
    if(inherits(try.element,'try-error')){
      stop(element," is not an available element to extract from an evorates fit: please specify one of the following: 'means', 'quantiles', 'chains', 'MAPs', or 'diagnostics'")
    }
    element<-try.element
    
    nchains<-fit$sampler.control$chains
    niter<-dim(fit$chains)[1]
    
    #for now, only takes chains of R array-->probably should generalize
    e<-nrow(fit$call$tree$edge)
    if(!is.null(R)){
      R<-.expand.element(R)
      if(any(dim(R)!=c(niter,e,nchains))){
        stop('mismatch between dimensions of provided R array and expected dimensions: please double-check input')
      }
      fit[c('quantiles','means','MAPs')]<-NULL
      fit$chains<-evorates:::.expand.element(fit$chains)
      if(grepl(dimnames('R_[1-9]',fit$chains))){
        fit$chains[,paste('R',1:e,sep='_'),]<-R
      }else{
        tmp.dims<-dim(fit$chains)
        tmp.dims[2]<-tmp.dims[2]+e
        tmp.dimnames<-dimnames(fit$chains)
        tmp.dimnames[[2]]<-c(tmp.dimnames[[2]],paste('R',1:e,sep='_'))
        fit$chains<-array(c(aperm(fit$chains,c(1,3,2)),aperm(R,c(1,3,2))),
                          tmp.dims[c(1,3,2)],
                          tmp.dimnames[c(1,3,2)])
        fit$chains<-aperm(fit$chains,c(1,3,2))
      }
      tmp.dimnames<-dimnames(fit$param.diags)
      fit$param.diags<-array(NA,c(4,dim(fit$chains)[2],nchains),
                             c(tmp.dimnames[1],dimnames(fit$chains)[2],tmp.dimnames[3]))
      fit$param.diags[2,,]<-apply(fit$chains,c(2,3),rstan::ess_bulk)
      fit$param.diags[3,,]<-apply(fit$chains,c(2,3),rstan::ess_tail)
      fit$param.diags[4,,]<-apply(fit$chains,c(2,3),rstan::Rhat)
    }
    
    if(include.post.probs){
      if(!is.null(post.probs)){
        post.probs<-.expand.element(post.probs)
        if(any(dim(post.probs)!=c(1,e,nchains))){
          stop('mismatch between dimensions of provided post.prob array and expected dimensions: please double-check input')
        }
        fit$post.probs<-post.probs[1,,]
      }
    }
    
    #deal with chain stuff-->remember chain can be of length 2!
    if(is.null(chain)&nchains>1){
      fit<-combine.chains(fit)
    }else if(!is.null(chain)){
      fit<-combine.chains(select.chains(fit,chain))
    }
    niter<-dim(fit$chains)[1]
    
    if(!is.null(select.extra)){
      if(length(select.extra)!=1){
        stop('Can only plot one set of edgewise quantities at a time: please make select.extra is of length 1.')
      }
    }else{
      if(element=='chains'){
        select.extra<-sample(niter,1)
      }else if(element=='quantiles'){
        select.extra<-0.5
      }else if(element=='diagnostics'){
        select.extra<-'bulk_ess'
      }
    }
    
    try.R<-try(.int.chains(fit,1),silent=TRUE)
    if(inherits(try.R,'try-error')){
      select<-0
    }else{
      select<-1:e
    }
    if(!is.null(select.extra)){
      select<-list(select,select.extra)
    }
    
    R<-do.call(paste0('.int.',element),list(fit=fit,select=select))
    if(inherits(try.R,'try-error')){
      tmp.dims<-dim(R)
      tmp.dims[2]<-e
      tmp.dimnames<-dimnames(R)
      tmp.dimnames[[2]]<-paste('R',1:e,sep='_')
      R<-array(aperm(R,c(1,3,2)),tmp.dims[c(1,3,2)],tmp.dimnames[c(1,3,2)])
      R<-aperm(R,c(1,3,2))
    }
    
    if(is.null(fit$post.probs)){
      fit$post.probs<-compare.params(params1=remove.trend(fit),
                                     params2=get.bg.rate(fit))$post.probs
    }
    
    
    
    #handle trait stuff! Maybe make a bit more flexible for showing uncertainty?
    if(element=='chains'){
      try.anc.recon<-try(match.arg(anc.recon,c('summary','sample')),silent=T)
      if(inherits(try.element,'try-error')){
        stop(anc.recon," is not an available method for reconstructing ancestral states given a single iteration: please specify either 'summary' or 'sample")
      }
      anc.recon<-try.anc.recon
      if(anc.recon=='summary'){
        X<-get.post.traits(in.fit,select.extra=select.extra,stochastic=FALSE)
        X<-X[1:(length(X)/2)]
      }else{
        X<-get.post.traits(in.fit,select.extra=select.extra)
      }
    }else{
      X<-apply(get.post.traits(in.fit),2,mean)
    }
    X<-as.matrix(X)
    colnames(X)<-colnames(fit$call$trait.data)
    rownames(X)<-gsub(paste0('^',colnames(X),'_'),'',rownames(X))
    
    sim<-list('tree'=fit$call$tree,'R'=as.vector(R),
              'X'=X)
    if(include.post.probs){
      sim$post.probs<-as.vector(fit$post.probs)
    }
    class(sim)<-'evorates'
  }
  
  #need to improve legend handling...but it works for now
  plot.args<-list(...)
  if(is.null(plot.args$lwd)){
    plot.args$lwd<-6
  }
  if(is.null(plot.args$legend.args$location)){
    plot.args$legend.args$location<-'bottomleft'
  }
  if(is.null(plot.args$legend.args$side)){
    if(grepl('left',plot.args$legend.args$location)){
      plot.args$legend.args$side<-1
    }else{
      plot.args$legend.args$side<-2
    }
  }
  if(is.null(plot.args$legend.args$main.side)){
    plot.args$legend.args$main.side<-3
  }
  if(is.null(plot.args$legend.args$main.srt)){
    plot.args$legend.args$main.srt<-90
  }
  if(is.null(plot.args$legend.args$main.adj)){
    plot.args$legend.args$main.adj<-c(0,0.5)
  }
  legend.coords<-do.call(evorates:::plot.evorates,
                         c(sim=list(sim),
                           color.element='R',
                           colvec=NULL,
                           plot.args[!(names(plot.args)%in%c('sim','color.element','colvec'))]))
  
  
  
  if(include.post.probs&!is.null(sim$post.probs)){
    if(is.null(post.probs.args$col)){
      post.probs.args$col<-c('gray80','gray20')
    }
    if(!('breaks'%in%names(post.probs.args))){
      post.probs.args$breaks<-c(0.05,0.1,0.9,0.95)
    }
    if(is.null(post.probs.args$lwd)){
      post.probs.args$lwd<-plot.args$lwd/3
    }
    post.probs.args<-c(post.probs.args,
                       plot.args[!(names(plot.args)%in%c('legend.args',names(post.probs.args)))])
    if(is.null(post.probs.args$legend.args$box.offset)){
      post.probs.args$legend.args$box.offset<-rep(NA,2)
      if(grepl('left',plot.args$legend.args$location)){
        post.probs.args$legend.args$box.offset[1]<-max(legend.coords$x)-par()$usr[1]
      }else{
        post.probs.args$legend.args$box.offset[1]<-par()$usr[2]-min(legend.coords$x)
      }
      if(grepl('top',plot.args$legend.args$location)){
        post.probs.args$legend.args$box.offset[2]<-par()$usr[4]-max(legend.coords$y)
      }else{
        post.probs.args$legend.args$box.offset[2]<-min(legend.coords$y)-par()$usr[3]
      }
    }
    if(is.null(post.probs.args$legend.args$side)){
      if(grepl('left',plot.args$legend.args$location)){
        post.probs.args$legend.args$side<-2
      }else{
        post.probs.args$legend.args$side<-1
      }
    }
    if(is.null(post.probs.args$legend.args$main)){
      post.probs.args$legend.args$main<-'post. prob.'
    }
    if(is.null(post.probs.args$legend.args$main.side)){
      post.probs.args$legend.args$main.side<-3
    }
    if(is.null(post.probs.args$legend.args$main.srt)){
      post.probs.args$legend.args$main.srt<-90
    }
    if(is.null(post.probs.args$legend.args$main.adj)){
      post.probs.args$legend.args$main.adj<-c(0,0.5)
    }
    if(is.null(post.probs.args$legend.args$select.levels)&!is.null(post.probs.args$breaks)){
      break.len<-length(post.probs.args$breaks)
      if(break.len%%2==0){
        post.probs.args$legend.args$select.levels<-(1:(break.len+1))[-(break.len/2+1)]
      }else{
        post.probs.args$legend.args$select.levels<-(1:(break.len+1))[-(break.len/2+c(0.5,1.5))]
      }
    }
    if(is.null(post.probs.args$alpha)){
      if(is.null(post.probs.args$breaks)){
        post.probs.args$alpha<-c(1,0,1)
      }else{
        break.len<-length(post.probs.args$breaks)
        post.probs.args$alpha<-rep(0,break.len)
        post.probs.args$alpha[post.probs.args$legend.args$select.levels]<-1
      }
    }
    post.probs.args$legend.args<-c(post.probs.args$legend.args,
                                   plot.args$legend.args[!(names(plot.args$legend.args)%in%names(post.probs.args$legend.args))])
    do.call(evorates:::plot.evorates,
            c(sim=list(sim),
              color.element='post.probs',
              colvec=NULL,
              add=TRUE,
              post.probs.args[!(names(plot.args)%in%c('sim','color.element','colvec','add'))]))
  }
  invisible(sim)
}