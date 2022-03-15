#later add support for plotting multiple samples from chains
#okay, can use select.extra to specify a particular iteration from chains too
#use chain to specify a specific chain using select.chains()
#if NULL, combine chains (with warning if trying to plot inits)
#god this function, could be SOOOO much nicer and cleaner-->I have to work on it more
#in particular, think about accomodating non-default rate and post.probs-->I think there's a better system!
#totally broken right now, but will fix soon
#' @export
plot.evorates_fit<-function(fit,chain=NULL,
                            type=c('quantiles','chains','means','diagnostics'),extra.select=NULL,
                            post.probs=TRUE,remove.trend=TRUE,geometric=TRUE,plot.Rdev=FALSE,
                            recon.type=NULL,recon.extra.select=NULL,
                            ...,
                            post.probs.args=NULL,
                            sim=NULL){
  if(is.null(sim)){ #if no evorates object is supplied, coerce to evorates_fit to evorates
    #check evorates_fit and combine chains, if necessary
    if(!inherits(fit,'evorates_fit')){
      stop("Fit must be an evorates_fit object")
    }
    nchains<-fit$sampler.control$chains
    niter<-dim(fit$chains)[1]
    e<-Nedge(fit)
    if(is.null(chain)&nchains>1){
      fit<-combine.chains(fit)
    }else if(!is.null(chain)){
      fit<-combine.chains(select.chains(fit,chain))
    }
    niter<-dim(fit$chains)[1]
    
    #process type stuff
    type<-.match.type(type,choices=c('quantiles','chains','means','diagnostics'))
    if(!is.null(extra.select)){
      extra.select<-unlist(extra.select,use.names=FALSE)[1]
    }else{
      extra.select<-switch(type,
                           chains=sample(niter,1),
                           quantiles=0.5,
                           diagnostics='bulk_ess')
    }
    
    #process post.prob stuff and get branchwise rates
    if(post.probs|plot.Rdev){
      R<-get.bg.rate(fit,simplify=FALSE,keep.R=TRUE,
                     remove.trend=remove.trend,geometric=geometric)
      R<-R$R-R$bg_rate
      if(post.probs){
        pp<-.call.op('means',list(chains=R>0,sampler.params=1),'.',FALSE)
      }
    }
    if(!plot.Rdev){
      R<-get.R(fit,type=type,extra.select=extra.select,simplify=FALSE)
    }else{
      R<-.call.op(type,list(chains=R,sampler.params=1),list('.',extra.select),FALSE)
    }
    
    #trait stuff
    if(is.null(recon.type)){
      if(type=='chains'){
        recon.type<-'chains'
      }else{
        recon.type<-'quantiles'
      }
    }
    recon.type<-.match.type(recon.type,c('quantiles','chains','means'))
    if(is.null(recon.extra.select)){
      recon.extra.select<-switch(recon.type,
                                 chains=if(type=='chains') extra.select else sample(niter,1),
                                 quantiles=0.5)
    }
    X<-get.post.traits(fit,type=recon.type,extra.select=recon.extra.select,trait.name='')
    X<-as.matrix(X)
    colnames(X)<-colnames(fit$call$trait.data)
    rownames(X)<-gsub('^_','',rownames(X))
    
    #final coercion
    sim<-list('tree'=fit$call$tree,'R'=as.vector(R),
              'X'=X)
    if(post.probs){
      sim$post.probs<-as.vector(pp)
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
  
  if(post.probs&!is.null(sim$post.probs)){
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