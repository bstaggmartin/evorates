#with new drop.tip functionality, might be nice to force using posterior probabilities in fit
#maybe just set to NA to use post probs in original fit?
#' @export
plot.evorates_fit<-function(fit,chain=NULL,
                            type=c('quantiles','chains','means','diagnostics'),extra.select=NULL,
                            post.probs=c("stored","recalculate","none"),devs=c("none","stored","recalculate"),
                            remove.trend=TRUE,geometric=TRUE,
                            recon.type=NULL,recon.extra.select=NULL,
                            ...,
                            post.probs.args=NULL,
                            sim=NULL){
  #for backwards compatibility
  if(is.logical(post.probs)){
    if(post.probs){
      post.probs<-"recalculate"
    }else{
      post.probs<-"none"
    }
  }
  #may want smarter default that checks for whether devs are stored in fit...
  post.probs<-.match.arg(post.probs,c("stored","recalculate","none"),"post.probs")
  #for backwards compatibility
  if(hasArg(plot.Rdev)){
    if(plot.Rdev){
      devs<-"recalculate"
    }else{
      devs<-"none"
    }
  }
  devs<-.match.arg(devs,c("none","stored","recalculate"),"devs")
  
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
    if(post.probs!="none"|devs!="none"){
      #check to see if stored even works here...
      if(post.probs=="stored"|devs=="stored"){
        if(!any(grepl('^Rdev_[1-9][0-9]*$',names(fit$chains)))){
          if(post.probs=="stored"){
            # warning() #warn that no stored rate deviations found
            post.probs<-"recalculate"
          }
          if(devs=="stored"){
            # warning() #warn that no stored rate deviations found
            devs<-"recalculate"
          }
        }
      }
      #recalculate devs if post.probs or devs need to be recalculated...
      if(post.probs=="recalculate"|devs=="recalculate"){
        R<-get.bg.rate(fit,simplify=FALSE,keep.R=TRUE,
                       remove.trend=remove.trend,geometric=geometric)
        R<-R$R-R$bg_rate
      }
      #recalculate or get post.probs as necessary...
      if(post.probs=="recalculate"){
        tmp<-R
        tmp[tmp==0]<-NA
        pp<-.call.op('means',list(chains=tmp>0,sampler.params=1),'.',FALSE)
        pp[is.na(pp)]<-0.5
      }else if(post.probs=="stored"){
        pp<-fit$post.probs
      }
    }
    R<-switch(devs,
              "none"=get.R(fit,type=type,extra.select=extra.select,simplify=FALSE),
              "stored"=.call.op(type,fit,select=list('^Rdev_[1-9][0-9]*$',extra.select),FALSE),
              "recalculate"=.call.op(type,list(chains=R,sampler.params=1),list('.',extra.select),FALSE))
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
    X<-get.post.traits(fit,type=recon.type,extra.select=recon.extra.select)
    X<-as.matrix(X)
    colnames(X)<-colnames(fit$call$trait.data)
    rownames(X)<-gsub(paste0('^',colnames(X),'_'),'',rownames(X))
    #final coercion
    sim<-list('tree'=fit$call$tree,'R'=as.vector(R),
              'X'=X)
    if(post.probs!="none"){
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
  
  if(post.probs!='none'&!is.null(sim$post.probs)){
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