#make it match args.integer with args.R[integer] -- check!
#make an option to graph everything in the same plot -- check!
#make it so that args.R[integer] options are passed to appropriate spot in args.R -- half check, easy work around by setting together to T, separate.R to F, and selecting the appropriate parameters
  #still might be nice for there to be 'cross-talk'; i.e., if args.R is around but separate.R = T, pass args.R to all rate params, if args.R[1] is around but separate.R = F, pass args.R[1] to 1st spot in vectors of parameters in args.R...
  #would only involve manipulating get.args.master function, I believe
#6/26: trace and profile plot are working well--I think the only thing left to do is add support for legends and POTENTIALLY
#parameter transformations (though this might involve altering the %% operators themselves...)
#' @export
trace.plot<-function(fit,select='.',separate.R=F,separate.X=F,separate.dev=F,together=F,
                     alpha=NA,exp=F,sqrt=F,legend=F,...){
  plot.args<-c(names(formals(plot.default)),
               names(formals(axis)),names(formals(box)),names(formals(plot.window)),names(formals(title)))
  plot.args<-plot.args[-which(plot.args=='...')]
  gen.args<-graphics:::.Pars
  matplot.args<-names(formals(matplot))
  matplot.args<-matplot.args[-which(matplot.args=='...')]
  chains<-chains.to.array(fit%chains%select)
  args.master<-do.call(get.args.master,
                       c(chains=list(chains),
                         separate.R=separate.R,
                         separate.X=separate.X,
                         separate.dev=separate.dev,
                         together=together,
                         alpha=alpha,
                         exp=exp,
                         sqrt=sqrt,
                         legend=legend,
                         list(...)))
  
  for(i in names(args.master)){
    param<-paste(strsplit(i,'\\.')[[1]][-1],collapse='.')
    incl.names<-numeric(0)
    if(param=='R'){
      incl.names<-grep('^R\\[\\d+\\]$',dimnames(chains)[[2]])
    }else if(param=='X'){
      incl.names<-grep('^R0$|^Rsig2$|^X0$|^bg.rate$|^R\\[\\d+\\]$|^R\\[\\d+\\] dev$|^Rmu$|^Ysig2$',
                       dimnames(chains)[[2]],invert=T)
    }else if(param=='dev'){
      incl.names<-grep('^R\\[\\d+\\] dev$',dimnames(chains)[[2]])
    }else{
      grep.param<-sub('\\[','\\\\\\[',param)
      grep.param<-sub('\\]','\\\\\\]',grep.param)
      incl.names<-grep(paste('^',grep.param,'$',sep=''),dimnames(chains)[[2]])
    }
    if(length(incl.names)==0){
      next
    }
    if(args.master[[i]]$exp){
      chains[,incl.names,]<-exp(chains[,incl.names,])
    }
    if(args.master[[i]]$sqrt){
      if(!all(chains[,incl.names,]>=0)){
        warning(paste('parameters matching with ',substr(i,6,nchar(i)),' include negative numbers: sqrt argument in ',i,' set to FALSE',sep=''))
        break
      }
      chains[,incl.names,]<-sqrt(chains[,incl.names,])
    }
  }
  
  if(together){
    if(is.null(list(...)$xlab)){
      xlab<-'iterations'
    }else{
      xlab<-list(...)$xlab
    }
    if(is.null(list(...)$ylab)){
      ylab<-'parameters'
    }else{
      ylab<-list(...)$ylab
    }
    if(is.null(list(...)$xlim)){
      xlim<-c(1,dim(chains)[1])
    }else{
      xlim<-list(...)$xlim
    }
    if(is.null(list(...)$ylim)){
      ylim<-range(chains,na.rm=T)
    }else{
      ylim<-list(...)$ylim
    }
    do.call(plot,c(x=0,
                   type='p',
                   pch=1,
                   col='white',
                   xlab=list(xlab),
                   ylab=list(ylab),
                   xlim=list(xlim),
                   ylim=list(ylim),
                   list(...)[!(names(list(...))%in%c('x','type','pch','col','xlab','ylab','xlim','ylim'))&
                               (names(list(...))%in%plot.args|names(list(...))%in%gen.args)]))
  }
  for(i in names(args.master)){
    param<-paste(strsplit(i,'\\.')[[1]][-1],collapse='.')
    if(param=='R'){
      incl.names<-grep('^R\\[\\d+\\]$',dimnames(chains)[[2]])
    }else if(param=='X'){
      incl.names<-grep('^R0$|^Rsig2$|^X0$|^bg.rate$|^R\\[\\d+\\]$|^R\\[\\d+\\] dev$|^Rmu$|^Ysig2$',
                       dimnames(chains)[[2]],invert=T)
    }else if(param=='dev'){
      incl.names<-grep('^R\\[\\d+\\] dev$',dimnames(chains)[[2]])
    }else{
      grep.param<-sub('\\[','\\\\\\[',param)
      grep.param<-sub('\\]','\\\\\\]',grep.param)
      incl.names<-grep(paste('^',grep.param,'$',sep=''),dimnames(chains)[[2]])
    }
    if(length(incl.names)==0){
      next
    }
    tmp.chains<-chains[,incl.names,]
    if(length(incl.names)==1){
      attr(tmp.chains,'parameters')<-param
    }
    tmp.chains<-chains.to.array(tmp.chains)
    if(!together){
      if(is.null(args.master[[i]]$xlab)){
        args.master[[i]]$xlab<-'iterations'
      }
      if(is.null(args.master[[i]]$ylab)){
        args.master[[i]]$ylab<-param
      }
      if(is.null(args.master[[i]]$xlim)){
        args.master[[i]]$xlim<-c(1,dim(chains)[1])
      }
      if(is.null(args.master[[i]]$ylim)){
        args.master[[i]]$ylim<-range(tmp.chains,na.rm=T)
      }
      do.call(plot,c(x=0,
                     type='p',
                     pch=1,
                     col='white',
                     args.master[[i]][!(names(args.master[[i]])%in%c('x','type','pch','col'))&
                                        (names(args.master[[i]])%in%plot.args|names(args.master[[i]])%in%gen.args)]))
    }
    if(is.null(args.master[[i]]$lty)){
      args.master[[i]]$lty<-rep(1:6,length.out=dim(chains)[3])
    }else{
      args.master[[i]]$lty<-rep(args.master[[i]]$lty,length.out=dim(chains)[3])
    }
    if(is.null(args.master[[i]]$col)){
      args.master[[i]]$col<-rep(alter.cols(palette(),alph=args.master[[i]]$alpha),
                                length.out=dim(tmp.chains)[2])
    }else{
      args.master[[i]]$col<-rep(alter.cols(args.master[[i]]$col,alph=args.master[[i]]$alpha),
                                length.out=dim(tmp.chains)[2])
    }
    for(j in 1:dim(chains)[3]){
      do.call(matplot,c(y=list(tmp.chains[,,j]),
                        type='l',
                        lty=args.master[[i]]$lty[j],
                        col=list(args.master[[i]]$col),
                        add=T,
                        args.master[[i]][!(names(args.master[[i]])%in%c('x','y','type','lty','col','add'))&
                                           (names(args.master[[i]])%in%matplot.args|
                                              names(args.master[[i]])%in%plot.args|
                                              names(args.master[[i]])%in%gen.args)]))
    }
  }
}
