#disabled adjust for plotting for now due its redudancy with calls to density and generally not being too helpful graphically...
#actually ADJUST is for density(),so hopefully argument sorting features will prevent it's redundancy...
#list of potential parameters...
# names(formals(matplot))
# names(formals(hist.default))
# names(formals(plot.default))
# graphics:::.Pars
# names(formals(density.default))
#check for names that match with these parameters should prevent most warnings and prevent you form having to call
#suppressWarnings()
deprecated.profile.plot<-function(fit,select='.',separate.R=F,separate.X=F,separate.dev=F,together=F,dens=F,
                       alpha=NA,exp=F,sqrt=F,...){
  plot.args<-c(names(formals(plot.default)),
               names(formals(axis)),names(formals(box)),names(formals(plot.window)),names(formals(title)))
  plot.args<-plot.args[-which(plot.args=='...')]
  gen.args<-graphics:::.Pars
  dens.args<-names(formals(density.default))
  dens.args<-dens.args[-which(dens.args=='...')]
  hist.args<-names(formals(hist.default))
  hist.args<-hist.args[-which(hist.args=='...')]
  poly.args<-names(formals(polygon))
  poly.args<-poly.args[-which(poly.args=='...')]
  chains<-chains.to.array(fit%chains%select)
  args.master<-do.call(get.args.master,
                       c(chains=list(chains),
                         separate.R=separate.R,
                         separate.X=separate.X,
                         separate.dev=separate.dev,
                         together=together,
                         alpha=alpha,
                         sqrt=sqrt,
                         exp=exp,
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
  
  #dealing with breaks/bandwidths for 'together' plots
  if(together){
    if(dens){
      def.bw<-bw.nrd0(chains)
      for(i in 1:length(args.master)){
        if(!('bw'%in%names(args.master[[i]]))){
          args.master[[i]]$bw<-def.bw
        }
      }
    }else{
      def.breaks<-pretty(chains,nclass.Sturges(chains))
      for(i in 1:length(args.master)){
        if(!('breaks'%in%names(args.master[[i]]))){
          args.master[[i]]$breaks<-def.breaks
        }else if(length(args.master[[i]]$breaks)==1&is.numeric(args.master[[i]]$breaks)){
          args.master[[i]]$breaks<-pretty(chains,args.master[[i]]$breaks)
        }
      }
    }
  }
  #loop to store the graph info for each histogram/kernel density function...
  graphs<-rep(list(vector(mode='list',length=dim(chains)[3])),dim(chains)[2])
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
    tmp.chains<-chains[,incl.names,]
    if(length(incl.names)==1){
      attr(tmp.chains,'parameters')<-param
    }
    tmp.chains<-chains.to.array(tmp.chains)
    if(!together){
      if(dens){
        if(!('bw'%in%names(args.master[[i]]))){
          args.master[[i]]$bw<-bw.nrd0(tmp.chains)
        }
      }else{
        if(!('breaks'%in%names(args.master[[i]]))){
          args.master[[i]]$breaks<-pretty(tmp.chains,nclass.Sturges(tmp.chains))
        }else if(length(args.master[[i]]$breaks)==1&is.numeric(args.master[[i]]$breaks)){
          args.master[[i]]$breaks<-pretty(tmp.chains,args.master[[i]]$breaks)
        }
      }
    }
    for(j in 1:dim(tmp.chains)[2]){
      for(k in 1:dim(chains)[3]){
        if(dens){
          graphs[[incl.names[j]]][[k]]<-
            do.call(density,
                    c(x=list(tmp.chains[,j,k]),
                      args.master[[i]][names(args.master[[i]])!='x'&
                                         names(args.master[[i]])%in%dens.args]))
        }else{
          graphs[[incl.names[j]]][[k]]<-
            do.call(hist,
                    c(x=list(tmp.chains[,j,k]),
                      warn.unused=F,
                      plot=F,
                      args.master[[i]][!(names(args.master[[i]])%in%c('x','warn.unused','plot'))&
                                         names(args.master[[i]])%in%hist.args]))
        }
      }
    }
  }
  if(together){
    if(is.null(list(...)$main)){
      main<-''
    }else{
      main<-list(...)$meain
    }
    if(is.null(list(...)$xlab)){
      xlab<-'parameters'
    }else{
      xlab<-list(...)$xlab
    }
    if(is.null(list(...)$ylab)){
      ylab<-'density'
    }else{
      ylab<-list(...)$ylab
    }
    if(is.null(list(...)$xlim)){
      if(dens){
        xlim<-range(unlist(sapply(graphs,function(ii) lapply(ii,function(jj) jj$x))))
      }else{
        xlim<-range(unlist(sapply(graphs,function(ii) lapply(ii,function(jj) jj$breaks))))
      }
    }else{
      xlim<-list(...)$xlim
    }
    if(is.null(list(...)$ylim)){
      if(dens){
        ylim<-c(0,max(unlist(sapply(graphs,function(ii) lapply(ii,function(jj) jj$y)))))
      }else{
        ylim<-c(0,max(unlist(sapply(graphs,function(ii) lapply(ii,function(jj) jj$density)))))
      }
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
                   list(...)[!(names(list(...))%in%c('x','y','xlab','ylab','xlim','ylim','type','pch','col'))&
                               (names(list(...))%in%gen.args|names(list(...))%in%plot.args)]))
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
        args.master[[i]]$xlab<-param
      }
      if(is.null(args.master[[i]]$ylab)){
        args.master[[i]]$ylab<-'density'
      }
      if(is.null(args.master[[i]]$xlim)){
        if(dens){
          args.master[[i]]$xlim<-range(unlist(sapply(graphs[incl.names],function(ii) lapply(ii,function(jj) jj$x))))
        }else{
          args.master[[i]]$xlim<-range(unlist(sapply(graphs[incl.names],function(ii) lapply(ii,function(jj) jj$breaks))))
        }
      }
      if(is.null(args.master[[i]]$ylim)){
        if(dens){
          args.master[[i]]$ylim<-c(0,max(unlist(sapply(graphs[incl.names],function(ii) lapply(ii,function(jj) jj$y)))))
        }else{
          args.master[[i]]$ylim<-c(0,max(unlist(sapply(graphs[incl.names],function(ii) lapply(ii,function(jj) jj$density)))))
        }
      }
      do.call(plot,c(x=0,
                     type='p',
                     pch=1,
                     col='white',
                     args.master[[i]][!(names(args.master[[i]])%in%c('x','y','type','pch','col'))&
                                        (names(args.master[[i]])%in%gen.args|names(args.master[[i]])%in%plot.args)]))
    }
    if(is.null(args.master[[i]]$density)){
      args.master[[i]]$density<-rep(c(-1,rep(20,4)),length.out=dim(chains)[3])
    }else{
      args.master[[i]]$density<-rep(args.master[[i]]$density,length.out=dim(chains)[3])
    }
    if(is.null(args.master[[i]]$angle)){
      args.master[[i]]$angle<-rep(c(0,36,72,108,144),length.out=dim(chains)[3])
    }else{
      args.master[[i]]$angle<-rep(args.master[[i]]$angle,length.out=dim(chains)[3])
    }
    if(is.null(args.master[[i]]$col)){
      args.master[[i]]$col<-rep(alter.cols(palette(),alph=args.master[[i]]$alpha),length.out=dim(tmp.chains)[2])
    }else{
      args.master[[i]]$col<-rep(alter.cols(args.master[[i]]$col,alph=args.master[[i]]$alpha),length.out=dim(tmp.chains)[2])
    }
    for(j in 1:dim(tmp.chains)[2]){
      for(k in 1:dim(chains)[3]){
        if(dens){
          do.call(polygon,c(x=list(graphs[[incl.names[j]]][[k]]$x),
                            y=list(graphs[[incl.names[j]]][[k]]$y),
                            density=args.master[[i]]$density[[k]],
                            angle=args.master[[i]]$angle[[k]],
                            col=args.master[[i]]$col[[j]],
                            args.master[[i]][!(names(args.master[[i]])%in%c('x','y','density','angle','col','add'))&
                                               (names(args.master[[i]])%in%gen.args|names(args.master[[i]])%in%poly.args)]))
        }else{
          do.call(plot,c(x=list(graphs[[incl.names[j]]][[k]]),
                         freq=F,
                         density=args.master[[i]]$density[[k]],
                         angle=args.master[[i]]$angle[[k]],
                         col=args.master[[i]]$col[[j]],
                         add=T,
                         args.master[[i]][!(names(args.master[[i]])%in%c('x','y','freq','density','angle','col','add'))&
                                            (names(args.master[[i]])%in%gen.args|names(args.master[[i]])%in%plot.args)]))
        }
      }
    }
  }
}
