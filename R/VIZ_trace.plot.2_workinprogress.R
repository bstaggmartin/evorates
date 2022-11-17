trace.plot<-function(x,p=0.05,n.windows=100,window.size=NULL,
                    lower.quant=NULL,upper.quant=NULL,
                    add.lines=NULL,
                    add=FALSE,make.legend=TRUE,...,
                    param.args=c('col','lines.col',
                                 'alpha','lines.alpha',
                                 'ribbon.col','ribbon.border',
                                 'ribbon.alpha','ribbon.border.alpha'),
                    chain.args=c('ribbon.angle','ribbon.density'),
                    inpar.args=c('lines.lty','lines.lwd')){
  
  ####getting arguments####
  
  .check.args(param.args,chain.args,inpar.args)
  .x<-.proc.x(x,list(...)[['fit']])
  if(ribbon&is.null(window.size)){
    window.size<-10*(dim(.x)[1]-1)/(n.windows-1)
  }
  nms<-dimnames(.x)
  param.names<-.get.param.names(nms[[2]],list(...)[['overwrite.param.names']])
  chain.names<-nms[[3]]
  n.params<-length(param.names)
  n.chains<-length(chain.names)
  n.tot<-n.params*n.chains
  n.inpar<-length(add.lines)
  i<-which(unlist(lapply(list(param.args,chain.args,inpar.args),
                         function(ii) 'angle'%in%ii)))
  if(!length(i)) i<-0
  foo<-function(n){
    if(length(n)){
      if(n==1){
        NULL
      }else{
        180/(n+1)*seq_len(n)
      }
    }
  }
  def.angle<-foo(switch(i,n.params,n.chains,n.inpar)) #replace with .get.def.angle function
  def.density<-if(!is.null(def.angle)) 20 else NULL #replace with .get.def.density function
  get.args<-.make.args.fxn(...,
                           param.args=param.args,
                           chain.args=chain.args,
                           inpar.args=inpar.args,
                           n.params=n.params,
                           n.chains=n.chains,
                           n.inpar=n.inpar,
                           inheritance=list(lines='gen',ribbon=c('gen','poly'),legend=NULL),
                           #make sure alpha arguments are properly inherited--did it with some edits to the default names
                           defaults=list(col=palette(),
                                         border=expression(col),
                                         lines.col=expression(border),
                                         alpha=1,
                                         border.alpha=expression(alpha),
                                         lines.alpha=expression(border.alpha),
                                         ribbon.col=expression(col),
                                         ribbon.alpha=expression(alpha/2), #good default I hope?
                                         angle=def.angle,
                                         density=def.density,
                                         lwd=4,
                                         ribbon.lwd=expression(lwd/2),
                                         lines.lwd=expression(ribbon.lwd),
                                         ribbon.lty=expression(lty+1),
                                         lines.lty=expression(ribbon.lty)))
  
  ####getting ribbons/lines####
  add.lines<-matrix(if(is.null(add.lines)) NA else add.lines,nrow=n.params*n.chains,ncol=n.inpar,byrow=TRUE)
  ribbons.and.lines<-.get.cuts(.x,p,
                               lower.quant,upper.quant,
                               lower.cut=NULL,upper.cut=NULL,
                               n.params,n.chains,
                               trace=TRUE,
                               add.lines=add.lines,
                               n.windows=n.windows,window.size=window.size)
  if(ribbon){
  }
  
  ####plotting window####
  
  if(!add){
    tmp.args<-get.args(~(plot|gen)&!recyc)
    do.call(.make.param.plot,
            c(x=list(.x),
              param.names=list(param.names),n.params=list(n.params),
              def.xlab=list(deparse(substitute(x))),
              tmp.args))
  }
  
  ####plotting profiles####
  
  .x<-asplit(.x,-1)
  for(i in seq_len(n.tot)){
    ii<-if(n.inpar) n.inpar*(i-1)+1 else i
    #maybe split .x for easier plotting?
    tmp.args<-get.args(~(poly|gen)&!recyc)
    tmp.recyc.args<-get.args(~(poly|gen)&recyc,ii)
    # do.call(polygon,
    #         )
    lines(.x[,(i-1)%%n.params+1,(i-1)%/%n.params+1])
    #need x coords for ribbon.quants...
    # lines(ribbon.quants[[1]][,1]~)
    prof<-tmp[[i]]
    full.xy<-.trim.prof(prof$x,prof$y,smooth,
                        range(.x[,(i-1)%%n.params+1,(i-1)%/%n.params+1]),
                        yran)
    part.xy<-.cut.prof(full.xy[[1]],cuts[[i]],smooth)
    tmp.args<-get.args(~(poly|gen)&nonborder&!recyc)
    tmp.recyc.args<-get.args(~(poly|gen)&nonborder&recyc,ii)
    do.call(polygon,
            c(x=part.xy[1],
              y=part.xy[2],
              border=NA,
              tmp.args,
              tmp.recyc.args))
    tmp.args<-get.args(~(poly|gen)&noncol&!recyc)
    tmp.recyc.args<-get.args(~(poly|gen)&noncol&recyc,ii)
    do.call(polygon,
            c(x=full.xy[[2]][1],
              y=full.xy[[2]][2],
              col=list(rgb(0,0,0,0)),
              tmp.args,
              tmp.recyc.args))
    ##NEED TO ADD SUPPORT FOR LINES OUT OF BOUNDARIES!
    if(n.inpar){
      ii<-seq.int(ii,ii+n.inpar-1)
      x0<-lines[[i]]
      y0<-0
      line.pt<-findInterval(x0,full.xy[[1]][[1]])
      y1<-.get.y(line.pt+1,full.xy[[1]][[1]],full.xy[[1]][[2]],x0,smooth)
      tmp.args<-get.args(~lines&!recyc)
      tmp.recyc.args<-get.args(~lines&recyc,ii)
      do.call(segments,
              c(x0=list(x0),
                y0=list(y0),
                y1=list(y1),
                tmp.args,
                tmp.recyc.args))
    }
  }
}