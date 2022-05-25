#notes

##p, quant, and cut can all be vectors
###all are recycled to proper length
###non-NA cut entries will override quant entries, and non-NA quant entries will override p entries in turn

##"inheritance" system--makes the whole thing a bit more cumbersome, but more flexible in the long run
###woot! All seems to work. Only hiccup is you might want to "hardwire" the param.args, chain.args, and inpar.args in
###it could be really easy to break things otherwise! but maybe keep it in for advanced users...

##add.lines corresponds to x-values if quant.lines is FALSE, and quantiles otherwise
###specifying NA for add.lines defaults to mean, regardless of quant.lines


#(low priority) future things:

##make an argument for where to cap off legend/plot title for a given number of parameters

##extra argument to combine chains, if desired?

##it might be good to include lines in legend, but I don't see a simple way to do this if different numbers of lines are allowed per
###profile

##allow lines to vary across parameters, like cuts? but how to do this efficiently? matrices/lists?
###maybe lists to allow varying numbers of lines... (see below)--will complicate recycling and indexing, but not by much
####marking divergent transitions using add.lines feature?
###or not--will be same number of divergences regardless of parameter! Keeping same number of lines per profile will be much easier

##on that note, allow quant.lines to be a vector recycled so lines corresponding to quantiles and cuts can be efficiently added?
###wait, same number per parameter, but NOT per chain, so it will have to allow different numbers...

#' Plot profiles of posterior samples from fitted 
#'
#'
#' This function processes tree and trait data, runs a Stan-based Hamiltonian Monte Carlo (HMC) sampler to fit
#' these data to an evorates model, and returns the output of this sampler in a (relatively) user-friendly format.
#'
#'
#' @param tree An object of class "\code{phylo}"
#' 
#' 
#' @return 
#' 
#' 
#' @details 
#' PARAMETER DEFINITIONS HERE
#' 
#' 
#' @family evorates plotting functions
#' 
#' 
#' @examples
#' 
#' 
#' @export
prof.plot<-function(x,smooth=FALSE,p=0.05,
                    lower.quant=NULL,upper.quant=NULL,
                    lower.cut=NULL,upper.cut=NULL,
                    add.lines=NULL,quant.lines=TRUE,
                    add=FALSE,make.legend=TRUE,...,
                    param.args=c('col','border','lines.col',
                                 'alpha','border.alpha','lines.alpha'),
                    chain.args=c('angle','density'),
                    inpar.args=c('lines.lty','lines.lwd')){
  
  ####getting arguments####
  
  .check.args(param.args,chain.args,inpar.args)
  .x<-.proc.x(x,list(...)[['fit']])
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
  def.angle<-foo(switch(i,n.params,n.chains,n.inpar))
  def.density<-if(!is.null(def.angle)) 20 else NULL
  get.args<-.make.args.fxn(...,
                           param.args=param.args,
                           chain.args=chain.args,
                           inpar.args=inpar.args,
                           n.params=n.params,
                           n.chains=n.chains,
                           n.inpar=n.inpar,
                           inheritance=list(lines='gen',legend=NULL),
                           defaults=list(col=palette(),
                                         border=expression(col),
                                         lines.col=expression(border),
                                         alpha=NA,
                                         border.alpha=NA,
                                         lines.alpha=expression(border.alpha),
                                         angle=def.angle,
                                         density=def.density))
  
  ####getting profiles####
  
  if(smooth){
    foo<-density.default
    nm<-'bw'
    tmp.args<-get.args(~dens&!recyc)
  }else{
    foo<-hist
    nm<-'breaks'
    tmp.args<-c(plot=FALSE,
                warn.unused=FALSE,
                get.args(~hist&!recyc))
  }
  tmp<-do.call(foo,c(x=list(.x),tmp.args))[[nm]]
  tmp.args[[nm]]<-tmp
  tmp<-apply(.x,c(2,3),function(ii)
    do.call(foo,c(x=list(ii),tmp.args)))
  if(!smooth){
    for(i in seq_len(n.tot)){
      names(tmp[[i]])<-gsub('counts','y',gsub('breaks','x',names(tmp[[i]])))
    }
  }
  
  ####plotting window####
  
  if(!add){
    tmp.args<-get.args(~(plot|gen)&!recyc)
    do.call(.make.param.plot,
            c(x=list(tmp),
              param.names=list(param.names),n.params=list(n.params),
              def.xlab=list(deparse(substitute(x))),
              smooth=list(smooth),
              tmp.args))
  }
  
  ####getting cuts/lines####
  
  lines<-.get.lines(.x,add.lines,quant.lines,n.tot) #not "vectorized" like cuts
  cuts<-.get.cuts(.x,p,
                  lower.quant,upper.quant,
                  lower.cut,upper.cut,
                  n.params,n.chains)
  
  ####plotting profiles####
  
  yran<-diff(par()$usr[c(3,4)])
  for(i in seq_len(n.tot)){
    ii<-if(n.inpar) n.inpar*(i-1)+1 else i
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
    ##NEED TO ADD SUPPORT FOR LINES OUT OF BOUNDARIES!--done by modifying .get.y, I believe
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
  
  ####plotting legend####
  
  if(make.legend){
    .make.param.plot.legend(param.names,chain.names,
                            n.params,n.chains,
                            get.args)
  }
}
