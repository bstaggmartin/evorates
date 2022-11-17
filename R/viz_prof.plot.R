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

#' Plot profiles of posterior samples from a fitted evorates model
#'
#'
#' This function plots histograms/densities ("profiles") of posterior samples for particular parameters from an \code{evorates_fit}
#' object or \code{param_block} array.
#'
#'
#' @param x The parameters to be plotted. Typically, this is a \code{param_block} array extracted from an \code{evorates_fit} object,
#' though it can also be a character or numeric vector specifying parameters to extract from an \code{evorates_fit} object passed to
#' the \code{fit} argument (see \code{...}). For details on how parameters are extracted, see \link{grapes-chains-grapes}. Multiple
#' \code{param_block} arrays and character/numeric selections can be combined as lists.
#' @param smooth \code{TRUE}or \code{FALSE}: should profiles be plotted as histograms or "smoothed" density plots? By default, set to
#' \code{FALSE} such that profiles are plotted as histograms. Densities are estimated via R's built-in \code{density()} function.
#' @param p A numeric vector controlling the widths of the shaded portions of profiles for each parameter. Specifically, the shaded
#' region is the \code{(1-p)}\% credible interval (so 95\% credible interval by default). The vector is recycled as necessary,
#' \code{NA}'s result in shading in the
#' entire profile (technically the 100\% credible interval), and entries below 0 or above 1 are rounded up and down to 0 and 1,
#' respectively. Overwritten by \code{lower.quant}, \code{upper.quant}, \code{lower.cut}, and \code{upper.cut}, which may be used to
#' provide finer control over the bounds of the shaded regions.
#' @param lower.quant,upper.quant Numeric vectors specifying the boundaries of the shaded portions of profiles for each
#' parameter based on posterior quantiles. Set to \code{NULL} (the default) to just use all boundaries specified by \code{p} (see above).
#' The vectors are recycled as necessary, \code{NA}'s also default to boundaries specified by \code{p}, and entries
#' below 0 or above 1 are rounded up and down to 0 and 1, respectively. Overwritten by \code{lower.cut} and \code{upper.cut}, which may
#' be used to provide finer control over the bounds of the shaded regions.
#' @param lower.cut,upper.cut Numeric vectors specifying the boundaries of the shaded portions of profiles for each
#' parameter (on the x-axis scale, rather than quantile-based). Notably, these can be used plot credible intervals based on the
#' highest posterior density regions, rather than quantiles as is done by default (see examples below). Set to \code{NULL}
#' (the default) to just use all boundaries specified by \code{p}, \code{lower.quant}, and/or \code{upper.quant} (see above).
#' The vectors are recycled as necessary and \code{NA}'s also default to boundaries specified by \code{p}, \code{lower.quant},
#' and/or \code{upper.quant}.
#' @param add.lines A numeric vector specifying where to plot lines over profiles. This is not yet "vectorized", and will plot the same lines
#' for all profiles (though I hope to change this in the future)! Can be specified via posterior quantiles if \code{quant.lines} is
#' \code{TRUE}, or via actual positions along x-axis otherwise. \code{NA}'s result in plotting lines at posterior means.
#' @param quant.lines \code{TRUE}or \code{FALSE}: should line positions (see above) correspond to posterior quantiles or positions along
#' the x-axis? Defaults to \code{TRUE}, meaning that \code{add.lines} corresponds to posterior quantiles; in this case entries of
#' \code{add.lines} below 0 or above 1 are rounded up and down to 0 and 1, respectively. This does not affect \code{NA} entries of
#' \code{add.lines}, which always result in plotting lines at posterior means.
#' @param add \code{TRUE}or \code{FALSE}: should profiles be plotted in an open plotting window, or should a new plot window be opened?
#' Defaults to \code{FALSE}, meaning that a new plot window is started.
#' @param make.legend \code{TRUE}or \code{FALSE}: should a legend be automatically generated for the plotted profiles? Defaults to
#' \code{TRUE}, but automatically suppressed when only a single profile is plotted.
#' @param ... Other arguments, such as graphical parameters or arguments to pass to \code{hist()}/\code{density()} for profile
#' construction. Here are the some commonly-used ones:
#' \itemize{
#' \item \code{fit} to specify an \code{evorates_fit} object from which to extract parameters (if \code{x} includes character or numeric
#' selections).
#' \item \code{overwrite.param.names} to specify alternate names for parameters in the automatically generated x-axis title and legend.
#' This is helpful if you want to include otherwise forbidden characters in parameter names, like greek symbols.
#' \item \code{breaks} to specify the number of breaks in histograms if \code{smooth} is \code{FALSE}.
#' \item \code{bw} and/or \code{adjust} to specify the width of kernels used to estimate density if \code{smooth} is \code{TRUE}.
#' \item Most base R graphical parameters should work, along with a few extras:
#' \itemize{
#' \item{\code{col} by default specifies the color of the profiles for each parameter. Defaults to \code{palette()}.}
#' \item{\code{border} by default specifies the color of the profile borders for each parameter. Inherits from \code{col} if not
#'       specified.}
#' \item{\code{alpha} by default controls the transparency of profile colors for each parameter, with 0, 1, and \code{NA} corresponding to
#'       completely transparent, opaque, and no modification, respectively. Defaults to \code{NA}.}
#' \item{\code{lwd},\code{lty} by default specify the line width and style for all profiles (you could mess with \code{param.args},
#'       etc. to change this--see below).}
#' \item{\code{angle} by default controls the angle of the lines used to shade profiles for each chain if \code{density} is a positive number.
#'       Defaults to a range of angles between 0 and 180 if profiles for multiple chains are plotted.}
#' \item{\code{density} by default controls the density of lines used to shade profiles for each chain, with \code{NULL} and \code{NA}
#'       resulting in solid shading (the default if profiles are plotted for a single chain). Defaults to 20 if profiles for multiple
#'       chains are plotted.}
#' \item{Some arguments like \code{col}, \code{alpha}, \code{lwd}, etc. can be modified with "\code{lines.}" or "\code{border.}" at the
#'       beginning to make them only affect certain elements of the plot. Perhaps most importantly, "\code{legend.}" can be used to 
#'       make arguments affect the automatically generated legend if \code{make.legend} is \code{TRUE}.}
#' }}
#' @param param.args,chain.args,inpar.args Generally, these should not be altered, but are exposed here for the curious to play around
#' with. These are vectors of argument names for graphical parameters that control whether the graphical parameters vary by the
#' parameter the profile corresponds to (\code{param.args}), the chain the profile corresponds to (\code{chain.args}), or within
#' profiles (\code{inpar.args}). Notably, default graphical parameters are tailored to the defaults for these arguments, and messing
#' with these will often result in ugly plots without a lot of extra customization!
#' 
#' 
#' @return Nothing as of now--this function just makes a plot!
#' 
#' 
#' @family evorates plotting functions
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #plot some parameters
#' prof.plot(cet_fit %chains% c("R_0", "R_mu"), alpha = 0.5)
#' #the above is equivalent to:
#' prof.plot(c("R_0", "R_mu"), alpha = 0.5, fit = cet_fit)
#' #could also something like this:
#' par <- get.bg.rate(fit = cet_fit,
#'                    node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                           c('Mesoplodon','Orca','Globicephalinae')),
#'                    )
#' prof.plot(par, alpha = 0.5)
#' #or even:
#' prof.plot(list(par, "R_0"), alpha = 0.5, fit = cet_fit)
#' 
#' #some ways the plot style might be tweaked:
#' prof.plot(list(par, "R_0"), fit = cet_fit,
#'           overwrite.param.names = expression("ln"~sigma["Meso"]^2,
#'                                              "ln"~sigma["Orca"]^2,
#'                                              "ln"~sigma["Glob"]^2,
#'                                              "ln"~sigma["Root"]^2),
#'           breaks = 100,
#'           col = c('blue','red','green','black'), alpha = 0.1, border.alpha = 0.25,
#'           add.lines = c(0.25, NA, 0.75), lines.lwd = c(2, 3, 2), lines.lty = c(3, 2, 3),
#'           bty='n', legend.bty = 'n', legend.inset = c(0.1, 0), lwd = 2)
#' #may want to "smooth out" profiles!
#' prof.plot(list(par, "R_0"), fit = cet_fit,
#'           overwrite.param.names = expression("ln"~sigma["Meso"]^2,
#'                                              "ln"~sigma["Orca"]^2,
#'                                              "ln"~sigma["Glob"]^2,
#'                                              "ln"~sigma["Root"]^2),
#'           smooth = TRUE, bw = 0.1,
#'           col = c('blue','red','green','black'), alpha = 0.1, border.alpha = 0.25,
#'           add.lines = c(0.25, NA, 0.75), lines.lwd = c(2, 3, 2), lines.lty = c(3, 2, 3),
#'           bty='n', legend.bty = 'n', legend.inset = c(0.1, 0), lwd = 2)
#' #try some narrower credible intervals?
#' prof.plot(list(par, "R_0"), fit = cet_fit, p = 0.5,
#'           overwrite.param.names = expression("ln"~sigma["Meso"]^2,
#'                                              "ln"~sigma["Orca"]^2,
#'                                              "ln"~sigma["Glob"]^2,
#'                                              "ln"~sigma["Root"]^2),
#'           smooth = TRUE, bw = 0.1,
#'           col = c('blue','red','green','black'), alpha = 0.1, border.alpha = 0.25,
#'           add.lines = c(0.25, NA, 0.75), lines.lwd = c(2, 3, 2), lines.lty = c(3, 2, 3),
#'           bty='n', legend.bty = 'n', legend.inset = c(0.1, 0), lwd = 2)
#' #or intervals based on highest posterior density regions
#' #(have to combine chains since different cuts for the same parameters aren't yet allowed)
#' #(requires HDInterval package)
#' fit <- combine.chains(cet_fit)
#' combined.par <- get.bg.rate(fit = fit,
#'                 node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                        c('Mesoplodon','Orca','Globicephalinae')),
#'                 )
#' tmp <- par.c(combined.par, "R_0", fit = fit)
#' tmp <- apply(tmp, 2, HDInterval::hdi)
#' prof.plot(list(combined.par, "R_0"), fit = fit, lower.cut = tmp[1,], upper.cut = tmp[2,],
#'           overwrite.param.names = expression("ln"~sigma["Meso"]^2,
#'                                              "ln"~sigma["Orca"]^2,
#'                                              "ln"~sigma["Glob"]^2,
#'                                              "ln"~sigma["Root"]^2),
#'           smooth = TRUE, bw = 0.1,
#'           col = c('blue','red','green','black'), alpha = 0.1, border.alpha = 0.25,
#'           add.lines = c(0.25, NA, 0.75), lines.lwd = c(2, 3, 2), lines.lty = c(3, 2, 3),
#'           bty='n', legend.bty = 'n', legend.inset = c(0.1, 0), lwd = 2)
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
                            n.params,n.chains,n.inpar,
                            get.args)
  }
}
