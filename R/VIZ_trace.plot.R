#THINGS TO FIX: LEGEND HANDLING COULD BE IMPROVED
#               CHAIN LEGEND POINTLESS IF RIBBON PLOTTING SUPPRESSED
#                 added option to suppress chain legend, but could be made more automatic
#               BETTER LINES--MODIFY LEGEND FUNCTION TO DO LITTLE SQUIGGLES?

#' Plot traces of posterior samples from a fitted evorates model
#'
#'
#' This function plots posterior samples as a function of iterations ("traces") for particular parameters from an \code{evorates_fit}
#' object or \code{param_block} array.
#'
#'
#' @param x The parameters to be plotted. Typically, this is a \code{param_block} array extracted from an \code{evorates_fit} object,
#' though it can also be a character or numeric vector specifying parameters to extract from an \code{evorates_fit} object passed to
#' the \code{fit} argument (see \code{...}). For details on how parameters are extracted, see \link{grapes-chains-grapes}. Multiple
#' \code{param_block} arrays and character/numeric selections can be combined as lists.
#' @param ribbon \code{TRUE}or \code{FALSE}: should ribbons of "rolling" posterior quantiles be plotted? Here, rolling means that
#' multiple quantiles per trace are calculated within sliding windows of posterior samples (the number of windows and their size can
#' be tweaked via the arguments \code{n.windows} and \code{window.size}; see below). The ribbons help give a sense of
#' how "stationary" the posterior samples are and thus provide insight into he convergence/mixing of chains. Set to \code{TRUE} by default.
#' @param ribbon.first \code{TRUE}or \code{FALSE}: should ribbons be plotted first, "underneath" their main trace lines, or vice versa?
#' Defaults to \code{TRUE}, such that ribbons are plotted second, "on top of" their main trace lines.
#' @param n.windows The number of windows used to calculate "rolling" quantiles for ribbons. Windows are centered such that chains
#' are broken into approximately equal intervals (i.e., specfying 10 windows for 1000 iterations will result in windows centered at
#' iterations 1, 112, 223, 334, ..., 1000). Also specifies number of windows used for lines defined by \code{add.lines} (see below).
#' @param window.size The size of windows used to calculate "rolling" quantiles for ribbons. If \code{NULL} (the default), window sizes are
#' set to about 10 times the interval between window centers. Larger window sizes result in smoother-looking ribbons. Also specifies the
#' window size used for lines defined by \code{add.lines} (see below).
#' @param p A numeric vector controlling the width of ribbons for each parameter. Specifically, the shaded
#' region is the "rolling" \code{(1-p)}\% credible interval (so 95\% credible interval by default). The vector is recycled as necessary,
#' \code{NA}'s suppress ribbon plotting, and entries below 0 or above 1 are rounded up and down to 0 and 1,
#' respectively. Overwritten by \code{lower.quant} and \code{upper.quant}, which may be used to
#' provide finer control over the bounds of the shaded regions.
#' @param lower.quant,upper.quant Numeric vectors specifying the boundaries of ribbons for each
#' parameter based on posterior quantiles. Set to \code{NULL} (the default) to just use all boundaries specified by \code{p} (see above).
#' The vectors are recycled as necessary, \code{NA}'s also default to boundaries specified by \code{p}, and entries
#' below 0 or above 1 are rounded up and down to 0 and 1, respectively.
#' @param add.lines A numeric vector specifying additional "rolling" quantiles to plot over traces. This is not yet "vectorized",
#' and will plot the same lines for all traces (though I hope to change this in the future)! \code{NA}'s result in plotting rolling
#' means.
#' @param add \code{TRUE} or \code{FALSE}: should traces be plotted in an open plotting window, or should a new plot window be opened?
#' Defaults to \code{FALSE}, meaning that a new plot window is started.
#' @param make.legend \code{TRUE}or \code{FALSE}: should a legend be automatically generated for the plotted traces? Defaults to
#' \code{TRUE}, but automatically suppressed when only a single trace is plotted.
#' @param include.chain.legend \code{TRUE} or \code{FALSE}: should the automatically-generated legend include a legend for the different
#' chains? \code{TRUE} by default, but automatically suppressed when traces from only a single chain are plotted. Particularly helpful here
#' since chains are indistinguishable by default when \code{ribbon} is \code{FALSE} (you could mess with \code{param.args}, etc. to change
#' this--see below).
#' @param ... Other arguments, such as graphical parameters. Here are the some commonly-used ones:
#' \itemize{
#' \item \code{fit} to specify an \code{evorates_fit} object from which to extract parameters (if \code{x} includes character or numeric
#' selections).
#' \item \code{overwrite.param.names} to specify alternate names for parameters in the automatically generated x-axis title and legend.
#' This is helpful if you want to include otherwise forbidden characters in parameter names, like greek symbols.
#' \item Most base R graphical parameters should work, along with a few extras:
#' \itemize{
#' \item{\code{col} by default specifies the color of main trace line for each parameter. Defaults to \code{palette()}.}
#' \item{\code{alpha} by default controls the transparency of main trace line colors for each parameter, with 0, 1, and \code{NA}
#'       corresponding to completely transparent, opaque, and no modification, respectively. Defaults to \code{NA}.}
#' \item{\code{ribbon.col} by default specifies the color of the ribbon fill for each parameter. Inherits from \code{col}
#'       if not specified.}
#' \item{\code{ribbon.border}/\code{border} by default specifies the color of the ribbon borders for each parameter.
#'       Inherits from \code{ribbon.col} if not specified.}
#' \item{\code{ribbon.angle}/\code{angle} by default controls the angle of the lines used to fill ribbons for each chain if
#'       \code{ribbon.density}/\code{density} is a positive number. Defaults to a range of angles between 0 and 180 if profiles
#'       for multiple chains are plotted.}
#' \item{\code{ribbon.density}/\code{density} by default controls the density of lines used to shade ribbons for each chain, with
#'       \code{NULL} and \code{NA} resulting in solid shading (the default if traces are plotted for a single chain). Defaults to 20
#'       if traces for multiple chains are plotted.}
#' \item{\code{lwd},\code{lty} by default specify the line width and style for all main trace lines (you could mess with \code{param.args},
#'       etc. to change this--see below).}
#' \item{Some arguments like \code{col}, \code{alpha}, \code{lwd}, etc. can be modified with "\code{lines.}", "\code{ribbon.}", or
#'       "\code{ribbon.border.}" at the beginning to make them only affect certain elements of the plot. Perhaps most importantly,
#'       "\code{legend.}" can be used to  make arguments affect the automatically generated legend if \code{make.legend} is \code{TRUE}.}
#' }}
#' @param param.args,chain.args,inpar.args Generally, these should not be altered, but are exposed here for the curious to play around
#' with. These are vectors of argument names for graphical parameters that control whether the graphical parameters vary by the
#' parameter the profile corresponds to (\code{param.args}), the chain the profile corresponds to (\code{chain.args}), or within
#' traces (\code{inpar.args}). Notably, default graphical parameters are tailored to the defaults for these arguments, and messing
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
#' trace.plot(cet_fit %chains% c("R_sig2", "R_mu"))
#' #the above is equivalent to:
#' trace.plot(c("R_sig2", "R_mu"), fit = cet_fit)
#' #could also do something like this:
#' par <- get.bg.rate(fit = cet_fit,
#'                    node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                           c('Mesoplodon','Orca','Globicephalinae')),
#' )
#' trace.plot(par)
#' #or even:
#' trace.plot(list(par, "R_0"), fit = cet_fit)
#' 
#' #some ways the plot style might be tweaked:
#' trace.plot(list(par, "R_0"), fit = cet_fit,
#'            overwrite.param.names = expression("ln"~sigma["Meso"]^2,
#'                                               "ln"~sigma["Orca"]^2,
#'                                               "ln"~sigma["Glob"]^2,
#'                                               "ln"~sigma["Root"]^2),
#'            col = c('blue','red','green','black'), alpha = 0.2, lty=1,
#'            ribbon = FALSE,
#'            add.lines = c(0.025, 0.25, NA, 0.75, 0.975),
#'            lines.alpha = 0.3, lines.lwd = c(1, 2, 3, 2, 1), lines.lty = c(3, 2, 1, 2, 3),
#'            bty='n', legend.bty = 'n', legend.x = "topleft", legend.horiz = TRUE, include.chain.legend = FALSE)
#' #may want rougher, more finely-grained rolling quantiles?
#' trace.plot(list(par, "R_0"), fit = cet_fit,
#'            overwrite.param.names = expression("ln"~sigma["Meso"]^2,
#'                                               "ln"~sigma["Orca"]^2,
#'                                               "ln"~sigma["Glob"]^2,
#'                                               "ln"~sigma["Root"]^2),
#'            col = c('blue','red','green','black'), alpha = 0.2, lty=1,
#'            ribbon = FALSE, n.windows = 100, window.size = 50,
#'            add.lines = c(0.025, 0.25, NA, 0.75, 0.975),
#'            lines.alpha = 0.3, lines.lwd = c(1, 2, 3, 2, 1), lines.lty = c(3, 2, 1, 2, 3),
#'            bty='n', legend.bty = 'n', legend.x = "topleft", legend.horiz = TRUE, include.chain.legend = FALSE)
#' 
#' 
#' @export
trace.plot<-function(x,
                     ribbon=TRUE,ribbon.first=FALSE,
                     n.windows=50,window.size=NULL,
                     p=0.05,lower.quant=NULL,upper.quant=NULL,
                     add.lines=NULL,
                     add=FALSE,make.legend=TRUE,include.chain.legend=TRUE,
                     ...,
                     param.args=c('col','lines.col',
                                  'alpha','lines.alpha',
                                  'ribbon.col','ribbon.border',
                                  'ribbon.alpha','ribbon.border.alpha'),
                     chain.args=c('ribbon.angle','ribbon.density'),
                     inpar.args=c('lines.lty','lines.lwd')){
  
  ####getting arguments####
  
  .check.args(param.args,chain.args,inpar.args)
  .x<-.proc.x(x,list(...)[['fit']])
  if((ribbon|!is.null(add.lines))&is.null(window.size)){
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
                         function(ii) 'ribbon.angle'%in%ii)))
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
                                         alpha=expression(ifelse(col2rgb(col,alpha=TRUE)[4,]<255,NA,0.3)),
                                         lwd=1,
                                         lty=1,
                                         ribbon.col=expression(col),
                                         ribbon.border=expression(ribbon.col),
                                         lines.col=expression(col),
                                         ribbon.alpha=expression(ifelse(col2rgb(ribbon.col,alpha=TRUE)[4,]<255,NA,0.2)),
                                         ribbon.border.alpha=expression(alpha),
                                         lines.alpha=expression(ribbon.border.alpha),
                                         ribbon.angle=def.angle,
                                         ribbon.density=def.density,
                                         ribbon.lwd=expression(2*lwd),
                                         lines.lwd=expression(max(1,ribbon.lwd/2)),
                                         ribbon.lty=1,
                                         lines.lty=expression(lty+1)),
                           trace=TRUE)
  
  ####getting ribbons/lines####
  
  if(ribbon|!is.null(add.lines)){
    add.lines<-matrix(if(is.null(add.lines)) NA else add.lines,nrow=n.params*n.chains,ncol=n.inpar,byrow=TRUE)
    ribbons.and.lines<-.get.cuts(.x,p,
                                 lower.quant,upper.quant,
                                 lower.cut=NULL,upper.cut=NULL,
                                 n.params,n.chains,
                                 trace=TRUE,ribbon=ribbon,
                                 add.lines=add.lines,
                                 n.windows=n.windows,window.size=window.size)
    lines.xx<-ribbons.and.lines[[1]][[1]]
    if(ribbon){
      ribbons<-lapply(ribbons.and.lines[[2]],function(ii) ii[,seq_len(2),drop=FALSE])
      lines<-lapply(ribbons.and.lines[[2]],function(ii) ii[,2+seq_len(n.inpar),drop=FALSE])
    }else{
      lines<-ribbons.and.lines[[2]]
    }
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
  
  ####plotting traces####
  
  .x<-asplit(.x,-1)
  xx<-seq_len(dim(.x[[1]])[1])
  for(i in seq_len(n.tot)){
    ii<-if(n.inpar) n.inpar*(i-1)+1 else i
    if(!ribbon.first){
      tmp.args<-get.args(~gen&!recyc)
      tmp.recyc.args<-get.args(~gen&recyc,ii)
      do.call(graphics::lines,
              c(x=list(xx),
                y=list(.x[[i]]),
                tmp.args,
                tmp.recyc.args))
    }
    if(ribbon){
      tmp.args<-get.args(~ribbon&nonborder&!recyc)
      tmp.recyc.args<-get.args(~ribbon&nonborder&recyc,ii)
      do.call(polygon,
              c(x=list(c(lines.xx,rev(lines.xx))),
                y=list(c(ribbons[[i]][,1],rev(ribbons[[i]][,2]))),
                border=NA,
                tmp.args,
                tmp.recyc.args))
      tmp.args<-get.args(~ribbon&noncol&!recyc)
      tmp.recyc.args<-get.args(~ribbon&noncol&recyc,ii)
      do.call(polygon,
              c(x=list(c(lines.xx,rev(lines.xx))),
                y=list(c(ribbons[[i]][,1],rev(ribbons[[i]][,2]))),
                col=list(rgb(0,0,0,0)),
                tmp.args,
                tmp.recyc.args))
    }
    if(ribbon.first){
      tmp.args<-get.args(~gen&!recyc)
      tmp.recyc.args<-get.args(~gen&recyc,ii)
      do.call(graphics::lines,
              c(x=list(xx),
                y=list(.x[[i]]),
                tmp.args,
                tmp.recyc.args))
    }
    if(n.inpar){
      ii<-seq.int(ii,ii+n.inpar-1)
      tmp.args<-get.args(~lines&!recyc)
      tmp.recyc.args<-get.args(~lines&recyc,ii)
      do.call(matplot,
              c(x=list(lines.xx),
                y=list(lines[[i]]),
                add=TRUE,
                type="l",
                tmp.args,
                tmp.recyc.args))
    }
  }
  
  ####plotting legend####
  
  if(make.legend){
    .make.param.plot.legend(param.names,chain.names,
                            n.params,n.chains,n.inpar,
                            get.args,include.chain.legend,
                            trace=TRUE,ribbon=ribbon,ribbon.first=ribbon.first)
  }
}