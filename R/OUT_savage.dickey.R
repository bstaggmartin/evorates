#' @export
get.sd<-function(fit){
  #interesting--logspline internally changes names of evorates element
  #this is potentially indicative of a more general undesirable side effect
  if(fit$call$constrain.Rsig2){
    stop("Fit needs unconstrained R_sig2 parameter to calculate Savage-Dickey ratio!")
  }
  0.5*dlogspline(0,logspline(.strip.par.class(.call.select(fit$chains,'R_sig2')),lbound=0))/
    dcauchy(0,0,fit$call$Rsig2_prior_sig)
}