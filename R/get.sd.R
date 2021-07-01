#' @export
get.sd<-function(fit){
  0.5*logspline::dlogspline(0,logspline::logspline(fit%chains%R_sig2,lbound=0))/
    dcauchy(0,0,fit$call$Rsig2_prior_sig)
}