#' @export
get.cov.mat<-function(fit,traits=colnames(fit$call$trait.data),cor=F,
                      element=c('chains','quantiles','means','MAPs','diagnostics'),
                      type=c('evocov','intracov'),
                      ret=c('complete','diagonal','lower.triangle'),
                      select.extra=NULL,simplify=T){
  if(!inherits(fit,'corateBM_fit')){
    stop("fit must be a fitted correlated rates BM fit (class 'corateBM_fit')")
  }
  try.type<-try(match.arg(type,c('evocov','intracov')),silent=T)
  if(inherits(try.type,'try-error')){
    stop(type," is not an available option for covariance type: please specifiy either 'evocov' for evolutionary covariance or 'intracov' for intraspecific covariance")
  }
  type<-try.type
  if(ncol(fit$call$trait.data)==1&type=='evocov'){
    stop('fit appears to be for univariate trait dataset: there are no evolutionary (co)variance parameters')
  }
  try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics')),silent=T)
  if(inherits(try.element,'try-error')){
    stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', or 'diagnostics'")
  }
  element<-try.element
  if(type=='intracov'&sum(grepl('intracov',dimnames(fit$chains)[['parameters']]))==0){
    warning("covariance type 'intracov' selected, but no intraspecific variance modeled in corateBM_fit: defaulted to 'evocov'")
    type<-'evocov'
  }
  try.ret<-try(match.arg(ret,c('complete','diagonal','lower.triangle')),silent=T)
  if(inherits(try.ret,'try-error')){
    stop(ret," is not an available option for return type: please specifiy one of the following: 'complete', 'diagonal', or 'lower.triangle'")
  }
  ret<-try.ret
  if(is.numeric(traits)){
    traits<-colnames(fit$call$trait.data)[traits]
  }
  #check if any trait names not available
  traits.exist<-traits%in%colnames(fit$call$trait.data)
  if(all(!traits.exist)){
    stop('none of the specified traits found')
  }
  if(any(!traits.exist)){
    warning(paste(traits[which(!traits.exist)],collapse=', '),' not found')
    traits<-traits[which(traits.exist)]
  }
  if(ret=='lower.triangle'&length(traits)==1){
    warning("since only a single trait is specified, the resulting covariance matrix will have no lower triangle: set ret to 'complete'")
  }
  if(cor){
    element<-'chains'
  }
  k<-length(traits)
  new.traits<-gsub(',','\\\\,',traits)
  new.traits<-gsub('_','\\\\_',new.traits)
  select<-paste(rep(new.traits,each=k),',',
                rep(new.traits,k),paste('_',type,sep=''),
                sep='')[lower.tri(matrix(NA,k,k),diag=T)]
  if(element=='quantiles'|element=='diagnostics'&!is.null(select.extra)){
    select<-list(select,select.extra)
  }
  tmp<-do.call(paste('.int.',element,sep=''),list(fit=fit,select=select))
  out<-array(NA,c(k,k,dim(tmp)[1],dim(tmp)[3]),c(rep(list(parameters=traits),2),dimnames(tmp)[c(1,3)]))
  for(i in 1:k){
    for(j in 1:i){
      tmp.name<-grep(paste(paste(traits[c(i,j)],',',traits[c(j,i)],'_',type,sep=''),collapse='|'),
                     dimnames(tmp)[[2]])
      out[i,j,,]<-tmp[,tmp.name,]
      if(i!=j){
        out[j,i,,]<-tmp[,tmp.name,]
      }
    }
  }
  if(cor){
    out[,,,]<-unlist(lapply(asplit(out,c(3,4)),cov2cor))
    element<-try.element
    if(element!='chains'){
      first.cor<-matrix(NA,length(traits),length(traits))
      first.cor[lower.tri(first.cor,diag=T)]<-fit%diagnostics%list(select,'inits')
      first.cor[upper.tri(first.cor)]<-t(first.cor)[upper.tri(first.cor)]
      first.cor<-suppressWarnings(cov2cor(first.cor))
      tmp<-list(chains=.index.element(.coerce.to.3D(.expand.element(out)),which(lower.tri(first.cor,diag=T)),2),
                sampler.params=fit$sampler.params)
      class(tmp)<-'corateBM_fit'
      tmp$param.diags<-.index.element(tmp$chains,1:4,1)
      dimnames(tmp$param.diags)<-c(diagnostics=list(c('inits','bulk_ess','tail_ess','Rhat')),
                                   dimnames(tmp$chains)[-1])
      tmp$param.diags[1,,]<-first.cor[lower.tri(first.cor,diag=T)]
      tmp$param.diags[2,,]<-apply(tmp$chains,c(2,3),rstan::ess_bulk)
      tmp$param.diags[3,,]<-apply(tmp$chains,c(2,3),rstan::ess_tail)
      tmp$param.diags[4,,]<-apply(tmp$chains,c(2,3),rstan::Rhat)
      select<-gsub('_evocov|_intracov','',select)
      if(element=='quantiles'|element=='diagnostics'&!is.null(select.extra)){
        select<-list(select,select.extra)
      }
      tmp<-do.call(paste('.int.',element,sep=''),list(fit=tmp,select=select))
      out<-array(NA,c(k,k,dim(tmp)[1],dim(tmp)[3]),c(rep(list(parameters=traits),2),dimnames(tmp)[c(1,3)]))
      for(i in 1:k){
        for(j in 1:i){
          tmp.name<-grep(paste(paste(traits[c(i,j)],traits[c(j,i)],sep=','),collapse='|'),
                         dimnames(tmp)[[2]])
          out[i,j,,]<-tmp[,tmp.name,]
          if(i!=j){
            out[j,i,,]<-tmp[,tmp.name,]
          }
        }
      }
    }
  }
  if(ret=='complete'){
    if(simplify){
      new.dims<-c(dim(out)[1:2],ifelse(dim(out)[3:4]==1,NA,dim(out)[3:4]))
      new.dimnames<-dimnames(out)[!is.na(new.dims)]
      new.out<-array(out,new.dims[!is.na(new.dims)],new.dimnames)
      for(i in which(is.na(new.dims))){
        attr(new.out,names(dimnames(out))[i])<-dimnames(out)[[i]]
      }
      out<-new.out
    }
  }else{
    if(ret=='diagonal'){
      new.out<-array(NA,dim(out)[-1],dimnames(out)[-1])
      new.out[,,]<-apply(out,c(3,4),diag)
      out<-aperm(new.out,c(2,1,3))
    }else if(ret=='lower.triangle'){
      out<-.index.element(.coerce.to.3D(.expand.element(out)),which(lower.tri(out[,,1,1])),2)
    }
    if(simplify){
      out<-.simplify.element(out)
    }
  }
  out
}

#corateBM_fit names are technically backwards if they're supposed to be the lower triangle--the column names come first!