#' @export
get.cov.mat<-function(fit,traits=colnames(fit$call$trait.data),
                      element=c('chains','quantiles','means','MAPs','diagnostics'),
                      type=c('evocov','intracov'),
                      output.list=F,select.extra=NULL,
                      simplify=T){
  if(!inherits(fit,'corateBM_fit')){
    stop("fit must be a fitted correlated rates BM fit (class 'corateBM_fit')")
  }
  if(ncol(fit$call$trait.data)==1){
    stop('fit appears to be for univariate trait dataset: there are no trait covariance parameters')
  }
  try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics')),silent=T)
  if(inherits(try.element,'try-error')){
    stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', or 'diagnostics'")
  }
  element<-try.element
  try.type<-try(match.arg(type,c('evocov','intracov')),silent=T)
  if(inherits(try.type,'try-error')){
    stop(type," is not an available option for covariance type: please specifiy either 'evocov' for evolutionary covariance or 'intracov' for intraspecific covariance")
  }
  type<-try.type
  if(type=='intracov'&sum(grepl('intracov',dimnames(fit$chains)[['parameters']]))==0){
    warning("covariance type 'intracov' selected, but no intraspecific variance modeled in corateBM_fit: defaulted to 'evocov'")
    type<-'evocov'
  }
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
  if(output.list){
    out<-asplit(out,4)
    out<-lapply(out,function(ii) if(dim(ii)[3]==1) ii else asplit(ii,3))
    if(length(out)==1){
      attr(out[[1]],'chains')<-names(out)
      out<-out[[1]]
    }
  }else if(simplify){
    new.dims<-c(dim(out)[1:2],ifelse(dim(out)[3:4]==1,NA,dim(out)[3:4]))
    new.dimnames<-dimnames(out)[!is.na(new.dims)]
    new.out<-array(out,new.dims[!is.na(new.dims)],new.dimnames)
    for(i in which(is.na(new.dims))){
      attr(new.out,names(dimnames(out))[i])<-dimnames(out)[[i]]
    }
    out<-new.out
  }
  out
}