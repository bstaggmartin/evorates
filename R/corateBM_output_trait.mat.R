#' @export
get.trait.mat<-function(fit,tips=fit$call$tree$tip.label,traits=colnames(fit$call$trait.data),
                      element=c('chains','quantiles','means','MAPs','diagnostics'),
                      include.dat=T,
                      select.extra=NULL,simplify=T){
  if(!inherits(fit,'corateBM_fit')){
    stop("fit must be a fitted correlated rates BM fit (class 'corateBM_fit')")
  }
  try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics')),silent=T)
  if(inherits(try.element,'try-error')){
    stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', or 'diagnostics'")
  }
  element<-try.element
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
  if(is.numeric(tips)){
    tips<-fit$call$tree$tip.label[tips]
  }
  #check if any trait names not available
  tips.exist<-tips%in%fit$call$tree$tip.label
  if(all(!tips.exist)){
    stop('none of the specified tips found')
  }
  if(any(!tips.exist)){
    warning(paste(tips[which(!tips.exist)],collapse=', '),' not found')
    tips<-tips[which(tips.exist)]
  }
  k<-length(traits)
  n<-length(tips)
  new.traits<-gsub(',','\\\\,',traits)
  new.traits<-gsub('_','\\\\_',new.traits)
  new.tips<-gsub(',','\\\\,',tips)
  new.tips<-gsub('_','\\\\_',new.tips)
  select<-paste(rep(new.traits,each=n),rep(tips,k),sep='_')
  if(element=='quantiles'|element=='diagnostics'&!is.null(select.extra)){
    select<-list(select,select.extra)
  }
  tmp<-suppressWarnings(try(do.call(paste('.int.',element,sep=''),list(fit=fit,select=select)),silent=T))
  if(inherits(tmp,'try-error')){
    if(!include.dat){
      stop('no tip mean parameters available: set include.dat to TRUE to return trait data')
    }else if(element=='diagnostics'){
      stop('no tip mean parameters available: no associated parameter diagnostics to extract')
    }else{
      warning('no tip mean parameters available: returning trait data')
      out<-array(NA,c(n,k,1,1),list(parameters=tips,parameters=traits,iterations=NULL,chains=NULL))
    }
  }else{
    out<-array(NA,c(n,k,dim(tmp)[1],dim(tmp)[3]),c(parameters=list(tips),parameters=list(traits),dimnames(tmp)[c(1,3)]))
  }
  for(i in tips){
    for(j in traits){
      tmp.name<-paste(j,i,sep='_')
      if(sum(grepl(paste(tmp.name,'$',sep=''),dimnames(tmp)[[2]]))>0){
        out[i,j,,]<-tmp[,tmp.name,]
      }else if(element%in%c('chains','quantiles','means','MAPs')&include.dat){
        out[i,j,,]<-fit$call$trait.data[i,j]
      }
    }
  }
  if(simplify){
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
