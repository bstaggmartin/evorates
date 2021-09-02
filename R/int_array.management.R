#newest version of check.n.proc and coerce.to.array--expands any given element to a 3D array (if chains,
#quantiles, sampler, or parameter diagnostics) or a 2D matrix (if means or MAPs)
#finally got it to a point where it should work with just about anything...
#I made it just always return 3D arrays, and label means and MAPs projects with 'iterations' along
#1st dimension; just makes things easier for other functions to deal with...
#rendered function more robust in the face of 4D arrays outputted from get.cov.mat and get.trait.mat: side-benefit of now
#also reorganizing any arrays back into iterations/quantiles/diagnostics, then parameters (however many there are), then
#chains order
.expand.element<-function(arr,simplify=FALSE){
  element.type<-.get.element.type(arr)
  if(length(dim(arr))==0){
    new.arr<-matrix(arr,length(arr))
    rownames(new.arr)<-names(.strip.ele(arr))
    new.dimnames<-dimnames(new.arr)
    if(is.null(dimnames(new.arr)[[1]])){
      new.dimnames<-rep(list(NULL),2)
      names(new.dimnames)[1]<-'iterations'
    }else if(sum(grepl('%$',dimnames(new.arr)[[1]]))!=0&is.null(attr(arr,'quantiles'))){
      names(new.dimnames)[1]<-'quantiles'
    }else if(sum(grepl('^inits$|^bulk_ess$|^tail_ess$|^Rhat$',dimnames(new.arr)[[1]]))!=0&is.null(attr(arr,'diagnostics'))){
      names(new.dimnames)[1]<-'diagnostics'
    }else if(sum(grepl('chain',dimnames(new.arr)[[1]]))!=0&is.null(attr(arr,'chains'))){
      names(new.dimnames)[1]<-'chains'
    }else{
      names(new.dimnames)[1]<-'parameters'
    }
    for(i in c('quantiles','diagnostics','parameters','chains')){
      if(!is.null(attr(arr,i))){
        attr(new.arr,i)<-attr(arr,i)
      }
    }
    class(new.arr)<-c(class(new.arr),'loose_element')
    dimnames(new.arr)<-new.dimnames
    arr<-new.arr
  }
  dims<-names(dimnames(arr))
  dim.map<-lapply(c('iterations|quantiles|diagnostics','parameters','chains'),function(ii) grep(ii,dims))
  if(any(lengths(dim.map)==0)){
    tmp<-dim.map
    tmp[lengths(tmp)==0]<-NA
    tmp<-unlist(tmp)
    out.dim<-dim(arr)[tmp]
    out.dim[is.na(tmp)]<-1
    out.dimnames<-dimnames(arr)[tmp]
    ndims<-length(tmp)
    if(lengths(dim.map)[1]==0){
      if('quantiles'%in%names(attributes(arr))){
        names(out.dimnames)[1]<-'quantiles'
        out.dimnames[1]<-list(attr(arr,'quantiles'))
      }else if('diagnostics'%in%names(attributes(arr))){
        names(out.dimnames)[1]<-'diagnostics'
        out.dimnames[1]<-list(attr(arr,'diagnostics'))
      }else{
        names(out.dimnames)[1]<-'iterations'
        out.dimnames[1]<-list(NULL)
      }
    }
    if(lengths(dim.map)[2]==0){
      names(out.dimnames)[-c(1,ndims)]<-'parameters'
      if('parameters'%in%names(attributes(arr))){
        out.dimnames[-c(1,ndims)]<-list(attr(arr,'parameters'))
      }else{
        out.dimnames[-c(1,ndims)]<-list(NULL)
      }
    }
    if(lengths(dim.map)[3]==0){
      names(out.dimnames)[ndims]<-'chains'
      if('chains'%in%names(attributes(arr))){
        out.dimnames[ndims]<-list(attr(arr,'chains'))
      }else{
        out.dimnames[ndims]<-list(NULL)
      }
    }
    if(any(is.na(names(dimnames(arr))))){
      arr<-array(arr,out.dim,out.dimnames)
    }else{
      arr<-array(aperm(arr,unlist(dim.map)),out.dim,out.dimnames)
    }
  }else{
    dim.map<-unlist(dim.map)
    if(any(dim.map!=1:length(dim.map))){
      arr<-array(aperm(arr,dim.map),dim(arr)[dim.map],dimnames(arr)[dim.map])
    }
  }
  attr(arr,'element')<-element.type
  arr<-.add.ele(arr)
  if(simplify){
    if(is.null(dimnames(arr)[[1]])&dim(arr)[1]==1){
      arr<-array(arr,dim(arr)[-1],dimnames(arr)[-1])
    }
  }
  arr
}

#opposite of expand.element--collapses dimensions of length 1 and adds them to parameters. In cases of
#length 1 output, prioritizes parameter names over everything else
##NEED TO UPDATE TO WORK WITH 2D ELEMENTS (like means and MAPs)--done 7/27
.simplify.element<-function(arr){
  old.dims<-dim(arr)
  if(is.null(old.dims)){
    arr
  }else{
    old.dimnames<-dimnames(arr)
    out<-do.call('[',c(list(arr),lapply(old.dims,function(ii) -(ii+1))))
    for(i in 1:length(old.dims)){
      if(old.dims[i]==1){
        if(!is.null(old.dimnames[[i]])){
          attr(out,names(old.dimnames)[i])<-old.dimnames[[i]]
        }
      }
    }
    if(length(out)==1){
      names(out)<-attr(out,names(old.dimnames)[2])
      attr(out,names(old.dimnames)[2])<-NULL
    }
    for(i in c('quantiles','diagnostics','parameters','chains','element')){
      if(!is.null(attr(arr,i))){
        attr(out,i)<-attr(arr,i)
      }
    }
    out<-.add.ele(out)
    out
  }
}

#newest, generalized version of reduce.array--selects indices from arbitrary arrays without
#simplification or destroying names
#tehcnically, this can be done with [,,...,drop=F], but I didn't know that at the time...
.index.element<-function(arr,inds,dims,invert=F,allow.reorder=F){
  inds.list<-vector('list',length(dim(arr)))
  if(is.numeric(inds)){
    inds<-list(inds)
  }
  inds.list[dims]<-inds
  if(invert){
    inds.list<-lapply(inds.list,function(ii) -1*ii)
  }else{
    inds.list<-lapply(1:length(inds.list),
                      function(ii) if(is.null(inds.list[[ii]])) 1:dim(arr)[ii] else inds.list[[ii]])
  }
  new.dims<-dim(arr)
  if(invert){
    new.dims<-new.dims-lengths(inds.list)
    inds.list[lengths(inds.list)==0]<- -dim(arr)[lengths(inds.list)==0]-1
  }else{
    new.dims<-lengths(inds.list)
  }
  if(!allow.reorder){
    inds.list<-lapply(inds.list,function(ii) unique(sort(ii)))
  }
  new.dimnames<-dimnames(arr)
  new.dimnames<-lapply(1:length(inds.list),
                       function(ii) new.dimnames[[ii]][inds.list[[ii]]])
  names(new.dimnames)<-names(dimnames(arr))
  out<-array(do.call('[',c(list(arr),inds.list)),new.dims,new.dimnames)
  for(i in c('quantiles','diagnostics','parameters','chains','element')){
    if(!is.null(attr(arr,i))){
      attr(out,i)<-attr(arr,i)
    }
  }
  out
}

.is.evorates.element<-function(params){
  check.vec<-c('iterations','quantiles','diagnostics','parameters','chains')
  any(names(dimnames(params))%in%check.vec)|any(names(attributes(params))%in%check.vec)
}

.get.element.type<-function(arr){
  attr(arr,'element')
}

.coerce.to.3D<-function(arr){
  if(length(dim(arr))<3){
    arr<-.expand.element(arr)
  }
  if(length(dim(arr))>3){
    new.dim<-dim(arr)
    new.dimnames<-dimnames(arr)
    param.dims<-which(names(dimnames(arr))=='parameters')
    param.names<-dimnames(arr)[param.dims]
    param.names<-do.call(paste,
                         c(lapply(1:length(param.names),function(ii)
                           rep(rep(param.names[[ii]],prod(lengths(param.names)[(1:length(param.names))[1:length(param.names)>ii]])),
                               each=prod(lengths(param.names)[(1:length(param.names))[1:length(param.names)<ii]]))),
                           sep=','))
    new.dim<-new.dim[-param.dims]
    new.dim<-append(new.dim,length(param.names),1)
    new.dimnames<-new.dimnames[-param.dims]
    new.dimnames<-append(new.dimnames,list(param.names),1)
    names(new.dimnames)<-c(names(new.dimnames)[1],'parameters',names(new.dimnames)[3])
    arr<-array(arr,new.dim,new.dimnames)
  }
  arr
}

.combine.elements<-function(in.params,fit=NULL,element=NULL,select.extra=NULL,simplify=T){
  params<-in.params
  if(.is.evorates.element(params)){
    params<-list(params)
  }
  if(!is.list(params)){
    params<-as.list(params)
  }
  params<-params[lengths(params)>0]
  types<-sapply(params,function(ii) if(.is.evorates.element(ii)) 'element' else 'select')
  if(all(types=='select')&is.null(fit)){
    stop(deparse(substitute(in.params)),' appears to consist of strings/numbers specifying parameters to extract out of a evorates_fit, but no evorates_fit is supplied')
  }
  if(any(types=='select')&is.null(fit)){
    warning(deparse(substitute(in.params)),' contains ',paste(params[which(types=='select')],collapse=', '),', which appear to be strings/numbers specifying parameters to extract out of a evorates_fit, but no evorates_fit is supplied: these strings/numbers were excluded')
    params<-params[-which(types=='select')]
    types<-types[-which(types=='select')]
  }
  #getting provided parameters
  params[types=='element']<-lapply(params[types=='element'],function(ii) .coerce.to.3D(ii))
  select.params<-NULL
  
  
  
  
  #getting parameters selected by characters/numbers
  if(sum(types=='select')>0){
    select<-unlist(params[types=='select'])
    #trying to find right element to use if none provided
    if(is.null(element)){
      if(sum(types=='element')>0){
        
        
        
        element.types<-sapply(params[types=='element'],.get.element.type)
        if(any(element.types%in%c('ambiguous','unrecognized'))|length(unique(element.types))>1){
          stop('element is unspecified, but element type based on provided loose elements is ambiguous: try specifying which element you wish to extract and double-check all loose elements are of the same type (chains, quantiles, means, etc.)')
        }else{
          element<-unique(element.types)
        }
      }else{
        element<-'chains'
      }
      
      
      
    }else{
      try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics','sampler')),silent=T)
      if(inherits(try.element,'try-error')){
        stop(element," is not an available element to extract from a evolving rates model fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', 'diagnostics', or 'sampler'")
      }
      element<-try.element
    }
    if(element=='quantiles'|element=='diagnostics'){
      if(is.null(select.extra)&sum(types=='select')>0){
        dim1.names<-lapply(params[types=='element'],function(ii) dimnames(ii)[[1]])
        if(length(unique(dim1.names))>1){
          stop('mismatching quantiles and/or parameter diagnostics in provided loose elements')
        }else if(element=='quantiles'){
          select.extra<-as.numeric(substr(dim1.names[[1]],0,nchar(dim1.names[[1]])-1))/100
        }else{
          select.extra<-dim1.names
        }
      }
    }
    if(!is.null(select.extra)){
      select<-list(select,select.extra)
    }
    select.params<-do.call(paste0('.int.',element),list(fit=fit,select=select))
  }
  
  
  
  if(is.null(params[types=='element'])){
    out<-select.params
  }else{
    if(!is.null(select.params)){
      out<-c(list(select.params),params[types=='element'])
    }else{
      out<-params[types=='element']
    }
    out.dim<-vector(mode='list',length=3)
    out.dimnames<-vector(mode='list',length=3)
    for(i in 1:3){
      if(i==2){
        out.dim[[i]]<-sum(sapply(out,function(ii) dim(ii)[i]))
        out.dimnames[[i]]<-unlist(lapply(out,function(ii) dimnames(ii)[[i]]))
      }else{
        out.dim[[i]]<-unique(sapply(out,function(ii) dim(ii)[i]))
        out.dimnames[[i]]<-unique(lapply(out,function(ii) dimnames(ii)[[i]]))
      }
    }
    if(all(lengths(c(out.dim[-2],out.dimnames[-2]))==1)){
      for(i in c(1,3)){
        out.dim[[i]]<-out.dim[[i]][[1]]
        out.dimnames[i]<-out.dimnames[[i]][1]
      }
      out.dim<-unlist(out.dim)
    }else{
      stop('dimensional mismatch in elements specified by ',deparse(substitute(in.params)),': did these all come from the same evorates_fit, and were they all extracted using the same parameters?')
    }
    names(out.dimnames)<-c(names(dimnames(out[[1]]))[1],'parameters','chains')
    out.arr<-array(NA,out.dim,out.dimnames)
    class(out.arr)<-c(class(out.arr),'loose_element')
    counter<-0
    for(i in 1:length(out)){
      tmp.dim<-dim(out[[i]])[2]
      out.arr[,counter+1:tmp.dim,]<-out[[i]]
      counter<-counter+tmp.dim
    }
    out<-out.arr
  }
  if(simplify){
    out<-.simplify.element(out)
  }
  out
}
#can get duplicates of the same parameter and probs not super efficient since it runs a separate call to %chains% each time it runs...
#could cannabalize this function to add parameters to your evorates_fit...
#yeah, I did a minor update to begin generalizing this, but ultimately I'd like it to be a combine.elements function...
#need to figure out what to do in cases where the element is ambiguous
#also, should 4-D arrays be coerced to 3-D with .expand.element?
#now generalized to simply combine all elements it's face with, including ones specified by select. Tries its best to match
#element type, but doesn't do anything like coercing provided elements to other types on the fly (considering doing this in the
#future...)

.strip.ele<-function(x){
  tmp<-class(x)
  class(x)<-tmp[!which(tmp=='loose_element')]
  x
}

.add.ele<-function(x){
  tmp<-class(x)
  class(x)<-c(tmp,'loose_element')
  x
}

.check.dims.compat<-function(l,r){
  l<-.expand.element(l)
  r<-.expand.element(r)
  ldims<-dim(l)
  rdims<-dim(r)
  lnames<-dimnames(l)
  rnames<-dimnames(r)
  if(ldims[1]==rdims[1]&ldims[length(ldims)]==rdims[length(rdims)]&
     all(lnames[[1]]==rnames[[1]])&all(lnames[[3]]==rnames[[3]])){
    list(l,r,ldims,rdims,lnames,rnames)
  }else{
    NULL
  }
}