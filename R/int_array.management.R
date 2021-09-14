#newest version of check.n.proc and coerce.to.array--expands any given element to a 3D array (if chains,
#quantiles, sampler, or parameter diagnostics) or a 2D matrix (if means or MAPs)
#finally got it to a point where it should work with just about anything...
#I made it just always return 3D arrays, and label means and MAPs projects with 'iterations' along
#1st dimension; just makes things easier for other functions to deal with...
#rendered function more robust in the face of 4D arrays outputted from get.cov.mat and get.trait.mat: side-benefit of now
#also reorganizing any arrays back into iterations/quantiles/diagnostics, then parameters (however many there are), then
#chains order

#less reliance on regular expressions might make this substantially faster...
#9/14/21--twice as fast! woo!
.expand.element<-function(arr,simplify=FALSE){
  element.type<-.get.element.type(arr)
  if(length(dim(arr))==0){
    new.arr<-matrix(arr,length(arr))
    check<-rownames(new.arr)<-names(.strip.ele(arr))
    new.dimnames<-dimnames(new.arr)
    if(is.null(check)){
      new.dimnames<-rep(list(NULL),2)
      names(new.dimnames)[1]<-'iterations'
    }else if(element.type=='quantiles'&is.null(attr(arr,'quantiles'))){
      names(new.dimnames)[1]<-'quantiles'
    }else if(element.type=='diagnostics'&is.null(attr(arr,'diagnostics'))){
      names(new.dimnames)[1]<-'diagnostics'
    }else if(is.null(attr(arr,'chains'))){
      names(new.dimnames)[1]<-'chains'
    }else{
      names(new.dimnames)[1]<-'parameters'
    }
    for(i in c('quantiles','diagnostics','parameters','chains')){
      if(!is.null(attr(arr,i))){
        attr(new.arr,i)<-attr(arr,i)
      }
    }
    dimnames(new.arr)<-new.dimnames
    arr<-.add.ele(new.arr)
  }
  dims<-names(dimnames(arr))
  dim.map<-vector('list',3)
  tmp<-c('quantiles','diagnostics','iterations')
  type.code<-tmp[match(element.type,tmp[-3],nomatch=3)]
  names(dim.map)<-c(type.code,'parameters','chains')
  for(i in names(dim.map)){
    dim.map[[i]]<-which(dims==i)
  }
  if(any(lengths(dim.map)==0)){
    tmp<-dim.map
    tmp[lengths(tmp)==0]<-NA
    lens<-lengths(tmp)
    tmp<-unlist(tmp)
    out.dim<-dim(arr)[tmp]
    out.dim[is.na(tmp)]<-1
    out.dimnames<-dimnames(arr)[tmp]
    names(out.dimnames)<-rep(names(dim.map),lens)
    ndims<-length(tmp)
    if(lengths(dim.map)[1]==0){
      out.dimnames[1]<-list(attr(arr,type.code))
    }
    if(lengths(dim.map)[2]==0){
      out.dimnames[-c(1,ndims)]<-list(attr(arr,'parameters'))
    }
    if(lengths(dim.map)[3]==0){
      out.dimnames[ndims]<-list(attr(arr,'chains'))
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
  }else if(any(old.dims==1)){
    old.dimnames<-dimnames(arr)
    inds<-old.dims>1
    out<-array(arr,old.dims[inds],old.dimnames[inds])
    if(length(dim(out))==1){
      out<-setNames(as.vector(out),dimnames(out)[[1]])
    }
    for(i in 1:length(old.dims)){
      if(!inds[i]){
        if(!is.null(old.dimnames[[i]])){
          attr(out,names(old.dimnames)[i])<-old.dimnames[[i]]
        }
      }
    }
    if(length(out)==1){
      names(out)<-attr(out,'parameters')
      attr(out,'parameters')<-NULL
    }
    for(i in c('quantiles','diagnostics','parameters','chains','element')){
      if(!is.null(attr(arr,i))){
        attr(out,i)<-attr(arr,i)
      }
    }
    out<-.add.ele(out)
    out
  }else{
    arr
  }
}

#newest, generalized version of reduce.array--selects indices from arbitrary arrays without
#simplification or destroying names
#technically, this can be done with [,,...,drop=F], but I didn't know that at the time...
#one of the ladder functions depends on this guys--I'll keep it for now 9/14/21
.index.element<-function(arr,inds,dims,invert=FALSE,allow.reorder=FALSE){
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

#make return one code if iterations are messed up, another if chains are messed up --> done
#returns a list of expanded elements, their dimensions, their dimension names, whether they match, and whether 1st/3rd dim match
.check.dims.compat<-function(...){
  exps<-lapply(list(...),.expand.element)
  dims<-lapply(exps,dim)
  nms<-lapply(exps,dimnames)
  tmp.dims<-dims[[1]]
  tmp.dims<-tmp.dims[c(1,length(tmp.dims))]
  tmp.nms<-nms[[1]]
  tmp.nms<-tmp.nms[c(1,length(tmp.nms))]
  foo<-function(ind){
    dims<-dims[[ind]]
    len<-length(dims)
    dims<-dims[c(1,len)]
    nms<-nms[[ind]][c(1,len)]
    out<-vector(length=2)
    for(i in 1:2){
      if(dims[i]==tmp.dims[i]){
        if(dims[i]==0){
          out[i]<-TRUE
        }else{
          out[i]<-all(nms[[i]]==tmp.nms[[i]])
        }
      }else{
        out[i]<-FALSE
      }
    }
    out
  }
  checks<-matrix(unlist(lapply(2:length(exps),foo)),
                 ncol=2,byrow=TRUE)
  dim1.check<-all(checks[,1])
  chain.check<-all(checks[,2])
  out<-list(exps,dims,nms,TRUE,c(dim1.check,chain.check))
  if(dim1.check&chain.check){
    out
  }else{
    out[[4]]<-FALSE
    out
  }
}

.is.evorates.element<-function(params){
  check.vec<-c('iterations','quantiles','diagnostics','parameters','chains')
  any(names(dimnames(params))%in%check.vec)|any(names(attributes(params))%in%check.vec)
}

.get.element.type<-function(arr){
  attr(arr,'element')
}

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
