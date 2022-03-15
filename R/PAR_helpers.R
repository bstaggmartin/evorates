#opposite of .simplify.par
#make param_blocks fully expanded, with dimension going in iterations > parameters > chains order
.expand.par<-function(x,simplify=FALSE){
  par.type<-.get.par.type(x)
  if(length(dim(x))==0){
    new.x<-matrix(x,length(x))
    check<-rownames(new.x)<-names(.strip.par.class(x))
    new.dimnames<-dimnames(new.x)
    if(is.null(check)){
      new.dimnames<-rep(list(NULL),2)
      names(new.dimnames)[1]<-'iterations'
    }else if(par.type=='quantiles'&is.null(attr(x,'quantiles'))){
      names(new.dimnames)[1]<-'quantiles'
    }else if(par.type=='diagnostics'&is.null(attr(x,'diagnostics'))){
      names(new.dimnames)[1]<-'diagnostics'
    }else if(is.null(attr(x,'chains'))){
      names(new.dimnames)[1]<-'chains'
    }else{
      names(new.dimnames)[1]<-'parameters'
    }
    for(i in c('quantiles','diagnostics','parameters','chains')){
      if(!is.null(attr(x,i))){
        attr(new.x,i)<-attr(x,i)
      }
    }
    dimnames(new.x)<-new.dimnames
    x<-.add.par.class(new.x)
  }
  dims<-names(dimnames(x))
  dim.map<-vector('list',3)
  tmp<-c('quantiles','diagnostics','iterations')
  type.code<-tmp[match(par.type,tmp[-3],nomatch=3)]
  names(dim.map)<-c(type.code,'parameters','chains')
  for(i in names(dim.map)){
    dim.map[[i]]<-which(dims==i)
  }
  if(any(lengths(dim.map)==0)){
    tmp<-dim.map
    tmp[lengths(tmp)==0]<-NA
    lens<-lengths(tmp)
    tmp<-unlist(tmp)
    out.dim<-dim(x)[tmp]
    out.dim[is.na(tmp)]<-1
    out.dimnames<-dimnames(x)[tmp]
    names(out.dimnames)<-rep(names(dim.map),lens)
    ndims<-length(tmp)
    if(lengths(dim.map)[1]==0){
      out.dimnames[1]<-list(attr(x,type.code))
    }
    if(lengths(dim.map)[2]==0){
      out.dimnames[-c(1,ndims)]<-list(attr(x,'parameters'))
    }
    if(lengths(dim.map)[3]==0){
      out.dimnames[ndims]<-list(attr(x,'chains'))
    }
    if(any(is.na(names(dimnames(x))))){
      x<-array(x,out.dim,out.dimnames)
    }else{
      x<-array(aperm(x,unlist(dim.map)),out.dim,out.dimnames)
    }
  }else{
    dim.map<-unlist(dim.map)
    if(any(dim.map!=1:length(dim.map))){
      x<-array(aperm(x,dim.map),dim(x)[dim.map],dimnames(x)[dim.map])
    }
  }
  attr(x,'param_type')<-par.type
  x<-.add.par.class(x)
  if(simplify){
    if(is.null(dimnames(x)[[1]])&dim(x)[1]==1){
      x<-array(x,dim(x)[-1],dimnames(x)[-1])
    }
  }
  x
}

#opposite of .expand.par
#collapse any dimensions of param_block of length 1, storing dimension name info as attributes
#preserves dimensions of length 0 now--though there's no guaruntee such blocks won't cause issues!!!
.simplify.par<-function(x){
  old.dims<-dim(x)
  if(is.null(old.dims)){
    x
  }else if(any(old.dims==1)){
    old.dimnames<-dimnames(x)
    inds<-old.dims!=1
    if(all(!inds)){
      out<-vector(mode(x),1)
      out[]<-x
    }else{
      out<-array(x,old.dims[inds],old.dimnames[inds])
    }
    if(length(dim(out))==1){
      out<-setNames(as.vector(out),dimnames(out)[[1]])
    }
    for(i in seq_along(old.dims)){
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
    for(i in c('quantiles','diagnostics','parameters','chains','param_type')){
      if(!is.null(attr(x,i))){
        attr(out,i)<-attr(x,i)
      }
    }
    out<-.add.par.class(out)
    out
  }else{
    x
  }
}

.make.par.3D<-function(x){
  if(length(dim(x))<3){
    x<-.expand.par(x)
  }
  if(length(dim(x))>3){
    new.dim<-dim(x)
    new.dimnames<-dimnames(x)
    param.dims<-which(names(dimnames(x))=='parameters')
    param.names<-dimnames(x)[param.dims]
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
    x<-array(x,new.dim,new.dimnames)
  }
  x
}

#make return one code if iterations are messed up, another if chains are messed up --> done
#returns a list of expanded param_blocks, their dimensions, their dimension names, whether they match, and whether 1st/3rd dim match
.compatible.dims.check<-function(...){
  exps<-lapply(list(...),.expand.par)
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

.is.par<-function(x){
  check.vec<-c('iterations','quantiles','diagnostics','parameters','chains')
  any(names(dimnames(x))%in%check.vec)|any(names(attributes(x))%in%check.vec)
}

.get.par.type<-function(x){
  attr(x,'param_type')
}

.strip.par.class<-function(x){
  tmp<-class(x)
  class(x)<-tmp[!which(tmp=='param_block')]
  x
}

.add.par.class<-function(x){
  tmp<-class(x)
  class(x)<-c('param_block',tmp)
  x
}
