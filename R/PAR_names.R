.param.code<-function(x){
  try1<-dimnames(x)
  if(length(try1)>0){
    try1<-which(names(try1)=='parameters')
  }
  try2<-attr(x,'parameters')
  if(length(try1)>0){
    list(1,c(try1))
  }else if(length(try2)>0){
    list(2)
  }else{
    list(3)
  }
}

#' @export
names.param_block<-function(x){
  code<-.param.code(x)
  out<-switch(code[[1]],
              dimnames(x)[code[[2]]],
              attr(x,'parameters'),
              names(.strip.par.class(x)))
  if(is.list(out)&length(out)==1){
    out<-out[[1]]
  }
  out
}

#' @export
`names<-.param_block`<-function(x,value){
  code<-.param.code(x)
  if(!is.list(value)){
    value<-list(value)
  }
  len<-length(value)
  lens<-lengths(value)
  if(code[[1]]==1){
    new.dimnames<-dimnames(x)
    dims<-code[[2]]
    dims.len<-length(dims)
    if(len!=dims.len){
      stop("Provided number of name vectors and parameters dimensions don't match")
    }
    for(i in seq_len(dims.len)){
      tmp.len<-length(new.dimnames[[dims[i]]])
      if(tmp.len){ #may want to return a warning otherwise...but fine for now
        if(lens[i]!=tmp.len){
          stop("Provided number of names and length for dimension ",
               dims[i],
               " don't match")
        }
        if(any(!nzchar(value[[i]]))){
          stop("Provided name of length 0 for dimension ",
               dims[i],
               ": this is not allowed")
        }
        new.dimnames[[dims[i]]]<-value[[i]]
      }
    }
    dimnames(x)<-new.dimnames
  }else if(code[[1]]==2){
    tmp<-attr(x,'parameters')
    tmp.len<-length(tmp)
    if(tmp.len){
      if(lens!=length(tmp)){
        stop("Provided number of names and length of parameter names don't match")
      }
      if(any(!nzchar(value[[1]]))){
        stop("Names of length 0 not allowed")
      }
      attr(x,'parameters')<-value[[1]]
    }
  }else{
    if(length(x)!=0){
      if(any(!nzchar(value[[1]]))){
        stop("Names of length 0 not allowed")
      }
      x<-.strip.par.class(x)
      names(x)<-value[[1]]
      x<-.add.par.class(x)
    }
  }
  x
}
