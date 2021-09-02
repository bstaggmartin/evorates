#' @export
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
names.loose_element<-function(x){
  code<-.param.code(x)
  out<-switch(code[[1]],
              dimnames(x)[code[[2]]],
              attr(x,'parameters'),
              names(.strip.ele(x)))
  if(is.list(out)&length(out)==1){
    out<-out[[1]]
  }
  out
}

#' @export
`names<-.loose_element`<-function(x,value){
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
    for(i in 1:dims.len){
      if(lens[i]!=length(new.dimnames[[dims[i]]])){
        stop("Provided number of names for dimension ",
             dims[i],
             " and number of parameter names for dimension ",
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
    dimnames(x)<-new.dimnames
  }else if(code[[1]]==2){
    tmp<-attr(x,'parameters')
    if(lens!=length(tmp)){
      stop("Provided number of names and length of parameter names don't match")
    }
    if(any(!nzchar(value[[1]]))){
      stop("Names of length 0 not allowed")
    }
    attr(x,'parameters')<-value[[1]]
  }else{
    if(any(!nzchar(value[[1]]))){
      stop("Names of length 0 not allowed")
    }
    x<-.strip.ele(x)
    names(x)<-value[[1]]
    x<-.add.ele(x)
  }
  x
}
