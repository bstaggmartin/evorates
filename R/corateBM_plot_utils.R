#a more bare-bones version of the coerce to array part of check.n.proc focused on the 'chains'
#element; for use in trace.plot and profile.plot
#also works for sampler, params, and probably quantiles...
#' @export
coerce.to.array<-function(arr){
  if(length(dim(arr))==0){
    new.arr<-as.matrix(arr)
    colnames(new.arr)<-attr(arr,'parameters')
    arr<-new.arr
  }
  if(length(dim(arr))==2){
    if(sum(grepl('chains',names(dimnames(arr))[2]))!=0){
      new.arr<-array(arr,c(dim(arr)[1],1,dim(arr)[2]))
      dimnames(new.arr)<-list(iterations=dimnames(arr)[[1]],
                              parameters=attr(arr,'parameters'),
                              chains=dimnames(arr)[[2]])
      arr<-new.arr
    }else{
      new.arr<-array(arr,c(dim(arr),1))
      dimnames(new.arr)<-list(iterations=dimnames(arr)[[1]],
                              parameters=dimnames(arr)[[2]],
                              chains='chain 1')
      arr<-new.arr
    }
  }
 arr
}

#function for getting the final list of arguments lists in a call to trace.plot or profile.plot
#' @export
get.args.master<-function(chains,separate.R,separate.X,separate.dev,together,...){
  args.master<-vector(mode='list',length=dim(chains)[2])
  names(args.master)<-paste('args.',dimnames(chains)[[2]],sep='')
  if(!separate.R){
    args.master<-args.master[!grepl('^args\\.R\\[\\d+\\]$',names(args.master))]
    if(is.null(args.master$args.R)){
      args.master<-c(args.master,'args.R'=list(NULL))
    }
  }
  if(!separate.X){
    args.master<-args.master[grepl(paste('^args.R0','Rsig2','X0','bg.rate','R\\[\\d+\\]',
                                         'R\\[\\d+\\] dev','Rmu','Ysig2','R','X','dev$',
                                         collapse='',sep='$|^args.'),
                                   names(args.master))]
    if(is.null(args.master$args.X)){
      args.master<-c(args.master,'args.X'=list(NULL))
    }
  }
  if(!separate.dev){
    args.master<-args.master[!grepl('^args\\.R\\[\\d+\\] dev$',names(args.master))]
    if(is.null(args.master$args.dev)){
      args.master<-c(args.master,'args.dev'=list(NULL))
    }
  }
  for(i in names(args.master)){
    if(i%in%names(list(...))){
      args.master[[i]]<-list(...)[[i]]
      if(grepl('^args\\.R\\[\\d+\\]$',i)){
        edge.ind<-paste('args.',
                        substr(i,regexpr('\\[',i)+1,regexpr('\\]',i)-1),
                        sep='')
        if(edge.ind%in%names(list(...))){
          warning(paste('multiple argument lists matched with parameter ',substr(i,6,nchar(i)),': argument list indicated by integer (',edge.ind,') not used',sep=''))
        }
      }
    }else if(grepl('^args\\.R\\[\\d+\\]$',i)){
      edge.ind<-paste('args.',
                      substr(i,regexpr('\\[',i)+1,regexpr('\\]',i)-1),
                      sep='')
      if(edge.ind%in%names(list(...))){
        args.master[[i]]<-list(...)[[edge.ind]]
      }
    }else{
      args.master[[i]]<-NULL
    }
    def.args<-list(...)
    #def.args<-def.args[!sapply(def.args,is.list)]
    def.args<-def.args[!grepl('args\\.',names(list(...)))]
    def.args<-def.args[!(names(def.args)%in%names(args.master[[i]]))]
    args.master[[i]]<-c(args.master[[i]],def.args)
  }
  args.master
}
