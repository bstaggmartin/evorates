#defunct covariance matrix fun
#need more stuff for handling commas in trait names, excluding all/partial NA rows/columns, etc.
#also separating intra/evocov, if necessary
#not a bad function, but actually, see a better alternative below...
# as.cov.mat<-function(arr,exclude.part.NAs=F,na.rm='none'){
#   arr<-.expand.element(arr)
#   dummy.fit<-list(chains=arr)
#   entries<-dimnames(arr)[[2]]
#   entries<-gsub('_evocov$','',entries)
#   entries<-gsub('_intracov$','',entries)
#   entries<-strsplit(entries,',')
#   parsed.entries<-unique(unlist(entries))
#   k<-length(parsed.entries)
#   param.names<-list(parameters=parsed.entries)
#   out<-array(NA,c(k,k,dim(arr)[1],dim(arr)[3]),c(param.names,param.names,dimnames(arr)[c(1,3)]))
#   for(i in 1:k){
#     for(j in 1:i){
#       tmp<-try(.int.chains(dummy.fit,paste(parsed.entries[i],parsed.entries[j],sep=',')),silent=T)
#       if(!inherits(tmp,'try-error')){
#         out[i,j,,]<-tmp
#         if(i!=j){
#           out[j,i,,]<-tmp
#         }
#       }
#     }
#   }
#   all.NAs<-which(sapply(1:k,function(ii) all(is.na(out[ii,,,]))))
#   if(length(all.NAs)==k){
#     stop('covariance matrices are all just NAs: double-check input')
#   }
#   na.rm<-match(na.rm,c('none','part','full'))
#   if(na.rm>1){
#     if(length(all.NAs)==0){
#       all.NAs<-k+1
#     }
#     out<-.index.element(out,rep(list(all.NAs),2),c(1,2),T)
#     k<-dim(out)[1]
#     if(na.rm>2){
#       any.NAs<-which(sapply(1:k,function(ii) any(is.na(out[ii,,,]))))
#       if(length(any.NAs)==k){
#         stop("all rows/columns of covariance matrices include NAs: double-check input or set na.rm to 'part'")
#       }
#       if(length(any.NAs)==0){
#         any.NAs<-k+1
#       }
#       out<-.index.element(out,rep(list(any.NAs),2),c(1,2),T)
#     }
#   }
#   if(output.list){
#     out<-asplit(out,4)
#     out<-lapply(out,asplit,MARGIN=3)
#     if(length(out)==1){
#       attr(out[[1]],'chains')<-names(out)
#       out<-out[[1]]
#     }
#   }else{
#     out<-.simplify.element(out)
#   }
#   out
# }