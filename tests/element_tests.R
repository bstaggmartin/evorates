library(evorates)
fit<-output.evorates('misc_R_tests/9.10.cet')
fit2<-combine.chains(fit)
eles<-setNames(rep(list(NULL),4),c('chains','quantiles','means','diagnostics'))
selects<-list(1:3,1,list(1:3,1),list(1,1))
for(i in names(eles)){
  for(j in list(fit,fit2)){
    for(k in selects){
      eles[[i]]<-c(eles[[i]],list(do.call(paste0('%',i,'%'),list(j,k))))
    }
  }
}

test<-list(eles[[1]][[1]],eles[[1]][[2]],eles[[2]][[2]],eles[[3]][[2]],eles[[4]][[2]])

fit.evorates(fit$call$tree,fit$call$trait.data,trend=TRUE,out.file='misc_R_tests/9.10.cet',
             iter=3000,warmup=1500,cores=4)
debug(fit.evorates)

er<-edge.ranges(fit)
Rmu<-fit%chains%R_mu
el<-fit$call$tree$edge.length
Rs<-fit%chains%1:Nedge(fit)
#something weird with order of operations--fixed it!
out<-Rs-(-log(abs(Rmu))-log(el)+log(abs(exp(Rmu*er[,2])-exp(Rmu*er[,1]))))
wgts<-el/sum(el)
tmp<-sum(wgts*out)
out<-out-tmp
names(out)<-paste0('Rdev_',1:Nedge(fit))

(log(abs(Rmu))-log(el)+log(abs(exp(Rmu*er[,2])-exp(Rmu*er[,1]))))
