library(phytools)
library(contSimmap)
setwd('/mnt/gs18/scratch/users/mart2121/prelim_evorates_sim_study')
set.seed(123)
reps<-10
Rmu.levs<-c(-4,0,4)
n.Rmu<-length(Rmu.levs)
Rsig2.levs<-c(0,3,6)
n.Rsig2<-length(Rsig2.levs)
size.levs<-c(200)
n.size<-length(size.levs)
dims<-c(reps,n.Rmu,n.Rsig2,n.size)
n<-prod(dims)
trees<-sapply(size.levs,function(ii) pbtree(n=ii,scale=1,nsim=n/n.size))
sims<-array(vector('list'),dims,
            dimnames=list(rep=NULL,R_mu=Rmu.levs,Rsig2=Rsig2.levs,size=size.levs))
args<-expand.grid(rep(Rmu.levs,each=reps),Rsig2.levs,size.levs)
sims[,,,]<-lapply(1:nrow(args),function(ii) sim.evorates(trees[[ii]],
                                                         Rmu=args[ii,1],
                                                         Rsig2=args[ii,2]))
# foo<-function(sims){
#   seed<-sample(length(sims),1)
#   plot(sims[[seed]],main=paste(paste(c('R_mu =','R_sig2 =','size ='),args[seed,]),collapse='; '))
# }
# foo(sims)
for(i in 1:reps){
  for(j in as.character(Rmu.levs)){
    for(k in as.character(Rsig2.levs)){
      for(l in as.character(size.levs)){
          fit.evorates(sims[[i,j,k,l]]$tree,
                       sims[[i,j,k,l]]$trait.data,
                       trait.se=0,
                       trend=TRUE,
                       cores=4,
                       warmup=700,
                       iter=1000,
                       return.as.obj=FALSE,
                       out.file=paste(paste0(c('rep','Rmu','Rsig','size'),c(i,j,k,l)),collapse='__'),
                       check.overwrite=FALSE)
      }
    }
  }
}
#I declare these priors fine!
#tested--the difference in bias between R0.prior.sig=3 and =5 is minimal for pretty crazy dataset with Rmu=-10 and
#Rsig2=6...so probably nothing to worry about
# fit<-fit.evorates(sims[[1,3,2,3]]$tree,
#                   sims[[1,3,2,3]]$trait.data,
#                   trait.se=0,
#                   trend=TRUE,
#                   cores=4,
#                   iter=1000)
# fit2<-fit.evorates(sims[[1,1,3,3]]$tree,
#                    sims[[1,1,3,3]]$trait.data,
#                    trait.se=0,
#                    trend=TRUE,
#                    cores=4,
#                    warmup=700,
#                    iter=1000,
#                    refresh=10)
#The above two look good! Probs good to go
