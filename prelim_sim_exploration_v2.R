library(phytools)
library(contSimmap)
setwd('~/../Desktop/prelim_evorates_sim_study')
set.seed(123)
reps<-10
Rmu.levs<-c(-4,0,4)
n.Rmu<-length(Rmu.levs)
Rsig2.levs<-c(0,3,6)
n.Rsig2<-length(Rsig2.levs)
size.levs<-c(50,100)
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
# res<-sims
# for(i in 1:reps){
#   for(j in as.character(Rmu.levs)){
#     for(k in as.character(Rsig2.levs)){
#       for(l in as.character(size.levs)){
#         res[[i,j,k,l]]<-output.evorates(paste(paste0(c('rep','Rmu','Rsig','size'),c(i,j,k,l)),collapse='__'))
#         cat(paste(paste0(c('rep','Rmu','Rsig','size'),c(i,j,k,l)),collapse='__'),'\n')
#       }
#     }
#   }
# }
# saveRDS(res,"prelim_res")
res<-readRDS("prelim_res")

check.conv<-function(x){
  if(x$sampler.control$chains==1){
    Rhats<-x%diagnostics%list('.','Rhat')
  }else{
    chains<-x%chains%.
    Rhats<-apply(chains,2,rstan::Rhat)
  }
  all(Rhats<1.05)
}

check.ess<-function(x,baseline=100){
  tmp<-x%diagnostics%list('.','ess')
  c(all(tmp[1,]>100),all(tmp[2,]>baseline))
}

#turn to data.frame
tmp<-dimnames(res)
tmp[[1]]<-1:reps
df<-do.call(expand.grid,tmp)
df[,]<-sapply(df,function(ii) as.numeric(as.character(ii)))
tmp<-vector('list',prod(dim(res)))
tmp[]<-res
df$res<-tmp
df$conv<-sapply(df$res,check.conv) #oooh, all converged!
df$res.com<-lapply(df$res,combine.chains)
ess.checks<-sapply(df$res.com,check.ess,baseline=200) #ess is wayyyy overshot, it seems! The vast majority are over 1,000
#now bulk is usually 1,000+, but tails are more like 200-300 (sufficient)

df$lower.Rsig2<-unlist(lapply(df$res.com,`%quantiles%`,list('R_sig2',0.025)))
df$upper.Rsig2<-unlist(lapply(df$res.com,`%quantiles%`,list('R_sig2',0.975)))
df$media.Rsig2<-unlist(lapply(df$res.com,`%quantiles%`,list('R_sig2',0.5)))

df$lower.Rmu<-unlist(lapply(df$res.com,`%quantiles%`,list('R_mu',0.025)))
df$upper.Rmu<-unlist(lapply(df$res.com,`%quantiles%`,list('R_mu',0.975)))
df$media.Rmu<-unlist(lapply(df$res.com,`%quantiles%`,list('R_mu',0.5)))

plot(0,col='white',xlim=range(Rsig2.levs)+c(-0.5,0.5),
     ylim=range(df[,c('lower.Rsig2','upper.Rsig2')]))
segments(x0=jitter(df$Rsig2[df$size==50]),y0=df$lower.Rsig2[df$size==50],y1=df$upper.Rsig2[df$size==50])
segments(x0=Rsig2.levs-0.75,x1=Rsig2.levs+0.75,y0=Rsig2.levs,lwd=4,col=2)

size.select<-100
pdf(paste0(size.select,'_initial_results.pdf'),width=10,height=7)
plot(0,col='white',xlim=range(Rsig2.levs)+c(-1,1),
     ylim=range(df[df$size==size.select,c('lower.Rsig2','upper.Rsig2')]),
     xlab=paste0('true R_sig2 (',size.select,' tip trees)'),ylab='95% cred. int.',
     xaxt='n')
axis(1,Rsig2.levs)
jit.seq1<-seq(-0.5,0.5,length.out=n.Rmu)
jit.seq2<-seq(-diff(jit.seq1)[1]/3,diff(jit.seq1)[1]/3,length.out=reps+2)[-c(1,reps+2)]
segments(x0=df$Rsig2[df$size==100]+
           jit.seq1[as.factor(df$R_mu[df$size==size.select])]+
           jit.seq2[as.factor(df$rep[df$size==size.select])],
         y0=df$lower.Rsig2[df$size==size.select],y1=df$upper.Rsig2[df$size==size.select],
         col=c('blue','black','green')[as.factor(df$R_mu[df$size==size.select])])
points(x=df$Rsig2[df$size==size.select]+
         jit.seq1[as.factor(df$R_mu[df$size==size.select])]+
         jit.seq2[as.factor(df$rep[df$size==size.select])],
       y=df$media.Rsig2[df$size==size.select],
       pch=16)
segments(x0=Rsig2.levs-1,x1=Rsig2.levs+1,y0=Rsig2.levs,lwd=4,col=rgb(1,0,0,0.5))
legend('topleft',lty=1,col=c('blue','black','green'),legend=c(paste('R_mu =',Rmu.levs)))

plot(0,col='white',xlim=range(Rmu.levs)+c(-1,1),
     ylim=range(df[df$size==size.select,c('lower.Rmu','upper.Rmu')]),
     xlab=paste0('true R_mu (',size.select,' tip trees)'),ylab='95% cred. int.',
     xaxt='n')
axis(1,Rmu.levs)
jit.seq1<-seq(-0.5,0.5,length.out=n.Rsig2)
jit.seq2<-seq(-diff(jit.seq1)[1]/3,diff(jit.seq1)[1]/3,length.out=reps+2)[-c(1,reps+2)]
segments(x0=df$R_mu[df$size==100]+
           jit.seq1[as.factor(df$Rsig2[df$size==size.select])]+
           jit.seq2[as.factor(df$rep[df$size==size.select])],
         y0=df$lower.Rmu[df$size==size.select],y1=df$upper.Rmu[df$size==size.select],
         col=c('blue','black','green')[as.factor(df$Rsig2[df$size==size.select])])
points(x=df$R_mu[df$size==size.select]+
         jit.seq1[as.factor(df$Rsig2[df$size==size.select])]+
         jit.seq2[as.factor(df$rep[df$size==size.select])],
       y=df$media.Rmu[df$size==size.select],
       pch=16)
segments(x0=Rmu.levs-1,x1=Rmu.levs+1,y0=Rmu.levs,lwd=4,col=rgb(1,0,0,0.5))
abline(h=0,lty=2)
legend('topleft',lty=1,col=c('blue','black','green'),legend=c(paste('R_sig2 =',Rsig2.levs)))

dev.off()

library(logspline)
foo<-function(x){
  2*dcauchy(0,0,x$call$Rsig2_prior_sig)/dlogspline(0,logspline(x%chains%R_sig2,lbound=0))
}
df$sd<-unlist(lapply(df$res.com,foo))
df$sd.cat<-cut(df$sd,breaks=c(0,1,3,10,100,max(df$sd)+1),labels=c('null','barely','substantial','strong','decisive'))
tmp<-tapply(df$sd.cat,list(df$sd.cat,df$Rsig2,df$size),length)
tmp[is.na(tmp)]<-0
pdf('savage_dickey.pdf',width=10,height=7)
barplot(tmp[,,1]/30,xlab='true R_sig2 (50 tips)',ylab='proportion',legend.text=TRUE)
barplot(tmp[,,2]/30,xlab='true R_sig2 (100 tips)',ylab='proportion',legend.text=TRUE)
dev.off()

#I think I should rerun with Rmu, R0 set to 10, Rsig2 set to 5
#Maybe just shut off tip error--seems like it's messing with things

check.divs<-function(x){
  sum(x%sampler%'div')/x$sampler.control$iter
}

#SHIT I think there's a mistake in the sim.evorates function!

#very few divergent transitions, though, much better

#export df
saveRDS(df,'sim.data.frame')
