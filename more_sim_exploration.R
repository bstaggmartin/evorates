rm(list=ls())
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
res<-readRDS('sim.data.frame')
res$true.Rs<-lapply(sims,'[[','R')
colnames(res)[colnames(res)=='R_mu']<-'Rmu'
res$Rsig2.ac<-res$media.Rsig2-res$Rsig2
res$Rmu.ac<-res$media.Rmu-res$Rmu
res$Rsig2.pr<-res$upper.Rsig2-res$lower.Rsig2
res$Rmu.pr<-res$upper.Rmu-res$lower.Rmu
res$Rsig2.co<-ifelse(res$Rsig2>res$lower.Rsig2&res$Rsig2<res$upper.Rsig2,TRUE,FALSE)
res$Rmu.co<-ifelse(res$Rmu>res$lower.Rmu&res$Rmu<res$upper.Rmu,TRUE,FALSE)
res$all<-factor(apply(res[,c('Rsig2','Rmu','size')],1,paste,collapse=','))
plot(Rsig2.ac~all,data=res)
abline(h=0)
plot(Rmu.ac~all,data=res)
abline(h=0)
plot(Rsig2.pr~all,data=res)
abline(h=0)
plot(Rmu.pr~all,data=res)
abline(h=0)
tapply(res$Rsig2.co,list(res$Rsig2,res$Rmu,res$size),sum)
tapply(res$Rmu.co,list(res$Rsig2,res$Rmu,res$size),sum)
#coverage of 80 to 90%, right around what we want
#mean's literally 94% in case of Rsig2 (excluding Rsig2 = 0) and 95% in case of Rmu
#need to beautify/summarize all this better

get.R.params<-function(x){
  est.R<-res[[x,'res.com']]%quantiles%list('R_[1-9]',c(0.025,0.5,0.975))
  sim.R<-res[[x,'true.Rs']]
  R.ac<-est.R[2,]-sim.R
  R.pr<-apply(est.R[-2,],2,diff)
  R.co<-ifelse(sim.R>est.R[1,]&sim.R<est.R[3,],TRUE,FALSE)
  R.cor<-cor(est.R[2,],sim.R)
  sapply(list(R.ac,R.pr,R.co,R.cor),mean)
}
tmp<-sapply(1:nrow(res),get.R.params)
rownames(tmp)<-c('R.ac','R.pr','R.co','R.cor')
res<-cbind(res,t(tmp))
plot(R.cor~Rmu,data=res)
for(i in which(res$R.cor<0.6)){
  plot(res[[i,'res.com']]%means%'R_[1-9]'~sims[[i]]$R,
       xlim=range(sims[[i]]$R),ylim=range(res[[i,'res.com']]%chains%'R_[1-9]'))
  segments(x0=sims[[i]]$R,
           y0=res[[i,'res.com']]%quantiles%list('R_[1-9]',0.025),
           y1=res[[i,'res.com']]%quantiles%list('R_[1-9]',0.975))
  abline(0,1)
}
#all cases of over-conservativism

#let's try something other than correlations
get.R.params<-function(x){
  est.R<-res[[x,'res.com']]%quantiles%list('R_[1-9]',c(0.025,0.5,0.975))
  sim.R<-res[[x,'true.Rs']]
  R.ac<-est.R[2,]-sim.R
  R.pr<-apply(est.R[-2,],2,diff)
  R.co<-ifelse(sim.R>est.R[1,]&sim.R<est.R[3,],TRUE,FALSE)
  tmp<-est.R[2,]-mean(sim.R)
  R.reg<-lm(tmp~scale(sim.R,scale=FALSE))
  sapply(list(R.ac,R.pr,R.co,coef(R.reg)[1],coef(R.reg)[2]),mean)
}
tmp<-sapply(1:nrow(res),get.R.params)
rownames(tmp)<-c('R.ac','R.pr','R.co','R.int','R.slope')
hist(tmp[5,],breaks=100)
hist(tmp[4,],breaks=100)
#hmmm...interesting need to think about centering more. Some slopes are pretty high
for(i in which(tmp[5,]>1.5)){
  plot(res[[i,'res.com']]%means%'R_[1-9]'~sims[[i]]$R,
       xlim=range(sims[[i]]$R),ylim=range(res[[i,'res.com']]%chains%'R_[1-9]'))
  segments(x0=sims[[i]]$R,
           y0=res[[i,'res.com']]%quantiles%list('R_[1-9]',0.025),
           y1=res[[i,'res.com']]%quantiles%list('R_[1-9]',0.975))
  abline(0,1)
}
for(i in which(tmp[5,]<0.5)){
  plot(res[[i,'res.com']]%means%'R_[1-9]'~sims[[i]]$R,
       xlim=range(sims[[i]]$R),ylim=range(res[[i,'res.com']]%chains%'R_[1-9]'))
  segments(x0=sims[[i]]$R,
           y0=res[[i,'res.com']]%quantiles%list('R_[1-9]',0.025),
           y1=res[[i,'res.com']]%quantiles%list('R_[1-9]',0.975))
  abline(0,1)
}
#all seem to be cases where magnitude of R_mu was overestimated

#fixed centering --> averages look fine

#check proportion conclusions R_sig2 and R_mu are significant

res$Rsig2.sig<-res$sd>3
tapply(res$Rsig2.sig,list(res$Rsig2,res$Rmu,res$size),sum)
res$Rmu.neg<-res$upper.Rmu<0
res$Rmu.pos<-res$lower.Rmu>0
tapply(res$Rmu.neg,list(res$Rsig2,res$Rmu,res$size),sum)
tapply(res$Rmu.pos,list(res$Rsig2,res$Rmu,res$size),sum)

twoway.plot<-function(x1,x2,y,x1lab=NULL,x2lab=NULL,ylab=NULL,col=palette(),jitter=0,
                      data=NULL,add=FALSE,
                      ...,
                      legend.args=list(),
                      density=FALSE){
  if(is.null(x1lab)) x1lab<-deparse(substitute(x1))
  if(is.null(x2lab)) x2lab<-deparse(substitute(x2))
  if(is.null(ylab)) ylab<-deparse(substitute(y))
  if(!is.null(data)){
    attach(data)
  }
  x1<-factor(x1)
  x2<-factor(x2)
  x1.levs<-levels(x1)
  x2.levs<-levels(x2)
  tmp<-matrix(1,length(x1.levs),length(x2.levs))
  tmp<-barplot(tmp,beside=TRUE,plot=FALSE,...)
  x12<-tmp[x1:x2]
  if(density){
    dens<-tapply(y,x1:x2,function(ii) do.call('density',
                                                    c(list(ii),...)))
    width<-diff(tmp[1:2,1])/4
    if(!add){
      plot(y~x12,xaxt='n',xlab=x1lab,ylab=ylab,col=rgb(0,0,0,0),...)
    }
    for(i in names(dens)){
      xx<-dens[[i]]$y
      xx<-c(0,xx/max(xx)*width,0)
      xx<-tmp[match(i,levels(x1:x2))]+rep(c(-1,1),each=length(xx))*c(xx,rev(xx))
      yy<-dens[[i]]$x
      yy<-c(yy[1],yy,yy[length(yy)])
      yy<-c(yy,rev(yy))
      polygon(yy~xx,col=col[match(strsplit(i,':')[[1]][2],levels(x2))],...)
    }
  }else{
    x12<-jitter(x12,factor=jitter)
    if(!add){
      plot(y~x12,xaxt='n',xlab=x1lab,ylab=ylab,col=col[x2],...)
    }else{
      points(y~x12,col=col[x2],...)
    }
  }
  tmp.dots<-list(...)[names(list(...))!='lwd']
  axis(1,colMeans(tmp),x1.levs,tmp.dots)
  if(is.null(legend.args$x)) legend.args$x<-'topleft'
  if(density){
    args.list<-setNames(vector('list',3),c('border','density','angle'))
    dots<-names(list(...))
    for(i in names(args.list)){
      if(!is.na(match(i,dots))){
        args.list[[i]]<-list(...)[[i]]
      }
    }
    do.call(legend,
            c(legend.args,args.list,
              fill=list(col[1:length(x2.levs)]),
              legend=list(x2.levs),
              title=list(x2lab)))
  }else{
    if(hasArg(pch)){
      pch<-list(...)$pch
    }else{
      pch<-par()$pch
    }
    do.call(legend,
            c(legend.args,
              pch=list(pch),
              col=list(col[1:length(x2.levs)]),
              legend=list(x2.levs),
              title=list(as.expression(x2lab))))
  }
  detach(data)
  tmp
}
twoway.plot(Rmu,Rsig2,media.Rmu,pch=16,data=res,jitter=1,space=c(0,0.5),
            legend.args=list(bty='n'),col=c('cornflowerblue','gray','indianred'),
            density=TRUE,lwd=2)

tmp.df<-data.frame(R_mu=rep(res$Rmu,each=1200),
                   R_sig2=rep(res$Rsig2,each=1200),
                   size=rep(res$size,each=1200),
                   R_sig2_est=unlist(lapply(res$res.com,function(ii) ii%chains%R_sig2)),
                   R_mu_est=unlist(lapply(res$res.com,function(ii) ii%chains%R_mu)))

tmp.df.100<-tmp.df[tmp.df$size==100,]
test<-twoway.plot(R_mu,R_sig2,R_mu_est,pch=16,data=tmp.df.100,space=c(0,0.5),
                  legend.args=list(bty='n'),col=c('cornflowerblue','gray','indianred'),
                  density=TRUE,border=NA)
segments(x0=colMeans(test)-1.2*diff(test[1:2,1]),x1=colMeans(test)+1.2*diff(test[1:2,1]),
         y0=unique(res$Rmu),
         lwd=4,lty=2,
         col=rgb(0,0,0,0.5))

#it may be better to use a normal prior for scale params...hmmm...
#definitely some overegularization, I think...

cols<-alter.cols(c('cornflowerblue','gray','indianred'),alpha=1)
test<-twoway.plot(Rsig2,Rmu,media.Rsig2,pch=16,data=res,space=c(0,0.5),
                  legend.args=list(bty='n',horiz=TRUE),col=rep('white',3),
                  density=FALSE,border=NA,bw=1,
                  ylim=c(0,20),
                  jitter=1,
                  x1lab=bquote('Simulated'*~'R'[sigma^2]),
                  x2lab=bquote('Simulated'*~'R'[mu]),
                  ylab=bquote('Estimated'*~'R'[sigma^2]))
segments(x0=test,
         y0=t(tapply(tmp.df$R_sig2_est,list(tmp.df$R_sig2,tmp.df$R_mu),function(ii) quantile(ii,0.25))),
         y1=t(tapply(tmp.df$R_sig2_est,list(tmp.df$R_sig2,tmp.df$R_mu),function(ii) quantile(ii,0.75))),
         lwd=8,
         col=cols)
segments(x0=test,
         y0=t(tapply(tmp.df$R_sig2_est,list(tmp.df$R_sig2,tmp.df$R_mu),function(ii) quantile(ii,0.025))),
         y1=t(tapply(tmp.df$R_sig2_est,list(tmp.df$R_sig2,tmp.df$R_mu),function(ii) quantile(ii,0.975))),
         lwd=2,
         col=cols)
test<-twoway.plot(Rsig2,Rmu,media.Rsig2,pch=16,data=res,space=c(0,0.5),
                  legend.args=list(bty='n',horiz=TRUE),col=alter.cols(cols,mod.val=-0.2),
                  density=FALSE,border=NA,bw=1,
                  ylim=c(0,20),add=TRUE,
                  jitter=1,
                  x1lab=bquote('Simulated'*~'R'[sigma^2]),
                  x2lab=bquote('Simulated'*~'R'[mu]))
segments(x0=colMeans(test)-1.2*diff(test[1:2,1]),x1=colMeans(test)+1.2*diff(test[1:2,1]),
         y0=unique(res$Rsig2),
         lwd=4,lty=2)

test<-twoway.plot(Rmu,Rsig2,media.Rmu,pch=16,data=res,space=c(0,0.5),
                  legend.args=list(bty='n',horiz=TRUE),col=rep('white',3),
                  density=FALSE,border=NA,bw=1,
                  ylim=c(-15,15),
                  jitter=1,
                  x2lab=bquote('Simulated'*~'R'[sigma^2]),
                  x1lab=bquote('Simulated'*~'R'[mu]),
                  ylab=bquote('Estimated'*~'R'[mu]))
segments(x0=test,
         y0=t(tapply(tmp.df$R_mu_est,list(tmp.df$R_mu,tmp.df$R_sig2),function(ii) quantile(ii,0.25))),
         y1=t(tapply(tmp.df$R_mu_est,list(tmp.df$R_mu,tmp.df$R_sig2),function(ii) quantile(ii,0.75))),
         lwd=8,
         col=c('cornflowerblue','gray','indianred'))
segments(x0=test,
         y0=t(tapply(tmp.df$R_mu_est,list(tmp.df$R_mu,tmp.df$R_sig2),function(ii) quantile(ii,0.025))),
         y1=t(tapply(tmp.df$R_mu_est,list(tmp.df$R_mu,tmp.df$R_sig2),function(ii) quantile(ii,0.975))),
         lwd=2,
         col=c('cornflowerblue','gray','indianred'))
test<-twoway.plot(Rmu,Rsig2,media.Rmu,pch=16,data=res,space=c(0,0.5),
                  legend.args=list(bty='n',horiz=TRUE),col=alter.cols(cols,mod.val=-0.2),
                  density=FALSE,border=NA,bw=1,
                  ylim=c(-15,15),
                  jitter=1,
                  add=TRUE,
                  x2lab=bquote('Simulated'*~'R'[sigma^2]),
                  x1lab=bquote('Simulated'*~'R'[mu]),
                  ylab=bquote('Estimated'*~'R'[mu]))
segments(x0=colMeans(test)-1.2*diff(test[1:2,1]),x1=colMeans(test)+1.2*diff(test[1:2,1]),
         y0=unique(res$Rmu),
         lwd=4,lty=2)

