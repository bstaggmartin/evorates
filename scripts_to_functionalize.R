library(contSimmap)
library(phytools)
tree<-pbtree(n=50,scale=5)
gen.eb<-function(R0,b,tree){
  mat<-mrca(tree)
  mat[]<-R0*(exp(b*node.depth.edgelength(tree)[mat])-1)/b
  mat
}
X<-mvnfast::rmvn(1,rep(0,50),gen.eb(1,1.5,tree))
names(X)<-tree$tip.label
X<-rnorm(length(X)*10,rep(X,each=10),3)
names(X)<-rep(tree$tip.label,each=10)
phenogram(tree,X,ftype='off')
test<-gen.corateBM(tree,Rsig2=0.4,intra.var=T,n_obs=rpois(length(tree$tip.label),10)+1,Xsig2=0.5,Rmu=-0.5)
#tree<-keep.tip(whales$phy,rownames(whales$dat))
#X<-whales$dat[,'lnlength']
plot(test,tree)
#test.fit<-fit.corateBM(tree,test$X_obs,intra.var=T,chains=1,iter=8000,control=list(adapt_delta=0.9),thin=10,warmup=2000)
test.fit<-fit.corateBM(tree,test$X_obs,intra.var=T,chains=2,trend=T)
#sample a random branch and plot the trace
branch<-sample(nrow(tree$edge),1)
plot((test.fit%chains%branch)[3:1000]~(test.fit%chains%branch)[1:998])
spacer<-3
hist(sapply(1:nrow(tree$edge),function(ii) cor((test.fit%chains%ii)[(1+spacer):1000],(test.fit%chains%ii)[1:(1000-spacer)])))
plot(test.fit%chains%branch,type='l',main=paste('ln(rate) for branch',branch))
abline(h=test$R[branch],col='red')

#look at overall correlation
plot(0,type='n',xlim=range(test$R),ylim=range(test.fit$chains),
     xlab='true ln(rate)',ylab='ln(rate) posterior distribution')
for(i in 1:nrow(tree$edge)){
  dens<-density(test.fit%chains%i,bw=0.2)
  yy<-c(dens$x,rev(dens$x))
  xx<-0.1*c(dens$y,-1*rev(dens$y))+test$R[i]
  polygon(yy~xx,border=NA,col='gray')
}
points(as.vector(apply(test.fit%chains%1:nrow(tree$edge),2,median))~test$R,pch=16)
abline(0,1,col='red')

plot(0,type='n',xlim=range(test$X),ylim=range(test.fit%chains%'^t'),
     xlab='true ln(rate)',ylab='ln(rate) posterior distribution')
for(i in tree$tip.label){
  dens<-density(test.fit%chains%paste(i,'$',sep=''),bw=0.2)
  yy<-c(dens$x,rev(dens$x))
  xx<-0.05*c(dens$y,-1*rev(dens$y))+test$X[i]
  polygon(yy~xx,border=NA,col='gray')
}
points(test.fit$MAPs[tree$tip.label,]~test$X,pch=16)
abline(0,1,col='red')

branch<-sample(nrow(tree$edge),2)
plot(test.fit%chains%branch[1]-test.fit%chains%branch[2],type='l',
     main=paste('ln(rate) diff for branch',branch[1],'and',branch[2]))
abline(h=test$R[branch[1]]-test$R[branch[2]],col='red')

hist(test.fit%chains%'Rsig2',breaks=40)

tol.el<-sum(tree$edge.length)
wgts<-tree$edge.length/tol.el
plot(0,type='n',xlim=c(1,nrow(tree$edge)),ylim=range(test.fit%chains%'dev'),
     xlab='branch',ylab='ln(rate) posterior distribution')
nrs<-dim(test.fit$quantiles)[1]
chain<-1
for(i in 1:nrow(tree$edge)){
  tmp<-grep(paste('\\[',i,'\\] dev',sep=''),dimnames(test.fit$chains)[[2]])
  segments(x0=i,
           y0=c(min(test.fit$chains[,tmp,chain]),test.fit$quantiles[c(1,nrs),tmp,chain]),
           y1=c(test.fit$quantiles[c(1,nrs),tmp,chain],max(test.fit$chains[,tmp,chain])),
           lty=c(2,1,2),col=c('gray','black','gray'))
}
t.bg.rate<-sum(test$R*wgts)
points(test$R-t.bg.rate,col='blue',pch=16)
abline(h=0,col='red',lty=2)
plot(tree,edge.width=10,cex=0.8)
tmp<-list(R=test.fit%means%'R[');class(tmp)<-'corateBM'
plot(tmp,tree,phylogram=F,edge.width=10,show.tip.label=F,val.range=range(test$R))
plot(tmp,tree,phylogram=F,edge.width=10)
plot(test,tree,phylogram=F,edge.width=10,show.tip.label=F,val.range=range(test$R))
coords<-get("last_plot.phylo",envir=.PlotPhyloEnv)
#trans<-abs(test.fit$post.probs[,1]-0.5)/0.5
trans<-test.fit$post.probs
segments(y0=coords$yy[as.vector(t(tree$edge))],
         y1=coords$yy[rep(tree$edge[,2],each=2)],
         x0=coords$xx[rep(tree$edge[,1],each=2)],
         x1=coords$xx[as.vector(t(tree$edge))],
         col=rep(rgb(trans,trans,trans),each=2),
         lwd=3)
plot(test,tree)
tmp$X<-tapply(X,names(X),mean)
plot(tmp,tree)
plot(test$R~apply(matrix(node.depth.edgelength(tree)[tree$edge],ncol=2),1,mean))
points(tmp$R~apply(matrix(node.depth.edgelength(tree)[tree$edge],ncol=2),1,mean),pch=16)

mat<-matrix(nrow=dim(rates)[3],ncol=dim(rates)[3])
simmat<-mat
for(i in 2:dim(rates)[3]){
  for(j in 1:(i-1)){
    mat[i,j]<-sum((rates[,1,j]-rates[,1,i])>0)/dim(rates)[1]
    simmat[i,j]<-test$R[j]-test$R[i]
  }
}
mat[upper.tri(mat)]<- 1-t(mat)[upper.tri(mat)]
simmat[upper.tri(simmat)]<- -t(simmat)[upper.tri(simmat)]
image(simmat)
image(mat)
plot(as.vector(mat)~as.vector(simmat),ylim=c(0,1),col=c('black','red')[(as.vector(mat)<0.025|as.vector(mat)>0.975)+1])
abline(h=c(0.025,0.975),col='red',lty=2)
sum(as.vector(mat)<0.025&as.vector(simmat)>0,na.rm=T)/length(!is.na(mat))
#~5% type I error!


par(mfrow=c(1,2))
plot(test,tree,phylogram=F,edge.width=3,val.range=range(test.fit$chains[,grep(paste('R\\[\\d+\\]$',sep=''),rownames(test.fit$MAPs)),1]))
plot(tmp,tree,phylogram=F,edge.width=3,val.range=range(test.fit$chains[,grep(paste('R\\[\\d+\\]$',sep=''),rownames(test.fit$MAPs)),1]))


hmmm<-apply(test.fit$quantiles,2,function(ii) abs(sum(diff(ii))))
par(mfrow=c(1,1))
plot(hmmm[grep(paste('R\\[\\d+\\]$',sep=''),rownames(test.fit$MAPs))]~tree$edge.length,col=rgb(1-wgts/max(wgts),1-wgts/max(wgts),1-wgts/max(wgts)),pch=16)
#shortest branches towards tips have the most precise estimates, but honestly the diff isn't great
#it seems you might actually estimate deep branches with more precision since they have more info
#'flowing' into them from the tips...





tol.el<-sum(tree$edge.length)
wgts<-tree$edge.length/tol.el
hist(apply(rates[,1,],1,function(ii) sum(ii*wgts)))
#quick function to find all edges in a monophyletic clade
find.subtree<-function(tree,mrca){
  involved.edges<-which(tree$edge[,1]==mrca)
  indices<-1:length(involved.edges)
  while(length(indices)>0){
    start<-length(involved.edges)+1
    involved.edges<-c(involved.edges,unlist(sapply(tree$edge[involved.edges[indices],2],function(ii) which(tree$edge[,1]==ii))))
    if(start==length(involved.edges)+1){
      break
    }
    indices<-start:length(involved.edges)
  }
  involved.edges
}
sub.tol.el<-sum(tree$edge.length[find.subtree(tree,189)])
sub.wgts<-tree$edge.length[find.subtree(tree,189)]/sub.tol.el
hist(apply(rates[,1,find.subtree(tree,189)],1,function(ii) sum(ii*sub.wgts)))


#comparing edge rates to the background rate to look for significant differences is also possible...
bg.rate<-apply(rates[,1,],1,function(ii) sum(ii*wgts))
pp<-apply(rates[,1,],2,function(ii) sum(((ii)-bg.rate)>0)/dim(rates)[1])
plot(pp,ylim=c(0,1),
     xlab='branches',xaxt='n',ylab='post prob branch rate < avg rate',
     col=c('black','red')[(pp<0.025|pp>0.975)+1])
abline(h=c(0.025,0.975),col='red',lty=2)

plot(test.fit$post.probs,ylim=c(0,1),
     xlab='branches',xaxt='n',ylab='post prob branch rate < avg rate',
     col=c('black','red')[(test.fit$post.probs<0.025|test.fit$post.probs>0.975)+1])
abline(h=c(0.025,0.975),col='red',lty=2)
plot(test.fit$post.probs~test.fit%MAPs%'dev',col=c('black','red')[(test.fit$post.probs<0.025|test.fit$post.probs>0.975)+1])
abline(h=c(0.025,0.975),col='red',lty=2)
sum((test.fit$post.probs<0.1&(test.fit%MAPs%'dev')>0)|(test.fit$post.probs>0.9&(test.fit%MAPs%'dev')<0),na.rm=T)/length(test.fit$post.probs)
#20% alpha results in type II of 0! (for 500 tip test)

sig.mat<-matrix(NA,nrow(tree$edge),nrow(tree$edge))
emp.mat<-sig.mat
len<-length(test.fit%chains%1)
for(i in 1:nrow(sig.mat)){
  for(j in i:ncol(sig.mat)){
    if(i==j){
      next
    }else{
      sig.mat[i,j]<-sum((test.fit%chains%i - test.fit%chains%j)>0)/len
      emp.mat[i,j]<-test$R[i]-test$R[j]
    }
    cat(i,', ',j,'\n',sep='')
  }
}
sig.mat[lower.tri(sig.mat)]<-NA
plot(as.vector(sig.mat)~as.vector(emp.mat),col=c('black','red')[(as.vector(sig.mat)<0.1|as.vector(sig.mat)>0.9)+1])
abline(h=c(0.1,0.9),col='red',lty=2)
sum((as.vector(sig.mat)<0.1&as.vector(emp.mat)>0)|(as.vector(sig.mat)>0.9&as.vector(emp.mat)<0),na.rm=T)/sum(!is.na(sig.mat))
cbind(which(sig.mat<0.1)%/%nrow(tree$edge),which(sig.mat<0.1)%%nrow(tree$edge))
#type II error is only 0.05% for alpha=5%...
#jumps to ~0.3% for alpha=10% and ~1% for alpha=20% (for 500 tip test)

tol.el<-sum(tree$edge.length)
wgts<-tree$edge.length/tol.el
plot(0,type='n',xlim=c(1,nrow(tree$edge)),ylim=range(test.fit$chains[,grep('dev',dimnames(test.fit$chains)[[2]]),1]),xlab='branch',ylab='ln(rate) posterior distribution')
nrs<-dim(test.fit$quantiles)[1]
chain<-1
for(i in 1:nrow(tree$edge)){
  tmp<-grep(paste('\\[',i,'\\] dev',sep=''),dimnames(test.fit$chains)[[2]])
  segments(x0=i,
           y0=c(min(test.fit$chains[,tmp,chain]),test.fit$quantiles[c(1,nrs),tmp,chain]),
           y1=c(test.fit$quantiles[c(1,nrs),tmp,chain],max(test.fit$chains[,tmp,chain])),
           lty=c(2,1,2),col=c('gray','black','gray'))
}
t.bg.rate<-sum(test$R*wgts)
points(test$R-t.bg.rate,col='blue',pch=16)
abline(h=0,col='red',lty=2)

#DUDE! 5 minutes to do 4 chains for a 100 tip tree

#5/25: weird phenomenon where if you have accelarating rates and low intraspecific sampling, modeling
#intraspecific variance creates big problems: some tip trait values start to 'jump' around range of
#extant variation! Best way to solve problem is increasing intraspecific sampling; narrowing the prior
#did little to stop it. Should see if similar thing happen with declining rates, but it seems fine, from
#what I've seen...

#5/26: hard to separately estimate Rmu and Rsig2 well...trends are detectable, it seems, but it's hard
#to detect their actual magnitude
#difficulty of reconstructing ancestral states...

#savage-dickey calc
dens<-density(test.fit%chains%'Rsig2')
dens<-density(test.fit%chains%'Rmu')
dcauchy(0,0,20)/dens$y[which.min(abs(dens$x-0))]
