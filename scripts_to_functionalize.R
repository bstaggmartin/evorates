library(phytools)
tree<-pbtree(n=50,scale=5)
test<-gen.corateBM(tree,Rsig2=0.8)
tree<-keep.tip(whales$phy,rownames(whales$dat))
X<-whales$dat[,'lnlength']
plot(test,tree)
test.fit<-fit.corateBM(tree,X,optimize = T)
#sample a random branch and plot the trace
branch<-sample(nrow(tree$edge),1)
plot(test.fit%chains%branch,type='l',main=paste('ln(rate) for branch',branch))
abline(h=test$R[branch],col='red')

#look at overall correlation
plot(0,type='n',xlim=range(test$R),ylim=range(test.fit$chains),
     xlab='true ln(rate)',ylab='ln(rate) posterior distribution')
for(i in 1:nrow(tree$edge)){
  dens<-density(test.fit%chains%i,bw=0.2)
  yy<-c(dens$x,rev(dens$x))
  xx<-0.05*c(dens$y,-1*rev(dens$y))+test$R[i]
  polygon(yy~xx,border=NA,col='gray')
}
points(as.vector(apply(test.fit%chains%1:nrow(tree$edge),2,median))~test$R,pch=16)
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
tmp<-list(R=apply(test.fit%chains%'R[',2,mean));class(tmp)<-'corateBM'
plot(tmp,tree,phylogram=F,edge.width=10,cex=0.7)
coords<-get("last_plot.phylo",envir=.PlotPhyloEnv)
#trans<-abs(test.fit$post.probs[,1]-0.5)/0.5
trans<-test.fit$post.probs
segments(y0=coords$yy[as.vector(t(tree$edge))],
         y1=coords$yy[rep(tree$edge[,2],each=2)],
         x0=coords$xx[rep(tree$edge[,1],each=2)],
         x1=coords$xx[as.vector(t(tree$edge))],
         col=rep(rgb(trans,trans,trans),each=2),
         lwd=4)

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