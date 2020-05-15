library(phytools)
tree<-pbtree(n=50,scale=5)
test<-gen.corateBM(tree,Rsig2=10)
plot(test,tree)
test.fit<-fit.corateBM(tree,test$X,chains=2)

rates<-test.fit$R

#sample a random branch and plot the trace
branch<-sample(dim(rates)[3],1)
plot(rates[,1,branch],type='l',main=paste('ln(rate) for branch',branch))
abline(h=test$R[branch],col='red')

#look at overall correlation
plot(0,type='n',xlim=range(test$R),ylim=range(test.fit$R),xlab='true ln(rate)',ylab='ln(rate) posterior distribution')
for(i in 1:(dim(test.fit$R)[3])){
  dens<-density(test.fit$R[,1,i],bw=0.2)
  yy<-c(dens$x,rev(dens$x))
  xx<-0.05*c(dens$y,-1*rev(dens$y))+test$R[i]
  polygon(yy~xx,border=NA,col='gray')
}
points(as.vector(apply(test.fit$R,c(2,3),median))~test$R,pch=16)
abline(0,1,col='red')

branch<-sample(dim(rates)[3],2)
plot(rates[,1,branch[1]]-rates[,1,branch[2]],type='l',
     main=paste('ln(rate) diff for branch',branch[1],'and',branch[2]))
abline(h=test$R[branch[1]]-test$R[branch[2]],col='red')
hist(test.fit$Rsig2,breaks=40)


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

tmp<-list(R=as.vector(apply(test.fit$R,c(2,3),median)));class(tmp)<-'corateBM'
plot(test,tree,phylogram=F,edge.width=3,val.range=range(tmp$R))
plot(tmp,tree,phylogram=F,edge.width=3)










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

plot(0,type='n',xlim=c(1,nrow(tree$edge)),ylim=range(test.fit$chains$rate.devs),xlab='branch',ylab='ln(rate) posterior distribution')
for(i in 1:nrow(tree$edge)){
  segments(x0=i,y0=c(min(test.fit$chains$rate.devs[,1,i]),test.fit$quantiles$rate.devs[c(1,dim(test.fit$quantiles$rate.dev)[1]),1,i]),
           y1=c(test.fit$quantiles$rate.devs[c(1,dim(test.fit$quantiles$rate.dev)[1]),1,i],max(test.fit$chains$rate.devs[,1,i])),lty=c(2,1,2),col=c('gray','black','gray'))
}
t.bg.rate<-sum(test$R*wgts)
points(test$R-t.bg.rate,col='blue',pch=16)
abline(h=0,col='red',lty=2)
