t1<-0;t2<-2;p1<-0;p2<-0;s2<-2;res<-100
ts<-seq(t1,t2,length.out=res+1)
n.bridges<-10000
mat<-matrix(NA,nrow=res+1,ncol=n.bridges)
mat[1,]<-p1;mat[res+1,]<-p2
for(t.ind in sample(2:res)){
  last.t.ind<-max(which(!is.na(mat[,1])&ts<ts[t.ind]))
  next.t.ind<-min(which(!is.na(mat[,1])&ts>ts[t.ind]))
  mat[t.ind,]<-rnorm(n.bridges,
                     (ts[t.ind]-ts[last.t.ind])/(ts[next.t.ind]-ts[last.t.ind])*(mat[next.t.ind,]-mat[last.t.ind,])+mat[last.t.ind,],
                     s2*(ts[t.ind]-ts[last.t.ind])*(ts[next.t.ind]-ts[t.ind])/(ts[next.t.ind]-ts[last.t.ind]))
}
for(n in 1:n.bridges){
  ord<-sample(2:res)
  for(t.ind in ord){
    last.t.ind<-max(which(!is.na(mat[,n])&ts<ts[t.ind]))
    next.t.ind<-min(which(!is.na(mat[,n])&ts>ts[t.ind]))
    mat[t.ind,n]<-rnorm(1,
                        (ts[t.ind]-ts[last.t.ind])/(ts[next.t.ind]-ts[last.t.ind])*(mat[next.t.ind,n]-mat[last.t.ind,n])+mat[last.t.ind,n],
                        s2*(ts[t.ind]-ts[last.t.ind])*(ts[next.t.ind]-ts[t.ind])/(ts[next.t.ind]-ts[last.t.ind]))
  }
}
ord<-sample(2:res,length(2:res))
for(t.ind in ord){
  last.t.ind<-max(which(!is.na(mat[,1])&ts<ts[t.ind]))
  next.t.ind<-min(which(!is.na(mat[,1])&ts>ts[t.ind]))
  mat[t.ind,]<-rnorm(n.bridges,
                      (ts[t.ind]-ts[last.t.ind])/(ts[next.t.ind]-ts[last.t.ind])*(mat[next.t.ind,]-mat[last.t.ind,])+mat[last.t.ind,],
                      s2*(ts[t.ind]-ts[last.t.ind])*(ts[next.t.ind]-ts[t.ind])/(ts[next.t.ind]-ts[last.t.ind]))
}
hist(apply(mat,2,mean));sd(apply(mat,2,mean))

#For a BB with var=1 and time interval=1, sd(mean) seems to be ~0.10-0.12; it seems to linearly increase with product of var and time
#interval.
matplot(mat[,sample(10000,1000)],type='l',col=rgb(0,0,0,0.1),lty=1)
yy<-s2*(ts-t1)*(t2-ts)/(t2-t1)
plot(as.vector(apply(mat,1,quantile,probs=c(0.16,0.84)))~rep(ts,each=2),ylim=c(-1,1))
lines(1*yy~ts,type='l')
lines(-1*yy~ts,type='l')


#Use cbind and sapply(x, sample(2:res)) and sapply(x,function(ii) max(which(!is.na(mat[,ii])&ts<ts[t.ind]))) to simultaneously
#sample BB with different sampling orders; should improve the pathologies you noticed in this exercise