library(phytools)
library(mvnfast)

tree<-pbtree(n=40)
C<-vcv(tree)


#generating random vcv matrixs


n<-10
#ang<-matrix(runif(n*(n-1),0,2*pi),ncol=n-1)
#coords<-sapply(1:nrow(ang),function(ii) c(cos(ang[ii,1]),sapply(ang[ii,],sin)))
coords<-matrix(rnorm(n^2),ncol=n)
#qr decomposes the matrix assuming it was orthogonal; qr.Q then reconstructs the original, orthognalized matrix that qr was based on
coords<-qr.Q(qr(coords))
vars<-diag(rexp(n,1))
R<-t(coords)%*%vars%*%coords
#inverse of sq matrix A multiplied by scalar c --> 1/c*inv(A)
#determinant of sq matrix A multiplied by scalar c --> c^nrow(A)*det(A)


M<-kronecker(R,C)
x<-matrix(rmvn(1,rep(0,nrow(M)),M),ncol=n)
rownames(x)<-tree$tip.label
pics<-sapply(1:n,function(ii) pic(x[,ii],tree))
cov(pics)
zlim<-c(min(cov(pics),R),max(cov(pics),R))
image(cov(pics),zlim=zlim)
image(R,zlim=zlim)
#nice