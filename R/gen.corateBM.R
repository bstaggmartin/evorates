#simulate trait and rate data under an autocorrelated Brownian motion model
#' @export
gen.corateBM<-function(tree,R0=0,Rsig2=1,X0=0,anc.states=F){
  n<-length(tree$tip.label)
  n_e<-nrow(tree$edge)
  X<-rep(NA,n+tree$Nnode)
  n_R<-X
  R<-rep(NA,n_e)
  names(X)<-c(tree$tip.label,1:tree$Nnode)
  X[n+1]<-X0;n_R[n+1]<-R0
  for(i in 1:n_e){
    n_R[tree$edge[i,2]]<-rnorm(1,n_R[tree$edge[i,1]],sqrt(tree$edge.length[i]*Rsig2))
    R[i]<-rnorm(1,mean(n_R[tree$edge[i,]]),sqrt(Rsig2*tree$edge.length[i]/12))
    X[tree$edge[i,2]]<-rnorm(1,X[tree$edge[i,1]],sqrt(tree$edge.length[i]*exp(R[i])))
  }
  if(anc.states){
    out<-list('R0'=R0,'Rsig2'=Rsig2,'X0'=X0,'R'=R,'X'=X)
  }else{
    out<-list('R0'=R0,'Rsig2'=Rsig2,'X0'=X0,'R'=R,'X'=X[1:n])
  }
  class(out)<-'corateBM'
  out
}