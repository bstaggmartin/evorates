gen.corateBM<-function(tree,r0=0,rtrend=0,rvar=1,x0=0,xtrend=0,internal=F){
  
  ##HELPER FUNCTIONS##
  trans.C<-function(C,rates){
    return(apply(array(sapply(1:length(rates),function(e) C[,,e]*rates[e]),dim=dim(C)),c(1,2),sum))
  }
  
  
  ##GENERATE VCV MATRIX FOR PHYLOGENY##
  C<-array(0,dim=c(rep(tree$Nnode*2+1,2),nrow(tree$edge)))
  for(e in 1:nrow(tree$edge)){
    D<-unique(c(tree$edge[e,2],getDescendants(tree,tree$edge[e,2])))
    C[as.matrix(cbind(expand.grid(D,D),rep(e,length(D)^2)))]<-tree$edge.length[e]
  }
  
  
  node.rates<-as.vector(rmvnorm(1,
                                mean=r0+rtrend*node.depth.edgelength(tree),
                                sigma=trans.C(C,rep(rvar,dim(C)[3]))))
  node.rates[length(tree$tip.label)+1]<-r0
  edge.rates<-exp(apply(cbind(node.rates[tree$edge[,1]],node.rates[tree$edge[,2]]),1,mean))
  
  
  if(internal){
    traits<-as.vector(rmvnorm(1,
                              mean=x0+xtrend*diag(trans.C(C,edge.rates)),
                              sigma=trans.C(C,edge.rates)))
    traits[length(tree$tip.label)+1]<-x0
    names(traits)<-c(tree$tip.label,(length(tree$tip.label)+1):(2*tree$Nnode+1))
  }else{
    traits<-as.vector(rmvnorm(1,
                              mean=x0+xtrend*diag(trans.C(C[1:length(tree$tip.label),1:length(tree$tip.label),],edge.rates)),
                              sigma=trans.C(C[1:length(tree$tip.label),1:length(tree$tip.label),],edge.rates)))
    names(traits)<-tree$tip.label
  }
  
  
  out<-list(traits=traits,edge.rates=edge.rates,ln.node.rates=node.rates,
            params=list(r0=r0,rtrend=rtrend,rvar=rvar,x0=x0,xtrend=xtrend,internal=internal))
  class(out)<-"corateBM"
  return(out)
}


plot.corateBM<-function(corateBM,tree,cols=heat.colors(100),log=F,
                        norm.lb=NULL,norm.ub=NULL,...){
  if(!is.list(corateBM)){
    corateBM<-list(edge.rates=corateBM)
  }
  if(is.null(norm.lb)){
    if(log){
      norm.lb<-min(log(corateBM$edge.rates))
    }else{
      norm.lb<-min(corateBM$edge.rates)
    }
  }
  if(is.null(norm.ub)){
    if(log){
      norm.ub<-max(log(corateBM$edge.rates))
    }else{
      norm.ub<-max(corateBM$edge.rates)
    }
  }
  if(log){
    rates.norm<-(log(corateBM$edge.rates)-norm.lb)/(norm.ub-norm.lb)
  }else{
    rates.norm<-(corateBM$edge.rates-norm.lb)/(norm.ub-norm.lb)
  }
  rates.norm<-rates.norm*(length(cols)-1)+1
  col.vec<-cols[round(rates.norm)]
  plot(tree,edge.color=col.vec,...)
}


get.edge.rates<-function(tree,par.mat){
  return(exp(sapply(1:ncol(par.mat),function(ii) apply(cbind(par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,1]],
                                                             par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,2]]),
                                                       1,mean))))
}




start<-proc.time()
par.mat<-relaxed.clock.BM(tree,corateBM$traits,n.iter=1e5,thin=25,tune.period=1e4,report.every=50,win.update=100,block.size=10)
elapsed<-proc.time()-start

thin<-par.mat[,seq(1200,ncol(par.mat),2)]
edge.rates<-get.edge.rates(tree,thin)
pdf("promising initial results.pdf",height=10,width=10)
par(mfrow=c(2,2))
plot(log(apply(edge.rates,1,median))~log(corateBM$edge.rates),pch=19,xlab="true ln(rate)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=rgb(1,0,0,0.5))
segments(x0=log(corateBM$edge.rates),x1=log(corateBM$edge.rates),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
abline(0,1,lwd=4)
plot(corateBM,tree,norm.lb=-4,norm.ub=4,
     cols=colorRampPalette(c("skyblue2","cyan2","chartreuse2","goldenrod","red"))(20),log=T,edge.width=5,main="true ln(rates)")
frame()
med.edge.rates<-apply(edge.rates,1,median)
plot.corateBM(med.edge.rates,tree,norm.lb=-4,norm.ub=4,
              cols=colorRampPalette(c("skyblue2","cyan2","chartreuse2","goldenrod","red"))(20),log=T,edge.width=5,
              main="median ln(rates) from posterior sample")
dev.off()