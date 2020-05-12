#' @export
gen.corateBM<-function(tree,r0=0,rtrend=0,rvar=1,x0=0,xtrend=0,internal=F){
  
  node.rates<-rep(NA,tree$Nnode*2+1);node.rates[length(tree$tip.label)+1]<-r0
  for(e in 1:nrow(tree$edge)){
    node.rates[tree$edge[e,2]]<-node.rates[tree$edge[e,1]]+rnorm(1,rtrend*tree$edge.length[e],rvar*tree$edge.length[e])
  }
  edge.rates<-exp(apply(cbind(node.rates[tree$edge[,1]],node.rates[tree$edge[,2]]),1,mean))
  
  traits<-rep(NA,tree$Nnode*2+1);traits[length(tree$tip.label)+1]<-x0
  for(e in 1:nrow(tree$edge)){
    traits[tree$edge[e,2]]<-traits[tree$edge[e,1]]+rnorm(1,xtrend*tree$edge.length[e],edge.rates[e]*tree$edge.length[e])
  }
  
  if(internal){
    names(traits)<-c(tree$tip.label,(length(tree$tip.label)+1):(2*tree$Nnode+1))
  }else{
    traits<-traits[1:length(tree$tip.label)]
    names(traits)<-tree$tip.label
  }
  
  out<-list(traits=traits,edge.rates=edge.rates,ln.node.rates=node.rates,
            params=list(r0=r0,rtrend=rtrend,rvar=rvar,x0=x0,xtrend=xtrend,internal=internal))
  class(out)<-"corateBM"
  return(out)
}

#' @export
plot.corateBM<-function(corateBM,tree,cols=colorRampPalette(c("skyblue2","cyan","chartreuse2","goldenrod","red"))(100),log=F,
                        norm.lb=NULL,norm.ub=NULL,phenogram=F,xlab="time",ylab="trait",lwd=1,...){
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
  rates.norm[rates.norm>1]<-1
  rates.norm[rates.norm<0]<-0
  rates.norm<-rates.norm*(length(cols)-1)+1
  col.vec<-cols[round(rates.norm)]
  if(!phenogram){
    plot(tree,edge.color=col.vec,edge.width=lwd,...)
  }else{
    if(!corateBM$params$internal){
      trait.vec<-c(corateBM$traits,fastAnc(tree,corateBM$traits))
    }else{
      trait.vec<-corateBM$traits
    }
    plot(1,type='n',xlab=xlab,ylab=ylab,
         xlim=c(0,node.depth.edgelength(tree)[1]),ylim=c(min(corateBM$traits),max(corateBM$traits)))
    segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
             y0=trait.vec[tree$edge[,1]],y1=trait.vec[tree$edge[,2]],col=col.vec,lwd=lwd,...)
  }
}
