#plot an autocorrelated Brownian motion simulation
#' @export
plot.corateBM<-function(sim,tree,cols=c('deepskyblue','darkgray','brown1'),
                        phenogram=T,val.range=range(sim$R),res=100,
                        xlab='time',ylab='trait value',...){
  colramp<-colorRampPalette(cols)(res)
  if((val.range[2]-val.range[1])==0){
    colvec<-colramp[round((res+1)/2)]
  }else{
    inds<-round((sim$R-val.range[1])/(val.range[2]-val.range[1])*(res-1))+1
    inds[inds<1]<-1;inds[inds>res]<-res
    colvec<-colramp[inds]
  }
  if(phenogram){
    n<-length(tree$tip.label)
    if(length(sim$X)==n){
      tmp<-fastAnc(tree,sim$X)
      plot(c(sim$X[tree$tip.label],tmp)~node.depth.edgelength(tree),col='white',
           xlab=xlab,ylab=ylab,...)
      segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
               y0=c(sim$X,tmp)[tree$edge[,1]],y1=c(sim$X,tmp)[tree$edge[,2]],col=colvec,...)
    }else{
      plot(sim$X[c(tree$tip.label,1:tree$Nnode)]~node.depth.edgelength(tree),col='white',
           xlab=xlab,ylab=ylab,...)
      segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
               y0=sim$X[tree$edge[,1]],y1=sim$X[tree$edge[,2]],col=colvec,...)
    }
  }else{
    plot(tree,edge.color=colvec,...)
  }
}