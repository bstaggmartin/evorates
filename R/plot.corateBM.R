#plot an autocorrelated Brownian motion simulation
#' @export
plot.corateBM<-function(sim,tree,cols=c('deepskyblue','darkgray','brown1'),phylogram=T,val.range=range(sim$R),res=100,...){
  colramp<-colorRampPalette(cols)(res)
  if((val.range[2]-val.range[1])==0){
    colvec<-colramp[round((res+1)/2)]
  }else{
    inds<-round((sim$R-val.range[1])/(val.range[2]-val.range[1])*(res-1))+1
    inds[inds<1]<-1;inds[inds>res]<-res
    colvec<-colramp[inds]
  }
  if(phylogram){
    n<-length(tree$tip.label)
    if(length(sim$X)==n){
      tmp<-fastAnc(tree,sim$X)
      plot(c(sim$X,tmp)~node.depth.edgelength(tree),col='white',...)
      segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
               y0=c(sim$X,tmp)[tree$edge[,1]],y1=c(sim$X,tmp)[tree$edge[,2]],col=colvec,...)
    }else{
      plot(sim$X~node.depth.edgelength(tree),col='white',...)
      segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
               y0=sim$X[tree$edge[,1]],y1=sim$X[tree$edge[,2]],col=colvec,...)
    }
  }else{
    plot(tree,edge.color=colvec,...)
  }
}