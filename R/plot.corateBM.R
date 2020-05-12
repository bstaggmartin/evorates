#plot an autocorrelated Brownian motion simulation
#' @export
plot.corateBM<-function(sim,tree,cols=c('deepskyblue','darkgray','brown1'),phylogram=T,...){
  colramp<-colorRampPalette(cols)(100)
  if((max(sim$R)-min(sim$R))==0){
    colvec<-colramp[51]
  }else{
    colvec<-colramp[round((sim$R-min(sim$R))/(max(sim$R)-min(sim$R))*99)+1]
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