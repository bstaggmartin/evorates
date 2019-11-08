#' @title Plots the simulated trait output from simBM onto a tree
#' @name plot.simBM
#' @author Bruce
#' 
#' @description 
#' 
#' @param sim returned list from sim.BM containing simulated trait values
#' @param tree tree of class 'phylo'
#' @param n.sample numeric value for number of simulations to run
#' @param lwd numeric value for line width
#' @param zoom.window 
#' @param add logical: permits (T) or denies (F) addition of plot to previous plot
#' 
#' @return plot of simulated traits over time
#' 
#' @example 
#' plot.simBM(sim, tree)

plot.simBM<-function(sim,tree,n.sample=10,
                     lwd=1,zoom.window=F,add=F,colmap="black",alph=NA,
                     xlab="time",ylab="trait value",rev.time=F,
                     xlim=c(0,max(sim$ts[length(sim$ts)],node.depth.edgelength(tree),na.rm=T)),
                     ylim=c(min(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],],
                                sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],],na.rm=T),
                            max(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],],
                                sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],],na.rm=T))){

    #Simple Error Checks
  if(!inherits(tree,"phylo")) stop ("tree must be of class phylo")
  
  samp<-sample(1:ncol(sim$nodes),n.sample,replace=T)
  edges<-tree$edge
  if(length(colmap)<nrow(edges)){
    colmap<-rep(colmap,length.out=nrow(edges))
  }
  if(length(alph)<nrow(edges)){
    alph<-rep(alph,length.out=nrow(edges))
  }
  colmap<-alter.cols(colmap,alph=alph)
  if(zoom.window){
    ylim=c(min(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],samp],
               sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],samp],na.rm=T),
           max(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],samp],
               sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],samp],na.rm=T))
  }
  if(!add){
    plot(0,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,type='n',xaxt='n')
    if(rev.time){
      axis(1,at=(max(node.depth.edgelength(tree))-ceiling(max(node.depth.edgelength(tree)))):
             ceiling(max(node.depth.edgelength(tree))),
           labels=ceiling(max(node.depth.edgelength(tree))):0,tick=F,padj=-1)
    }else{
      axis(1,at=round(xlim[1]):round(xlim[2]),labels=round(xlim[1]):round(xlim[2]),tick=F,padj=-1)
    }
  }
  
  for(i in samp){
    for(e in 1:nrow(edges)){
      lines(c(sim$nodes[edges[e,1],i],sim$edges[e,,i][!is.na(sim$edges[e,,i])],sim$nodes[edges[e,2],i])~
              c(node.depth.edgelength(tree)[edges[e,1]],sim$ts[!is.na(sim$edges[e,,i])],node.depth.edgelength(tree)[edges[e,2]]),
            col=colmap[e],lwd=lwd)
    }
  }
}
