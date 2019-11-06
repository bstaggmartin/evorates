plot.simBM<-function(sim,tree,n.sample=10,
                     lwd=1,zoom.window=F,add=F,colmap="black",alph=NA,
                     xlab="time",ylab="trait value",rev.time=F,
                     xlim=c(0,max(sim$ts[length(sim$ts)],node.depth.edgelength(tree),na.rm=T)),
                     ylim=c(min(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],],
                                sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],],na.rm=T),
                            max(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],],
                                sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],],na.rm=T))){
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
