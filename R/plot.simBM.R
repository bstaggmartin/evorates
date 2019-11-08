#' @export
plot.simBM<-function(sim,tree,n.sample=10,choose=NULL,lines=T,pts=F,polyg=F,
                     lwd=1,lty=1,cex=1,pch=16,border=NA,zoom.window=F,add=F,colmap="black",alph=NA,
                     xlab="time",ylab="trait value",rev.time=F,
                     xlim=c(0,max(sim$ts[length(sim$ts)],node.depth.edgelength(tree),na.rm=T)),
                     ylim=c(min(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],],
                                sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],],na.rm=T),
                            max(sim$edges[,sim$ts>=xlim[1]&sim$ts<=xlim[2],],
                                sim$nodes[node.depth.edgelength(tree)>=xlim[1]&node.depth.edgelength(tree)<=xlim[2],],na.rm=T))){
  
  ##BASIC ERROR CHECKS AND FIXES##
  if(!inherits(tree,'phylo')){
    stop('tree should be object of class \"phylo\"')
  }
  if(!inherits(sim,'simBM')){
    stop('sim should be object of class \"simBM\"')
  }
  n.sample<-round(n.sample)
  if(n.sample<1){
    stop('n.sample is less than 1')
  }
  ####
  
  if(n.sample>ncol(sim$nodes)){
    samp<-sample(1:ncol(sim$nodes),n.sample,replace=T)
  }else{
    samp<-sample(1:ncol(sim$nodes),n.sample)
  }
  if(!is.null(choose)){
    samp<-choose
  }
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
  if(lines){
    for(i in samp){
      for(e in 1:nrow(edges)){
        lines(c(sim$nodes[edges[e,1],i],sim$edges[e,,i][!is.na(sim$edges[e,,i])],sim$nodes[edges[e,2],i])~
                c(node.depth.edgelength(tree)[edges[e,1]],sim$ts[!is.na(sim$edges[e,,i])],node.depth.edgelength(tree)[edges[e,2]]),
              col=colmap[e],lwd=lwd,lty=lty)
      }
    }
  }
  if(pts){
    node.colmap<-get.node.colmap(colmap,tree)
    node.vec<-as.vector(sim$nodes[,samp])
    edge.vec<-as.vector(sim$edges[,,samp])
    node.t.vec<-rep(node.depth.edgelength(tree),length(samp))
    edge.t.vec<-rep(rep(sim$ts,each=nrow(edges)),length(samp))
    node.cols<-rep(node.colmap,length(samp))
    edge.cols<-rep(colmap,length(sim$ts)*length(samp))
    edge.cols<-edge.cols[!is.na(edge.vec)]
    edge.t.vec<-edge.t.vec[!is.na(edge.vec)]
    edge.vec<-edge.vec[!is.na(edge.vec)]
    points(c(node.vec,edge.vec)~c(node.t.vec,edge.t.vec),col=c(node.cols,edge.cols),cex=cex,pch=pch)
  }
  if(polyg){
    for(e in 1:nrow(edges)){
      if(length(which(!is.na(sim$edges[e,,1])))==0){
        upy<-apply(cbind(sim$nodes[edges[e,1],samp],
                         sim$nodes[edges[e,2],samp]),2,max)
        dwny<-rev(apply(cbind(sim$nodes[edges[e,1],samp],
                              sim$nodes[edges[e,2],samp]),2,min))
      }else if(length(which(!is.na(sim$edges[e,,1])))==1){
        upy<-apply(cbind(sim$nodes[edges[e,1],samp],
                         sim$edges[e,,samp][which(!is.na(sim$edges[e,,1])),],
                         sim$nodes[edges[e,2],samp]),2,max)
        dwny<-rev(apply(cbind(sim$nodes[edges[e,1],samp],
                              sim$edges[e,,samp][which(!is.na(sim$edges[e,,1])),],
                              sim$nodes[edges[e,2],samp]),2,min))
      }else{
        upy<-apply(cbind(sim$nodes[edges[e,1],samp],
                         t(sim$edges[e,,samp][which(!is.na(sim$edges[e,,1])),]),
                         sim$nodes[edges[e,2],samp]),2,max)
        dwny<-rev(apply(cbind(sim$nodes[edges[e,1],samp],
                              t(sim$edges[e,,samp][which(!is.na(sim$edges[e,,1])),]),
                              sim$nodes[edges[e,2],samp]),2,min))
      }
      x<-c(node.depth.edgelength(tree)[edges[e,1]],sim$ts[which(!is.na(sim$edges[e,,1]))],node.depth.edgelength(tree)[edges[e,2]])
      polygon(c(upy,dwny)~c(x,rev(x)),col=colmap[e],border=border,lwd=lwd,lty=lty)
    }
  }
}
