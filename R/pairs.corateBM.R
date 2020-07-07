#improve label handling
#' @export
pairs.corateBM<-function(sim,tree,trait=1:ncol(sim$X),lwd=1,col=c('deepskyblue','darkgray','brown1'),alpha=NA,
                         val.range=if(is.null(sim$R)) c(0,0) else range(sim$R),res=100,...){
  n<-length(tree$tip.label)
  if(nrow(sim$X)==n){
    scaled.tree<-tree
    if(!is.null(sim$R)){
      scaled.tree$edge.length<-tree$edge.length*exp(sim$R)
    }
    anc.states<-matrix(NA,tree$Nnode,ncol(sim$X))
    rownames(anc.states)<-n+1:tree$Nnode
    for(i in trait){
      anc.states[,i]<-fastAnc(scaled.tree,sim$X[,i])
    }
    sim$X<-rbind(sim$X,anc.states)
  }
  sim$X<-as.matrix(sim$X[c(tree$tip.label,n+1:tree$Nnode),])
  old.par<-par(no.readonly=T)
  par(mfrow=c(length(trait),length(trait)),mar=c(0,0,0,0),oma=c(5.1,4.1,0,0),xpd=T)
  for(i in 1:length(trait)){
    for(j in 1:length(trait)){
      if(j==1){
        yaxt=NULL
        args.y.mtext<-list(text=paste('trait',trait[i]),
                         side=2,
                         line=3,
                         cex=0.75)
      }else{
        yaxt='n'
        args.y.mtext<-list(NULL)
      }
      if(i==length(trait)){
        xaxt=NULL
        args.x.mtext<-list(text=paste('trait',trait[j]),
                         side=1,
                         line=3,
                         cex=0.75)
      }else{
        xaxt='n'
        args.x.mtext<-list(NULL)
      }
      if(i==j){
        if(i==length(trait)){
          plot(sim,tree,trait=c(trait[i],trait[j]),alpha=0,
               xaxt=xaxt,yaxt=yaxt,...)
          new.range<-range(sim$X[,i])
          node.depths<-node.depth.edgelength(tree)
          node.depths<-node.depths/max(node.depths)*diff(new.range)+min(new.range)
          plot(sim,tree,trait=i,lwd=lwd,col=col,alpha=alpha,val.range=val.range,res=res,
               xaxt='n',yaxt=yaxt,add=T,node.depths=node.depths,...)
        }else{
          plot(sim,tree,trait=trait[i],lwd=lwd,col=col,alpha=alpha,val.range=val.range,res=res,
               xaxt=xaxt,yaxt=yaxt,...)
        }
        do.call(mtext,args.x.mtext)
        do.call(mtext,args.y.mtext)
      }else{
        plot(sim,tree,trait=c(trait[j],trait[i]),lwd=lwd,col=col,alpha=alpha,val.range=val.range,res=res,
             xaxt=xaxt,yaxt=yaxt,...)
        do.call(mtext,args.x.mtext)
        do.call(mtext,args.y.mtext)
      }
    }
  }
  par(old.par)
}
