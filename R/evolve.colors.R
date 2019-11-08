#' @description Works by creating a 1 to 3-dimensional color space and simulating BM evolution within that space; to make sure lineages
#' smoothly grade into one another, internal node values are estimated using ancestral state reconstruction, 
#' and edge-wise colors are averaged across adjacent nodes
#' 
#'@title Produces an edge-wise color map for a given phylogenetic tree that exhibits phylogenetic signal.
#'@name evolve.colors
#'@author Bruce
#'
#'@param tree tree of class phylo
#'@param col.space.res 
#'@param res 
#'@param col.space.d
#'@param rate
#'@param return.nodes
#'@param ... (other arguments are passed)
#'

evolve.colors<-function(tree,col.space.res=100,col.space.d=2,rate=0.1,return.nodes=F,
                        hlim=c(0,1),s=1,v=1,slim=c(0.4,1),vlim=c(0.4,1),circular.h=ifelse(all(hlim==c(0,1)),T,F),alph=1,
                        plot=F,...){
  if(col.space.d!=1){
    sv.vals<-switch(col.space.d-1,
                    cbind(seq(slim[1],slim[2],length.out=col.space.res),seq(vlim[1],vlim[2],length.out=col.space.res)),
                    expand.grid(seq(slim[1],slim[2],length.out=col.space.res),seq(vlim[1],vlim[2],length.out=col.space.res)))
    col.space<-array(rainbow(col.space.res,rep(sv.vals[,1],each=col.space.res),rep(sv.vals[,2],each=col.space.res),
                             hlim[1],hlim[2],alpha=alph),
                     dim=rep(col.space.res,col.space.d))
    if(is.vector(rate)){
      rate<-diag(rate,nrow=col.space.d,ncol=col.space.d)
    }
    node.cols<-matrix(rmvn(1,rep(0,col.space.d*length(tree$tip.label)),kronecker(rate,vcv(tree))),ncol=col.space.d)
    rownames(node.cols)<-tree$tip.label
    node.cols<-rbind(node.cols,anc.recon(node.cols,tree))
    edge.cols<-sapply(1:col.space.d,function(ii) apply(matrix(node.cols[as.vector(tree$edge),ii],ncol=2),1,mean))
    if(return.nodes){
      edge.cols<-rbind(edge.cols,node.cols)
    }
    if(circular.h){
      problem.indices<-which(edge.cols[,1]>hlim[2])
      while(!all(edge.cols[,1]<=hlim[2])){
        edge.cols[problem.indices,1]<-edge.cols[problem.indices,1]-(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols[,1]>hlim[2])
      }
      problem.indices<-which(edge.cols[,1]<hlim[1])
      while(!all(edge.cols[,1]>=hlim[1])){
        edge.cols[problem.indices,1]<-edge.cols[problem.indices,1]+(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols[,1]<hlim[1])
      }
      edge.cols[,1]<-round(edge.cols[,1]*(col.space.res-1)+1)
    }else{
      edge.cols[,1]<-round((edge.cols[,1]-min(edge.cols[,1]))/(max(edge.cols[,1])-min(edge.cols[,1]))*(col.space.res-1)+1)
    }
    for(i in 2:col.space.d){
      edge.cols[,i]<-round((edge.cols[,i]-min(edge.cols[,i]))/(max(edge.cols[,i])-min(edge.cols[,i]))*(col.space.res-1)+1)
    }
    edge.cols<-col.space[edge.cols]
    if(plot){
      plot(tree,edge.color=edge.cols,edge.width=4,...)
    }
    if(return.nodes){
      list(colmap=edge.cols[1:nrow(tree$edge)],node.cols=edge.cols[(nrow(tree$edge)+1):length(edge.cols)])
    }else{
      edge.cols
    }
  }else{
    col.space<-rainbow(col.space.res,s,v,hlim[1],hlim[2])
    node.cols<-t(rmvn(1,rep(0,length(tree$tip.label)),rate*vcv(tree)))
    rownames(node.cols)<-tree$tip.label
    node.cols<-c(node.cols,anc.recon(node.cols,tree))
    edge.cols<-apply(matrix(node.cols[as.vector(tree$edge)],ncol=2),1,mean)
    if(return.nodes){
      edge.cols<-c(edge.cols,node.cols)
    }
    if(circular.h){
      problem.indices<-which(edge.cols>hlim[2])
      while(!all(edge.cols<=hlim[2])){
        edge.cols[problem.indices]<-edge.cols[problem.indices]-(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols>hlim[2])
      }
      problem.indices<-which(edge.cols<hlim[1])
      while(!all(edge.cols>=hlim[1])){
        edge.cols[problem.indices]<-edge.cols[problem.indices]+(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols<hlim[1])
      }
      edge.cols<-round(edge.cols*(col.space.res-1)+1)
    }else{
      edge.cols<-round((edge.cols-min(edge.cols))/(max(edge.cols)-min(edge.cols))*(col.space.res-1)+1)
    }
    edge.cols<-col.space[edge.cols]
    if(plot){
      plot(tree,edge.color=edge.cols,edge.width=4,...)
    }
    if(return.nodes){
      list(colmap=edge.cols[1:nrow(tree$edge)],node.cols=edge.cols[(nrow(tree$edge)+1):length(edge.cols)])
    }else{
      edge.cols
    }
  }
}

#Function to modify colors of edges in specified lineages/clades (can be inverted)
##Helpful for modifying existing colmaps or masking everything except for certain focal lineages/clades
alter.colmap<-function(tree,colmap,MRCAs=NULL,lins=NULL,alph=NA,mod.val=0,col=NULL,invert=F){
  if(!is.null(lins)){
    if(is.character(lins)){
      lins<-match(lins,tree$tip.label)
    }
    get.lin.edges<-function(tree,lin){
      prev.len<-0
      pp<-which(tree$edge[,2]==lin)
      while((length(pp)-prev.len)>0){
        prev.len<-length(pp)
        pp<-c(pp,which(tree$edge[,2]==tree$edge[pp[length(pp)],1]))
      }
      pp
    }
    lin.edges<-lapply(lins,get.lin.edges,tree=tree)
  }else{
    lin.edges<-NULL
  }
  if(!is.null(MRCAs)){
    clade.edges<-lapply(MRCAs,function(nn) which(tree$edge[,2]%in%phytools::getDescendants(tree,nn)))
  }else{
    clade.edges<-NULL
  }
  edges<-unique(c(unlist(clade.edges),unlist(lin.edges)))
  if(invert){
    edges<-(1:nrow(tree$edge))[!(1:nrow(tree$edge)%in%edges)]
  }
  if(!is.null(col)){
    cols<-rep(col,length.out=length(edges))
  }else{
    cols<-colmap[edges]
  }
  colmap[edges]<-alter.cols(cols,alph=alph,mod.val=mod.val)
  colmap
}