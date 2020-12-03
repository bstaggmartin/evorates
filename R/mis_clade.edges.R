#finds most recent common ancestral INTERNAL node of specified nodes (so, if tip, finds node ancestral
#to tip rather than tip itself) and then gets numbers of all edge descending from that node
#' @export
get.clade.edges<-function(tree,node){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  if(length(node)>1){
    try.node<-try(getMRCA(tree,node),silent=T)
    if(inherits(try.node,'try-error')){
      stop('one or more of the specified nodes (',paste(node,collapse=', '),') do not exist in tree')
    }
    node<-try.node
  }
  if(is.character(node)){
    node<-which(tree$tip.label==node)
    if(length(node)==0){
      stop('specified node ',node,' does not exist in tree')
    }
  }
  if(!(node%in%tree$edge)){
    stop('node ',node,' does not exist in tree')
  }
  out<-which(tree$edge[,1]==node)
  if(length(out)==0){
    out<-which(tree$edge[,2]==node)
  }else{
    new.edges<-out
    while(length(new.edges)>0){
      new.edges<-unlist(lapply(tree$edge[new.edges,2],function(ii) which(tree$edge[,1]==ii)))
      out<-c(out,new.edges)
    }
  }
  sort(out)
}

#function to take 2 clades(/edge groups) and cut smaller edge.group out of larger one
#' @export
exclude.clade<-function(tree,node1=NULL,node2=NULL,return.both=F,
                        edge.group1=NULL,edge.group2=NULL){
  if(is.null(edge.group1)){
    if(!is.null(node1)){
      edge.group1<-get.clade.edges(tree,node1)
    }else{
      stop('must specify first set of edges of interest by either setting node1 equal to a set of nodes defining a monophyletic clade OR setting edge.group1 equal to a vector of edge indices')
    }
  }
  if(is.null(edge.group2)){
    if(!is.null(node2)){
      edge.group2<-get.clade.edges(tree,node2)
    }else{
      stop('must specify second set of edges of interest by either setting node2 equal to a set of nodes defining a monophyletic clade OR setting edge.group2 equal to a vector of edge indices')
    }
  }
  out<-list(edge.group1,edge.group2)
  out<-out[order(lengths(out),decreasing=T)]
  out[[1]]<-out[[1]][!(out[[1]]%in%out[[2]])]
  if(!return.both){
    out[[1]]
  }else{
    out
  }
}

# colvec<-rep('black',nrow(tree$edge))
# clade<-sample(tree$tip.label,4)
# clade<-cut.out.clade(fit$call$tree,clade[1:2],clade[3:4],return.both=T)
# colvec[clade[[1]]]<-'red'
# colvec[clade[[2]]]<-'blue'
# plot(tree,edge.color=colvec)
# compare.params(fit,get.bg.rate(fit,edge.group=clade[[2]]),get.bg.rate(fit,edge.group=clade[[1]]))
