#works for the most part, except if new node has created new edge and position is specified to be past old
#edge...
#works in those cases now--simple sorting issue
#wrapper for bind.tip to allow adding multiple tips to a phylogeny at once
#' @export
multi.bind.tip<-function(tree,names,edge.lengths=NULL,nodes,positions=0){
  if(is.null(edge.lengths)){
    edge.lengths<-NA
  }
  args.list<-list(names=names,edge.lengths=edge.lengths,nodes=nodes,positions=positions)
  max.len<-max(sapply(args.list,length))
  for(i in 1:length(args.list)){
    args.list[[i]]<-rep(args.list[[i]],length.out=max.len)
  }
  nodes<-args.list$nodes
  positions<-args.list$positions
  names(nodes)<-ifelse(nodes<=length(tree$tip.label),nodes+tree$Nnode,nodes-length(tree$tip.label))
  args.list<-lapply(args.list,function(ii) ii[order(as.numeric(names(nodes)),-positions)])
  nodes<-nodes[order(as.numeric(names(nodes)))]
  int.nodes<-nodes>length(tree$tip.label)
  tips<-nodes<=length(tree$tip.label)
  for(i in 1:max.len){
    try.tree<-try(phytools::bind.tip(tree,args.list$names[i],
                                     if(is.na(args.list$edge.lengths[i])) NULL else args.list$edge.lengths[i],
                                     args.list$nodes[i],
                                     args.list$positions[i]),silent=T)
    if(inherits(try.tree,'try-error')){
      warning("failed to bind tip '",args.list$names[i],"' to node ",nodes[i],' due to following error:\n',
              try.tree)
    }else{
      tree<-try.tree
      args.list$nodes[int.nodes]<-args.list$nodes[int.nodes]+if(args.list$positions[i]<=0) 1 else 2
      tmp<-which(tree$tip.label==args.list$names[i])
      tmp<-args.list$nodes>=tmp
      args.list$nodes[tips&tmp]<-args.list$nodes[tips&tmp]+1
    }
  }
  tree
}