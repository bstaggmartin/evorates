#'@title Same function as above, but for lineages rather than clades
#'@name color.lineages
#'@author Bruce
#'
#'@description In the case of overlapping lineages, the color is resolved by averaging rgb values
#'
#'@param tree tree of class phylo
#'@param lins character vector of tip labels or numeric vector of node ids
#'@param alph numeric, 0-1 with 0 being transparent and 1 being opaque
#'
#'

color.lineages<-function(tree, lins, cols=palette()[2:(length(lins)+1)], alph = NA, base.cols=palette()[1], base.alph=NA){
  
  if(length(cols) < length(lins)){
    cols<-rep(cols,length.out=length(lins))
  }
  if(length(alph)<length(lins)){
    alph<-rep(alph,length.out=length(lins))
  }
  cols<-alter.cols(cols,alph=alph)
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
  lin.scores<-rep(1:length(lin.edges),lengths(lin.edges))
  tmp.cols<-t(col2rgb(cols[lin.scores],alpha=T))/255
  edge.wise.lin.scores<-split(tmp.cols,unlist(lin.edges))
  lin.scores<-as.numeric(names(edge.wise.lin.scores))
  cols<-rgb(t(sapply(edge.wise.lin.scores,function(ii) apply(matrix(ii,ncol=4),2,mean))))
  if(length(base.cols)<nrow(tree$edge)){
    base.cols<-rep(base.cols,length.out=nrow(tree$edge))
  }
  if(length(base.alph)<nrow(tree$edge)){
    base.alph<-rep(base.alph,length.out=nrow(tree$edge))
  }
  colmap<-alter.cols(base.cols,base.alph)
  colmap[lin.scores]<-cols
  colmap
}