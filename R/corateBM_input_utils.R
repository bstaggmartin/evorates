#' @export
prep.tree<-function(tree,tol=0){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladwise')
  tree<-di2multi(tree,tol=tol)
}