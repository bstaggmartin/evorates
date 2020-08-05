#' @export
prep.tree<-function(tree,tol=0){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=tol)
}