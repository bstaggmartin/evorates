#' @title Gets coordinates of nodes from phylogenetic plotting environment
#' @name get.node.coords
#' @author Bruce
#' 
#' @return coordinates of nodes based on previous plot
#' 
#' @example plot.simBM(tree)+ get.node.coords()

get.node.coords<-function(){
  tmp<-get("last_plot.phylo",envir=.PlotPhyloEnv)[c('xx','yy')]
  cbind(tmp[[1]],tmp[[2]])
}
