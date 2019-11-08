#' @description I just find myself using this over and over again. It's not as 'safe' as Liam's function and doesn't do
#' tip labels, but it's a nice in how simple and flexible it is. Also works with the coloring system I created. Update: did not seem to make the bridges more noisy...but keeping the function anyways; I commented out a line of code in the
#' SimBM' function so I can easily switch between random sampling and this type of sampling when simulating BM bridges
#' 
#' @title fast phenogram function
#' @name fastpgram
#' @author Bruce
#' 
#' @param x trait data
#' @param tree tree of type 'phylo'
#' @param add boolean which permits or denies addition of plot to previous plot
#' @param ... additional arguments
#' 
fastpgram<-function(x,tree,add=F,node.col='black',edge.col='black',...){
  if(length(x)==length(tree$tip.label)){
    x<-c(x,fastAnc(tree,x))
  }
  if(add){
    points(x~node.depth.edgelength(tree),col=node.col,...)
  }else{
    plot(x~node.depth.edgelength(tree),col=node.col,...)
  }
  segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
           y0=x[tree$edge[,1]],y1=x[tree$edge[,2]],col=edge.col,...)
}
