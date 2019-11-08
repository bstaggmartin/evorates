#'
#'@title Extract edge-wise rates from the output of a relaxed.clock.BM MCMC
#'@author Bruce
#'@name  get.edge.rates
#'
#'@description Helpful since the autocorrelated rate model output gives (log) node-wise rate parameters, rather than the edge-wise rates (which
#'can be approximated by averaging across adjacent nodes)
#'
#'@param tree tree of class phylo
#'@param par.mat matrix of 
#'@return Edge-wise rates from the rownian-Motion relaxed clock MCMC

get.edge.rates<-function(tree, par.mat){
  exp(sapply(1:ncol(par.mat),function(ii) {
    apply(cbind(par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,1]],par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,2]]),1,mean)
  }))
}
