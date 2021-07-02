#' evorates: Fitting continuously variable rate models to comparative data on continuous traits
#'
#' evorates more detailed description
#' 
#' @section function category 1:
#' functions...
#'
#' @docType package
#' @name evorates
#' @import Rcpp methods
#' @importFrom ape plot.phylo drop.tip multi2di di2multi getMRCA node.depth.edgelength
#' @importFrom phytools bind.tip
#' @importFrom graphics plot
#' @importFrom rstan sampling extract Rhat ess_bulk ess_tail
#' @importFrom logspline logspline dlogspline
#' @useDynLib evorates, .registration=TRUE
NULL