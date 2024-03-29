#' evorates: Fitting continuously variable rate models to comparative data on continuous traits
#'
#' This package implements a Bayesian method for fitting "relaxed" models of continuous trait evolution
#' to comparative data (i.e., a phylogeny and associated trait data), whereby the rate of trait evolution
#' itself gradually changes over time and across lineages. More details on the model and method can
#' be found in the
#' \link[=https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syac068/6830631]{associated manuscript}.
#' The package provides additional tools for simulating data, as well as manipulating, analyzing, and visualizing
#' estimated model parameters. Below is a broad overview of the available functions in this package:
#' 
#' @section Basic functions:
#' \itemize{
#'    \item{To simulate data, use \code{\link{sim.evorates}()}.}
#'    \item{To fit data to models, use \code{\link{fit.evorates}()} (which is a wrapper for the functions
#'    \code{\link{input.evorates}()}, \code{\link{run.evorates}()}, and \code{\link{output.evorates}()},
#'    which may be used to provide finer control over certain aspects of the model fitting process).}
#'    \item{To check that fitted model converged and sampled posterior distributions adequately, use
#'    \code{\link{check.mix}()} and \code{\link{check.ess}()}.}
#'    \item{To subset, combine, or thin Hamiltonian Monte Carlo chains in a fitted model, use
#'    \code{\link{select.chains}()}, \code{\link{combine.chains}()}, \code{\link{exclude.warmup}()}, and
#'    \code{\link{thin.chains}()}}.
#' }
#' 
#' @section Analysis functions:
#' \itemize{
#'    \item{To calculate Savage-Dickey raitos and see if a fitted model yields "substantial" evidence for
#'    rate heterogeneity, use \code{\link{get.sd}()}.}
#'    \item{For extraction, summarization, and manipulation of posterior samples (including 
#'    calculating posterior probabilities!), see documentation on
#'    the \code{\link[=param_block-class]{param_block}} class and associated documentation on the
#'    \code{param_block} operators \code{\link[=grapes-chains-grapes]{\%chains\%}()},
#'    \code{\link[=grapes-quantiles-grapes]{\%quantiles\%}()},
#'    \code{\link[=grapes-means-grapes]{\%means\%}()},
#'    \code{\link[=grapes-diagnostics-grapes]{\%diagnostics\%}()},
#'    and \code{\link[=grapes-select-grapes]{\%select\%}()}, as well as the functions
#'    \code{\link{par.c}()} (for combining parameter blocks), \code{\link{rnorm.par}()}
#'    (for generating parameter blocks of normal random samples), and \code{\link{pwc}()}
#'    (for doing pairwise comparisons among parameter blocks).}
#'    \item{For convenient extraction of average rates along each branch, use \code{\link{get.R}()}}.
#'    \item{To "adjust" rates for trends, use \code{\link{remove.trend}()}}.
#'    \item{To calculate "background rates" and summarize rates over different parts of a phylogeny,
#'    use \code{\link{get.bg.rate}()}.}
#'    \item{To sample posterior distributions of trait values at nodes (i.e., ancestral state estimation)
#'    in a phylogeny given a fitted model, use \code{\link{get.post.traits}()}.}
#' }
#' 
#' @section Plotting functions:
#' \itemize{
#'    \item{There are plotting methods for simulated data and fitted models; see
#'    \code{\link{plot.evorates}()} and \code{\link{plot.evorates_fit}()} (NOTE:
#'    documentation still under construction).}
#'    \item{You WILL be able to plot plots of posterior samples vs. chain iterations
#'    ("traces") very soon using \code{\link{trace.plot}()}!}
#'    \item{You can plot histograms/density plots ("profiles") of posterior samples
#'    using \code{\link{prof.plot}()}.}
#'    \item{A lot of these functions rely on a nice little helper function
#'    \code{\link{alter.cols}()}, which can be used to mix a vector of colors with
#'    other colors or modify transparency.}
#'    \item{I plan on developing another function some day, perhaps called \code{level.plot()} or
#'    something like that, which will plot 2D histograms/density plots similarly to 
#'    \code{\link[graphics]{smoothScatter}()} and \code{\link[graphics]{contour}()}. The idea
#'    is to be able to look at posterior correlations among parameter estimates.}
#' }
#' 
#' @section Miscellaneous functions (NOTE: documentation still largely under construction):
#' \itemize{
#'    \item{Use \code{\link{get.clade.edges}()} to extract the edge indices in a phylogeny
#'    associated with a particular clade (defined by either its most recent common ancestor or
#'    a group of tip labels). You can also use \code{\link{exclude.clade}()} to take edges
#'    associated with some nested subclade out of a larger clade.}
#'    \item{The function \code{\link{edge.vcv}()} calculates the "edgewise" variance-covariance
#'    matrix of a phylogeny. Normally, phylogenetic variance-covariance matrices describe
#'    the covariance structure of trait values \emph{at nodes} expected under a Brownian Motion
#'    model. Edgewise variance-covariance matrices instead describe the covariance structure
#'    of \emph{average trait values along each branch} expected under a Brownian Motion model.}
#'    \item{Some functions rely a nice helper function, \code{\link{multi.bind.tip}()}, which
#'    wraps and generalizes the \pkg{phytools} function \code{\link[phytools]{bind.tip}()} to handle
#'    binding multiple tips to a phylogeny at once.}
#'    \item{There are \code{Ntip()}, \code{Nedge()}, and \code{Nnode()} methods for simulated
#'    data and fitted models to quickly extract these quantities as necessary.}
#'    \item{Further, there are some new convenient methods for extracting topological information
#'    from phylogenies, simulated data, and fitted models:
#'    \itemize{
#'       \item{Use \code{\link{edge.ranges}()} to get start and end times of each edge in a phylogeny
#'       in matrix form (with a row for each edge and two columns for start and end times).}
#'       \item{There are new edgewise
#'       "tree-walking" functions for getting edge indices corresponding to the ancestor,
#'       descendants, or sisters for each edge in a phylogeny: \code{\link{anc.edges}()},
#'       \code{\link{des.edges}()}, and \code{\link{sis.edges}()}.}
#'       \item{You can use \code{\link{tip.edges}()} to get the edge indices for each tip
#'       in a phylogeny and \code{\link{root.edges}()} to get indices of edges descending from
#'       the root node.}
#'       \item{Lastly, you can use \code{\link{ladder}()} to "ladderize" edges in a phylogeny.
#'       In the context of simulated data and fitted evorates models, this also rearranges
#'       edgewise information like branchwise rates accordingly.}
#'    }}
#' }
#'
#' @docType package
#' @name evorates
#' @import Rcpp methods graphics ape
#' @importFrom phytools bind.tip
#' @importFrom rstan sampling extract Rhat ess_bulk ess_tail
#' @importFrom logspline logspline dlogspline
#' @useDynLib evorates, .registration=TRUE
NULL