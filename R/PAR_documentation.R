#' Class \code{param_block}: parameter block arrays
#' 
#' 
#' @description The class "\code{param_block}" is an S3 class I created for the \code{evorates}
#' package that is intended to make extraction, manipulation, and summarization of parameter
#' distributions as convenient and transparent as possible. In a nutshell, \code{param_block}
#' (short for "parameter block") objects
#' are arrays of numeric data that behave a bit differently than normal array objects and
#' come with some extra rules. It seemed particularly important to develop this class
#' because evorates models often include many parameters that users will want to manipulate in
#' various ways (e.g., calculating average rates for clades of interest, comparing rates among
#' clades, etc.).
#' 
#' 
#' @format
#' \code{param_block} arrays are typically 3D arrays of numeric data with a
#' variable first dimension/rows, a second dimension/columns corresponding to different
#' parameters, and a third dimension/slices corresponding to different chains (note that this
#' is distinct from the \pkg{rstan} convention of having columns and slices correspond to
#' chains and parameters, respectively). These arrays come in four different types, as
#' indicated by their \code{param_type} attribute:
#' \describe{
#'    \item{"\code{chains}"}{Arrays of posterior samples where the rows correspond to
#'          iterations of Hamiltonian Monte Carlo chains. These typically represent the entire
#'          posterior distribution of particular parameters. Currently, specific iterations
#'          are unlabeled, but I may change this in a future update.}
#'    \item{"\code{quantiles}"}{Arrays of posterior quantiles where the rows correspond to
#'          different quantiles labeled by their associated percentile. These are typically
#'          used to summarize posterior distributions in terms of point estimates (e.g., 50\%
#'          quantiles/medians) or credible intervals (e.g., 2.5\% and 97.5\% quantiles, the
#'          lower and upper bounds of the 95\% equal-tail interval).}
#'    \item{"\code{means}"}{Arrays of posterior means which always consist of 1 row. The row
#'          dimension is sometimes labeled "iterations" for compatibility, but individual rows
#'          are always unlabeled. These are typically used to summarize posterior distributions
#'          as point estimates or calculate posterior probabilities (see 
#'          \code{\link[=grapes-means-grapes]{\%chains\%}()} for more details on this).}
#'    \item{"\code{diagnostics}"}{Arrays of posterior diagnostics where the rows correspond to
#'          one of four different diagnostic statistics: initial values ("\code{inits}"), bulk
#'          effective sample sizes ("\code{bulk_ess}"), tail effective sample sizes
#'          ("\code{tail_ess}"), or Rhat diagnostics ("\code{Rhat}"). These are typically used
#'          to check if the posterior distribution for a particular parameter was sampled well.
#'          You can look at \code{\link[=grapes-diagnostics-grapes]{\%diagnostics\%}()}, \code{\link{check.mix}()}, \code{\link{check.ess}()},
#'          and \link[rstan]{Rhat} for more information on these quantities.}
#' }
#' Notably, \code{param_block} arrays of \code{param_type} "\code{chains}" can be converted to
#' any other type using \code{param_block} operators: \code{\link[=grapes-quantiles-grapes]{\%quantiles\%}()},
#' \code{\link[=grapes-means-grapes]{\%means\%}()}, and \code{\link[=grapes-diagnostics-grapes]{\%diagnostics\%}()}.
#' \code{param_block} arrays can also always be converted back to normal arrays via the \code{unclass()} function.
#' 
#' \subsection{Automatic Parameter Renaming}{
#'    For increased transparency and clarity, parameters names in \code{param_block} arrays are
#'    automatically updated to reflect parameter manipulations.
#'    For example, adding two parameters named "\code{par1}" and "\code{par2}" will result in a new parameter
#'    named "\code{par1\%+\%par2}". Similarly, applying a function like \code{exp()} to \code{par1} will result
#'    in the name "\code{exp(par1)}", and converting from a \code{param_type} of "\code{chains}" will rename
#'    parameters to "\code{quantiles(par1)}", "\code{means(par1)}", etc. Of course, this system sometimes
#'    results in really long, unnecessarily confusing names! The parameter names of
#'    \code{param_block} arrays can always be accessed and replaced using the \code{names()} function
#'    (see examples below).
#' }
#' 
#' \subsection{Dimension Collapsing}{
#'    As with base R arrays, dimensions of length 1 in \code{param_block} arrays are sometimes "collapsed".
#'    For comparison, think about when you take a single row or column out of a normal matrix--in most
#'    cases, the result is automatically converted to a vector (the only way to prevent this in base R is to
#'    specify \code{drop = FALSE}, e.g., \code{matrix[i, , drop = FALSE]}). \code{param_block}
#'    arrays work much the same way: for example, a 1000 (iteration) x 1 (parameter) x 4 (chain) 
#'    \code{param_block} array is typically converted to a 1000 (iteration) x 4 (chain) matrix. 
#'    However, to prevent loss of information, this matrix will \emph{also} have a "\code{parameters}"
#'    attribute storing the associated parameter name. For another example, consider a 1 x 4 x 1 array
#'    consisting of medians for 4 different parameters: such an array is typically converted to a vector 
#'    of length 4 with a "\code{quantiles}" attribute specifying the numbers correspond to "\code{50\%}"
#'    quantiles or medians, as well as a "\code{chains}" attribute with the associated chain name.
#' }
#' 
#' @section Creation and Extraction:
#' 
#' Since \code{param_block} arrays require information that is often impossible to generate on
#' the fly, I have not yet created functions for converting arbitrary objects to \code{param_block}
#' arrays. However, fitted evorates models (class "\code{evorates_fit}") are mostly just lists of
#' \code{param_block} arrays (see \code{\link{fit.evorates}()} for more
#' on the format of fitted evorates models) which can be accessed
#' via normal list subsetting (e.g., \code{cet_fit$chains}). Parts of these arrays
#' can also be directly extracted from fitted evorates models using the \code{param_block} operators: 
#' \code{\link[=grapes-chains-grapes]{\%chains\%}()}, \code{\link[=grapes-quantiles-grapes]{\%quantiles\%}()},
#' \code{\link[=grapes-means-grapes]{\%means\%}()}, and
#' \code{\link[=grapes-diagnostics-grapes]{\%diagnostics\%}()}. Some functions also output \code{param_block}
#' arrays, most
#' notably \code{\link{get.R}()}, \code{\link{get.post.traits}()}, \code{\link{get.bg.rate}()}, and
#' \code{\link{remove.trend}()}.
#' 
#' There is one special function that allows users to create \code{param_block} arrays almost "from
#' scratch", \code{\link{rnorm.par}()}, which creates \code{param_block} arrays of standard normal random samples
#' (i.e., mean 0 and standard deviation 1). This function can be helpful for
#' certain parameter manipulations, like simulating the distribution of rate changes over some time
#' interval or sampling ancestral trait values. I may try to generalize this function later, but there
#' hasn't really been a need for any other "\code{param_block} constructor" functions as of yet.
#' 
#' @section Subsetting and Appending:
#' 
#' \code{param_block} arrays are subsetted via operators: 
#' \code{\link[=grapes-chains-grapes]{\%chains\%}()},
#' \code{\link[=grapes-quantiles-grapes]{\%quantiles\%}()},
#' \code{\link[=grapes-means-grapes]{\%means\%}()},
#' \code{\link[=grapes-diagnostics-grapes]{\%diagnostics\%}()},
#' and \code{\link[=grapes-select-grapes]{\%select\%}()} (see examples below for more details
#' on usage). These functions allow you to easily
#' subset the rows and columns, but not slices,
#' of \code{param_block} arrays (though chains may
#' be manipulated for entire fitted evorates models using \code{\link{select.chains}()} and
#' \code{\link{combine.chains}()}). You can obviously use normal R indexing too select particular slices
#' (e.g., \code{array[, , k]}), but this converts \code{param_block} arrays to normal arrays. This is helpful for
#' compatibility of \code{param_block} arrays with other functions, but can obviously be frustrating.
#' I'm still working on ways to improve this particular quirk of \code{param_block} array behavior.
#' 
#' Multiple \code{param_block} arrays, provided they have compatible rows and slices (i.e., same number
#' of chains, same names), can be combined into larger arrays using the \code{c()} function.
#' \code{param_block} arrays with a \code{param_type} of "\code{chains}" have rows compatible with any other
#' \code{param_block} array (they can be coerced to another \code{param_type} on the fly), with the exception
#' of other \code{chains} \code{param_block} arrays that
#' have a different number of rows. Otherwise, \code{param_block} arrays generally have compatible rows with other arrays
#' of the same \code{param_type} provided they share least one row (i.e., quantile or diagnostic) in common. See
#' examples below for more details, as well as the relevant function \code{\link{par.c}()}.
#' 
#' @section Manipulation:
#' 
#' \code{param_block} arrays behave a bit differently than normal arrays when you do math with them.
#' In addition to automatic parameter renaming, \code{param_block} arrays enforce that math operations
#' combine act on \emph{parameters}. For example, \code{param_block \%select\% c("par1", "par2") + c(1, 2)} will result
#' in adding 1 to the parameter \code{par1} and 2 to \code{par2}, rather than recycling \code{c(1, 2)} to be the
#' same length as \code{param_block \%select\% c("par1", "par2")} and adding this to
#' \code{param_block \%select\% c("par1", "par2")} as a whole. Similarly,
#' \code{param_block \%select\% c("par1", "par2", "par3") + param_block \%select\% c("par4", "par5")} typically results
#' in adding \code{par1} to \code{par4}, \code{par2} to \code{par4}, and recycling the \emph{parameters} of 
#' \code{param_block \%select\% c("par4", "par5")} to add \code{par3} to \code{par4}. Overall, this makes parameter
#' manipulations more convenient and readable. Notably, \code{param_block} arrays must have strictly matching rows
#' and slices (i.e., same number and names) to be combined using math operations, and \code{param_block} arrays
#' cannot yet be combined with normal arrays via math operations.
#' See examples below for more details, as well as the relevant function
#' \code{\link{pwc}()}.
#' 
#' @section Visualization:
#' 
#' \code{param_block} arrays may be passed to the functions \code{\link{trace.plot}()} and \code{\link{prof.plot}()} to
#' visualize their contents as traces (i.e., parameter values on y-axis, iterations on x-axis) or histogram/density 
#' style plots (i.e., parameter values on x-axis, relative frequencies on y-axis). In the future, I hope make another
#' function for visualizing joint distributions of two parameters. In any case, since \code{param_block} objects are
#' really just arrays, they should be compatible with many other plotting functions as long as you subset
#' them to just 2 parameters and 1 chain.
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #get entire param_blocks for evorates fit
#' cet_fit$chains
#' cet_fit$quantiles
#' cet_fit$means
#' cet_fit$diagnostics
#' 
#' #create param_blocks using some functions
#' get.R(cet_fit)
#' rnorm.par(1500, 2, 4)
#' 
#' #accessing/changing names
#' parblock <- rnorm.par(1500, 2, 4)
#' names(parblock)
#' names(parblock) <- c("par1", "par2")
#' 
#' #extract particular parameters using operators
#' #see param_block operator documentation for more details
#' cet_fit %chains% c("R_mu", "R_sig2")
#' #notice dimension collapsing here
#' cet_fit %chains% "R_mu"
#' select.chains(cet_fit, 1) %chains% "R_mu"
#' cet_fit %means% "R_mu"
#' select.chains(cet_fit, 1) %means% "R_mu"
#' #numeric index selection generally corresponds to edge indices
#' cet_fit %chains% c(1, 2)
#' #these can also be used on param_blocks themselves, both for subsetting or conversion
#' parblock <- cet_fit %chains% 1:10
#' parblock %chains% 2
#' parblock %quantiles% "." #the . indicates all parameters in parblock
#' 
#' #can also extract particular rows using list selections
#' #first element of list corresponds to parameters, second to rows
#' #see param_block operator documentation for more details
#' cet_fit %chains% list(c("R_mu", "R_sig2"), c(1, 2, 47))
#' cet_fit %quantiles% list(c("R_mu", "R_sig2"), c(0.32, 0.64))
#' cet_fit %diagnostics% list(c("R_mu", "R_sig2"), "ess")
#' #also works with param_blocks themselves
#' parblock %quantiles% list(".", 0.1)
#' #%select% is a more general-purpose subsetting function built for any param_type
#' parblock %select% list(c(1, 2), 3)
#' 
#' #combining param_blocks--how are average rates for some focal clades affected by R_sig2 estimate?
#' parblock <- get.bg.rate(fit = cet_fit,
#'                         node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                                c('Mesoplodon','Orca','Globicephalinae')),
#'                        )
#' parblock<-c(cet_fit %chains% "R_sig2", parblock)
#' plot(parblock %select% 2 ~ parblock %select% 1)
#' #note compatibility rules
#' q.parblock1 <- cet_fit %quantiles% "R_mu"
#' q.parblock2 <- cet_fit %quantiles% list("R_mu", 0.5)
#' q.parblock3 <- cet_fit %quantiles% list("R_mu", 0.64)
#' c(parblock, q.parblock1) #chains are converted to quantiles
#' c(parblock, q.parblock1, q.parblock2) #only quantile in common is 50%/median
#' \dontrun{c(parblock, q.parblock1, q.parblock2, q.parblock3) #returns error since no shared quantiles}
#' 
#' #manipulating param_blocks
#' #calculate expected changes in average (rather than median) rates
#' Rdel <- cet_fit %chains% "R_mu" + cet_fit %chains% "R_sig2" / 2
#' #notice new name
#' names(Rdel)
#' #but we can of course still change this
#' names(Rdel) <- "R_del"
#' #recycling behavior with vectors
#' names(cet_fit %chains% 1:2 + c(1, 2))
#' #with other param_blocks
#' names(cet_fit %chains% 1:3 + cet_fit %chains% 4:5)
#' #the below, however, result in errors due to incompatible dimensions
#' \dontrun{cet_fit %chains% 1:2 + q.parblock1}
#' \dontrun{cet_fit %chains% 1:2 + array(0, c(1500, 1, 2))}
#' #special behavior for "summary" functions all, any, sum, prod, min, max, and range
#' #basically applied across rows of param_blocks for consistency
#' #"effective evolutionary time" of lineage from root to blue whale
#' sum(exp(cet_fit %chains% 1:5) * cet_fit$call$tree$edge.length[1:5])
#' #admittedly, naming could be improved here--it's weird having such long numbers in the names!
#' #could of course always do this
#' setNames(sum(exp(cet_fit %chains% 1:5) * cet_fit$call$tree$edge.length[1:5]), "blue_whale_evo_time")
#' #could also look at range of rates
#' range(cet_fit %chains% 1:5)
#' #nifty trick: calculating posterior probabilites with boolean logic and %means%
#' #are rates for edges 2 to 5 greater than that for edge 1?
#' (cet_fit %chains% 2:5 > cet_fit %chains% 1) %means% "."
#' 
#' 
#' @name param_block-class
NULL