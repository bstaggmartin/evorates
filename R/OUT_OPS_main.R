####MAIN FUNCTIONS####

#' Extract posterior distribution samples from a fitted evorates model
#'
#'
#' This operator extracts samples from the posterior distributions (i.e., iterations from the chains)
#' for particular parameters from an \code{evorates_fit} object or \code{param_block} array.
#'
#'
#' @param fit An object of class "\code{evorates_fit}" or "\code{param_block}". If \code{fit} is a
#' \code{param_block} array, it must have a \code{param_type} of "\code{chains}".
#' @param select A list with two elements (2nd element is optional):
#' \itemize{
#' \item{A character or numeric vector for selecting parameters. If a character vector, entries are
#' matched to parameter names using regular expressions, with one key exception: if any entries are
#' \emph{exact} matches to sampling parameter names, these will be used to select parameters from
#' the \code{sampler.params} element of \code{fit} instead, if it exists (see details). If a 
#' numeric vector, entries are matched to edge indices to select branchwise rate parameters;
#' these can be negative to instead \emph{exclude} branchwise rates. If no branchwise rate parameters
#' are found, a number \code{i} will instead select the \code{i}th parameter in \code{fit}.}
#' \item{A numeric vector for selecting samples from the distribution (i.e., iterations
#' from the chains). If unsupplied, all samples are selected.}
#' }
#' 
#' 
#' @return An array of class "\code{param_block}" with a \code{param_type} of "\code{chains}".
#' The dimension of these arrays will generally go in the order of iterations, then parameters,
#' then chains. Any dimensions of length 1 are collapsed and stored as attributes.
#' 
#' 
#' @details In the case that a numeric vector is provided to select parameters and no parameters
#' with names following the pattern "\code{R_i}" are found, the function then looks for the pattern
#' "\code{Rdev_i}", then "\code{uncent_Rdev_i}". If neither of these are found, then it finally
#' defaults to selecting the \code{i}th parameters. If a single parameter name involves multiple
#' "\code{R_i}" patterns, the first pattern is always used (e.g., \code{R_1\%-\%R_2} would correspond
#' to \code{1}). Similarly, if multiple parameters match to a single number, the first parameter is
#' always used (similar to behavior of \link{match}).
#' 
#' The \code{sampler.params} element of \code{fit} always includes 9 parameters, named:
#' "\code{accept_stat__}", "\code{treedepth__}", "\code{stepsize__}", "\code{divergent__}",
#' "\code{n_leapfrog__}", "\code{energy__}", "\code{prior__}", "\code{lik__}", and "\code{post__}".
#' This tends to include both warmup and non-warmup samples, but warmup samples are
#' automatically excluded if both sampler parameters and normal parameters are selected by this
#' function. Most of these parameters are rather technical quantities used to tune the
#' Hamiltonian Monte Carlo sampler run by Stan (see Stan manual for further details on what they
#' mean). Generally, users will only want to look at the last 3 parameters, which give the
#' (log) prior probability, likelihood, and posterior probability, respectively, of sampled
#' parameters. Note that these are on the sampling scale and will differ from those for the
#' originally-scaled data by a constant. Also, the posterior probability will be affected
#' by what \code{lik.power} was set to for fitting the model. In some cases, users may also
#' wish to look at what parameter values are associated with divergent transitions (i.e.,
#' iterations where \code{divergent__ = 1}), which indicate regions of parameter space where
#' the sampler got "stuck", yielding potentially misleading posterior distribution estimates.
#' 
#' 
#' @family param_block operators
#'  
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #extracting directly from evorates fit
#' cet_fit %chains% "R_mu"
#' #regular expressions
#' cet_fit %chains% "R"
#' #using . is a quick way to extract ALL parameters!
#' cet_fit %chains% "."
#' #numeric index-based selection
#' cet_fit %chains% 1
#' cet_fit %chains% -1
#' #select particular samples
#' cet_fit %chains% list("R_mu", 1)
#' #getting sampler parameters
#' cet_fit %chains% "lik__"
#' #note warmup samples automatically excluded from "lik__" if combined with "R_mu"
#' cet_fit %chains% c("R_mu", "lik__")
#' 
#' #extracting from a param_block array
#' par <- get.bg.rate(fit = cet_fit,
#'                    node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                           c('Mesoplodon','Orca','Globicephalinae')),
#'                    )
#' par %chains% list("Mesoplodon", 36)
#' #note change in numeric index behavior
#' par %chains% 1
#' 
#' 
#' @export
`%chains%`<-function(fit,select){
  .proc.op('chains',fit,select,deparsed.select=deparse(substitute(select)))
}

#' Extract posterior distribution quantiles from a fitted evorates model
#'
#'
#' This operator extracts quantiles from the posterior distributions for particular parameters from
#' an \code{evorates_fit} object or \code{param_block} array. This can be used to, for example, get
#' posterior medians or credible intervals.
#'
#'
#' @param fit An object of class "\code{evorates_fit}" or "\code{param_block}". If \code{fit} is a
#' \code{param_block} array, it must have a \code{param_type} of "\code{chains}" or
#' "\code{quantiles}".
#' @param select A list with two elements (2nd element is optional):
#' \itemize{
#' \item{A character or numeric vector for selecting parameters. If a character vector, entries are
#' matched to parameter names using regular expressions, with one key exception: if any entries are
#' \emph{exact} matches to sampling parameter names, these will be used to select parameters from
#' the \code{sampler.params} element of \code{fit} instead, if it exists (see details). If a 
#' numeric vector, entries are matched to edge indices to select branchwise rate parameters;
#' these can be negative to instead \emph{exclude} branchwise rates. If no branchwise rate parameters
#' are found, a number \code{i} will instead select the \code{i}th parameter in \code{fit}.}
#' \item{A numeric, integer, or character vector for specifying quantiles. If unsupplied, default quantiles
#' are selected. Default quantiles are are set to the quantiles already stored in \code{fit}, if they exist,
#' otherwise they are set to the lower bound of the 95\% credible interval, median, and upper bound (i.e., 2.5\%,
#' 50\%, and 97.5\%). If a numeric vector, this will simply be the quantile to extract (i.e., 0.1 would be
#' the 10\% quantile). If an integer vector, an integer \code{i} will extract the \code{i}th default quantile.
#' If a character vector, numbers as taken as percents (i.e., 10\% or even 10 would correspond to 0.1).
#' Quantiles greater than 100\% or less than 0\% are rounded to 100\% and 0\%, respectively.}
#' }
#' 
#' 
#' @return An array of class "\code{param_block}" with a \code{param_type} of "\code{quantiles}".
#' The dimension of these arrays will generally go in the order of quantiles, then parameters,
#' then chains. Any dimensions of length 1 are collapsed and stored as attributes. If \code{fit}
#' is a \code{chains param_block} array, parameters are automatically renamed
#' "\code{quantiles(<parameter name>)}" to help keep track of parameter manipulations.
#' 
#' 
#' @details In the case that a numeric vector is provided to select parameters and no parameters
#' with names following the pattern "\code{R_i}" are found, the function then looks for the pattern
#' "\code{Rdev_i}", then "\code{uncent_Rdev_i}". If neither of these are found, then it finally
#' defaults to selecting the \code{i}th parameters. If a single parameter name involves multiple
#' "\code{R_i}" patterns, the first pattern is always used (e.g., \code{R_1\%-\%R_2} would correspond
#' to \code{1}). Similarly, if multiple parameters match to a single number, the first parameter is
#' always used (similar to behavior of \link{match}).
#' 
#' The \code{sampler.params} element of \code{fit} always includes 9 parameters, named:
#' "\code{accept_stat__}", "\code{treedepth__}", "\code{stepsize__}", "\code{divergent__}",
#' "\code{n_leapfrog__}", "\code{energy__}", "\code{prior__}", "\code{lik__}", and "\code{post__}".
#' This tends to include both warmup and non-warmup samples, but warmup samples are
#' automatically excluded if both sampler parameters and normal parameters are selected by this
#' function. Most of these parameters are rather technical quantities used to tune the
#' Hamiltonian Monte Carlo sampler run by Stan (see Stan manual for further details on what they
#' mean). Generally, users will only want to look at the last 3 parameters, which give the
#' (log) prior probability, likelihood, and posterior probability, respectively, of sampled
#' parameters. Note that these are on the sampling scale and will differ from those for the
#' originally-scaled data by a constant. Also, the posterior probability will be affected
#' by what \code{lik.power} was set to for fitting the model. In some cases, users may also
#' wish to look at what parameter values are associated with divergent transitions (i.e.,
#' iterations where \code{divergent__ = 1}), which indicate regions of parameter space where
#' the sampler got "stuck", yielding potentially misleading posterior distribution estimates.
#' 
#' This function uses a custom quantiles estimation function that yields the same results as
#' \link{quantile} with \code{type = 7}. \code{NA}s are always ignored when estimating
#' quantiles, unless all samples for a parameter are \code{NA}, in which case all
#' quantile estimates will also be \code{NA}.
#' 
#' 
#' @family param_block operators
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #extracting directly from evorates fit
#' cet_fit %quantiles% "R_mu"
#' #regular expressions
#' cet_fit %quantiles% "R"
#' #using . is a quick way to extract ALL parameters!
#' cet_fit %quantiles% "."
#' #numeric index-based selection
#' cet_fit %quantiles% 1
#' cet_fit %quantiles% -1
#' #select particular quantiles
#' #e.g., medians
#' cet_fit %quantiles% list("R_mu", 0.5)
#' cet_fit %quantiles% list("R_mu", 2L)
#' cet_fit %quantiles% list("R_mu", "50%")
#' #e.g., 95% credible intervals
#' cet_fit %quantiles% list("R_mu", c(0.025, 0.975))
#' cet_fit %quantiles% list("R_mu", c(1L, 3L))
#' cet_fit %quantiles% list("R_mu", c("2.5%", "97.5%"))
#' #getting sampler parameters
#' cet_fit %quantiles% "lik__"
#' #note warmup samples automatically excluded from "lik__" if combined with "R_mu"
#' cet_fit %quantiles% c("R_mu", "lik__")
#' 
#' #extracting from a param_block array
#' par <- get.bg.rate(fit = cet_fit,
#'                    node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                           c('Mesoplodon','Orca','Globicephalinae')),
#'                    )
#' par %quantiles% list("Mesoplodon", 0.5)
#' #note change in numeric index behavior
#' par %quantiles% 1
#' 
#' 
#' @export
`%quantiles%`<-function(fit,select){
  .proc.op('quantiles',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','quantiles'))
}

#' Extract posterior distribution means from a fitted evorates model
#'
#'
#' This operator extracts means from the posterior distributions for particular parameters from
#' an \code{evorates_fit} object or \code{param_block} array. This can be used to, for example, get
#' a point estimate for a parameter or posterior probabilities.
#'
#'
#' @param fit An object of class "\code{evorates_fit}" or "\code{param_block}". If \code{fit} is a
#' \code{param_block} array, it must have a \code{param_type} of "\code{chains}" or
#' "\code{means}".
#' @param select A character or numeric vector for selecting parameters. If a character vector, entries are
#' matched to parameter names using regular expressions, with one key exception: if any entries are
#' \emph{exact} matches to sampling parameter names, these will be used to select parameters from
#' the \code{sampler.params} element of \code{fit} instead, if it exists (see details). If a 
#' numeric vector, entries are matched to edge indices to select branchwise rate parameters;
#' these can be negative to instead \emph{exclude} branchwise rates. If no branchwise rate parameters
#' are found, a number \code{i} will instead select the \code{i}th parameter in \code{fit}.
#' 
#' 
#' @return An array of class "\code{param_block}" with a \code{param_type} of "\code{means}".
#' The dimension of these arrays will generally go in the order of means, then parameters,
#' then chains. Any dimensions of length 1 are collapsed and stored as attributes. While
#' the means dimension is always of length 1 and will be collapsed in most cases, it is
#' sometimes kept and labelled as an iterations dimension for compatibility reasons.
#' If \code{fit} is a \code{chains param_block} array, parameters are automatically renamed
#' "\code{means(<parameter name>)}" to help keep track of parameter manipulations.
#' 
#' 
#' @details In the case that a numeric vector is provided to select parameters and no parameters
#' with names following the pattern "\code{R_i}" are found, the function then looks for the pattern
#' "\code{Rdev_i}", then "\code{uncent_Rdev_i}". If neither of these are found, then it finally
#' defaults to selecting the \code{i}th parameters. If a single parameter name involves multiple
#' "\code{R_i}" patterns, the first pattern is always used (e.g., \code{R_1\%-\%R_2} would correspond
#' to \code{1}). Similarly, if multiple parameters match to a single number, the first parameter is
#' always used (similar to behavior of \link{match}).
#' 
#' The \code{sampler.params} element of \code{fit} always includes 9 parameters, named:
#' "\code{accept_stat__}", "\code{treedepth__}", "\code{stepsize__}", "\code{divergent__}",
#' "\code{n_leapfrog__}", "\code{energy__}", "\code{prior__}", "\code{lik__}", and "\code{post__}".
#' This tends to include both warmup and non-warmup samples, but warmup samples are
#' automatically excluded if both sampler parameters and normal parameters are selected by this
#' function. Most of these parameters are rather technical quantities used to tune the
#' Hamiltonian Monte Carlo sampler run by Stan (see Stan manual for further details on what they
#' mean). Generally, users will only want to look at the last 3 parameters, which give the
#' (log) prior probability, likelihood, and posterior probability, respectively, of sampled
#' parameters. Note that these are on the sampling scale and will differ from those for the
#' originally-scaled data by a constant. Also, the posterior probability will be affected
#' by what \code{lik.power} was set to for fitting the model. In some cases, users may also
#' wish to look at what parameter values are associated with divergent transitions (i.e.,
#' iterations where \code{divergent__ = 1}), which indicate regions of parameter space where
#' the sampler got "stuck", yielding potentially misleading posterior distribution estimates.
#' 
#' \code{NA}s are always ignored when estimating means, unless all samples for a parameter
#' are \code{NA}, in which case all mean estimates will also be \code{NA}.
#' 
#' Technically, \code{select} can be a list of two elements for compatibility reasons, but
#' the second element is always ignored.
#' 
#' 
#' @family param_block operators
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #extracting directly from evorates fit
#' cet_fit %means% "R_mu"
#' #regular expressions
#' cet_fit %means% "R"
#' #using . is a quick way to extract ALL parameters!
#' cet_fit %means% "."
#' #numeric index-based selection
#' cet_fit %means% 1
#' cet_fit %means% -1
#' #getting sampler parameters
#' cet_fit %means% "lik__"
#' #note warmup samples automatically excluded from "lik__" if combined with "R_mu"
#' cet_fit %means% c("R_mu", "lik__")
#' 
#' #extracting from a param_block array
#' par <- get.bg.rate(fit = cet_fit,
#'                    node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                           c('Mesoplodon','Orca','Globicephalinae')),
#'                    )
#' par %means% "Mesopldon"
#' #note change in numeric index behavior
#' par %means% 1
#' 
#' #getting posterior probabilities
#' #basically, we can transform a parameter to consist of 0s or 1s, based on some logic statement
#' #for example, we can recode R_mu to be 1 if it's above 0 (positive trend/late burst) and 0 otherwise (negative trend/early burst)
#' Rmu <- cet_fit %chains% "R_mu" > 0
#' #note that 1s are indicated by TRUE and 0s by FALSE
#' Rmu
#' #now, if we take the average of these 0s and 1s, we actually get the proportion of 1s in this transformed parameter
#' #in other words, this is a way of calculating the probability that R_mu is greater than 0 in the posterior distribution
#' Rmu %means% "."
#' 
#' 
#' @export
`%means%`<-function(fit,select){
  .proc.op('means',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','means'))
}

#' Extract posterior distribution diagnostics from a fitted evorates model
#'
#'
#' This operator extracts diagnostic summaries of the posterior distributions for particular
#' parameters from an \code{evorates_fit} object or \code{param_block} array. This can be used to
#' check that a posterior distribution was sampled thoroughly.
#'
#'
#' @param fit An object of class "\code{evorates_fit}" or "\code{param_block}". If \code{fit} is a
#' \code{param_block} array, it must have a \code{param_type} of "\code{chains}" or
#' "\code{diagnostics}".
#' @param select A list with two elements (2nd element is optional):
#' \itemize{
#' \item{A character or numeric vector for selecting parameters. If a character vector, entries are
#' matched to parameter names using regular expressions, with one key exception: if any entries are
#' \emph{exact} matches to sampling parameter names, these will be used to select parameters from
#' the \code{sampler.params} element of \code{fit} instead, if it exists (see details). If a 
#' numeric vector, entries are matched to edge indices to select branchwise rate parameters;
#' these can be negative to instead \emph{exclude} branchwise rates. If no branchwise rate parameters
#' are found, a number \code{i} will instead select the \code{i}th parameter in \code{fit}.}
#' \item{A character or numeric vector for specifying diagnostics. If unsupplied, default diagnostics
#' are selected. Default diagnostics are are set to the diagnostics already stored in \code{fit}, if they exist,
#' otherwise they are set to initial values ("\code{inits}"), bulk effective sample sizes ("\code{bulk_ess}"),
#' tail effective samples sizes ("\code{tail_ess}"), and Rhat's ("\code{Rhats}") (see details for what these
#' quantities mean). If a character vector, entries are matched to diagnostic names using regular expressions.
#' If a numeric vector, a number \code{i} extract the \code{i}th default diagnostic.}
#' }
#' 
#' 
#' @return An array of class "\code{param_block}" with a \code{param_type} of "\code{diagnostics}".
#' The dimension of these arrays will generally go in the order of diagnostics, then parameters,
#' then chains. Any dimensions of length 1 are collapsed and stored as attributes. If \code{fit}
#' is a \code{chains param_block} array, parameters are automatically renamed
#' "\code{diagnostics(<parameter name>)}" to help keep track of parameter manipulations.
#' 
#' 
#' @details In the case that a numeric vector is provided to select parameters and no parameters
#' with names following the pattern "\code{R_i}" are found, the function then looks for the pattern
#' "\code{Rdev_i}", then "\code{uncent_Rdev_i}". If neither of these are found, then it finally
#' defaults to selecting the \code{i}th parameters. If a single parameter name involves multiple
#' "\code{R_i}" patterns, the first pattern is always used (e.g., \code{R_1\%-\%R_2} would correspond
#' to \code{1}). Similarly, if multiple parameters match to a single number, the first parameter is
#' always used (similar to behavior of \link{match}).
#' 
#' The \code{sampler.params} element of \code{fit} always includes 9 parameters, named:
#' "\code{accept_stat__}", "\code{treedepth__}", "\code{stepsize__}", "\code{divergent__}",
#' "\code{n_leapfrog__}", "\code{energy__}", "\code{prior__}", "\code{lik__}", and "\code{post__}".
#' This tends to include both warmup and non-warmup samples, but warmup samples are
#' automatically excluded if both sampler parameters and normal parameters are selected by this
#' function. Most of these parameters are rather technical quantities used to tune the
#' Hamiltonian Monte Carlo sampler run by Stan (see Stan manual for further details on what they
#' mean). Generally, users will only want to look at the last 3 parameters, which give the
#' (log) prior probability, likelihood, and posterior probability, respectively, of sampled
#' parameters. Note that these are on the sampling scale and will differ from those for the
#' originally-scaled data by a constant. Also, the posterior probability will be affected
#' by what \code{lik.power} was set to for fitting the model. In some cases, users may also
#' wish to look at what parameter values are associated with divergent transitions (i.e.,
#' iterations where \code{divergent__ = 1}), which indicate regions of parameter space where
#' the sampler got "stuck", yielding potentially misleading posterior distribution estimates.
#' 
#' Initial values are simply the starting values for each parameter. If fit is a 
#' \code{chains param_block} array, the initial values will always be \code{NA} since it
#' is generally unknown if a \code{chains param_block} array contains all warmup
#' iterations (and they generally do not, in any case). Bulk and tail effective sample
#' sizes (ESS) measure how well the centers and tails of posterior
#' distributions are sampled, respectively. A high enough bulk ESS (~100 or more, but ideally
#' 100 times the number of chains) indicates that posterior median/mean estimates are reliable,
#' while a high enough tail ESS (again, ~100 or 100 times the number of chains) indicates that
#' posterior credible interval estimates are reliable. Rhat measures how well a chain
#' samples "all" parts of a posterior distribution, so to speak, acting as a
#' measure of "mixing" (i.e., the sampler converging on the actual posterior distribution).
#' In this case, lower Rhat's are better, with ~1 being the "best" (a good cutoff is 
#' 1.05 or lower, with 1.01 being ideal). This function calculates ESS and Rhat using
#' functions from the \code{rstan} package, and more details about these quantities can
#' be found in the Stan manual.
#' 
#' Notably, this function calculates ESS and Rhat's afor \emph{each chain}. This
#' is helpful for diagnosing problematic behavior in individual chains/parameters, but may not
#' reflect the "healthiness" of the entire fit. To check the ESS and Rhat's for all chains
#' combined, you can use \link{check.ess} and \link{check.mix}. You should consider your fit
#' "good to go" if the \emph{highest} combined Rhat doesn't exceed 1.01/1.05 and the
#' \emph{lowest} combined ESS doesn't drop below 100/100 times the number of chains. Also,
#' \link{combine.chains} should only be used if chains properly converged, as indicated
#' by a low enough maximum Rhat.
#' 
#' 
#' @family param_block operators
#'
#' 
#' @seealso See \link{check.ess} and \link{check.mix} for functions tailored toward
#' diagnosing the "healthiness" of a fit as a whole. For the \code{rstan} functions
#' this function relies on, see \link[rstan]{Rhat}.
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #extracting directly from evorates fit
#' cet_fit %diagnostics% "R_mu"
#' #regular expressions
#' cet_fit %diagnostics% "R"
#' #using . is a quick way to extract ALL parameters!
#' cet_fit %diagnostics% "."
#' #numeric index-based selection
#' cet_fit %diagnostics% 1
#' cet_fit %diagnostics% -1
#' #select particular diagnostics
#' cet_fit %diagnostics% list("R_mu", "ess")
#' cet_fit %diagnostics% list("R_mu", 1)
#' cet_fit %diagnostics% list("R_mu", -1)
#' #getting sampler parameters
#' cet_fit %diagnostics% "lik__"
#' #warmup samples are still automatically excluded from "lik__" if combined with "R_mu"
#' #BUT this only affects bulk_ess, tail_ess, and Rhat
#' #the function DOES extract the actual initial value of the likelihood!
#' cet_fit %diagnostics% c("R_mu", "lik__")
#' 
#' #extracting from a param_block array
#' par <- get.bg.rate(fit = cet_fit,
#'                    node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                           c('Mesoplodon','Orca','Globicephalinae')),
#'                    )
#' par %diagnostics% list("Mesoplodon", -1)
#' #note that inits returns NA
#' par %diagnostics% list("Mesoplodon", "inits")
#' #note change in numeric index behavior
#' par %diagnostics% 1
#' 
#' 
#' @export
`%diagnostics%`<-function(fit,select){
  .proc.op('diagnostics',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','diagnostics'))
}

####SIMPLIFIED FUNCTION NAMES####

#' @rdname grapes-chains-grapes
#' 
#' 
#' @export
`%c%`<-function(fit,select){
  .proc.op('chains',fit,select,deparsed.select=deparse(substitute(select)))
}

#' @rdname grapes-quantiles-grapes
#' 
#' 
#' @export
`%q%`<-function(fit,select){
  .proc.op('quantiles',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','quantiles'))
}

#' @rdname grapes-means-grapes
#' 
#' 
#' @export
`%m%`<-function(fit,select){
  .proc.op('means',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','means'))
}

#' @rdname grapes-diagnostics-grapes
#' 
#' 
#' @export
`%d%`<-function(fit,select){
  .proc.op('diagnostics',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','diagnostics'))
}

####SELECT####

#' Extract subsets from a parameter block
#' 
#' 
#' This operator extracts particular parameter quantities from a \code{param_block} array.
#' 
#' 
#' @param x An object of class "\code{param_block}" with any \code{param_type}.
#' @param select A list with two elements (2nd element is optional):
#' \itemize{
#' \item{A character or numeric vector for selecting parameters. If a character vector, entries are
#' matched to parameter names using partial matching. If a numeric vector, a number \code{i} will
#' select the \code{i}th parameter in \code{x}.}
#' \item{A numeric or character vector for specifying samples, quantiles, or diagnostics
#' based on the \code{param_type} of \code{x}. This is ignored if \code{x} is of \code{param_type}
#' "\code{means}". If unsupplied, defaults to including all samples/quantiles/diagnostics. If a
#' numeric vector, a number \code{i} will select the \code{i}th sample/quantile/diagnostic in
#' \code{x}. If a character vector, entries are matched to sample/quantile/diagnostic names using
#' partial matching. Note that samples don't have names in the case \code{chains param_block}
#' arrays--in this case, character vectors are automatically converted to numeric vectors.}
#' }
#' 
#' 
#' @return An array of class "\code{param_block}" with the same \code{param_type} as \code{x}.
#' The dimension of these arrays will generally go in the order of samples/quantiles/diagnostics,
#' then parameters, then chains. Any dimensions of length 1 are collapsed and stored as attributes.
#' 
#' 
#' @details This function can be used instead of the other \code{param_block} operators to emulate
#' a more "traditional" means of subsetting arrays comparable to \code{data.frame} objects. Since it
#' relies on partial matching, rather than regular expressions, any name will extract either 0 or 1
#' parameters, samples, quantiles, or diagnostics, depending on whether the name had a clear match
#' or not. Furthermore, any number \code{i} will extract that \code{i}th parameter, sample,
#' quantile, or diagnostic, rather than using "special" rules (branchwise rate searching, etc.).
#' Note that out of bound numbers are ignored, rather than returning \code{NA}'s. Other
#' \code{param_block} conventions are still enforced: collapsed dimensions have their information
#' stored as attributes and specific chains cannot be selected. You still have to use
#' \link{select.chains} on \emph{entire} \code{evorates_fit} objects for selecting particular chains.
#' 
#' 
#' @family param_block operators
#' 
#' 
#' @seealso This function relies on \link{pmatch} for partial matching.
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #get a chains param_block array
#' par <- cet_fit %chains% "R"
#' #note partial match vs. regular expressions behavior
#' par %select% "R_1"
#' par %select% list(NULL)
#' #note numeric index selection behavior change too
#' par %select% 1
#' 
#' #get a quantiles param_block array
#' par <- cet_fit %quantiles% "R"
#' #note numeric index selection behavior change for quantiles
#' par %select% list(1, 1)
#' 
#' #get a means param_block array
#' par <- cet_fit %means% "R"
#' #note 2nd element of list is ignored
#' par %select% list(1, 2)
#' 
#' #get a diangostics param_block array
#' par <- cet_fit %diagnostics% "R"
#' #note partial vs. regular expressions behavior for selecting diagnostics
#' par %select% list(1, "bulk")
#' par %select% list(1, "ess") #returns empty array with warning
#' 
#' 
#' @export
`%select%`<-function(x,select){
  .proc.select(x,select,deparsed.select=deparse(substitute(select)))
}

#' @rdname grapes-select-grapes
#' 
#' 
#' @export
`%s%`<-function(x,select){
  .proc.select(x,select,deparsed.select=deparse(substitute(select)))
}