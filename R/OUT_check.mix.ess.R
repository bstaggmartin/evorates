#these could use some cleaning up...eventually

#' Check whether chains of a fitted evorates model adequately converged
#' 
#' 
#' This function calculates and uses the Rhat diagnostic to check for adequate convergence, or
#' mixing, of Hamiltonian Monte Carlo chains given a fitted evorates model. Rhat's should generally
#' only be slightly above 1 for all parameters. Notably, this function takes \emph{all} chains into
#' account to calculate Rhat's--diagnostics param_block arrays instead calculate Rhat's for each
#' individual chain.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}".
#' @param printlen Sets the maximum number of parameter names to report for each category of Rhat
#' values: "good", "okay", and "bad". For no maximum, set this to 0, a negative number, or infinity.
#' 
#' 
#' @return Prints to console how many parameters exhibit "good" (Rhat <= 1.01), "okay" (1.01 < 
#' Rhat <= 1.05), and "bad" (Rhat > 1.05) mixing, along with some general advice depending on
#' the worst observed Rhat. If the number of parameters in a particular category is less than or
#' equal to \code{print.len}, the function will also print the names of the parameter in that
#' category.
#' 
#' Invisibly, the function also returns a list of all calculated Rhat's, named by
#' their corresponding parameters and split into good, okay, and bad categories.
#' 
#' 
#' @details Rhat's are calculated using \code{rstan::Rhat()}. More details on this diagnostic can be found
#' in the associated documentation: \link[rstan]{Rhat}. Empirically-calculated Rhat's can
#' sometimes be slightly below 1 when chains mix well. Even though 1.01 or below is ideal, if
#' Rhat's for all parameters are around 1.05 or below, the inferred posterior distribution should
#' be pretty reliable in my experience. In any case, users should make sure all the Rhat's are okay
#' to good before running \code{combine.chains(fit)} or trusting the estimated posterior distribution!
#' 
#' In the case of high model misspecification (i.e., an evorates model does not fit the data well),
#' one could theoretically get high Rhat's and running longer chains could result in virtually
#' no benefit. I have yet to see this situation, but it is certainly possible.
#' 
#' This function complements Rhat calculation via \link{grapes-diagnostics-grapes} since it calculates
#' Rhat's for \emph{all} parameters using \emph{all} chains. It thus is tailored for diagnosing the
#' "healthiness" of a fit as a whole, rather than the healthiness of chains for individual
#' parameters/chains.
#' 
#' 
#' @family evorates diagnostic functions
#' 
#' 
#' @seealso \link[rstan]{Rhat} for the function used to calculate Rhat's. Also,
#' \link{grapes-diagnostics-grapes} for a function tailored towards diagnosing individual chains for
#' specific parameters.
#' 
#' 
#' @examples 
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #check how well the chains converged
#' Rhats <- check.mix(cet_fit) #yay!
#' Rhats #invisibly returned object
#' 
#' #note the differences between...
#' Rhats; cet_fit %diagnostics% list(".", "Rhat")
#' 
#' #generally you should make sure the Rhat's are all okay before running combine.chains()
#' 
#' #fit a new model (poorly--please never set your warmup this low in practice)
#' #this is just to see what happens with bad mixing
#' #this will take just a little bit to run...
#' fit <- fit.evorates(cet_fit$call$tree, cet_fit$call$trait.data, 
#'                     constrain.Rsig2 = TRUE,
#'                     warmup = 10, iter = 1010)
#' Rhats <- check.mix(fit) #not yay
#' 
#' 
#' @export
check.mix<-function(fit,printlen=6){
  printlen<-printlen[1]
  if(printlen<=0) printlen<-Inf
  if(fit$sampler.control$chains>1){
    tmp<-apply(fit$chains,2,rstan::Rhat)
  }else{
    tmp<-.strip.par.class(fit%diagnostics%list('.','Rhat'))
  }
  len.check<-length(tmp)
  rhats<-list(names(which(tmp<=1.01)),
              names(which(tmp<=1.05&tmp>1.01)),
              names(which(tmp>1.05)))
  lens<-lengths(rhats)
  check<-which(lens==len.check)
  if(length(check)>0){
    cat(' - all Rhats look ',c('good','okay','bad')[check],'\n')
  }else{
    for(i in seq_len(3)){
      if(lens[i]>0){
        cat(' - ',
            if(lens[i]<=printlen)'the following ',
            lens[i],
            ' parameters have ',
            c('good','okay','bad')[i],
            ' Rhats',
            if(lens[i]<=printlen)': ',
            if(lens[i]<=printlen) paste(rhats[[i]],collapse=', '),
            '\n',
            sep='')
      }
    }
  }
  if(lens[3]>0){
    cat("chains haven't mixed! Try running chain(s) with a longer warmup period or for longer overall")
  }else if(lens[2]>0){
    cat('consider running chain(s) with a longer warmup period or for longer overall')
  }
  invisible(list('good'=tmp[rhats[[1]]],'okay'=tmp[rhats[[2]]],'bad'=tmp[rhats[[3]]]))
}

#' Check whether chains adequately sampled the posterior distribution
#' 
#' 
#' This function calculates and uses the bulk/tail effective sample size (ESS) diagnostics to check for
#' adequate sampling of the posterior distribution given a fitted evorates model. Generally, ESS's should
#' at least be greater than 100 and preferably greater than 100 times the number of chains. Notably,
#' this function takes \emph{all} chains into account to calculate ESS's--diagnostics param_block arrays
#' instead calculate ESS's for each individual chain.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}".
#' @param printlen Sets the maximum number of parameter names to report for each category of bulk/tail
#' ESS values: "good", "okay", and "bad". For no maximum, set this to 0, a negative number, or infinity.
#' 
#' 
#' @return Prints to console how many parameters exhibit "good" (ESS >= 100 times number of chains),
#' "okay" (100 <= ESS  < 100 times number of chains), and "bad" (ESS < 100) sampling of both the bulks
#' (i.e., central tendency) and tails (i.e., "bounds" for lack of a better word) of the posterior
#' distribution, along with some general advice depending on the worst observed bulk/tail ESS.
#' If the number of
#' parameters in a particular category is less than or equal to \code{print.len}, the function will
#' also print the names of the parameter in that category.
#' 
#' Invisibly, the function also returns a list of all calculated bulk/tail ESS's, named by
#' their corresponding parameters and split into good, okay, and bad categories.
#' 
#' 
#' @details ESS's are calculated using \code{rstan::bulk_ess()} and \code{rstan::tail_ess()}. More
#' details on this diagnostic can be found in the associated documentation: \link[rstan]{Rhat}.
#' Even though 100 times the number of chains or above is ideal, if bulk/tail
#' ESS's for all parameters are around 200 or more, the inferred posterior distribution should
#' be pretty reliable in my experience. In any case, users should ensure make sure all the bulk
#' ESS's are okay to good before trusting parameter point estimates (e.g., means or medians),
#' and that all the tail ESS's are okay to good before trusting parameter uncertainty estimates
#' (e.g., 95\% credible intervals).
#' 
#' In the case of high model misspecification (i.e., an evorates model does not fit the data well),
#' one could theoretically get low ess's due to extremely high autocorrelation of posterior samples,
#' and running longer chains would thus yield little benefit. I have yet to see this situation, but it
#' is certainly possible.
#' 
#' This function complements ESS calculation via \link{grapes-diagnostics-grapes} since it calculates
#' ESS's for \emph{all} parameters using \emph{all} chains. It thus is tailored for diagnosing the
#' "healthiness" of a fit as a whole, rather than the healthiness of chains for individual
#' parameters/chains.
#' 
#' 
#' @family evorates diagnostic functions
#' 
#' 
#' @seealso \link[rstan]{Rhat} for the functions used to calculate bulk/tail ESS's. Also,
#' \link{grapes-diagnostics-grapes} for a function tailored towards diagnosing individual chains for
#' specific parameters.
#' 
#' 
#' @examples 
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #check how well the chains sampled the posterior
#' ess <- check.ess(cet_fit) #yay!
#' ess #invisibly returned object
#' 
#' #note the differences between...
#' ess; cet_fit %diagnostics% list(".", "ess")
#' 
#' #generally you should make sure all the ESS's are okay before trusting your parameter estimates
#' #especially in the case of bulk ESS!
#' 
#' #fit a new model (poorly--please never set your iterations this low in practice)
#' #this is just to see what happens with bad sampling
#' #this will take just a little bit to run...
#' fit <- fit.evorates(cet_fit$call$tree, cet_fit$call$trait.data,
#'                     constrain.Rsig2 = TRUE,
#'                     warmup = 1000, iter = 1040)
#' ess <- check.ess(fit) #not yay
#' 
#' 
#' @export
check.ess<-function(fit,printlen=6){
  printlen<-printlen[1]
  if(printlen<=0) printlen<-Inf
  nchains<-fit$sampler.control$chains
  if(nchains>1){
    bulk.ess<-apply(fit$chains,2,rstan::ess_bulk)
    tail.ess<-apply(fit$chains,2,rstan::ess_tail)
  }else{
    bulk.ess<-.strip.par.class(fit%diagnostics%list('.','bulk_ess'))
    tail.ess<-.strip.par.class(fit%diagnostics%list('.','tail_ess'))
  }
  len.check<-length(bulk.ess)
  bulks<-list(names(which(bulk.ess>=100*nchains)),
              names(which(bulk.ess>=100&bulk.ess<100*nchains)),
              names(which(bulk.ess<100)))
  lens<-lengths(bulks)
  check<-which(lens==len.check)
  if(length(check)>0){
    cat(' - all bulk effective sample sizes look ',c('good','okay','bad')[check],'\n')
  }else{
    for(i in seq_len(3)){
      if(lens[i]>0){
        cat(' - ',
            if(lens[i]<=printlen) 'the following ',
            lens[i],
            ' parameters have ',
            c('good','okay','bad')[i],
            ' bulk effective sample sizes',
            if(lens[i]<=printlen)': ',
            if(lens[i]<=printlen) paste(bulks[[i]],collapse=', '),
            '\n')
      }
    }
  }
  if(lens[3]>0){
    cat('parameter estimates unreliable! Run chain(s) for longer\n')
  }else if(lens[2]>0){
    cat('consider running chain(s) for longer for better parameter estimates\n')
  }
  tails<-list(names(which(tail.ess>=100*nchains)),
              names(which(tail.ess>=100&tail.ess<100*nchains)),
              names(which(tail.ess<100)))
  lens<-lengths(tails)
  check<-which(lens==len.check)
  if(length(check)>0){
    cat(' - all tail effective sample sizes look ',c('good','okay','bad')[check],'\n')
  }else{
    for(i in seq_len(3)){
      if(lens[i]>0){
        cat(' - ',
            if(lens[i]<=printlen) 'the following ',
            lens[i],
            ' parameters have ',
            c('good','okay','bad')[i],
            ' tail effective sample sizes',
            if(lens[i]<=printlen)': ',
            if(lens[i]<=printlen) paste(tails[[i]],collapse=', '),
            '\n')
      }
    }
  }
  if(lens[3]>0){
    cat('estimates of parameter uncertainty unreliable! Run chain(s) for longer')
  }else if(lens[2]>0){
    cat('consider running chain(s) for longer for better estimates of parameter uncertainty')
  }
  invisible(list('good_bulk'=bulk.ess[bulks[[1]]],'good_tail'=tail.ess[tails[[1]]],
                 'okay_bulk'=bulk.ess[bulks[[2]]],'okay_tail'=tail.ess[tails[[2]]],
                 'bad_bulk'=bulk.ess[bulks[[3]]],'bad_tail'=tail.ess[tails[[3]]]))
}
