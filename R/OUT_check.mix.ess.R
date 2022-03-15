#these could use some cleaning up...eventually

#' @export
check.mix<-function(fit){
  if(fit$sampler.control$chains>1){
    rhats<-apply(fit$chains,2,rstan::Rhat)
  }else{
    rhats<-fit%diagnostics%list('.','Rhat')
  }
  len.check<-length(rhats)
  rhats<-list(names(which(rhats<=1.01)),
              names(which(rhats<=1.05&rhats>1.01)),
              names(which(rhats>1.05)))
  lens<-lengths(rhats)
  check<-which(lens==len.check)
  if(length(check)>0){
    cat(' - all Rhats look ',c('good','okay','bad')[check],'\n')
  }else{
    for(i in seq_len(3)){
      if(lens[i]>0){
        cat(' - ',
            if(lens[i]<=5)'the following ',
            lens[i],
            ' parameters have ',
            c('good','okay','bad')[i],
            ' Rhats',
            if(lens[i]<=5)': ',
            if(lens[i]<=5) paste(rhats[[i]],collapse=', '),
            '\n')
      }
    }
  }
  if(lens[3]>0){
    cat("...\nChains haven't mixed! Try running the chains with a longer warmup period or for longer overall")
  }else if(lens[2]>0){
    cat('...\nConsider running the chains with a longer warmup period or for longer overall')
  }else{
    cat('\n\n')
  }
}

#' @export
check.ess<-function(fit){
  if(fit$sampler.control$chains>1){
    bulk.ess<-apply(fit$chains,2,rstan::ess_bulk)
    tail.ess<-apply(fit$chains,2,rstan::ess_tail)
  }else{
    bulk.ess<-fit%diagnostics%list('.','bulk_ess')
    tail.ess<-fit%diagnostics%list('.','tail_ess')
  }
  len.check<-length(bulk.ess)
  bulks<-list(names(which(bulk.ess>=400)),
              names(which(bulk.ess>=100&bulk.ess<400)),
              names(which(bulk.ess<100)))
  lens<-lengths(bulks)
  check<-which(lens==len.check)
  if(length(check)>0){
    cat(' - all bulk effective sample sizes look ',c('good','okay','bad')[check],'\n')
  }else{
    for(i in seq_len(3)){
      if(lens[i]>0){
        cat(' - ',
            if(lens[i]<=5) 'the following ',
            lens[i],
            ' parameters have ',
            c('good','okay','bad')[i],
            ' bulk effective sample sizes',
            if(lens[i]<=5)': ',
            if(lens[i]<=5) paste(bulks[[i]],collapse=', '),
            '\n')
      }
    }
  }
  if(lens[3]>0){
    cat('...\nParameter estimates unreliable! Run the chain for longer\n')
  }else if(lens[2]>0){
    cat('...\nConsider running the chain for longer for better parameter estimates\n')
  }else{
    cat('\n\n')
  }
  tails<-list(names(which(tail.ess>=400)),
              names(which(tail.ess>=100&tail.ess<400)),
              names(which(tail.ess<100)))
  lens<-lengths(tails)
  check<-which(lens==len.check)
  if(length(check)>0){
    cat(' - all tail effective sample sizes look ',c('good','okay','bad')[check],'\n')
  }else{
    for(i in seq_len(3)){
      if(lens[i]>0){
        cat(' - ',
            if(lens[i]<=5) 'the following ',
            lens[i],
            ' parameters have ',
            c('good','okay','bad')[i],
            ' tail effective sample sizes',
            if(lens[i]<=5)': ',
            if(lens[i]<=5) paste(tails[[i]],collapse=', '),
            '\n')
      }
    }
  }
  if(lens[3]>0){
    cat('...\nEstimates of parameter uncertainty unreliable! Run the chain for longer')
  }else if(lens[2]>0){
    cat('...\nConsider running the chain for longer for better estimates of parameter uncertainty')
  }else{
    cat('\n\n')
  }
}
