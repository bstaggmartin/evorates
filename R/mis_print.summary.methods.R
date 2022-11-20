.print.info<-function(fit,printlen){
  nobs<-nrow(fit$call$trait.data)
  ntips<-Ntip(fit)
  nnodes<-Nnode(fit)
  hgt<-max(node.depth.edgelength(fit$call$tree))
  if(ntips<=printlen){
    tip.nms<-paste0(fit$call$tree$tip.label,collapse=", ")
  }else{
    tip.nms<-paste0(c(fit$call$tree$tip.label[1:printlen],"..."),collapse=", ")
  }
  node.flag<-FALSE
  if(!is.null(fit$call$tree$node.label)){
    node.flag<-TRUE
    if(nnodes<=printlen){
      node.nms<-paste0(fit$call$tree$node.label,collapse=", ")
    }else{
      node.nms<-paste0(c(fit$call$tree$node.label[1:printlen],"..."),collapse=", ")
    }
  }
  cat("An evorates model fitted to comparative data consisting of ",
      nobs," observations of ",
      ntips," tips/",nnodes," internal nodes spanning ",format(hgt,digits=3)," unit(s) of time",".\n\n",
      "Tip labels:\n\t",tip.nms,
      if(node.flag) "\n\nNode labels:\n\t",if(node.flag) node.nms,
      "\n\n",sep="")
  nchains<-fit$sampler.control$chains
  niter<-fit$sampler.control$iter
  warmup<-fit$sampler.control$warmup
  warmup.flag<-nrow(fit$chains)==fit$sampler.control$iter-fit$sampler.control$warmup
  thin<-fit$sampler.control$thin
  cat(" - currently consists of ",
      nchains," ",niter,
      "-iteration chain(s) with ",
      warmup,
      " warmup iterations",
      if(warmup.flag) " (warmup excluded from parameter posterior samples)",
      " and thinning rate of ",
      thin,
      "\n",
      sep="")
  nchains<-fit$sampler.control$call.chains
  niter<-fit$sampler.control$call.iter
  warmup<-fit$sampler.control$call.warmup
  thin<-fit$sampler.control$call.thin
  cat(" - originally fitted via ",
      nchains," ",niter,
      "-iteration chain(s) with ",
      warmup,
      " warmup iterations and thinning rate of ",
      thin,
      "\n",
      sep="")
}

.print.main.params<-function(fit,quants,max.params,max.chains){
  nchains<-fit$sampler.control$chains
  chain.flag<-param.flag<-FALSE
  if(nchains>max.chains){
    fit<-select.chains(fit,1:max.chains)
    chain.flag<-TRUE
  }
  params<-names(fit$chains)
  Rdev.flag<-any(grepl("^Rdev_\\d+$",params))
  nparams<-length(params)
  if(nparams>max.params){
    params<-params[1:max.params]
    param.flag<-TRUE
  }
  params<-paste0("^",params,"$")
  params<-fit%quantiles%list(params,quants)
  cat("\n")
  print(format(params,digits=3),quote=FALSE)
  attrs<-unlist(setNames(lapply(c("quantiles","parameters","chains"),function(ii) attr(params,ii)),
                         c("quantiles","parameters","chains")))
  if(length(attrs)){
    cat("\n",
        paste0(substr(names(attrs),1,nchar(names(attrs))-1),": ",attrs,collapse="\n"),
        "\n\n",
        sep="")
  }
  if(chain.flag|param.flag){
    cat("   [ omitted ",
        if(chain.flag) paste0(nchains-max.chains," chains"),
        if(chain.flag&param.flag) " and ",
        if(param.flag) paste0(nparams - max.params," parameters"),
        if(Rdev.flag) paste0(", including rate deviations"),
        " ]",
        sep="")
  }
}

#' @export
print.evorates_fit<-function(fit,printlen=6,quants=c(0.025,0.5,0.975),
                             max.params=printlen,max.chains=2,...){
  .print.info(fit,printlen)
  .print.main.params(fit,quants,max.params,max.chains)
}

#' @export
summary.evorates_fit<-function(fit,printlen=6,quants=c(0.025,0.5,0.975),
                               max.params=printlen,max.chains=2,
                               remove.trend=TRUE,geometric=TRUE,max.edges=printlen,...){
  .print.info(fit,printlen)
  check.mix(fit)
  check.ess(fit)
  tmp<-c("decisive","strong","substantial")
  tmp<-c(tmp,"none",rev(tmp))
  if(fit$call$trend){
    post.prob<-mean(as.vector((fit%chains%"^R_mu$">0)%means%"."))
    code<-findInterval(post.prob,c(0.001,0.025,0.05,0.95,0.975,0.999))+1
    if(code!=4){
      cat(" - posterior probability of increasing trend is ",
          format(post.prob,digits=3),
          ", indicating ",
          tmp[code],
          " evidence for ",
          if(code<4) "a decreasing" else "an increasing",
          " trend in rates over time\n",
          sep="")
    }
  }
  if(!fit$call$constrain.Rsig2){
    sd<-get.sd(fit)
    code<-findInterval(sd,c(1/100,1/10,1/3,3,10,100))+1
    if(code!=4){
      cat(" - Savage-Dickey ratio for no rate heterogeneity is ",
          format(sd,digits=3),
          ", indicating ",
          tmp[code],
          " evidence ",
          if(code<4) "for" else "against",
          " rate heterogeneity\n",
          sep="")
    }
  }
  if(fit$call$trend|!fit$call$constrain.Rsig2){
    bg<-get.bg.rate(fit,remove.trend=remove.trend,geometric=geometric,keep.R=TRUE)
    post.probs<-(bg$R-bg$bg_rate>0)%means%"."
    if(length(dim(post.probs))==2){
      post.probs<-rowMeans(post.probs)
    }
    sigs.flag<-FALSE
    sigs<-post.probs<0.05
    if(any(sigs)){
      sigs.flag<-TRUE
      nsigs<-sum(sigs)
      sig.nms<-unname(which(sigs))
      if(nsigs<=max.edges){
        sig.nms<-paste0(sig.nms,collapse=", ")
      }else{
        sig.nms<-paste0(c(sig.nms[1:max.edges],"..."),collapse=", ")
      }
      cat(" - ",nsigs," edge(s) exhibit(s) anomalously low average rate(s) (posterior probability > 0.95):\n\t",
          sig.nms,"\n",
          sep="")
    }
    sigs<-post.probs>0.95
    if(any(sigs)){
      sigs.flag<-TRUE
      nsigs<-sum(sigs)
      sig.nms<-unname(which(sigs))
      if(nsigs<=max.edges){
        sig.nms<-paste0(sig.nms,collapse=", ")
      }else{
        sig.nms<-paste0(c(sig.nms[1:max.edges],"..."),collapse=", ")
      }
      cat(" - ",nsigs," edge(s) exhibit(s) anomalously high average rate(s) (posterior probability > 0.95):\n\t",
          sig.nms,"\n",
          sep="")
    }
    if(sigs.flag){
      cat("   [ rates ",
          if(remove.trend) "adjusted" else "unadjusted",
          " for trends in rates over time and compared to ",
          if(geometric) "geometrically-" else "arithmetically-",
          "averaged background rate ]\n",
          " - plot the fitted model to visualize where these edges are phylogenetically located!\n",
          sep="")
    }
  }
  #indicate number of branches with significant rate deviations
  .print.main.params(fit,quants,max.params,max.chains)
}

