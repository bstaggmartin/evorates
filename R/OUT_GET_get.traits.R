#' Sample posterior trait values at nodes given a fitted evorates model
#' 
#' 
#' This function samples (or extracts previously-sampled) trait values at tips
#' or internal nodes (a form of stochastic ancestral state reconstruction) given
#' an \code{evorates_fit} object.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}".
#' @param select A numeric vector specifying node indices in \code{fit$call$tree} for  which
#' to extract trait values. These can be negative to exclude nodes as well. As with all
#' objects of class "\code{phylo}", the first \code{n} nodes correspond to tips (where \code{n}
#' is the number of tips in the tree) and the rest correspond to internal nodes. All nodes (tip
#' and internal alike) are extracted by default.
#' @param type A string specifying whether to extract posterior samples ("\code{chains}"),
#' quantiles ("\code{quantiles}"), means ("\code{means}"), or diagnostics ("\code{diagnostics}").
#' Can be any unambiguous abbreviation of these strings as well. Defaults to extracting
#' posterior samples.
#' @param extra.select A numeric, integer, or character vector specifying the specific samples/
#' quantiles/diagnostics to extract, depending on \code{type}. Defaults to \code{NULL},
#' which generally extracts all available samples/quantiles/diagnostics. See documentation on
#' param_block operators for details (\link{grapes-chains-grapes}, \link{grapes-quantiles-grapes},
#' \link{grapes-means-grapes}, or \link{grapes-diagnostics-grapes}). 
#' @param simplify \code{TRUE} or \code{FALSE}: should the resulting \code{param_block} array be 
#' simplified? If \code{TRUE} (the default), dimensions of length 1 in the result are
#' automatically collapsed, with corresponding information stored as attributes (this is the default
#' behavior of param_block operators).
#' @param store.in.fit (new in 0.1.1) \code{TRUE} or \code{FALSE}: should the resulting \code{param_block}
#' array be added to the \code{param_block} arrays in \code{fit}? Often makes the \code{fit} object quite a
#' bit larger, but makes later computation of trait values faster and more reproducible. Defaults to
#' \code{FALSE}. Note that \code{select}, \code{type}, and \code{extra.select} are ignored if this is
#' \code{TRUE}.
#' @param force.resample (new in 0.1.1) \code{TRUE} or \code{FALSE}: should trait values be resampled
#' even if already stored in \code{fit}? Defaults to \code{FALSE}.
#' 
#' 
#' @return If \code{store.in.fit} is \code{FALSE}: an array of class "\code{param_block}" with a
#' \code{param_type} set to whatever
#' \code{type} is. The dimension of these arrays will generally go in the order of
#' iterations/quantiles/diagnostics, then parameters, then chains. Any dimensions of length 1 are
#' collapsed and stored as attributes if \code{simplify} is \code{TRUE}. 
#' If \code{store.in.fit} is \code{TRUE}: an \code{evorates_fit} object including extra parameters corresponding
#' to trait values. Parameters go by the
#' name \code{<trait.name>_i}, where \code{i} is the label or index of the corresponding tip or internal
#' node, respectively, in \code{fit$call$tree}. The trait name
#' corresponds to the column name of \code{fit$call$trait.data}, which is "X1" if not provided in the
#' original input to \code{fit.evorates()} or \code{input.evorates()}.
#' 
#' 
#' @details The tip and node indices of \code{fit$call$tree} can be viewed by running:
#' \code{plot(fit$call$tree); tiplabels(); nodelabels()}.
#' 
#' Note that this function does random sampling internally, so you will get different results
#' every time you run it unless you use a consistent seed (e.g., using \code{set.seed()}) or
#' store sampled trait values.
#' 
#' 
#' @family parameter extraction functions
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #get posterior samples of all trait values
#' x <- get.post.traits(cet_fit)
#' #let's say we want to get some trait values towards the root
#' plot(cet_fit$call$tree); nodelabels()
#' #select particular edges
#' x <- get.post.traits(cet_fit, select = c(89, 90, 94, 103, 104, 105))
#' #maybe specific samples too?
#' x <- get.post.traits(cet_fit, select = c(89, 90, 94, 103, 104, 105),
#'                      extra.select=c(1, 23, 47))
#' 
#' #this may be a more common way to get node indices; say we want to look at the genus Mesoplodon
#' edges <- get.clade.edges(cet_fit$call$tree, "Mesoplodon")
#' nodes <- unique(as.vector(cet_fit$call$tree$edge[edges,]))
#' Meso.x <- get.post.traits(cet_fit, select = nodes)
#' #or at everything EXCEPT Mesoplodon
#' notMeso.x <- get.post.traits(cet_fit, select = -nodes)
#' #you could also rely on string matching if you're interested in the tips only
#' Meso.tips.x <- get.post.traits(cet_fit, select = "Mesoplodon")
#' 
#' #could also look at quantiles, means, or diagnostics
#' med.x <- get.post.traits(cet_fit, select = nodes,
#'                      type = "quantiles",
#'                      extra.select = 0.5)
#' mean.x <- get.post.traits(cet_fit, select = nodes,
#'                           type = "means")
#' init.x <- get.post.traits(cet_fit, select = nodes,
#'                           type = "diagnostics",
#'                           extra.select = c("inits", "ess"))
#'                  
#' #here's an example of what happens when you don't simplify the result
#' med.x <- get.R(cet_fit, select = nodes,
#'                type = "quantiles",
#'                extra.select = 0.5,
#'                simplify = FALSE)
#'                  
#' #note that, depending on your goals, it may be more efficient to do things like this instead
#' x <- get.post.traits(cet_fit)
#' med.Meso.x <- x %quantiles% list(nodes, 0.5)
#' 
#' 
#' @export
get.post.traits<-function(fit,
                          select=seq_len(Ntip(fit)+Nnode(fit)),
                          type=c('chains','quantiles','means','diagnostics'),
                          extra.select=NULL,
                          simplify=TRUE,
                          store.in.fit=FALSE,force.resample=FALSE){
  #helper function which I may want to generalize to other things...
  .make.par<-function(x){
    x<-.add.par.class(x)
    attr(x,'param_type')<-'chains'
    x
  }
  #helper for fast row sums
  .fast.sum<-function(x,n){
    eval(str2lang(paste0('x[,',seq_len(n),',,drop=FALSE]',collapse='+')))
  }
  
  #create function to check validity of node indices...
  type<-.match.type(type)
  
  #check for stored trait data
  tree<-fit$call$tree
  ntips<-Ntip(fit)
  nnodes<-Nnode(fit)
  trait.name<-colnames(fit$call$trait.data)
  param.names<-paste0(trait.name,'_',c(tree$tip.label,seq_len(nnodes)+ntips))
  has.stored<-all(param.names%in%names(fit$chains))
  
  if(!has.stored|force.resample){
    #need to resample traits
    
    #default settings if looking to store within fit
    if(store.in.fit){
      type<-"chains"
      select<-seq_len(Ntip(fit)+Nnode(fit))
      extra.select<-NULL
    }
    
    #tmp.extra.select --> initial selection to calculate whatever extra.select specifies
    if(type=="chains"&!is.null(extra.select)){
      tmp.extra.select<-extra.select
      extra.select<-NULL
    }else{
      tmp.extra.select<-NULL
    }
    
    R<-get.R(fit,type='chains',extra.select=tmp.extra.select,simplify=FALSE)
    niter<-dim(R)[1]
    #add in inits
    if(type=='diagnostics'|store.in.fit){
      inds<-c(1,seq_len(niter))
      R<-.call.select(R,list(NULL,inds))
      R[1,,]<-get.R(fit,type='diagnostics',extra.select='inits',simplify=FALSE)
      niter<-niter+1
    }
    #convert NAs in R to 0 for any remaining "pesudo branches"
    R[is.na(R)]<-0
    R<-exp(R)*tree$edge.length
    
    #trait.data code, mainly taken from FIT_main--might be good to consolidate...
    trait.data<-fit$call$trait.data
    parsed.trait.data<-split(trait.data,rownames(trait.data))
    #return NA is every trait val is NA, otherwise ignore NA entries
    foo<-function(x){
      if(all(is.na(x))){
        NA
      }else{
        mean(x,na.rm=TRUE)
      }
    }
    X<-sapply(parsed.trait.data,foo)[tree$tip.label]
    n_obs<-lengths(parsed.trait.data)[tree$tip.label]
    X[is.na(X)]<-n_obs[is.na(X)]<-0
    trait.se<-as.vector(fit$call$trait.se)
    if(!is.null(fit$call$Ysig2_prior_sig)){
      Ysig2<-.call.op('chains',fit,list('^Y_sig2$',tmp.extra.select),FALSE)
      if(type=='diagnostics'|store.in.fit){
        Ysig2<-.call.select(Ysig2,list(NULL,inds))
        Ysig2[1,,]<-.call.select(fit$diagnostics,list('Y_sig2','inits'))
      }
      Ysig2<-Ysig2/n_obs
      names(Ysig2)<-tree$tip.label
      Ysig2<-.strip.par.class(Ysig2)
    }else{
      dimnms<-dimnames(R)
      dims<-dim(R)
      dimnms[[2]]<-tree$tip.label
      dims[2]<-ntips
      Ysig2<-array(0,dims,dimnms)
    }
    prior.se<-!is.na(trait.se)
    if(any(prior.se)){
      Ysig2[,which(prior.se),]<-rep(trait.se[prior.se]^2,each=niter)
    }
    
    #get topology info
    des.e<-c(list(root.edges(fit$call$tree)),des.edges(fit))
    names(des.e)<-seq_along(des.e)-1
    des.e<-rev(des.e)
    ndes<-lengths(des.e)
    non.tip<-ndes>0
    des.e<-des.e[non.tip]
    ndes<-ndes[non.tip]
    prune.seq<-seq_along(des.e)
    des.n<-split(tree$edge[unlist(des.e,use.names=FALSE),2],
                 rep.int(prune.seq,ndes))
    foc.e<-as.numeric(names(des.e))
    foc.n<-c(tree$edge[foc.e,2],ntips+1)
    
    #initializing X and P arrays...
    dimnms<-dimnames(R)
    dims<-dim(R)
    dimnms[[2]]<-param.names
    if(is.numeric(select)){
      select<-paste0('^',dimnms[[2]][select],'$')
    }
    dims[2]<-ntips+nnodes
    PP<-XX<-.make.par(array(dim=dims,dimnames=dimnms))
    XX[,seq_len(ntips),]<-rep(X,each=niter)
    PP[,seq_len(ntips),]<-1/(Ysig2+R[,tip.edges(fit),,drop=FALSE])
    
    #preorder traversal
    for(i in prune.seq){
      tmp.foc.e<-foc.e[i]
      tmp.ndes<-ndes[i]
      tmp.des.n<-des.n[[i]]
      tmp.foc.n<-foc.n[i]
      des.P<-PP[,tmp.des.n,,drop=FALSE]
      tmp.P<-sum.P<-.fast.sum(des.P,tmp.ndes)
      if(tmp.foc.e){ #adjusting for ancestral edge if this is non-root node
        inf.prec<-is.infinite(sum.P)
        tmp.R<-R[,tmp.foc.e,,drop=FALSE]
        tmp.P[!inf.prec]<-sum.P[!inf.prec]/(1+tmp.R[!inf.prec]*sum.P[!inf.prec])
        tmp.P[inf.prec]<-1/tmp.R[inf.prec]
      }
      PP[,tmp.foc.n,]<-tmp.P
      des.P<-des.P/sum.P[,rep.int(1,tmp.ndes),,drop=FALSE]
      #assumes only 1 descendant has infinite precision!
      des.P[is.nan(des.P)]<-1
      des.P[is.infinite(des.P)]<-0
      des.X<-XX[,tmp.des.n,,drop=FALSE]
      XX[,tmp.foc.n,]<-.fast.sum(des.X*des.P,tmp.ndes)
    }
    
    #postorder traversal
    R<-R[,order(tree$edge[,2]),,drop=FALSE]
    PL<-PP[,-(ntips+1),,drop=FALSE]*R
    R<-1/R
    SS<-PP
    SS[,ntips+1,]<-1/SS[,ntips+1,,drop=FALSE]
    SS[,-(ntips+1),]<-1/(SS[,-(ntips)+1,,drop=FALSE]/(1-PL)+R)
    SS[is.nan(SS)]<-0 #should take care of infinities...
    SS<-sqrt(SS)
    seed<-rnorm.par(n=dim(XX)[2],fit=XX)*SS
    XX[,ntips+1,]<-XX[,ntips+1,,drop=FALSE]+seed[,ntips+1,,drop=FALSE]
    for(i in rev(prune.seq)){
      tmp.ndes<-ndes[i]
      PL.inds<-tmp.des.n<-des.n[[i]]
      tmp.inds<-PL.inds>ntips
      PL.inds[tmp.inds]<-PL.inds[tmp.inds]-1
      tmp.foc.n<-foc.n[i]
      tmp.X<-XX[,rep.int(tmp.foc.n,tmp.ndes),,drop=FALSE]
      XX[,tmp.des.n,]<-PL[,PL.inds,,drop=FALSE]*(XX[,tmp.des.n,,drop=FALSE]-tmp.X)+tmp.X+seed[,tmp.des.n,,drop=FALSE]
    }
    
    #store or output stuff
    if(type=='diagnostics'|store.in.fit){
      inits.save<-XX[1,,,drop=FALSE]
      XX<-.make.par(XX[-1,,,drop=FALSE])
    }
    out<-list(sampler.control=1,chains=XX)
    if(store.in.fit){
      for(i in c('quantiles','means','diagnostics')){
        if(!is.null(fit[[i]])){
          extra.select<-dimnames(.make.par.3D(fit[[i]]))[[1]]
          if(i=='diagnostics'){
            tmp.inds<-extra.select=='inits'
          }
          out[[i]]<-.call.op(i,out,list(".",extra.select),FALSE)
          if(i=='diagnostics'){
            tmp.n<-sum(tmp.inds)
            out[[i]][tmp.inds,,]<-rep(inits.save,each=tmp.n)
          }
        }
      }
      out<-out[-1]
      for(i in names(out)){
        if(has.stored){
          fit[[i]][,param.names,]<-out[[i]]
        }else{
          fit[[i]]<-.combine.par(list(fit[[i]],out[[i]]))
        }
        if(simplify){
          fit[[i]]<-.simplify.par(fit[[i]])
        }
      }
      for(i in c("sampler.params","post.probs")){
        if(!is.null(fit[[i]])){
          if(simplify){
            fit[[i]]<-.simplify.par(fit[[i]])
          }else{
            fit[[i]]<-.expand.par(fit[[i]])
          }
        }
      } 
      out<-fit
    }else{
      out<-.call.op(type,out,list(select,extra.select),FALSE)
      if(type=='diagnostics'){
        inds<-dimnames(out)[[1]]=='inits'
        out[inds,,]<-rep(inits.save,each=sum(inds))
      }
      if(simplify){
        out<-.simplify.par(out)
      }
    }
  }else{
    if(store.in.fit){
      # warning("store.in.fit is TRUE, but force.resample is FALSE--trait data in provided not resampled")
      par.inds<-which(names(fit)!='call'&names(fit)!='sampler.control')
      for(i in par.inds){
        if(simplify){
          fit[[i]]<-.simplify.par(fit[[i]])
        }else{
          fit[[i]]<-.expand.par(fit[[i]])
        }
      }
      out<-fit
    }else{
      #just select from fit
      #no (re)storage necessary--just output
      #note here that you might want to keep better track of original input for better error messages
      if(is.character(select)){
        select<-unlist(lapply(select,grep,x=param.names,value=TRUE),use.names=FALSE)
      }
      if(is.numeric(select)){
        select<-param.names[select]
      }
      select<-paste0('^',select,'$')
      
      out<-.call.op(type,fit,list(select,extra.select),FALSE)
      if(simplify){
        out<-.simplify.par(out)
      }
    }
  }
  out
}

#helper function
#drops any stored trait data from fit, also expands all stored param_blocks as side effect
.drop.stored.traits<-function(fit){
  tree<-fit$call$tree
  ntips<-Ntip(fit)
  nnodes<-Nnode(fit)
  trait.name<-colnames(fit$call$trait.data)
  param.names<-paste0(trait.name,'_',c(tree$tip.label,seq_len(nnodes)+ntips))
  has.stored<-all(param.names%in%names(fit$chains))
  pos<-match(param.names,names(fit$chains))
  for(i in c("chains","quantiles","means","diagnostics")){
    if(!is.null(fit[[i]])){
      fit[[i]]<-.expand.par(fit[[i]])
      if(has.stored){
        fit[[i]]<-.add.par.class(fit[[i]][,-pos,,drop=FALSE])
        attr(fit[[i]],'param_type')<-i
      }
    }
  }
  for(i in c("sampler.params","post.probs")){
    if(!is.null(fit[[i]])){
      fit[[i]]<-.expand.par(fit[[i]])
    }
  }
  fit
}