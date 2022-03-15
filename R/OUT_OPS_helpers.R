#one last thing that might be nice is taking "shortcuts" when one desires everything along dimension 1/2
#would require some potentially complicated switches, though

####GENERAL STUFF####
#new simpler, more generalizable code
##.proc.op process the left-hand input and simplifies the output
##.call.op is the main "switchboard" that calls .int.type as appropriate
###now with check.sampler TRUE or FALSE based on whether to check if sampler params are being extracted

.proc.op<-function(type,fit,select,deparsed.select,
                   choices='chains'){
  flag<-FALSE
  par<-.is.par(fit)
  if(par){
    in.type<-.get.par.type(fit)
    if(all(in.type!=choices)){
      flag<-TRUE
    }else{
      fit<-setNames(list(fit),in.type)
    }
  }else if(!inherits(fit,'evorates_fit')){
    flag<-TRUE
  }
  if(flag){
    stop(paste0("the %",type,"% operator only accepts evorates_fit or param_block (with param_type ",
                paste0(choices,collapse=" or "),
                ") objects on left hand side")
    )
  }
  if(inherits(try(is.character(select),silent=TRUE),'try-error')){
    select<-deparsed.select
  }
  out<-.call.op(type,fit,select,check.sampler=!par)
  .simplify.par(out)
}

.call.op<-function(type,fit,select,check.sampler=TRUE){
  if(check.sampler&!is.null(fit[['sampler.params']])){
    sampler.tmp<-.add.sampler.select(type,fit,select)
    sampler.out<-sampler.tmp[[1]]
    select<-sampler.tmp[[2]]
  }else{
    sampler.out<-NULL
  }
  if(!is.null(fit[[type]])){
    fit[[type]]<-.make.par.3D(fit[[type]])
  }
  out<-do.call(paste0('.int.',type),list(fit=fit,select=select))
  out<-.add.par.class(out)
  attr(out,'param_type')<-type
  if(type!='chains'){
    if(is.null(fit[['sampler.control']])&is.null(fit[[type]])&length(select)){
      nms<-names(out)
      pre.parens<-grepl('^\\(',nms)&grepl('\\)$',nms)
      nms[pre.parens]<-paste0(type,nms[pre.parens])
      nms[!pre.parens]<-paste0(type,'(',nms[!pre.parens],')')
      names(out)<-nms
    }
  }
  if(!is.null(sampler.out)){
    out<-.combine.par(list(out,sampler.out))
  }
  out
}

####SAMPLER STUFF####

#might be nice, for consistency's sake, to figure out a way to rearrange param_block to reflect order of inputs better...
#right now, sampler and non-sampler are split, and sampler is grafted onto the end of non-sampler portion...
.add.sampler.select<-function(type,fit,select){
  avail<-paste0(c('accept_stat','stepsize','treedepth','n_leapfrog','divergent','energy','prior','lik','post'),'__')
  sampler.params<-.make.par.3D(fit[['sampler.params']])
  if(is.list(select)){
    extra.select<-select[[2]]
    select<-select[[1]]
  }else{
    extra.select<-NULL
  }
  sampler.matches<-match(select,avail)
  did.match<-!is.na(sampler.matches)
  if(any(did.match)){
    sampler.select<-sampler.matches[did.match]
    if(is.null(extra.select)&(type=='quantiles'|type=='diagnostics')&!is.null(fit[[type]])){
      extra.select<-dimnames(.make.par.3D(fit[[type]]))[[1]]
    }
    select<-list(select[!did.match],extra.select)
    if(any(!did.match)){
      diags.len<-dim(sampler.params)[1]
      target.len<-dim(.make.par.3D(fit[['chains']]))[1] #should never get here with param_block on left hand side...
      tmp<-diags.len-target.len
      if(tmp){
        if(type=='diagnostics'){
          inits.save<-sampler.params[1,sampler.select,,drop=FALSE]
        }
        sampler.params<-sampler.params[-seq_len(diags.len-target.len),sampler.select,,drop=FALSE]
      }else if(type=='diagnostics'){
        inits.save<-NULL
      }
    }else{
      sampler.params<-sampler.params[,sampler.select,,drop=FALSE]
    }
    attr(sampler.params,'param_type')<-'chains'
    if(type!='chains'){
      sampler.params<-.call.op(type,list(chains=sampler.params,sampler.control=1),list('.',extra.select),FALSE) #"cheating" to have non-null sampler.params
    }
    if(type=='diagnostics'){
      if(!is.null(inits.save)){
        inds<-dimmnames(out)[[1]]=='inits'
        sampler.params[inds,,]<-rep(inits.save,each=sum(inds))
      }
    }
  }else{
    select<-list(select,extra.select)
    sampler.params<-NULL
  }
  list(sampler.params,select)
}

####WORKHORSES####

#now handles subsetting iterations, quantiles, or diagnostics (for use in %s%, not %q%/%d%)
.subset.dim1<-function(type,x,select){
  if(is.list(select)){
    out.type<-switch(type,
                     chains='iteration',
                     quantiles='quantile',
                     diagnostics='diagnostic',
                     NULL)
    if(is.null(out.type)){
      extra.select<-NULL
    }else{
      extra.select<-select[[2]]
      select<-select[[1]]
    }
    if(!is.null(extra.select)){
      if(length(extra.select)){
        probs<-is.na(extra.select)|is.infinite(extra.select)
        if(any(probs)){
          extra.select<-extra.select[!probs]
          warning('Some ',out.type,' selections were NA/NaN/Inf and ignored')
        }
        probs<-FALSE
        if(is.character(extra.select)){
          extra.select<-match(extra.select,dimnames(x)[[1]])
          probs<-is.na(matches)
        }
        if(!is.numeric(extra.select)){
          extra.select<-NULL
          warning("The format of ",out.type,' selections was not recognized and ignored: this should be a numeric ',
                  if(out.type=='iteration') '' else 'or character ',
                  'vector')
        }else{
          probs<-(extra.select>dim(x)[1])|probs #only remaining NAs will have probs=TRUE, so this will capture them
          if(any(probs)){
            extra.select<-extra.select[!probs]
            warning('Some ',out.type,' selections were out of bounds and ignored')
          }
        }
      }
    }
  }else{
    extra.select<-NULL
  }
  list(extra.select,select)
}

#better behavior with numeric indices, but now only selects the FIRST match based on its FIRST pattern recognition!
#i.e., 1 will match up to R_1%-%R_2, but NOT R_2%-%R_1; the latter would match with 2
#if multiple columns match up to 1, etc., only the first column will be extracted
#first checks for R_ params, then Rdev_ params, then just uses numeric index selection
#s.flag can be set to TRUE to default to numeric index selection
#now s.flag causes this to exhibit "normal R" conventions of numeric indexing and no regexpr-based matching
#the various %m%, %c%, etc. operators can be called to do regexpr stuff
.subset.dim2<-function(type,x,select,s.flag=FALSE){
  select<-unlist(select,use.names=TRUE)
  out<-numeric(0)
  if(length(select)){
    probs<-is.na(select)|is.infinite(select)
    if(any(probs)){
      select<-select[!probs]
      warning('Some parameter selections were NA/NaN/Inf and ignored')
    }
    if(length(select)){
      select<-as.vector(select)
      param.names<-dimnames(x)[[2]]
      if(is.numeric(select)){
        if(any(select!=0)){ #make sure at least one element of select is not 0
          if(!s.flag){
            matches<-regexpr('^R_[1-9][0-9]*$|^R_[1-9][0-9]*[^0-9]|[%\\(]R_[1-9][0-9]*$|[%\\(]R_[1-9][0-9]*[^0-9]',param.names)
            match.offset<-2
            if(all(matches)==-1){
              matches<-regexpr('^Rdev_[1-9][0-9]*$|^Rdev_[1-9][0-9]*[^0-9]|[%\\(]Rdev_[1-9][0-9]*$|[%\\(]Rdev_[1-9][0-9]*[^0-9]',param.names)
              match.offset<-4
            }
            if(all(matches)==-1){
              matches<-regexpr('^uncent_Rdev_[1-9][0-9]*$|^uncent_Rdev_[1-9][0-9]*[^0-9]|[%\\(]uncent_Rdev_[1-9][0-9]*$|[%\\(]uncent_Rdev_[1-9][0-9]*[^0-9]',param.names)
              match.offset<-11
            }
          }else{
            matches<- -1
          }
          did.match<-matches> -1
          if(any(did.match)){
            matches.lengths<-attr(matches,'match.length')[did.match]
            matches<-matches[did.match]
            inds<-substr(param.names[did.match],matches+match.offset,matches-1+matches.lengths)
            inds<-as.numeric(gsub('^_|[^0-9]$','',inds))
            dummy.vec<-rep(NA,max(inds))
            #I believe this should work for only grabbing first columns...
            not.dups<-!duplicated(inds)
            dummy.vec[inds[not.dups]]<-which(did.match)[not.dups]
            out<-dummy.vec[select]
            probs<-is.na(out)
            out<-out[!probs]
            if(any(probs)&all(select>=0)){ #only want this message if index inclusions rather than exclusions were specified
              warning('Some parameter selections were out of bounds and ignored')
            }
          }else{
            probs<-(select>dim(x)[2])
            if(any(probs)){
              select<-select[!probs]
              warning('Some parameter selections were out of bounds and ignored')
            }
            out<-select
          }
        }
      }else if(is.character(select)){
        if(s.flag){
          out<-match(select,param.names)
          probs<-is.na(out)
          if(any(probs)){
            out<-out[!probs]
            warning('Some parameter selections were out of bounds and ignored')
          }
        }else{
          out<-lapply(select,grep,x=param.names)
          probs<-lengths(out)==0
          if(any(probs)){
            warning("Couldn't find parameters corresponding to ",
                    unique(paste(select[probs],collapse=', ')))
          }
          out<-unlist(out,use.names=FALSE)
        }
      }else{
        stop("The format of parameter selections was not recognized: this should be a numeric or character vector")
      }
    }
  }
  out
}

####SPECIFIC STUFF####

#9/1/21: could improve the behavior of these by not forcing 4D arrays to be 3D...but that can wait until multivariate update
#would be pretty simple in case of .int.chains and .int.means, I think, but more complex for other things...
#could get weird with sampler additions, but that shouldn't be a concern since all 4D arrays should be param_blocks

.int.chains<-function(fit,select){
  out<-fit[['chains']]
  tmp<-.subset.dim1('chains',out,select)
  extra.select<-tmp[[1]]
  select<-.subset.dim2('chains',out,tmp[[2]])
  if(!is.null(extra.select)){
    out[extra.select,select,,drop=FALSE]
  }else{
    out[,select,,drop=FALSE]
  }
}

.int.quantiles<-function(fit,select){
  #get default quantiles
  if(is.null(fit[['quantiles']])){
    def.quantiles<-c(0.025,0.5,0.975)
  }else{
    def.quantiles<-dimnames(fit[['quantiles']])[[1]]
    ends<-nchar(def.quantiles)-1
    def.quantiles<-as.numeric(substr(def.quantiles,0,ends))/100
  }
  
  #process extra.select
  #no need to worry about less than 0s/greater than 1s--these are handled by .int.quant now
  if(is.list(select)){
    extra.select<-select[[2]]
    select<-select[[1]]
    if(!is.null(extra.select)){
      if(length(extra.select)){
        if(is.character(extra.select)){
          extra.select<-as.numeric(gsub('%','',extra.select))/100
        }
        if(!is.numeric(extra.select)){
          extra.select<-NULL
          warning("The format of quantile selections was not recognized and ignored: this should be a numeric vector of either numbers between 0 and 1 or integers, or alternatively a character vector of numbers between 0 and 100%")
        }else{
          probs<-is.na(extra.select)|is.infinite(extra.select)
          if(any(probs)){
            extra.select<-extra.select[!probs]
            warning('Some quantile selections were NA/NaN/Inf and ignored')
          }
        }
      }
    }
  }else{
    extra.select<-NULL
  }
  if(is.null(extra.select)){
    extra.select<-def.quantiles
  }
  
  #match extra.select to available quantiles
  #extra.select basically becomes the labels of the quantiles
  #match denotes which row to pull from in quantiles param_block
  #probs denotes which matches are unavailable (NA)
  if(!is.null(fit[['quantiles']])){
    if(is.integer(extra.select)){
      matches<-extra.select
      probs<-matches>length(def.quantiles)
      if(any(probs)){
        matches<-matches[!probs]
        warning('Some quantile selections were out of bounds and ignored')
      }
      extra.select<-extra.select[matches]
      probs<-rep(FALSE,length(extra.select))
    }else{
      matches<-match(extra.select,def.quantiles)
      probs<-is.na(matches)
      if(any(probs)){
        if(is.null(fit[['chains']])){
          extra.select<-extra.select[!probs]
          matches<-matches[!probs]
          probs<-rep(FALSE,length(extra.select))
          warning('Some quantile selections were out of bounds and ignored')
        }else{
          fit[['chains']]<-.make.par.3D(fit[['chains']])
        }
      }
    }
    dimnms<-dimnames(fit[['quantiles']])
    select<-.subset.dim2('quantiles',fit[['quantiles']],select)
  }else{
    matches<-rep(NA,length(extra.select))
    probs<-rep(TRUE,length(extra.select))
    fit[['chains']]<-.make.par.3D(fit[['chains']])
    dimnms<-dimnames(fit[['chains']])
    select<-.subset.dim2('quantiles',fit[['chains']],select)
  }
  
  #initialize output array
  if(length(extra.select)){
    dimnms[[1]]<-paste0(extra.select*100,'%')
  }
  dimnms[[2]]<-dimnms[[2]][select]
  names(dimnms)[1]<-'quantiles'
  dim3<-length(dimnms[[3]])
  out<-array(NA,dim=c(length(extra.select),length(select),dim3),dimnms)
  
  #fill output array, first with pre-stored quantiles...
  if(!is.null(fit[['quantiles']])){
    if(any(!probs)){
      out[!probs,,]<-fit[['quantiles']][matches[!probs],select,,drop=FALSE]
    }
  }
  #then with stuff coerced from chains...
  if(!is.null(fit[['chains']])){
    if(any(probs)){
      out[probs,,]<-apply(fit[['chains']][,select,,drop=FALSE],
                          -1,
                          .int.quant,
                          n=dim(fit[['chains']])[1],
                          p=extra.select[probs],
                          sorted=FALSE)
    }
  }
  out
}

.int.means<-function(fit,select){
  if(is.list(select)){
    select<-select[[1]]
  }
  if(is.null(fit[['means']])){
    fit[['chains']]<-.make.par.3D(fit[['chains']])
    dimnms<-dimnames(fit[['chains']])
    niter<-dim(fit[['chains']])[1]
    select<-.subset.dim2('means',fit[['chains']],select)
    dimnms[[2]]<-dimnms[[2]][select]
    names(dimnms)[1]<-'iterations'
    dims<-lengths(dimnms)
    dims[1]<-1
    out<-array(apply(fit[['chains']][,select,,drop=FALSE],3,.colMeans,m=niter,n=dims[2],na.rm=TRUE),dims,dimnms)
  }else{
    select<-.subset.dim2('means',fit[['means']],select)
    out<-fit[['means']][,select,,drop=FALSE]
  }
  out
}

#use grep here, but still use match in .subset.dim1 --> makes %s% behavior complementary to %d% behavior
#(just like you can use %s% to go back to normal numeric index selection, you can use %s% to go back to normal R behavior for selection along dim 1 too)
.int.diagnostics<-function(fit,select){
  #get default diagnostics
  avail<-c('inits','bulk_ess','tail_ess','Rhat')
  if(is.null(fit[['diagnostics']])){
    def.diags<-seq_len(4)
  }else{
    def.diags<-match(dimnames(fit[['diagnostics']])[[1]],avail)
  }
  
  #process extra.select
  if(is.list(select)){
    extra.select<-select[[2]]
    select<-select[[1]]
    if(!is.null(extra.select)){
      if(length(extra.select)){
        probs<-is.na(extra.select)|is.infinite(extra.select)
        if(any(probs)){
          extra.select<-extra.select[!probs]
          warning('Some diagnostic selections were NA/NaN/Inf and ignored')
        }
        if(length(extra.select)){
          if(is.character(extra.select)&length(extra.select)){
            tmp<-lapply(extra.select,grep,x=avail)
            probs<-lengths(tmp)==0
            if(any(probs)){
              warning("Couldn't find diagnostics corresponding to ",
                      unique(paste(extra.select[probs],collapse=', ')))
            }
            extra.select<-unlist(tmp,use.names=FALSE)
          }else if(is.numeric(extra.select)){
            extra.select<-match(def.diags[extra.select],avail)
            probs<-is.na(extra.select)
            if(any(probs)){
              extra.select<-extra.select[!probs]
              warning('Some diagnostic selections were out of bounds and ignored')
            }
          }else{
            extra.select<-NULL
            warning("The format of diagnostic selections was not recognized and ignored: this should be a numeric or character vector")
          }
        }
      }
    }
  }else{
    extra.select<-NULL
  }
  if(is.null(extra.select)){
    extra.select<-def.diags
  }
  
  #match extra.select to available diagnostics
  #extra.select basically becomes the labels of the diagnostics
  #match denotes which row to pull from in diagnostics param_block
  #probs denotes which matches are unavailable (NA)
  if(!is.null(fit[['diagnostics']])){
    matches<-match(extra.select,def.diags)
    probs<-is.na(matches)
    if(any(probs)){
      if(is.null(fit[['chains']])){
        extra.select<-extra.select[!probs]
        matches<-matches[!probs]
        probs<-rep(FALSE,length(extra.select))
        warning('Some diagnostics selections were out of bounds and ignored')
      }else{
        fit[['chains']]<-.make.par.3D(fit[['chains']])
      }
    }
    dimnms<-dimnames(fit[['diagnostics']])
    select<-.subset.dim2('diagnostics',fit[['diagnostics']],select)
  }else{
    matches<-rep(NA,length(extra.select))
    probs<-rep(TRUE,length(extra.select))
    fit[['chains']]<-.make.par.3D(fit[['chains']])
    dimnms<-dimnames(fit[['chains']])
    select<-.subset.dim2('diagnostics',fit[['chains']],select)
  }
  
  #initialize output array
  dimnms[[1]]<-avail[extra.select]
  dimnms[[2]]<-dimnms[[2]][select]
  names(dimnms)[1]<-'diagnostics'
  dim3<-length(dimnms[[3]])
  out<-array(NA,dim=c(length(extra.select),length(select),dim3),dimnms)
  
  #fill output array, first with pre-stored diagnostics...
  if(!is.null(fit[['diagnostics']])){
    if(any(!probs)){
      out[!probs,,]<-fit[['diagnostics']][matches[!probs],select,,drop=FALSE]
    }
  }
  #then with stuff coerced from chains...
  if(!is.null(fit[['chains']])){
    if(any(probs)){
      fit[['chains']]<-fit[['chains']][,select,,drop=FALSE]
      #most convenient, but will give wrong answer most of the time
      inits.foo<-function(x){
        x[1]
      }
      funs<-list(inits.foo,rstan::ess_bulk,rstan::ess_tail,rstan::Rhat)
      for(i in seq_len(4)){
        inds<-(extra.select==i)&probs
        if(any(inds)){
          out[inds,,]<-rep(apply(fit[['chains']],-1,funs[[i]]),each=sum(inds))
        }
      }
    }
  }
  out
}

####SELECT####

.proc.select<-function(x,select,deparsed.select){
  if(!.is.par(x)){
    stop("the %select% operator only accepts param_block objects on left hand side")
  }
  if(inherits(try(is.character(select),silent=TRUE),'try-error')){
    select<-deparsed.select
  }
  out<-.call.select(x,select)
  .simplify.par(out)
}

.call.select<-function(x,select){
  x<-.make.par.3D(x)
  type<-.get.par.type(x)
  tmp<-.subset.dim1(type,x,select)
  extra.select<-tmp[[1]]
  null.extra.select<-is.null(extra.select)
  select<-tmp[[2]]
  null.select<-is.null(select)
  if(!null.select){
    select<-.subset.dim2(type,x,select,s.flag=TRUE)
  }
  if(!null.select){
    if(!null.extra.select){
      x<-x[extra.select,select,,drop=FALSE]
    }else{
      x<-x[,select,,drop=FALSE]
    }
  }else if(!null.extra.select){
    x<-x[extra.select,,,drop=FALSE]
  }
  x<-.add.par.class(x)
  attr(x,'param_type')<-type
  x
}