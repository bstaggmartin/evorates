#9/1/21: could improve the behavior of these by not forcing 4D arrays to be 3D...but that can wait until multivariate update
#would be pretty simple in case of .int.chains and .int.means, I think, but more complex for other things...
#could get weird with .add.sampler too, but that shouldn't be a concern since all 4D arrays should be loose elements

#maybe make everything skip .int.op when supplied '.' for a decent speed gain?

#also make it skip quantiles() when extracting from quantiles object (can this be done?)

#stripping double parentheses?

.int.chains<-function(fit,select){
  fit[['chains']]<-.coerce.to.3D(fit[['chains']])
  #exact matches to sampling parameters...
  if(!is.null(fit$sampler.params)){
    sampler.tmp<-.add.sampler.select(fit,select,'chains')
    sampler.out<-sampler.tmp[[1]]
    select<-sampler.tmp[[2]]
  }else{
    sampler.out<-NULL
  }
  tmp<-.select.iterations(fit[['chains']],select)
  fit[['chains']]<-tmp[[1]]
  select<-tmp[[2]]
  out<-.int.op(fit[['chains']],select)
  if(!is.null(sampler.out)){
    out<-.combine.elements(list(out,sampler.out))
  }
  attr(out,'element')<-'chains'
  out<-.add.ele(out)
  out
}

.int.quantiles<-function(fit,select){
  if(is.null(fit[['quantiles']])){
    def.report.quantiles<-c(0.025,0.5,0.975)
  }else{
    fit[['quantiles']]<-.coerce.to.3D(fit[['quantiles']])
    def.report.quantiles<-as.numeric(substr(dimnames(fit[['quantiles']])[[1]],0,
                                            nchar(dimnames(fit[['quantiles']])[[1]])-1))/100
  }
  if(is.list(select)){
    if(length(select)>1){
      if(is.character(select[[2]])){
        select[[2]]<-as.numeric(gsub('%','',select[[2]]))
        tmp<-select[[2]]>1
        tmp[is.na(tmp)]<-FALSE
        if(any(tmp)){
          select[[2]]<-select[[2]]/100
        }
      }else if(!is.numeric(select[[2]])){
        select[[2]]<-numeric(0)
      }
      probs<-is.na(select[[2]])
      if(any(probs)){
        select[[2]]<-select[[2]][!probs]
        if(all(probs)){
          warning('All specified quantiles are NA: resorted to default quantiles')
        }else{
          warning('Some specified quantiles are NA: these were ignored')
        }
      }
      if(length(select[[2]])==0){
        select[[2]]<-def.report.quantiles
      }else{
        lessthans<-select[[2]]<0
        greaterthans<-select[[2]]>1
        if(any(lessthans)){
          select[[2]][lessthans]<-0
          warning('Some specified quantiles are less than 0 (they should be between 0 and 1): these were set to 0.')
        }
        if(any(greaterthans)){
          select[[2]][greaterthans]<-1
          warning('Some specified quantiles are greater than 1 (they should be between 0 and 1): these were set to 1.')
        }
      }
    }else{
      select[[2]]<-def.report.quantiles
    }
  }else{
    select<-list(select,def.report.quantiles)
  }
  #exact matches to sampling parameters...
  if(!is.null(fit$sampler.params)){
    sampler.tmp<-.add.sampler.select(fit,select,'quantiles')
    sampler.out<-sampler.tmp[[1]]
    select<-sampler.tmp[[2]]
  }else{
    sampler.out<-NULL
  }
  if(!is.null(fit$quantiles)){
    matches<-match(select[[2]],def.report.quantiles)
    if(is.null(fit$chains)){
      if(all(is.na(matches))){
        warning('All specified quantiles not found in provided quantiles element: resorted to default quantiles')
        select[[2]]<-def.report.quantiles
        matches<-1:length(def.report.quantiles)
      }else if(any(is.na(matches))){
        warning('Some specified quantiles not found in provided quantiles element: these were ignored')
        select[[2]]<-select[[2]][!is.na(matches)]
        matches<-matches[!is.na(matches)]
      }
      tmp<-.int.op(fit[['quantiles']],select[[1]])
    }
  }
  if(!is.null(fit$chains)){
    tmp<-.int.chains(fit,select[[1]])
  }
  out<-array(NA,dim=c(length(select[[2]]),dim(tmp)[-1]))
  dimnames(out)<-c('quantiles'=list(paste(select[[2]]*100,'%',sep='')),
                   dimnames(tmp)[-1])
  if(!is.null(fit[['quantiles']])){
    out[(1:dim(out)[1])[!is.na(matches)],,]<-
      fit[['quantiles']][matches[!is.na(matches)],dimnames(tmp)[[2]],]
    tmp.inds<-is.na(matches)
  }else{
    tmp.inds<-rep(TRUE,length(select[[2]]))
  }
  if(length(select[[2]])>0){
    out[tmp.inds,,]<-apply(tmp,c(2,3),quantile,probs=select[[2]][tmp.inds],na.rm=TRUE)
  }
  out<-.add.ele(out)
  if(is.null(fit$sampler.control)&length(select)>0){
    names(out)<-paste0('quantiles(',names(out),')')
  }
  if(!is.null(sampler.out)){
    names(sampler.out)<-gsub('\\)$','',gsub('^quantiles\\(','',names(sampler.out)))
    out<-.combine.elements(list(out,sampler.out))
  }
  attr(out,'element')<-'quantiles'
  out
}

.int.means<-function(fit,select){
  if(is.list(select)){
    select<-select[[1]]
  }
  if(!is.null(fit$sampler.params)){
    sampler.tmp<-.add.sampler.select(fit,select,'means')
    sampler.out<-sampler.tmp[[1]]
    select<-sampler.tmp[[2]]
  }else{
    sampler.out<-NULL
  }
  if(is.null(fit[['means']])){
    tmp<-.int.chains(fit,select)
    out<-array(apply(tmp,c(2,3),mean,na.rm=TRUE),c(1,dim(tmp)[-1]),dimnames(tmp))
  }else{
    fit[['means']]<-.coerce.to.3D(fit[['means']])
    out<-.int.op(fit[['means']],select)
  }
  out<-.add.ele(out)
  if(is.null(fit$sampler.control)&length(select)>0){
    names(out)<-paste0('means(',names(out),')')
  }
  if(!is.null(sampler.out)){
    names(sampler.out)<-gsub('\\)$','',gsub('^means\\(','',names(sampler.out)))
    out<-.combine.elements(list(out,sampler.out))
  }
  attr(out,'element')<-'means'
  out
}

#MAPs and sampler no longer existent

#need to update error handling to better reflect input
.int.diagnostics<-function(fit,select){
  avail<-c('inits','bulk_ess','tail_ess','Rhat')
  if(is.null(fit[['param.diags']])){
    fit$chains<-.coerce.to.3D(fit$chains)
    def.diags<-avail
  }else{
    fit[['param.diags']]<-.coerce.to.3D(fit[['param.diags']])
    def.diags<-dimnames(fit$param.diags)[[1]]
  }
  if(is.list(select)){
    if(length(select)>1){
      if(is.numeric(select[[2]])){
        select[[2]]<-def.diags[select[[2]]]
      }else if(!is.character(select[[2]])){
        select[[2]]<-character(0)
      }
      probs<-is.na(select[[2]])
      if(any(probs)){
        select[[2]]<-select[[2]][!probs]
        if(all(probs)){
          warning('All specified diagnostics are NA: resorted to default diagnostics')
        }else{
          warning('Some specified diagnostics are NA: these were ignored')
        }
      }
      if(length(select[[2]])==0){
        select[[2]]<-def.diags
      }
    }else{
      select[[2]]<-def.diags
    }
  }else{
    select<-list(select,def.diags)
  }
  #exact matches to sampling parameters...
  if(!is.null(fit$sampler.params)){
    sampler.tmp<-.add.sampler.select(fit,select,'diagnostics')
    sampler.out<-sampler.tmp[[1]]
    select<-sampler.tmp[[2]]
  }else{
    sampler.out<-NULL
  }
  matches<-lapply(select[[2]],function(ii) grep(ii,avail))
  probs<-lengths(matches)==0
  matches<-unlist(matches)
  if(all(probs)){
    warning('All specified diagnostics did not match to available diagnostics: resorted to default diagnostics')
    matches<-1:length(def.diags)
    select[[2]]<-def.diags[matches]
  }else if(any(probs)){
    warning('Some specified diagnostics did not match to available diagnostics: these were ignored')
    select[[2]]<-avail[matches]
  }
  if(is.null(fit$chains)){
    matches<-lapply(select[[2]],function(ii) grep(ii,def.diags))
    probs<-lengths(matches)==0
    matches<-unlist(matches)
    if(all(probs)){
      warning('All specified diagnostics not found in provided diagnostics element: resorted to default diagnostics')
      matches<-1:length(def.diags)
      select[[2]]<-def.diags
    }else if(any(probs)){
      warning('Some specified diagnostics not found in provided diagnostics element: these were ignored')
      select[[2]]<-def.diags[matches]
    }
  }
  #no need to ever combine things unless I allow people to stick half-complete diagnostics back into fit objects...
  if(!is.null(fit$param.diags)){
    out<-.int.op(fit$param.diags,select[[1]])
    out<-out[matches,,,drop=FALSE]
  }else{
    tmp<-.int.op(fit$chains,select[[1]])
    out<-array(NA,dim=c(length(matches),dim(tmp)[-1]))
    dimnames(out)<-c('diagnostics'=select[2],
                     dimnames(tmp)[-1])
    #can never know if chain starts at beginning...
    inits.foo<-function(x){
      NA
    }
    funs<-setNames(list(inits.foo,rstan::ess_bulk,rstan::ess_tail,rstan::Rhat),
                   avail)
    for(i in 1:length(matches)){
      out[i,,]<-apply(tmp,c(2,3),funs[[select[[2]][i]]])
    }
  }
  out<-.add.ele(out)
  if(is.null(fit$sampler.control)&length(select)>0){
    names(out)<-paste0('diagnostics(',names(out),')')
  }
  if(!is.null(sampler.out)){
    names(sampler.out)<-gsub('\\)$','',gsub('^diagnostics\\(','',names(sampler.out)))
    out<-.combine.elements(list(out,sampler.out))
  }
  attr(out,'element')<-'diagnostics'
  out
}

#internal operator -- select parameters out of fit$element based on edge numbers or regular expressions
.int.op<-function(element,select){
  if(length(select)!=0){
    if(is.list(select)){
      stop('Only non-list parameter selections allowed')
    }
    problems<-is.na(select)|is.infinite(select)
    if(all(problems)){
      stop("NA or Inf selections not allowed")
    }else if(any(problems)){
      warning("Some entries of select were NA/Inf, these were removed")
      select<-select[!problems]
    }
  }else{
    return(element[,numeric(0),,drop=FALSE])
  }
  select<-as.vector(select)
  param.names<-dimnames(element)[[2]]
  if(is.numeric(select)){
    tmp<-paste0('^R_',select,'$|^R_',select,'[^0-9]',
                '|[%\\(]R_',select,'$|[%\\(]R',select,'[^0-9]')
    tmp<-lapply(tmp,function(ii) grep(ii,param.names,perl=T))
    connect.inds<-NULL
  }else if(is.character(select)){
    select<-gsub('\\\\,',paste(rep('~',17),collapse=''),select)
    select<-gsub('\\\\_',paste(rep('@',17),collapse=''),select)
    if(any(grepl(',',select))){
      tmp<-strsplit(select,split=',')
      if(any(lengths(tmp)==3)){
        warning('text before and after comma not swapped for ',
                paste(select[which(lengths(tmp)==3)],collapse=', '),
                " since more than 1 comma was found: for any commas that are part of trait/tip names, please use '\\,' instead")
      }
      if(any(lengths(tmp)==2)){
        connect.inds<-which(lengths(tmp)==2)
        connect.inds<-cbind(connect.inds,length(select)+connect.inds)
        new.select<-tmp[connect.inds[,1]]
        new.select.l<-sapply(new.select,'[',1)
        new.select.r<-sapply(new.select,'[',2)
        new.select.r<-strsplit(new.select.r,'_')
        new.select<-sapply(1:length(new.select),
                           function(ii) paste(paste(new.select.r[[ii]][1],new.select.l[[ii]][1],sep=','),
                                              paste(new.select.r[[ii]][-1],collapse='.'),
                                              sep=if(length(new.select.r[[ii]])>1) '.' else ''))
        new.select<-gsub('\\?([\\=\\!])','~~~~~~~~~~~~~~~~~~~\\1',new.select)
        new.select<-gsub('\\?<([\\=\\!])','@@@@@@@@@@@@@@@@@@@\\1',new.select)
        new.select<-gsub('~{19}([\\=\\!])','?<\\1',new.select)
        new.select<-gsub('@{19}([\\=\\!])','?\\1',new.select)
        select<-gsub('_','.',select)
        select<-c(select,new.select)
      }else{
        connect.inds<-NULL
      }
    }else{
      connect.inds<-NULL
    }
    select<-gsub('~{17}',',',select)
    select<-gsub('@{17}','_',select)
    tmp<-lapply(select,function(ii) grep(ii,param.names,perl=T))
  }else{
    stop("the %",element,"% operator only accepts numeric or character vectors on right hand side")
  }
  forbidden.inds<-grep('^R_[1-9][0-9]*_dev$',param.names)
  if(length(forbidden.inds)>0){
    for(i in grep('dev\\$*$',select,invert=TRUE)){
      tmp[[i]]<-tmp[[i]][!(tmp[[i]]%in%forbidden.inds)]
    }
  }
  if(any(lengths(tmp)==0)){
    problem.select<-which(lengths(tmp)==0)
    if(!is.null(connect.inds)){
      for(i in 1:nrow(connect.inds)){
        if(any(lengths(tmp[connect.inds[i,]])!=0)){
          problem.select<-problem.select[-which(problem.select%in%connect.inds[i,])]
        }
      }
    }
    if(length(problem.select)>0){
      warning("couldn't find parameters corresponding to: ",
              unique(paste(select[problem.select],collapse=', ')))
    }
  }
  if(!is.numeric(select)){
    inds<-sort(unique(unlist(tmp)))
  }else{
    inds<-unlist(tmp)
  }
  element[,inds,,drop=FALSE]
}
#flip flop any question mark dohickeys...

.select.iterations<-function(element,select){
  if(is.list(select)){
    if(length(select)>1){
      if(!is.null(select[[2]])){
        all.neg<-all(select[[2]]<=0)
        all.pos<-all(select[[2]]>=0)
        if(all.neg|all.pos){
          exceeds<-select[[2]]>dim(element)[1]
          if(all(exceeds)){
            warning('All specified iterations are out of bounds (i.e., above the number of available iterations in chains): defaulted to including all iterations')
          }else{
            if(any(exceeds)){
              select[[2]]<-select[[2]][!exceeds]
              warning('Some specified iterations are out of bounds (i.e., above the number of iterations in chains): these iterations were ignored')
            }
            element<-element[select[[2]],,,drop=F]
          }
        }else{
          warning('A mix of negative and positive iterations are specified: please specify only negative iterations (to exclude iterations) or positive iterations (to include iterations). Defaulted to including all iterations.')
        }
      }
    }
    select<-select[[1]]
  }
  list(element,select)
}

.add.sampler.select<-function(fit,select,type){
  sampler.params<-.coerce.to.3D(fit$sampler.params)
  if(is.list(select)){
    select.extra<-select[[2]]
    select<-select[[1]]
  }else{
    select.extra<-NULL
  }
  sampler.matches<-grep(.get.sampler.names(),select)
  if(length(sampler.matches)>0){
    sampler.select<-paste0('^',select[sampler.matches],'$')
    select<-select[-sampler.matches]
    if(!is.null(select.extra)){
      sampler.select<-list(sampler.select,select.extra)
    }
    if(type=='diagnostics'){
      inits.save<-sampler.params[1,,,drop=FALSE]
    }
    if(length(select)>0){
      diags.len<-dim(sampler.params)[1]
      target.len<-dim(.coerce.to.3D(fit$chains))[1] #should never get here with element on left hand side...
      sampler.params<-sampler.params[-(1:(diags.len-target.len)),,,drop=FALSE]
    }
    attr(sampler.params,'element')<-'chains'
    out<-do.call(paste0('.int.',type),
                 list(fit=list(chains=sampler.params),select=sampler.select))
    if(type=='diagnostics'){
      inits.inds<-which(dimnames(out)[[1]]=='inits')
      if(length(inits.inds)>0){
        tmp<-if(is.list(sampler.select)) sampler.select[[1]] else sampler.select
        tmp<-.int.op(inits.save,tmp)
        for(i in inits.inds){
          out[inits.inds,,]<-tmp
        }
      }
    }
  }else{
    out<-NULL
  }
  if(!is.null(select.extra)){
    select<-list(select,select.extra)
  }
  list(out,select)
}
