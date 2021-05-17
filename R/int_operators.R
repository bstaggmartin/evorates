.int.chains<-function(fit,select){
  fit[['chains']]<-.expand.element(fit[['chains']])
  tmp<-.select.iterations(fit[['chains']],select)
  fit[['chains']]<-tmp[[1]]
  select<-tmp[[2]]
  out<-.int.op(fit[['chains']],select)
  out
}

.int.quantiles<-function(fit,select){
  if(is.null(fit[['quantiles']])){
    def.report.quantiles<-c(0.025,0.5,0.975)
  }else{
    fit[['quantiles']]<-.expand.element(fit[['quantiles']])
    def.report.quantiles<-as.numeric(substr(dimnames(fit[['quantiles']])[[1]],0,
                                            nchar(dimnames(fit[['quantiles']])[[1]])-1))/100
  }
  if(is.list(select)){
    if(length(select)>1){
      if(is.null(select[[2]])){
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
  tmp<-.int.chains(fit,select[[1]])
  out<-array(NA,dim=c(length(select[[2]]),dim(tmp)[-1]))
  dimnames(out)<-c('quantiles'=list(paste(select[[2]]*100,'%',sep='')),
                   dimnames(tmp)[-1])
  if(!is.null(fit[['quantiles']])){
    matches<-match(select[[2]],def.report.quantiles)
    out[(1:dim(out)[1])[!is.na(matches)],,]<-
      fit[['quantiles']][matches[!is.na(matches)],dimnames(tmp)[[2]],]
    select[[2]]<-select[[2]][is.na(matches)]
  }
  if(length(select[[2]])>0){
    out[is.na(out[,1,1]),,]<-apply(tmp,c(2,3),quantile,probs=select[[2]],na.rm=TRUE)
  }
  out
}

.int.means<-function(fit,select){
  if(is.null(fit[['means']])){
    tmp<-.int.chains(fit,select)
    out<-array(apply(tmp,c(2,3),mean,na.rm=TRUE),c(1,dim(tmp)[-1]),dimnames(tmp))
  }else{
    fit[['means']]<-.expand.element(fit[['means']])
    out<-.int.op(fit[['means']],select)
  }
  out
}

.int.MAPs<-function(fit,select){
  if(is.null(fit[['MAPs']])){
    tmp<-.int.chains(fit,select)
    chains.len<-dim(tmp)[1]
    nchain<-dim(tmp)[3]
    post<-.int.sampler(fit,'post')
    diags.len<-dim(post)[1]
    if(diags.len-chains.len>0){
      post<-post[-(1:(diags.len-chains.len)),,,drop=FALSE]
    }
    MAP.inds<-sapply(1:nchain, function(ii) which.max(post[,,ii]))
    out<-array(sapply(1:nchain,function(ii) tmp[MAP.inds[ii],,ii]),
               c(1,dim(tmp)[-1]),dimnames(tmp))
  }else{
    fit[['MAPs']]<-.expand.element(fit[['MAPs']])
    out<-.int.op(fit[['MAPs']],select)
  }
  out
}

#allow it to select specific iterations like the new %chains% operator
.int.sampler<-function(fit,select){
  fit[['sampler.params']]<-.expand.element(fit[['sampler.params']])
  tmp<-.select.iterations(fit[['sampler.params']],select)
  fit[['sampler.params']]<-tmp[[1]]
  select<-tmp[[2]]
  if(is.numeric(select)){
    select<-c(paste(c('accept_stat','stepsize','treedepth','n_leapfrog','divergent','energy'),
                  '__',sep=''),c('prior','lik','post'))[select]
  }
  out<-.int.op(fit[['sampler.params']],select)
  out
}

#need to update error handling
.int.diagnostics<-function(fit,select){
  if(is.list(select)){
    if(length(select)<2|is.null(select[[2]])){
      select[[2]]<-1:4
    }
  }else{
    select<-list(select,1:4)
  }
  fit[['param.diags']]<-.expand.element(fit[['param.diags']])
  out<-.int.op(fit[['param.diags']],select[[1]])
  if(is.numeric(select[[2]])){
    select[[2]]<-dimnames(out)[[1]][select[[2]]]
    if(length(select[[2]])==0){
      select[[2]]<-NA
    }
  }
  tmp<-lapply(select[[2]],function(ii) grep(ii,dimnames(out)[[1]]))
  if(all(lengths(tmp)==0)){
    stop("couldn't find any corresponding parameter diagnostics")
  }
  if(any(lengths(tmp)==0)){
    warning("couldn't find parameter diagnostics corresponding to: ",
            paste(select[[2]][which(lengths(tmp)==0)],collapse=', '))
  }
  inds<-sort(unique(unlist(tmp)))
  out<-out[inds,,,drop=FALSE]
  out
}

#internal operator -- select parameters out of fit$element based on edge numbers or regular expressions
.int.op<-function(element,select){
  select<-as.vector(select)
  param.names<-dimnames(element)[[2]]
  if(is.numeric(select)){
    n.R<-suppressWarnings(
      max(as.numeric(substr(param.names,regexpr('R_',param.names)+2,nchar(param.names))),na.rm=T)
    )
    if(all(select>n.R|select<0)){
      stop("all numbers out of range for edge rate (R) parameters: are edge indices based on right tree?")
    }
    if(any(select>n.R|select<0)){
      warning(paste(select[select>n.R|select<0],collapse=', '),
              " out of range for edge rate (R) parameters: are edge indices based on right tree?")
      select<-select[select<=n.R&select>=0]
    }
    select<-paste('R_',sort(unique(select)),'$',sep='')
  }else if(!is.character(select)){
    stop("the %",element,"% operator only accepts numeric or character vectors on right hand side")
  }
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
  forbidden.inds<-grep('_dev$',param.names)
  for(i in grep('dev$',select,invert=T)){
    tmp[[i]]<-tmp[[i]][!(tmp[[i]]%in%forbidden.inds)]
  }
  if(all(lengths(tmp)==0)){
    stop("couldn't find any corresponding parameters")
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
  inds<-sort(unique(unlist(tmp)))
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
