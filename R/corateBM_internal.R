####OPERATORS####


.int.chains<-function(fit,select){
  fit[['chains']]<-.expand.element(fit[['chains']])
  out<-.int.op(fit,'chains',select)
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
    if(length(select)<2|is.null(select[[2]])){
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
    out[is.na(out[,1,1]),,]<-apply(tmp,c(2,3),quantile,probs=select[[2]])
  }
  out
}

.int.means<-function(fit,select){
  if(is.null(fit[['means']])){
    tmp<-.int.chains(fit,select)
    out<-array(apply(tmp,c(2,3),mean),c(1,dim(tmp)[-1]),dimnames(tmp))
  }else{
    fit[['means']]<-.expand.element(fit[['means']])
    out<-.int.op(fit,'means',select)
  }
  out
}

.int.MAPs<-function(fit,select){
  if(is.null(fit[['MAPs']])){
    tmp<-.int.chains(fit,select)
    chains.len<-dim(tmp)[1]
    nchain<-dim(tmp)[3]
    lp<-.int.sampler(fit,'lp')
    diags.len<-dim(lp)[1]
    if(diags.len-chains.len>0){
      lp<-.index.element(lp,1:(diags.len-chains.len),1,T)
    }
    MAP.inds<-sapply(1:nchain, function(ii) which.max(lp[,,ii]))
    out<-array(sapply(1:nchain,function(ii) tmp[MAP.inds[ii],,ii]),
               c(1,dim(tmp)[-1]),dimnames(tmp))
  }else{
    fit[['MAPs']]<-.expand.element(fit[['MAPs']])
    out<-.int.op(fit,'MAPs',select)
  }
  out
}

.int.sampler<-function(fit,select){
  fit[['sampler.params']]<-.expand.element(fit[['sampler.params']])
  if(is.numeric(select)){
    select<-paste(c('accept_stat','stepsize','treedepth','n_leapfrog','diverget','energy','lp'),
                  '__',sep='')[select]
  }
  out<-.int.op(fit,'sampler.params',select)
  out
}

.int.diagnostics<-function(fit,select){
  if(is.list(select)){
    if(length(select)<2){
      select[[2]]<-1:4
    }
  }else{
    select<-list(select,1:4)
  }
  fit[['param.diags']]<-.expand.element(fit[['param.diags']])
  out<-.int.op(fit,'param.diags',select[[1]])
  if(is.numeric(select[[2]])){
    select[[2]]<-dimnames(out)[[1]][select[[2]]]
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
  out<-.index.element(out,inds,1)
  out
}

#internal operator -- select parameters out of fit$element based on edge numbers or regular expressions
.int.op<-function(fit,element,select){
  select<-as.vector(select)
  param.names<-dimnames(fit[[element]])[[2]]
  if(is.numeric(select)){
    if(element=='sampler'){
      select<-param.names[select]
    }else{
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
    }
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
                                            paste(new.select.r[[ii]][-1],collapse='_'),
                                            sep=if(length(new.select.r[[ii]])>1) '_' else ''))
      select<-c(select,new.select)
    }else{
      connect.inds<-NULL
    }
  }else{
    connect.inds<-NULL
  }
  select<-gsub('~{17}',',',select)
  select<-gsub('@{17}','_',select)
  tmp<-lapply(select,function(ii) grep(ii,param.names))
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
  .index.element(fit[[element]],inds,2)
}


####ARRAY MANAGEMENT####


#newest version of check.n.proc and coerce.to.array--expands any given element to a 3D array (if chains,
#quantiles, sampler, or parameter diagnostics) or a 2D matrix (if means or MAPs)
#finally got it to a point where it should work with just about anything...
#I made it just always return 3D arrays, and label means and MAPs projects with 'iterations' along
#1st dimension; just makes things easier for other functions to deal with...
.expand.element<-function(arr,simplify=F){
  if(length(dim(arr))==0){
    new.arr<-as.matrix(arr)
    new.dimnames<-dimnames(new.arr)
    if(is.null(dimnames(new.arr)[[1]])){
      names(new.dimnames)[1]<-'iterations'
    }else if(sum(grepl('%',dimnames(new.arr)[[1]]))!=0){
      names(new.dimnames)[1]<-'quantiles'
    }else if(sum(grepl('^inits$|^bulk_ess$|^tail_ess$|^Rhat$',dimnames(new.arr)[[1]]))!=0){
      names(new.dimnames)[1]<-'diagnostics'
    }else if(sum(grepl('_',dimnames(new.arr)[[1]]))!=0){
      names(new.dimnames)[1]<-'parameters'
    }else if(sum(grepl('chain',dimnames(new.arr)[[1]]))!=0){
      names(new.dimnames)[1]<-'chains'
    }
    dimnames(new.arr)<-new.dimnames
    for(i in c('quantiles','diagnostics','parameters','chains')){
      if(!is.null(attr(arr,i))){
        attr(new.arr,i)<-attr(arr,i)
      }
    }
    arr<-new.arr
  }
  if(length(dim(arr))==2){
    if('quantiles'%in%names(dimnames(arr))|'quantiles'%in%names(attributes(arr))){
      first.dim<-'quantiles'
    }else if('diagnostics'%in%names(dimnames(arr))|'diagnostics'%in%names(attributes(arr))){
      first.dim<-'diagnostics'
    }else{
      first.dim<-'iterations'
    }
    dim.map<-lapply(c(first.dim,'parameters','chains'),
                    function(ii) which(names(dimnames(arr))==ii))
    out.dim<-rep(NA,3)
    out.dimnames<-setNames(vector('list',3),c(first.dim,'parameters','chains'))
    for(i in 1:3){
      if(length(dim.map[[i]])==0){
        out.dim[i]<-1
        out.dimnames[i]<-list(attr(arr,names(out.dimnames)[i]))
      }else{
        out.dim[i]<-dim(arr)[dim.map[[i]]]
        out.dimnames[i]<-list(dimnames(arr)[[dim.map[[i]]]])
      }
    }
    dim.map[lengths(dim.map)==0]<-NA
    perm<-unique(unlist(dim.map))
    if(length(perm)==2){
      perm[is.na(perm)]<-if(perm[!is.na(perm)]==1) 2 else 1
    }else{
      perm<-perm[!is.na(perm)]
    }
    arr<-array(aperm(arr,perm),out.dim,out.dimnames)
  }
  if(simplify){
    if(is.null(dimnames(arr)[[1]])&dim(arr)[1]==1){
      arr<-array(arr,dim(arr)[-1],dimnames(arr)[-1])
    }
  }
  arr
}

#opposite of expand.element--collapses dimensions of length 1 and adds them to parameters. In cases of
#length 1 output, prioritizes parameter names over everything else
##NEED TO UPDATE TO WORK WITH 2D ELEMENTS (like means and MAPs)--done 7/27
.simplify.element<-function(arr){
  old.dimnames<-dimnames(arr)
  old.dims<-dim(arr)
  out<-do.call('[',c(list(arr),lapply(old.dims,function(ii) 1:ii)))
  for(i in 1:length(old.dims)){
    if(old.dims[i]==1){
      if(!is.null(old.dimnames[[i]])){
        attr(out,names(old.dimnames)[i])<-old.dimnames[[i]]
      }
    }
  }
  if(length(out)==1){
    names(out)<-attr(out,names(old.dimnames)[2])
    attr(out,names(old.dimnames)[2])<-NULL
  }
  out
}

#newest, generalized version of reduce.array--selects indices from arbitrary arrays without
#simplification or destroying names
.index.element<-function(arr,inds,dims,invert=F){
  inds.list<-vector('list',length(dim(arr)))
  if(is.numeric(inds)){
    inds<-list(inds)
  }
  inds.list[dims]<-inds
  if(invert){
    inds.list<-lapply(inds.list,function(ii) -1*ii)
  }else{
    inds.list<-lapply(1:length(inds.list),
                      function(ii) if(is.null(inds.list[[ii]])) 1:dim(arr)[ii] else inds.list[[ii]])
  }
  new.dims<-dim(arr)
  if(invert){
    new.dims<-new.dims-lengths(inds.list)
    inds.list[lengths(inds.list)==0]<- -dim(arr)[lengths(inds.list)==0]-1
  }else{
    new.dims<-lengths(inds.list)
  }
  inds.list<-lapply(inds.list,sort)
  new.dimnames<-dimnames(arr)
  new.dimnames<-lapply(1:length(inds.list),
                       function(ii) new.dimnames[[ii]][inds.list[[ii]]])
  names(new.dimnames)<-names(dimnames(arr))
  array(do.call('[',c(list(arr),inds.list)),new.dims,new.dimnames)
}


####OTHER####


.coerce.to.cov.mat<-function(in.mat,vars){
  dims<-dim(in.mat)
  if(length(dims)>0){
    new.diag<-rep(diag(in.mat),length.out=vars)
    mat<-do.call(rbind,c(list(in.mat),rep(0,vars-nrow(in.mat))))
    if(ncol(mat)>vars){
      mat<-mat[,1:vars]
    }else{
      mat<-do.call(cbind,c(list(mat),rep(0,vars-ncol(mat))))
    }
    diag(mat)<-new.diag
    if(!isSymmetric(mat)){
      warning(paste(deparse(substitute(in.mat)),'is not symmetric: reflected lower triangle into upper triangle'))
      mat[upper.tri(mat)]<-t(mat)[upper.tri(mat)]
    }
    if(any(eigen(mat)$values<=0)){
      stop(paste('failed to create properly-formed covariance matrix from ',deparse(substitute(in.mat)),': make sure variance-covariance structure makes sense',sep=''))
    }
    mat
  }else{
    mat<-matrix(0,vars,vars)
    diag(mat)<-rep(in.mat,length.out=vars)
    mat
  }
}

.quick.recon<-function(X,tree){
  if(!is.null(names(X))){
    X<-X[tree$tip.label]
  }
  XX<-rep(NA,nrow(tree$edge)+1)
  names(XX)<-c(tree$tip.label,length(tree$tip.label)+1:tree$Nnode)
  PP<-rep(NA,nrow(tree$edge)+1)
  for(e in nrow(tree$edge):0){
    if(e==0){
      n<-length(tree$tip.label)+1
      d<-tree$edge[which(tree$edge[,1]==n),2]
      sum.P<-sum(PP[d])
      XX[n]<-sum(XX[d]*PP[d]/sum.P)
      PP[n]<-sum(sum.P)
      break
    }
    n<-tree$edge[e,2]
    if(length(which(tree$edge[,1]==n))==0){
      XX[n]<-X[n]
      PP[n]<-1/tree$edge.length[e]
    }
    if(length(which(tree$edge[,1]==n))>0){
      d<-tree$edge[which(tree$edge[,1]==n),2]
      sum.P<-sum(PP[d])
      XX[n]<-sum(XX[d]*PP[d]/sum.P)
      PP[n]<-sum.P/(1+tree$edge.length[e]*sum.P)
    }
  }
  for(e in 1:nrow(tree$edge)){
    n<-tree$edge[e,2]
    if(length(which(tree$edge[,1]==n))!=0){
      a<-tree$edge[e,1]
      XX[n]<-XX[n]*PP[n]*tree$edge.length[e]+XX[a]-XX[a]*PP[n]*tree$edge.length[e]
      #no need for variance if only goal is ancestral states
    }
  }
  XX[-(1:length(tree$tip.label))]
}

.to.chains.array<-function(fit,in.params){
  params<-in.params
  if(any(names(dimnames(params))=='parameters')|!is.null(attr(params,'parameters'))){
    params<-list(params)
  }
  if(!is.list(params)){
    params<-as.list(params)
  }
  types<-lapply(params,
                function(ii) if(any(names(dimnames(ii))=='parameters')|!is.null(attr(ii,'parameters'))) 'chains' else 'select')
  if(all(types=='select')&is.null(fit)){
    stop(deparse(substitute(in.params)),' appears to consist of strings/numbers specifying parameters to extract out of a corateBM_fit, but no corateBM_fit is supplied')
  }
  if(any(types=='select')&is.null(fit)){
    warning(deparse(substitute(in.params)),' contains ',paste(params[which(types=='select')],collapse=', '),', which appear to be strings/numbers specifying parameters to extract out of a corateBM_fit, but no corateBM_fit is supplied: these strings/numbers were excluded')
    params<-params[-which(types=='select')]
    types<-types[-which(types=='select')]
  }
  tmp<-lapply(1:length(params),function(ii) if(types[[ii]]=='select') .int.chains(fit,params[[ii]]) else .expand.element(params[[ii]]))
  iterations.dim<-unique(sapply(tmp,function(ii) dim(ii)[1]))
  if(length(iterations.dim)!=1){
    stop('differing number of iterations in chains specified by ',deparse(substitute(in.params)),': did these all come from the same corateBM_fit, and were they all extracted from the chains element?')
  }
  out.dimnames<-lapply(tmp,dimnames)
  chains.dim<-unique(sapply(out.dimnames,'[',3))
  if(length(chains.dim)!=1){
    stop('different chains specified by ',deparse(substitute(in.params)),': did these all come from the same corateBM_fit?')
  }
  chains.dim<-unlist(chains.dim)
  params.dim<-unlist(lapply(out.dimnames,'[',2))
  out<-array(NA,
             c(iterations.dim,length(params.dim),length(chains.dim)),
             dimnames=list(iterations=NULL,parameters=params.dim,chains=chains.dim))
  for(i in 1:length(tmp)){
    tmp.param.names<-dimnames(tmp[[i]])[[2]]
    out[,tmp.param.names,]<-tmp[[i]]
  }
  out
}
#can get duplicates of the same parameter and probs not super efficient since it runs a separate call to %chains% each time it runs...
#could cannabalize this function to add parameters to your corateBM_fit...
