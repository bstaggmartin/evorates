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
  .index.element(fit[[element]],inds,2)
}
#flip flop any question mark dohickeys...

####ARRAY MANAGEMENT####


#newest version of check.n.proc and coerce.to.array--expands any given element to a 3D array (if chains,
#quantiles, sampler, or parameter diagnostics) or a 2D matrix (if means or MAPs)
#finally got it to a point where it should work with just about anything...
#I made it just always return 3D arrays, and label means and MAPs projects with 'iterations' along
#1st dimension; just makes things easier for other functions to deal with...
#rendered function more robust in the face of 4D arrays outputted from get.cov.mat and get.trait.mat: side-benefit of now
#also reorganizing any arrays back into iterations/quantiles/diagnostics, then parameters (however many there are), then
#chains order
.expand.element<-function(arr,simplify=F){
  if(length(dim(arr))==0){
    new.arr<-as.matrix(arr)
    new.dimnames<-dimnames(new.arr)
    if(is.null(dimnames(new.arr)[[1]])){
      new.dimnames<-rep(list(NULL),2)
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
    for(i in c('quantiles','diagnostics','parameters','chains')){
      if(!is.null(attr(arr,i))){
        attr(new.arr,i)<-attr(arr,i)
      }
    }
    dimnames(new.arr)<-new.dimnames
    arr<-new.arr
  }
  dims<-names(dimnames(arr))
  dim.map<-lapply(c('iterations|quantiles|diagnostics','parameters','chains'),function(ii) grep(ii,dims))
  if(any(lengths(dim.map)==0)){
    tmp<-dim.map
    tmp[lengths(tmp)==0]<-NA
    tmp<-unlist(tmp)
    out.dim<-dim(arr)[tmp]
    out.dim[is.na(tmp)]<-1
    out.dimnames<-dimnames(arr)[tmp]
    ndims<-length(tmp)
    if(lengths(dim.map)[1]==0){
      if('quantiles'%in%names(attributes(arr))){
        names(out.dimnames)[1]<-'quantiles'
        out.dimnames[1]<-list(attr(arr,'quantiles'))
      }else if('diagnostics'%in%names(attributes(arr))){
        names(out.dimnames)[1]<-'diagnostics'
        out.dimnames[1]<-list(attr(arr,'diagnostics'))
      }else{
        names(out.dimnames)[1]<-'iterations'
        out.dimnames[1]<-list(NULL)
      }
    }
    if(lengths(dim.map)[2]==0){
      names(out.dimnames)[-c(1,ndims)]<-'parameters'
      if('parameters'%in%names(attributes(arr))){
        out.dimnames[-c(1,ndims)]<-list(attr(arr,'parameters'))
      }else{
        out.dimnames[-c(1,ndims)]<-list(NULL)
      }
    }
    if(lengths(dim.map)[3]==0){
      names(out.dimnames)[ndims]<-'chains'
      if('chains'%in%names(attributes(arr))){
        out.dimnames[ndims]<-list(attr(arr,'chains'))
      }else{
        out.dimnames[ndims]<-list(NULL)
      }
    }
    if(any(is.na(names(dimnames(arr))))){
      arr<-array(arr,out.dim,out.dimnames)
    }else{
      arr<-array(aperm(arr,unlist(dim.map)),out.dim,out.dimnames)
    }
  }else{
    dim.map<-unlist(dim.map)
    if(any(dim.map!=1:length(dim.map))){
      arr<-array(aperm(arr,dim.map),dim(arr)[dim.map],dimnames(arr)[dim.map])
    }
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
  for(i in c('quantiles','diagnostics','parameters','chains')){
    if(!is.null(attr(arr,i))){
      attr(out,i)<-attr(arr,i)
    }
  }
  out
}

#newest, generalized version of reduce.array--selects indices from arbitrary arrays without
#simplification or destroying names
.index.element<-function(arr,inds,dims,invert=F,allow.reorder=F){
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
  if(!allow.reorder){
    inds.list<-lapply(inds.list,function(ii) unique(sort(ii)))
  }
  new.dimnames<-dimnames(arr)
  new.dimnames<-lapply(1:length(inds.list),
                       function(ii) new.dimnames[[ii]][inds.list[[ii]]])
  names(new.dimnames)<-names(dimnames(arr))
  out<-array(do.call('[',c(list(arr),inds.list)),new.dims,new.dimnames)
  for(i in c('quantiles','diagnostics','parameters','chains')){
    if(!is.null(attr(arr,i))){
      attr(out,i)<-attr(arr,i)
    }
  }
  out
}

.is.corateBM.element<-function(params){
  check.vec<-c('iterations','quantiles','diagnostics','parameters','chains')
  any(names(dimnames(params))%in%check.vec)|any(names(attributes(params))%in%check.vec)
}

.get.element.type<-function(arr){
  dim1<-setNames(dim(arr)[1],names(dimnames(arr))[1])
  if(names(dim1)=='iterations'){
    if(dim1>1){
      if(grepl(paste(paste0(c('^accept_stat','stepsize','treedepth','n_leapfrog','divergent','energy','lp'),'__$'),collapse='|^'),dimnames(arr)[[2]])){
        out<-'sampler'
      }else{
        out<-'chains'
      }
    }else{
      out<-'ambiguous'
    }
  }else if(names(dim1)=='quantiles'){
    out<-'quantiles'
  }else if(names(dim1)=='diagnostics'){
    out<-'diagnostics'
  }else{
    out<-'unrecognized'
  }
  out
}

.coerce.to.3D<-function(arr){
  if(length(dim(arr))>3){
    new.dim<-dim(arr)
    new.dimnames<-dimnames(arr)
    param.dims<-which(names(dimnames(arr))=='parameters')
    param.names<-dimnames(arr)[param.dims]
    param.names<-do.call(paste,
                         c(lapply(1:length(param.names),function(ii)
                           rep(rep(param.names[[ii]],prod(lengths(param.names)[(1:length(param.names))[1:length(param.names)>ii]])),
                               each=prod(lengths(param.names)[(1:length(param.names))[1:length(param.names)<ii]]))),
                           sep=','))
    new.dim<-new.dim[-param.dims]
    new.dim<-append(new.dim,length(param.names),1)
    new.dimnames<-new.dimnames[-param.dims]
    new.dimnames<-append(new.dimnames,list(param.names),1)
    arr<-array(arr,new.dim,new.dimnames)
  }
  arr
}

.combine.elements<-function(in.params,fit,element=NULL,select.extra=NULL,simplify=T){
  params<-in.params
  if(.is.corateBM.element(params)){
    params<-list(params)
  }
  if(!is.list(params)){
    params<-as.list(params)
  }
  types<-sapply(params,function(ii) if(.is.corateBM.element(ii)) 'element' else 'select')
  if(all(types=='select')&is.null(fit)){
    stop(deparse(substitute(in.params)),' appears to consist of strings/numbers specifying parameters to extract out of a corateBM_fit, but no corateBM_fit is supplied')
  }
  if(any(types=='select')&is.null(fit)){
    warning(deparse(substitute(in.params)),' contains ',paste(params[which(types=='select')],collapse=', '),', which appear to be strings/numbers specifying parameters to extract out of a corateBM_fit, but no corateBM_fit is supplied: these strings/numbers were excluded')
    params<-params[-which(types=='select')]
    types<-types[-which(types=='select')]
  }
  params[types=='element']<-lapply(params[types=='element'],function(ii) .coerce.to.3D(.expand.element(ii)))
  select.params<-NULL
  if(sum(types=='select')>0){
    select<-unlist(params[types=='select'])
    if(is.null(element)){
      if(sum(types=='element')>0){
        element.types<-sapply(params[types=='element'],.get.element.type)
        if(any(element.types%in%c('ambiguous','unrecognized'))|length(unique(element.types))>1){
          stop('element is unspecified, but element type based on provided loose elements is ambiguous: try specifying which element you wish to extract and double-check all loose elements are of the same type (chains, quantiles, means, etc.)')
        }else{
          element<-unique(element.types)
        }
      }else{
        element<-'chains'
      }
    }else{
      try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics','sampler')),silent=T)
      if(inherits(try.element,'try-error')){
        stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', 'diagnostics', or 'sampler'")
      }
      element<-try.element
    }
    if(element=='quantiles'|element=='diagnostics'){
      if(is.null(select.extra)&sum(types=='select')>0){
        dim1.names<-lapply(params[types=='element'],function(ii) dimnames(ii)[[1]])
        if(length(unique(dim1.names))>1){
          stop('mismatching quantiles and/or parameter diagnostics in provided loose elements')
        }else if(element=='quantiles'){
          select.extra<-as.numeric(substr(dim1.names[[1]],0,nchar(dim1.names[[1]])-1))/100
        }else{
          select.extra<-dim1.names
        }
      }
    }
    if(!is.null(select.extra)){
      select<-list(select,select.extra)
    }
    select.params<-do.call(paste0('.int.',element),list(fit=fit,select=select))
  }
  if(is.null(params[types=='element'])){
    out<-select.params
  }else{
    if(!is.null(select.params)){
      out<-c(list(select.params),params[types=='element'])
    }else{
      out<-params[types=='element']
    }
    out.dim<-vector(mode='list',length=3)
    out.dimnames<-vector(mode='list',length=3)
    for(i in 1:3){
      if(i==2){
        out.dim[[i]]<-sum(sapply(out,function(ii) dim(ii)[i]))
        out.dimnames[[i]]<-unlist(lapply(out,function(ii) dimnames(ii)[[i]]))
      }else{
        out.dim[[i]]<-unique(sapply(out,function(ii) dim(ii)[i]))
        out.dimnames[[i]]<-unique(lapply(out,function(ii) dimnames(ii)[[i]]))
      }
    }
    if(all(lengths(c(out.dim[-2],out.dimnames[-2]))==1)){
      for(i in c(1,3)){
        out.dim[[i]]<-out.dim[[i]][[1]]
        out.dimnames[i]<-out.dimnames[[i]][1]
      }
      out.dim<-unlist(out.dim)
    }else{
      stop('dimensional mismatch in elements specified by ',deparse(substitute(in.params)),': did these all come from the same corateBM_fit, and were they all extracted using the same parameters?')
    }
    names(out.dimnames)<-c(names(dimnames(out[[1]]))[1],'parameters','chains')
    out.arr<-array(NA,out.dim,out.dimnames)
    for(i in 1:length(out)){
      tmp.param.names<-dimnames(out[[i]])[[2]]
      out.arr[,tmp.param.names,]<-out[[i]]
    }
    out<-out.arr
  }
  if(simplify){
    out<-.simplify.element(out)
  }
  out
}
#can get duplicates of the same parameter and probs not super efficient since it runs a separate call to %chains% each time it runs...
#could cannabalize this function to add parameters to your corateBM_fit...
#yeah, I did a minor update to begin generalizing this, but ultimately I'd like it to be a combine.elements function...
#need to figure out what to do in cases where the element is ambiguous
#also, should 4-D arrays be coerced to 3-D with .expand.element?
#now generalized to simply combine all elements it's face with, including ones specified by select. Tries its best to match
#element type, but doesn't do anything like coercing provided elements to other types on the fly (considering doing this in the
#future...)

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

.lin.interp<-function(x,length.out){
  xx<-seq(1,length(x),length.out=length.out)
  in.inds<-xx%in%(1:length(xx))
  out<-rep(NA,length.out)
  out[in.inds]<-x[xx[in.inds]]
  xx<-xx[!in.inds]
  out[!in.inds]<-(x[ceiling(xx)]-x[floor(xx)])/(ceiling(xx)-floor(xx))*(xx-floor(xx))+x[floor(xx)]
  out
}
