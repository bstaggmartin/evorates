#' @export
opt.fun<-function(tree,trait.data,R0.prior=10,Rsig2.prior=20,X0.prior=100,
                  intra.var=F,intrasig2.prior=50,
                  trend=F,Rmu.prior=Rsig2.prior,
                  evosig2.prior=1,evocor.prior=1,intracor.prior=1,
                  constrain.Rsig2=F,
                  report.quantiles=c(0.025,0.5,0.975),report.means=T,report.devs=T,report.MAPs=T,
                  return.stanfit=F,
                  include.warmup=F,
                  ...){
  #initial checks
  Y<-trait.data
  Ysig2.prior<-intrasig2.prior
  Xsig2.prior<-evosig2.prior
  Xcor.prior<-evocor.prior
  Ycor.prior<-intracor.prior
  if(is.data.frame(Y)){
    if(!is.null(Y$tip.label)){
      labels<-Y$tip.label
      labels.col<-which(colnames(Y)=='tip.label')
    }else if(!is.null(rownames(Y))){
      labels<-rownames(Y)
      warning("trait data is a data.frame with no 'tip.label' column: used rownames as labels instead",
              immediate.=T)
    }else{
      stop('trait data is unlabelled: please name each row with its corresponding tip label or add tip.label column')
    }
    non.num.cols<-sapply(1:ncol(Y),function(ii) !is.numeric(Y[,ii]))
    if(any(non.num.cols)){
      if(all(non.num.cols)){
        stop('no numeric trait data provided')
      }
      non.num.cols.inds<-which(non.num.cols)
      if(exists('labels.col')){
        non.num.cols.inds<-non.num.cols.inds[non.num.cols.inds!=labels.col]
        tmp.non.num.cols.inds<-c(non.num.cols.inds[non.num.cols.inds<labels.col],
                                 non.num.cols.inds[non.num.cols.inds>labels.col]-1)
        if(length(X0.prior)==ncol(Y)-1){
          X0.prior<-X0.prior[-tmp.non.num.cols.inds]
        }
        if(length(Xsig2.prior)==ncol(Y)-1){
          Xsig2.prior<-Xsig2.prior[-tmp.non.num.cols.inds]
        }
        if(length(Ysig2.prior)==ncol(Y)-1){
          Ysig2.prior<-Ysig2.prior[-tmp.non.num.cols.inds]
        }
        Y<-Y[,-c(labels.col,non.num.cols.inds)]
        if(length(non.num.cols.inds)>0){
          warning('trait data in column(s) ',
                  paste(non.num.cols.inds,collapse=', '),' is non-numeric: these data were removed',
                  immediate.=T)
        }
      }else{
        if(length(X0.prior)==ncol(Y)-1){
          X0.prior<-X0.prior[-non.num.cols.inds]
        }
        if(length(Xsig2.prior)==ncol(Y)-1){
          Xsig2.prior<-Xsig2.prior[-non.num.cols.inds]
        }
        if(length(Ysig2.prior)==ncol(Y)-1){
          Ysig2.prior<-Ysig2.prior[-non.num.cols.inds]
        }
        Y<-Y[,-non.num.cols.inds]
        warning('trait data in column(s) ',
                paste(non.num.cols.inds,collapse=', '),'is non-numeric: these data were removed',
                immediate.=T)
      }
    }
    Y<-as.matrix(Y)
    rownames(Y)<-labels
  }else if(is.numeric(Y)){
    is.mat<-is.matrix(Y)
    Y<-as.matrix(Y)
    if(is.null(rownames(Y))){
      stop('trait data is unlabelled: please name each ',
           if(is.mat) 'row' else 'element',' with its corresponding tip label')
    }
  }else{
    stop("trait data format isn't recognized: please format it as a numeric vector/matrix or data.frame: see help file for more details")
  }
  k<-ncol(Y)
  colname.lens<-nchar(colnames(Y))
  if(is.null(colnames(Y))){
    colnames(Y)<-paste('X',1:k,sep='')
  }else if(any(colname.lens==0)){
    colnames(Y)<-ifelse(colname.lens==0,paste('X',1:k,sep=''),colnames(Y))
  }
  unmatched.names<-!(rownames(Y)%in%tree$tip.label)
  if(sum(unmatched.names)>0){
    if(all(unmatched.names)){
      stop('no matches between tree tip and trait data labels')
    }
    ind.unmatched.names<-which(unmatched.names)
    warning('could not find tree tip labels matching with ',
            paste(unique(rownames(Y)[unmatched.names]),collapse=', '),': these data were removed',
            immediate.=T)
    Y<-as.matrix(Y[-ind.unmatched.names,])
  }
  if(constrain.Rsig2&!trend&report.devs){
    report.devs<-F
    warning('report.devs was set to FALSE: rate deviations are meaningless with no rate heterogeneity',
            immediate.=T)
  }
  if(hasArg(chains)){
    nchain<-list(...)$chains
  }else{
    nchain<-4
  }
  
  #process/format data
  n<-length(tree$tip.label)
  e<-nrow(tree$edge)
  eV<-edge.vcv(tree)
  if(intra.var){
    is.obs<-!as.matrix(apply(Y,1,is.na))
    if(k>1){
      is.obs<-t(is.obs)
    }
    codes<-apply(is.obs,1,function(ii) sum(ii*10^((length(ii)-1):0)))
    if(any(codes==0)){
      null.obs<-which(codes==0)
      tmp.codes<-codes
      Y<-Y[-null.obs,]
      is.obs<-is.obs[-null.obs,]
      codes<-codes[-null.obs]
      for(i in sort(ind.unmatched.names)){
        tmp.codes<-append(tmp.codes,-1,i)
      }
      report.null.obs<-which(tmp.codes==0)
      warning('found no trait data for observation(s) ',
              paste(report.null.obs,collapse=', '),': these data were removed')
    }
    X_id<-match(rownames(Y),tree$tip.label)
    code_sizes<-tapply(codes,codes,length)
    obs_code<-unlist(lapply(names(code_sizes),function(ii) which(codes==as.numeric(ii))))
    names(code_sizes)<-paste(lapply(k-nchar(names(code_sizes)),function(ii) paste(rep(0,ii),collapse='')),
                             names(code_sizes),sep='')
    n_code<-length(code_sizes)
    parsed.codes<-gregexpr('1',names(code_sizes))
    code_ks<-lengths(parsed.codes)
    code_key<-unlist(parsed.codes)
  }else{
    col.names<-colnames(Y)
    n.obs<-sapply(rownames(Y),function(ii) sum(rownames(Y)==ii))
    if(any(n.obs>1)){
      warning("multiple observations per tip found: observations were averaged, but it's strongly recommend to run function with intra.var set to TRUE",
              immediate.=T)
      Y<-as.matrix(sapply(1:k,function(ii) tapply(Y[,ii],rownames(Y),mean,na.rm=T)))
      Y[is.nan(Y)]<-NA
    }
    missing.Y<-names(which(sapply(tree$tip.label,function(ii) sum(rownames(Y)==ii))==0))
    Y<-do.call(rbind,c(list(Y),setNames(rep(list(NA),length(missing.Y)),missing.Y)))
    Y<-as.matrix(Y[tree$tip.label,])
    colnames(Y)<-col.names
    is.obs<-!as.matrix(apply(Y,1,is.na))
    if(k>1){
      is.obs<-t(is.obs)
    }
    which_mis<-lapply(asplit(Y,2),function(ii) which(is.na(ii)))
    if(length(which_mis)==0){
      k_mis<-rep(0,k)
    }else{
      k_mis<-lengths(which_mis)
      which_mis<-unlist(which_mis)
    }
  }
  untrans.Y<-Y
  o.hgt<-max(eV)
  tree$edge.length<-tree$edge.length/o.hgt
  eV<-eV/o.hgt
  o.Xsig2<-rep(NA,k)
  o.X0<-rep(NA,k)
  for(i in 1:k){
    tmp.Y<-Y[,i]
    missing.dat<-names(which(tapply(tmp.Y,names(tmp.Y),function(ii) sum(!is.na(ii)))==0))
    missing.dat<-c(missing.dat,tree$tip.label[!(tree$tip.label%in%names(tmp.Y))])
    if(length(missing.dat)>0){
      tmp.tree<-drop.tip(tree,missing.dat)
    }else{
      tmp.tree<-tree
    }
    tmp.X<-tapply(tmp.Y,names(tmp.Y),mean,na.rm=T)[tmp.tree$tip.label]
    o.Xsig2[i]<-sum(pic(tmp.X,multi2di(tmp.tree))^2)/n
    xx<-c(tmp.X,rep(NA,tmp.tree$Nnode))
    for(ee in nrow(tmp.tree$edge):0){
      if(ee==0){
        des<-which(tmp.tree$edge[,1]==n+1)
        xx[n+1]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
        break
      }
      nn<-tmp.tree$edge[ee,2]
      if(nn<=n){
        next
      }
      des<-which(tmp.tree$edge[,1]==nn)
      xx[nn]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
    }
    o.X0[i]<-xx[n+1]
    Y[,i]<-(Y[,i]-o.X0[i])/sqrt(o.Xsig2[i])
  }
  if(!is.binary(tree)){
    poly.nodes<-which(sapply(1:max(tree$edge),function(ii) length(which(ii==tree$edge[,1])))>2)
    d_poly<-lapply(poly.nodes,function(ii) which(ii==tree$edge[,1]))
    tmp<-sort(unlist(lapply(d_poly,function(ii) ii[seq(2,length(ii)-1)])))
    tmp<-tmp+0:(length(tmp)-1)
    real_e<-(1:(2*n-2))[-tmp]+1
    tree<-multi2di(tree,random=F)
  }else{
    real_e<-1:e+1
  }
  des_e<-sapply(1:nrow(tree$edge),function(ii) which(tree$edge[,1]==tree$edge[ii,2])+1)
  tip_e<-which(lengths(des_e)==0)
  tip_e<-tip_e[order(tree$edge[tip_e,2])]
  des_e[tip_e]<-list(rep(-1,2))
  des_e<-matrix(unlist(des_e),ncol=2,byrow=T)
  root.edges<-which(tree$edge[,1]==n+1)+1
  des_e<-rbind(root.edges,des_e)
  prune_T<-c(0,tree$edge.length)
  postorder<-((2*n-2):1)[-((2*n-2)-tip_e)]
  tip_e<-tip_e+1
  Y[!is.obs]<-0
  dat<-list('n'=n,'e'=e,'Y'=t(Y),'eV'=eV,
            'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,
            'R0_prior'=R0.prior,'Rsig2_prior'=Rsig2.prior,'Rmu_prior'=Rmu.prior,
            'X0_prior'=rep(X0.prior,length.out=k),
            'constr_Rsig2'=as.numeric(constrain.Rsig2),'constr_Rmu'=as.numeric(!trend))
  
  #run mcmc
  if(intra.var){
    dat$obs<-nrow(Y)
    datX_id<-X_id
    dat$Ysig2_prior<-rep(Ysig2.prior,length.out=k)
    if(k>1){
      dat$k<-k
      dat$obs_code<-obs_code
      dat$n_code<-n_code
      dat$code_sizes<-code_sizes
      dat$code_ks<-as.array(code_ks)
      dat$code_key<-code_key
      dat$Xsig2_prior=rep(Xsig2.prior,length.out=k)
      dat$Xcor_prior=Xcor.prior
      dat$Ycor_prior=Ycor.prior
      ret<-rstan::optimizing(object=stanmodels$intravar_multivar_corateBM,data=dat,...)
    }else{
      dat$Y<-as.vector(Y)
      ret<-rstan::optimizing(object=stanmodels$intravar_univar_corateBM,data=dat,...)
    }
  }else{
    dat$which_mis<-which_mis
    dat$postorder<-postorder
    if(k>1){
      dat$k<-k
      dat$k_mis<-k_mis
      dat$Xsig2_prior=rep(Xsig2.prior,length.out=k)
      dat$Xcor_prior=Xcor.prior
      ret<-rstan::optimizing(object=stanmodels$multivar_corateBM,data=dat,...)
    }else{
      dat$Y<-as.vector(Y)
      dat$mis<-length(which_mis)
      ret<-rstan::optimizing(object=stanmodels$univar_corateBM,data=dat,...)
    }
  }
  ret
}