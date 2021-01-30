#remove non-numeric/all NA rows, coerce to matrix with labelled rows/columns
#allow naming of prior components, and make sure to coerce them to have same number of columns as Y and
#such before exporting them from this function...
.format.trait.data<-function(Y,X0.prior.mu=NULL,X0.prior.sig=NULL,
                             X.prior.mu=NULL,X.prior.sig=NULL,scale.X.prior=T,
                             Xsig2.prior=NULL,
                             Ysig2.prior=NULL,tree){
  #if any of the priors have names, attempt to reorder according to trait.data
  for(i in c('X0.prior.mu','X0.prior.sig','X.prior.mu','X.prior.sig','scale.X.prior','Xsig2.prior','Ysig2.prior')){
    if(i=='Xsig2.prior'){
      j<-'evosig2.prior'
    }else if(i=='Ysig2.prior'){
      j<-'intrasig2.prior'
    }else{
      j<-i
    }
    tmp<-get(i)
    if(!is.null(tmp)){
      if(!is.null(names(tmp))){
        if(!is.null(colnames(Y))){
          prior.name.mismatch<-which(!(names(tmp)%in%colnames(Y)))
          if(all(prior.name.mismatch)){
            warning('none of names for ',j,' matched trait names supplied in trait.data: names were ignored and ',j,' was assumed to be in same order as columns of trait.data (missing entries set to default prior values)',
                    immediate.=T)
          }
          if(any(prior.name.mismatch)){
            warning('names ',names(tmp)[prior.name.mismatch],' in ',j,' did not match trait names supplied in trait.data: these entries in ',j,' were ignored',
                    immediate.=T)
          }
          tmp<-tmp[colnames(Y)]
          names(tmp)<-colnames(Y)
          label.col.test<-which(colnames(Y)=='tip.label')
          if(length(label.col.test)>0){
            tmp<-tmp[-label.col.test]
          }
          assign(i,tmp)
        }else{
          warning(j,' has names, but trait.data has no trait names: names were ignored and ',j,' was assumed to be in same order as columns of trait.data (missing entries set to default prior values)',
                  immediate.=T)
          assign(i,unname(tmp))
        }
      }
    }
  }
  if(is.data.frame(Y)){
    if(!is.null(Y$tip.label)){
      Y$tip.label<-as.character(Y$tip.label)
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
        if(length(non.num.cols.inds)>0){
          tmp.non.num.cols.inds<-c(non.num.cols.inds[non.num.cols.inds<labels.col],
                                   non.num.cols.inds[non.num.cols.inds>labels.col]-1)
          for(i in c('X0.prior.mu','X0.prior.sig','X.prior.mu','X.prior.sig','scale.X.prior','Xsig2.prior','Ysig2.prior')){
            tmp<-get(i)
            if(!is.null(tmp)){
              if(length(tmp)==ncol(Y)-1){
                assign(i,tmp[-tmp.non.num.cols.inds])
              }
            }
          }
          Y<-Y[,-c(labels.col,non.num.cols.inds),drop=F]
          warning('trait data in column(s) ',
                  paste(non.num.cols.inds,collapse=', '),' is non-numeric: these data were removed',
                  immediate.=T)
        }else{
          Y<-Y[,-labels.col,drop=F]
        }
      }else{
        for(i in c('X0.prior.mu','X0.prior.sig','X.prior.mu','X.prior.sig','scale.X.prior','Xsig2.prior','Ysig2.prior')){
          tmp<-get(i)
          if(!is.null(tmp)){
            if(length(tmp)==ncol(Y)){
              assign(i,tmp[-non.num.cols.inds])
            }
          }
        }
        Y<-Y[,-non.num.cols.inds]
        warning('trait data in column(s) ',
                paste(non.num.cols.inds,collapse=', '),' is non-numeric: these data were removed',
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
  colname.lens<-nchar(colnames(Y))
  if(is.null(colnames(Y))){
    colnames(Y)<-paste('X',1:ncol(Y),sep='')
  }else if(any(colname.lens==0)){
    colnames(Y)<-ifelse(colname.lens==0,paste('X',1:ncol(Y),sep=''),colnames(Y))
  }
  for(i in c('X0.prior.mu','X0.prior.sig','X.prior.mu','X.prior.sig','scale.X.prior','Xsig2.prior','Ysig2.prior')){
    if(i=='Xsig2.prior'){
      j<-'evosig2.prior'
    }else if(i=='Ysig2.prior'){
      j<-'intrasig2.prior'
    }else{
      j<-i
    }
    tmp<-get(i)
    if(!is.null(tmp)){
      len<-length(tmp)
      if(len>ncol(Y)){
        warning(j,' has more entries than there are columns in trait.data: excess entries trimmed from end',
                immediate.=T)
        assign(i,tmp[1:ncol(Y)])
      }else if(len<ncol(Y)){
        warning(j,' has less entries than there are columns in trait.data: missing entries set to defaults',
                immediate.=T)
        assign(i,c(tmp,rep(NA,ncol(Y)-len)))
      }
    }
  }
  list(Y=Y,X0.prior.mu=X0.prior.mu,X0.prior.sig=X0.prior.sig,Xsig2.prior=Xsig2.prior,Ysig2.prior=Ysig2.prior)
}

#collapses non-tip edges of length 0 to form polytomies, makes sure edge indices are preordered, and collapses
#any tips with a cophenetic distance of 0 (renaming a provided trait matrix accordingly)
.format.tree<-function(tree,Y=NULL){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  hgts<-node.depth.edgelength(tree)[1:length(tree$tip.label)]
  if(any(hgts==0)){
    warning("corateBM currently doesn't support tips at root of tree: these tips were removed; please specify information regarding trait values at the root using the X0 prior",
            immediate.=T)
    tree<-drop.tip(tree,tree$tip.label[hgts==0])
    #alternatively, could just prevent specifying a *fixed* value at the root...technically, unfixed values
    #ought to be fine (though having two perfectly correlated parameters might get odd?)
  }
  dist.mat<-cophenetic(tree)
  diag(dist.mat)<-NA
  if(any(dist.mat==0,na.rm=T)){
    problem.indices<-which(dist.mat==0,arr.ind=T)
    problem.indices<-t(apply(problem.indices,1,sort))
    problem.indices<-problem.indices[!duplicated(problem.indices),,drop=F]
    collapses<-list(NULL)
    for(i in 1:nrow(problem.indices)){
      test<-which(sapply(collapses,function(ii) problem.indices[i,1]%in%ii))
      if(length(test)>0){
        if(!(problem.indices[i,2]%in%collapses[[test]])){
          collapses[[test]]<-c(collapses[[test]],problem.indices[i,2])
        }
      }else{
        collapses[[length(collapses)+1]]<-problem.indices[i,]        
      }
    }
    collapses<-collapses[-1]
    label.pos<-sapply(collapses,'[[',1)
    collapses<-lapply(collapses,function(ii) tree$tip.label[ii])
    warning('detected group(s) of tips representing the same exact phylogenetic position(s): ',
            paste(sapply(collapses,paste,collapse=', '),collapse='; '),
            ': these tips were collapsed into single tips and any trait data was relabelled accordingly',
            immediate.=T)
    new.tip.labels<-sapply(collapses,function(ii) paste(unique(ii),collapse='&'))
    if(length(Y)>0){
      for(i in 1:length(collapses)){
        rownames(Y)<-gsub(paste(paste0('^',collapses[[i]],'(|_mean|_sd)$'),collapse='|'),
                          paste0(new.tip.labels[i],'\\1'),rownames(Y))
      }
    }
    tree$tip.label[label.pos]<-new.tip.labels
    tree<-drop.tip(tree,unlist(sapply(collapses,'[',-1)))
  }
  if(any(duplicated(tree$tip.label))){
    stop('two tips share the same name: please relabel the tips of the tree such that each tip has a unique label')
  }
  list(tree=tree,Y=Y)
}
#alright, so depending on how you treat the 0-length edges, you might need to modify the tree when plotting
#in non-phenogram style to drop tips of 0-length vefore plotting (but then how to treat collapsed edges???)
#I think the easiest solution is to just return 0-length edge params...iirc, ape doesn't support singleton
#nodes...
# node.sample<-sample(1:tree$Nnode+length(tree$tip.label),10,T)
# lens<-lapply(node.sample,function(ii) which(tree$edge[,2]==ii))
# lens[lengths(lens)==0]<-NA
# lens<-tree$edge.length[unlist(lens)]
# lens[is.na(lens)]<-0
# edge.length.sample<-runif(length(node.sample),0,lens)
# tree<-multi.bind.tip(tree,names=letters[1:10],edge.lengths=NULL,nodes=node.sample,positions=edge.length.sample)
# tree<-.format.tree(tree,Y=Y)
# plot(tree)
# nodelabels()
# edgelabels()

#should account for cases where Y may be empty...
#split tip priors from trait data
.split.tp<-function(Y,tree,intra.var){
  #merge rows with same name, stop if there are two entries for the same tip/trait combo
  .check.dups<-function(tp){
    dup.rows<-which(duplicated(rownames(tp)))
    if(length(dup.rows)>0){
      for(i in dup.rows){
        base<-match(rownames(tp)[i],rownames(tp))
        for(j in 1:ncol(tp)){
          if(!is.na(tp[i,j])){
            if(is.na(tp[base,j])){
              tp[base,j]<-tp[i,j]
            }else{
              stop('duplicate tip prior means specified for tip ',
                   rownames(tp)[i],
                   ' found: please only specify only one')
            }
          }
        }
      }
      tp<-tp[!dup.rows,,drop=F]
    }
    tp
  }
  #remove empty rows
  .rm.empty<-function(mat,row.tracker=NULL){
    null.rows<-apply(mat,1,function(ii) all(is.na(ii)))
    if(any(null.rows)){
      null.rows<-which(null.rows)
      report.rows<-rownames(mat)[null.rows]
      tmp<-deparse(substitute(mat))
      if(tmp=='Y'){
        report.name<-'trait data'
        if(!is.null(row.tracker)){
          report.rows<-row.tracker[null.rows]
        }
      }else if(tmp=='tp.mu'){
        report.name<-'tip prior means'
      }else{
        report.name<-'tip prior standard deviations'
      }
      warning('found no non-missing data for row(s) ',paste(report.rows,collapse=', '),
              ' in ',report.name,': these data were removed',
              immediate.=T)
      mat<-mat[-null.rows,,drop=F]
    }
    mat
  }
  tmp.rownames<-gsub('_mean$|_sd$','',rownames(Y))
  row.tracker<-1:length(tmp.rownames)
  unmatched.names<-!(tmp.rownames%in%tree$tip.label)
  if(sum(unmatched.names)>0){
    if(all(unmatched.names)){
      stop('no matches between tree tip and trait data labels')
    }
    ind.unmatched.names<-which(unmatched.names)
    warning('could not find tree tip labels matching with ',
            paste(unique(tmp.rownames[unmatched.names]),collapse=', '),': these data were removed',
            immediate.=T)
    Y<-Y[-ind.unmatched.names,,drop=F]
    row.tracker<-row.tracker[-ind.unmatched.names]
  }
  tp.mu<-.check.dups(Y[grep('_mean$',rownames(Y)),,drop=F])
  if(length(tp.mu)>0){
    rownames(tp.mu)<-gsub('_mean$','',rownames(tp.mu))
    tp.mu.drops<-grep('_mean$',rownames(Y))
    Y<-Y[-tp.mu.drops,,drop=F]
    row.tracker<-row.tracker[-tp.mu.drops]
  }
  tp.sig<-.check.dups(Y[grep('_sd$',rownames(Y)),,drop=F])
  #merge rows with same name, stop if there are two entries for the same tip/trait combo
  if(length(tp.sig)>0){
    rownames(tp.sig)<-gsub('_sd$','',rownames(tp.sig))
    tp.sig.drops<-grep('_sd$',rownames(Y))
    Y<-Y[-tp.sig.drops,,drop=F]
    row.tracker<-row.tracker[-tp.sig.drops]
  }
  if(length(tp.mu)>0&length(tp.sig)>0){
    missing.sig<-rownames(tp.mu)[!(rownames(tp.mu)%in%rownames(tp.sig))]
    missing.mu<-rownames(tp.sig)[!(rownames(tp.sig)%in%rownames(tp.mu))]
    if(length(missing.sig)>0){
      tp.sig<-do.call(rbind,c(list(tp.sig),rep(list(NA),length(missing.sig))))
      rownames(tp.sig)[-(length(missing.sig)-1):0+nrow(tp.sig)]<-missing.sig
    }
    tp.sig<-tp.sig[sort(rownames(tp.sig)),,drop=F]
    if(length(missing.mu)>0){
      tp.mu<-do.call(rbind,c(list(tp.mu),rep(list(NA),length(missing.mu))))
      rownames(tp.mu)[-(length(missing.mu)-1):0+nrow(tp.mu)]<-missing.mu
    }
    tp.mu<-tp.mu[sort(rownames(tp.mu)),,drop=F]
    problem.inds<-which(tp.sig<0,arr.ind=T)
    if(length(problem.inds)>0){
      warning('negative tip prior standard deviation(s) for the following:',
              paste(paste(rownames(tp.sig)[problem.inds[,1]],
                          colnames(tp.sig)[problem.inds[,2]],sep=','),
                    collapse='; '),
              ': these tip prior standard deviation(s) were removed',
              immediate.=T)
      tp.sig[problem.inds]<-NA
    }
    check.nas<-tp.mu+tp.sig
    check.nas[,]<-is.na(check.nas)
    if(all(as.logical(check.nas))){
      warning('incomplete specification of tip prior(s) (missing means and/or standard deviations): tip prior(s) ignored',
              immediate.=T)
      tp.mu<-NULL
      tp.sig<-NULL
    }else if(any(as.logical(check.nas))){
      problem.inds<-which(check.nas==1,arr.ind=T)
      warning('incomplete specification of tip prior(s) for the following: ',
              paste(paste(rownames(tp.mu)[problem.inds[,1]],
                          colnames(tp.mu)[problem.inds[,2]],sep=','),
                    collapse='; '),
              ': these tip priors were removed',
              immediate.=T)
      tp.mu[check.nas==1]<-NA
      tp.sig[check.nas==1]<-NA
      tp.mu<-.rm.empty(tp.mu)
      tp.sig<-.rm.empty(tp.sig)
    }
  }else if((length(tp.mu)==0&length(tp.sig)>0)|
           (length(tp.sig)==0&length(tp.mu)>0)){
    warning('tip prior ',
            if(length(tp.mu)>0) 'mean(s)' else 'standard deviation(s)',
            ' supplied, but not ',
            if(length(tp.mu)>0) 'standard deviation(s)' else 'mean(s)',
            ': tip prior(s) ignored',
            immediate.=T)
    tp.mu<-NULL
    tp.sig<-NULL
  }
  if(length(Y)>0&!intra.var){
    n.obs<-sapply(rownames(Y),function(ii) sum(rownames(Y)==ii))
    if(any(n.obs>1)){
      warning("multiple observations per tip found: observations were averaged, but it's strongly recommend to run function with intra.var set to TRUE",
              immediate.=T)
      Y<-apply(Y,2,function(ii) tapply(ii,rownames(Y),mean,na.rm=T))
      Y[is.nan(Y)]<-NA
    }
    if(length(tp.mu)>0){
      overlap<-Y[rownames(Y)%in%rownames(tp.mu),,drop=F]
      if(length(overlap)>0){
        check.nas<-overlap+tp.mu[rownames(overlap),,drop=F]
        check.nas[,]<-!is.na(check.nas)
        if(all(as.logical(check.nas))){
          warning('when intra.var is set to FALSE, tip priors may only be set for tip/traits without observed values, so all provided trait observations were removed',
                  immediate.=T)
          Y<-matrix(NA,0,ncol(tp.mu))
        }else if(any(as.logical(check.nas))){
          problem.inds<-which(check.nas==1,arr.ind=T)
          warning('when intra.var is set to FALSE, tip priors may only be set for tip/traits without observed values, so the following trait observations were removed: ',
                  paste(unique(paste(rownames(Y)[problem.inds[,1]],
                                     colnames(Y)[problem.inds[,2]],sep=',')),
                        collapse='; '),
                  immediate.=T)
          Y[rownames(overlap),][check.nas==1]<-NA
        }
      }
    }
  }
  if(length(Y)>0){
    Y<-.rm.empty(Y,if(intra.var) row.tracker else NULL)
  }
  list(Y=Y,tp.mu=tp.mu,tp.sig=tp.sig)
}

#get constants to transform data and tree to "unit scale"
#tree height of 1
#sig2 under BM is 1
.get.trans.const<-function(tree,Y,tp.mu,tp.sig,
                           X0.prior.mu=NULL,X0.prior.sig=NULL,
                           X.prior.mu=NULL,X.prior.sig=NULL,scale.X.prior=T,
                           R0.prior.mu=NULL,Rsig2.prior=NULL,
                           Ysig2.prior=NULL,
                           Rmu.prior.mu=NULL,Rmu.prior.sig=NULL){
  trans.const<-list(hgt=NULL,X_sig2=rep(NA,ncol(Y)),X_0=rep(NA,ncol(Y)))
  trans.const$hgt<-max(node.depth.edgelength(tree))
  tree$edge.length<-tree$edge.length/trans.const$hgt
  for(i in c('Rsig2.prior','Rmu.prior.mu','Rmu.prior.sig')){
    tmp<-get(i)
    if(!is.null(tmp)){
      assign(i,tmp*trans.const$hgt)
    }
  }
  for(i in 1:ncol(Y)){
    tmp.X<-c(tp.mu[,i],tapply(Y[,i],rownames(Y),mean,na.rm=T))
    tmp.X<-tmp.X[!is.na(tmp.X)&!is.nan(tmp.X)]
    tmp.X<-tmp.X[!duplicated(names(tmp.X))]
    if(length(tmp.X)<length(tree$tip.label)){
      tmp.tree<-keep.tip(tree,names(tmp.X))
    }else{
      tmp.tree<-tree
    }
    trans.const$X_sig2[i]<-mean(pic(tmp.X,multi2di(tmp.tree))^2)
    trans.const$X_0[i]<-.quick.recon(tmp.X,tmp.tree,just.root=T)
    if(length(Y)>0){
      Y[,i]<-(Y[,i]-trans.const$X_0[i])/sqrt(trans.const$X_sig2[i])
    }
    if(length(tp.mu)>0){
      tp.mu[,i]<-(tp.mu[,i]-trans.const$X_0[i])/sqrt(trans.const$X_sig2[i])
      tp.sig[,i]<-tp.sig[,i]/sqrt(trans.const$X_sig2[i])
    }
  }
  if(!is.null(X0.prior.mu)){
    X0.prior.mu<-ifelse(is.na(X0.prior.mu),NA,(X0.prior.mu-trans.const$X_0)/sqrt(trans.const$X_sig2))
  }
  if(!is.null(X.prior.mu)){
    X.prior.mu<-ifelse(is.na(X.prior.mu),NA,
                       ifelse(scale.X.prior,(X.prior.mu-trans.const$X_0)/sqrt(trans.const$X_sig2),X.prior.mu))
  }
  if(!is.null(X.prior.sig)){
    X.prior.sig<-ifelse(is.na(X.prior.sig),NA,
                        ifelse(scale.X.prior,X.prior.sig/sqrt(trans.const$X_sig2),X.prior.sig))
  }
  for(i in c('X0.prior.sig','Ysig2.prior')){
    tmp<-get(i)
    if(!is.null(tmp)){
      tmp<-ifelse(is.na(tmp),NA,tmp/sqrt(trans.const$X_sig2))
      assign(i,tmp)
    }
  }
  if(!is.null(R0.prior.mu)){
    R0.prior.mu<-R0.prior.mu+log(trans.const$hgt)-log(mean(trans.const$X_sig2))
  }
  list(tree=tree,Y=Y,tp.mu=tp.mu,tp.sig=tp.sig,trans.const=trans.const,
       X0.prior.mu=X0.prior.mu,X0.prior.sig=X0.prior.sig,
       X.prior.mu=X.prior.mu,X.prior.sig=X.prior.sig,
       R0.prior.mu=R0.prior.mu,Rsig2.prior=Rsig2.prior,
       Ysig2.prior=Ysig2.prior,
       Rmu.prior.mu=Rmu.prior.mu,Rmu.prior.sig=Rmu.prior.sig)
}

.prep.run<-function(corateBM.form.dat,return.as.obj=T,out.file=NULL,check.overwrite=T,nchain,
                    constrain.Rsig2=NULL,trend=NULL,lik.power=NULL,
                    X0.prior.mu=NULL,X0.prior.sig=NULL,
                    evosig2.prior=NULL,evocor.prior=NULL,
                    R0.prior.mu=NULL,R0.prior.sig=NULL,Rsig2.prior=NULL,
                    intrasig2.prior=NULL,intracor.prior=NULL,
                    Rmu.prior.mu=NULL,Rmu.prior.sig=NULL){
  ##HANDLING FILE OUTPUTS##
  if(!return.as.obj|!is.null(out.file)){
    if(is.null(out.file)){
      message("no file output name specified, but model results aren't set to be returned in R working environment:\n picked file output name automatically based off current date and time")
      time<-Sys.time()
      out.file<-paste0('corateBM_',gsub(' |-|:','_',time))
    }
    directory<-dirname(out.file)
    if(!file.exists(directory)){
      directory.check<-dir.create(directory)
      if(!directory.check){
        stop('the directory for the output file could not be created: make sure you have permission to write in the specified directory and that there are no illegal characters or trailing spaces/periods')
      }
    }
    file.strings<-paste0(gsub('\\.csv$','',out.file),c(paste0('_',1:nchain,'.csv'),'_info'))
    if(check.overwrite){
      existing.files<-list.files(directory)
      file.strings.check<-sapply(basename(file.strings),'%in%',existing.files)
      if(any(file.strings.check)){
        readline('WARNING: running this function will overwrite existing files! Press Enter to continue or Esc to cancel')
      }
    }
    file.check<-file.create(file.strings)
    if(any(!file.check)){
      stop('could not create test output files: check for any illegal characters or trailing spaces/periods')
    }else{
      file.remove(file.strings)
    }
    if(nchain==1){
      out.file<-paste0(gsub('\\.csv$','',out.file),'_1.csv')
    }
    if(!grepl('\\.csv$',out.file)){
      out.file<-paste0(out.file,'.csv')
    }
  }
  
  ##FINAL PRIOR/CONSTRAINT SETTINGS##
  dat<-corateBM.form.dat$dat
  #recheck for any inappropriately negative priors
  for(i in c('X0.prior.sig',
             'evosig2.prior','evocor.prior',
             'R0.prior.sig','Rsig2.prior',
             'intrasig2.prior','intracor.prior',
             'Rmu.prior.sig')){
    tmp<-get(i)
    if(any(tmp<0)){
      warning(i,' must consist of postive numbers: entries that were negative or 0 set to defaults',
              immediate.=T)
      if(length(tmp)==1){
        tmp<-NULL
      }else{
        tmp[tmp<=0]<-NA
      }
      assign(i,tmp)
    }
  }
  for(i in c('X0.prior.mu','X0.prior.sig',
             'evosig2.prior','evocor.prior',
             'R0.prior.mu','R0.prior.sig','Rsig2.prior',
             'intrasig2.prior','intracor.prior',
             'Rmu.prior.mu','Rmu.prior.sig')){
    tmp<-get(i)
    if(!is.null(tmp)){
      if(!all(is.na(tmp))){
        dat[[gsub('\\.','_',gsub('Ysig2','intrasig2',gsub('Xsig2','evosig2',i)))]]<-tmp
      }
    }
  }
  def.priors<-setNames(list(rep(0,dat$k),rep(20,dat$k),rep(1,dat$k),1,0,10,10,rep(10,dat$k),1,0,10),
                       c('X0_prior_mu','X0_prior_sig','Xsig2_prior','Xcor_prior',
                         'R0_prior_mu','R0_prior_sig','Rsig2_prior',
                         'Ysig2_prior','Ycor_prior',
                         'Rmu_prior_mu','Rmu_prior_sig'))
  for(i in grep('prior',names(dat),value=T)){
    if(is.null(dat[[i]])){
      dat[[i]]<-def.priors[[i]]
    }else if(any(is.na(dat[[i]]))){
      dat[[i]]<-ifelse(is.na(dat[[i]]),def.priors[[i]],dat[[i]])
    }
    if(dat$k>1&i%in%c('X0_prior_mu','X0_prior_sig','Xsig2_prior','Ysig2_prior')){
      dat[[i]]<-array(dat[[i]])
    }
  }
  if(!is.null(constrain.Rsig2)){
    dat$constr_Rsig2<-as.numeric(constrain.Rsig2)
  }
  if(!is.null(trend)){
    dat$constr_Rmu<-as.numeric(!trend)
  }
  if(!is.null(lik.power)){
    if(lik.power<0){
      warning('likelihood powers must be 0 or positive: likelihood power set to default of 1 (i.e., sampling from normal posterior)',
              immediate.=T)
      lik.power<-1
    }
    dat$lik_power<-lik.power
  }
  intra.var<-any(grepl('X_id',names(dat)))
  corateBM.form.dat$dat<-dat
  
  ##RUN MODEL##
  exclude.pars<-c('std_R0','std_X0','std_Rsig2','std_Rmu','raw_R','trans_tp')
  if(dat$constr_Rsig2){
    exclude.pars<-c(exclude.pars,'Rsig2')
  }
  if(dat$constr_Rmu){
    exclude.pars<-c(exclude.pars,'Rmu')
  }
  if(dat$constr_Rsig2&dat$constr_Rmu){
    exclude.pars<-c(exclude.pars,'R')
  }
  if(intra.var){
    exclude.pars<-c(exclude.pars,'std_Ysig2','raw_X','cent_Y')
    if(dat$k>1){
      exclude.pars<-c(exclude.pars,'std_Xsig2','Xsig2','Ysig2','Ycov')
      stanobj<-'intravar_multivar_corateBM'
    }else{
      dat<-dat[-which(names(dat)%in%c('Xsig2_prior','Xcor_prior','Ycor_prior'))]
      dat$Y<-as.vector(dat$Y)
      stanobj<-'intravar_univar_corateBM'
    }
  }else{
    dat<-dat[-which(names(dat)%in%c('Ysig2_prior','Ycor_prior'))]
    if(dat$k>1){
      exclude.pars<-c(exclude.pars,'std_Xsig2','Xsig2','Xcov')
      stanobj<-'multivar_corateBM'
    }else{
      dat<-dat[-which(names(dat)%in%c('Xsig2_prior','Xcor_prior'))]
      dat$Y<-as.vector(dat$Y)
      stanobj<-'univar_corateBM'
    }
  }
  
  out<-list(return.as.obj=return.as.obj,stanobj=stanobj,dat=dat,exclude.pars=exclude.pars,out.file=out.file,
            corateBM.form.dat=corateBM.form.dat)
  if(!is.null(out.file)){
    out<-c(out,file.strings=list(file.strings),directory=directory)
  }
  out
}

#shit, I would need a way to get the trait matrix...ugh!
#I think...just run it at the end of the function. Stanfit should still be in memory, can just re-extract
#less efficient but clean way to do it
#I think the below should work
.get.lik<-function(fit,stanfit,trans.const,dat){
  trait.mat<-get.trait.mat(.coerce.fit(fit,stanfit,trans.const,dat),simplify=F)
  dims<-dim(trait.mat)
  n<-dims[1];k<-dims[2];iter<-dims[3];chains<-dims[4]
  for(i in 1:k){
    trait.mat[,i,,]<-(trait.mat[,i,,]-trans.const$X_0[i])/sqrt(trans.const$X_sig2)[i]
  }
  out<-matrix(NA,iter,chains)
  #if intra.var, just a bunch multivariate normal liks
  if(fit$call$intra.var){
    foo<-function(iter,chain,cent.trait.data,inv.chol.intracov,obs,codes,n_code,log.dets){
      trans.trait.data<-cent.trait.data[,,iter,chain]
      for(i in 1:dat$n_code){
        tmp<-as.matrix(cent.trait.data[obs[[i]],codes[[i]],iter,chain])
        if(length(obs[[i]])==1){
          tmp<-t(tmp)
        }
        trans.trait.data[obs[[i]],codes[[i]]]<-
          tmp%*%
          inv.chol.intracov[[i]][,,iter,chain]
      }
      sum(dnorm(trans.trait.data,log=T),na.rm=T)+sum(lengths(obs)*unlist(sapply(log.dets,'[[',iter,chain)))
    }
    chol.intracov<-extract(stanfit,"chol_Ycov",permute=F,inc_warmup=T)
    #this is all  wrong, unfortunately--need to rework array management
    #I think the below should work
    inv.chol.intracov<-array(0,c(k,k,iter,chains))
    for(i in 1:k){
      for(j in 1:i){
        inv.chol.intracov[j,i,,]<-as.vector(chol.intracov[,,paste('chol_Ycov[',j,',',i,']',sep='')])
      }
    }
    inv.chol.intracov<-rep(list(inv.chol.intracov),dat$n_code)
    log.dets<-vector('list',dat$n_code)
    codes<-split(dat$code_key,rep(1:dat$n_code,dat$code_ks))
    obs<-split(dat$obs_code,rep(1:dat$n_code,dat$code_sizes))
    for(i in 1:dat$n_code){
      inv.chol.intracov[[i]]<-inv.chol.intracov[[i]][codes[[i]],codes[[i]],,,drop=F]
      inv.chol.intracov[[i]][,,,]<-unlist(lapply(asplit(inv.chol.intracov[[i]],c(3,4)),solve))
      log.dets[[i]]<-apply(inv.chol.intracov[[i]],c(3,4),function(ii) log(prod(diag(ii))))
    }
    trait.data<-dat$Y
    if(k>1){
      trait.dat<-t(trait.data)
    }
    exp.trait.mat<-trait.mat[rownames(trait.data),,,,drop=F]
    cent.trait.data<-exp.trait.mat
    cent.trait.data[,,,]<-trait.data
    cent.trait.data<-cent.trait.data-exp.trait.mat
    for(i in 1:chains){
      out[,i]<-unlist(lapply(1:iter,foo,i,cent.trait.data,inv.chol.intracov,obs,codes,dat$n_code,log.dets))
    }
  #if no intra.var, use pruning algorithm
  }else{
    out[,]<-.pruning.alg(trait.mat,fit,stanfit,dat,n,k,iter,chains)
  }
  out
}

#or, you could integrate the prior into the get.lik function, and return them both as a list?
.get.prior<-function(fit,stanfit,trans.const,dat){
  trait.mat<-get.trait.mat(.coerce.fit(fit,stanfit,trans.const,dat),simplify=F)
  dims<-dim(trait.mat)
  n<-dims[1];k<-dims[2];iter<-dims[3];chains<-dims[4]
  for(i in 1:k){
    trait.mat[,i,,]<-(trait.mat[,i,,]-trans.const$X_0[i])/sqrt(trans.const$X_sig2)[i]
  }
  out<-matrix(0,iter,chains)
  #high level
  R0<-extract(stanfit,"R0",permute=F,inc_warmup=T)
  out[,]<-out+as.vector(dcauchy(R0,dat$R0_prior_mu,dat$R0_prior_sig,log=T))
  X0<-extract(stanfit,"X0",permute=F,inc_warmup=T)
  X0<-aperm(X0,c(1,3,2))
  tmp<-X0
  tmp[,,]<-dcauchy(X0,rep(dat$X0_prior_mu,each=iter),rep(dat$X0_prior_sig,each=iter),log=T)
  out[,]<-out+as.vector(apply(tmp,c(1,3),sum))
  if(fit$call$trend){
    Rmu<-extract(stanfit,"Rmu",permute=F,inc_warmup=T)
    out[,]<-out+as.vector(dcauchy(Rsig2,dat$Rmu_prior_mu,dat$Rmu_prior_sig,log=T))
  }
  if(!fit$call$constrain.Rsig2){
    Rsig2<-extract(stanfit,"Rsig2",permute=F,inc_warmup=T)
    out[,]<-out+as.vector(dcauchy(Rsig2,0,dat$Rsig2_prior,log=T))+log(2)
    #R prior
    invchol.eV<-solve(chol(dat$eV))
    log.det<-sum(log(diag(invchol.eV)))
    R<-extract(stanfit,"R",permute=F,inc_warmup=T)
    mus<-sigs<-R
    mus[,,]<-R0
    sigs[,,]<-Rsig2
    if(fit$call$trend){
      Tl<-dat$prune_T[dat$real_e]
      Tpts<-diag(dat$eV)
      T1<-rep(Tpts-Tl/3,each=iter)
      T2<-rep(Tpts+2*Tl/3,each=iter)
      Tl<-rep(Tl,each=iter)
      mus[,,]<- -log(abs(Rmu))-log(Tl)+log(abs(exp(Rmu*T2)-exp(Rmu*T1)))
    }
    mus<-aperm(mus,c(1,3,2))
    sigs<-aperm(sigs,c(1,3,2))
    R<-aperm(R,c(1,3,2))
    tmp<-R
    tmp[,,]<-unlist(lapply(1:chains,function(ii) ((R[,,ii]-mus[,,ii])/sqrt(sigs[,,ii]))%*%invchol.eV))
    out[,]<-out+as.vector(apply(dnorm(tmp,log=T)-0.5*log(sigs),c(1,3),sum))+log.det
  }
  if(fit$call$intra.var){
    #pruning alg if intravar
    dat$postorder<-((2*n-2):1)[-((2*n-2)-dat$tip_e)]
    out[,]<-out+as.vector(.pruning.alg(trait.mat,fit,stanfit,dat,n,k,iter,chains))
    #intraspecific variance(s) prior
    Ysig2<-extract(stanfit,"Ysig2",permute=F,inc_warmup=T)
    Ysig2<-aperm(Ysig2,c(1,3,2))
    tmp<-Ysig2
    tmp[,,]<-dcauchy(X0,0,rep(dat$Ysig2_prior,each=iter),log=T)+log(2)
    out[,]<-out+as.vector(apply(tmp,c(1,3),sum))
    #intraspecific correlation prior
    if(k>1){
      chol.Ycor<-extract(stanfit,"chol_Ycor",permute=F,inc_warmup=T)
      Ycor<-array(0,c(k,k,iter,chains))
      for(i in 1:k){
        for(j in 1:i){
          Ycor[i,j,,]<-as.vector(chol.Ycor[,,paste('chol_Ycor[',i,',',j,']',sep='')])
        }
      }
      Ycor[,,,]<-unlist(lapply(asplit(Ycor,c(3,4)),function(ii) ii%*%t(ii)))
      out[,]<-out+unlist(lapply(1:chains,function(ii) dlkj(Ycor[,,,ii],dat$Ycor_prior,log=T)))
    }
  }
  #if multivariate, evolutionary variances and correlation prior
  if(k>1){
    Xsig2<-extract(stanfit,"Xsig2",permute=F,inc_warmup=T)
    Xsig2<-aperm(Xsig2,c(1,3,2))
    Xsig2<-Xsig2/k
    out[,]<-out+unlist(lapply(1:chains,function(ii) ddiri(Xsig2[,,ii],dat$Xsig2_prior,log=T)))
    chol.Xcor<-extract(stanfit,"chol_Xcor",permute=F,inc_warmup=T)
    Xcor<-array(0,c(k,k,iter,chains))
    for(i in 1:k){
      for(j in 1:i){
        Xcor[i,j,,]<-as.vector(chol.Xcor[,,paste('chol_Xcor[',i,',',j,']',sep='')])
      }
    }
    Xcor[,,,]<-unlist(lapply(asplit(Xcor,c(3,4)),function(ii) ii%*%t(ii)))
    out[,]<-out+unlist(lapply(1:chains,function(ii) dlkj(Xcor[,,,ii],dat$Xcor_prior,log=T)))
  }
  #tip priors
  if(length(fit$call$tip.prior.means)>0){
    if(fit$call$intra.var){
      #technically inefficient, but just to make sure all indices line up...
      if(k>1){
        X<-extract(stanfit,"X",permute=F,inc_warmup=T)
        new.X<-array(0,c(n,k,iter,chains))
        for(i in 1:n){
          for(j in 1:k){
            new.X[i,j,,]<-as.vector(X[,,paste('X[',j,',',i,']',sep='')])
          }
        }
        fac<-rep(1:k,n_tp)
        which_tps<-split(dat$which_tp,fac)
        tp_mus<-split(dat$tp_mu,fac)
        tp_sigs<-split(dat$tp_sig,fac)
        tmp<-array(0,c(iter,n,chains))
        for(i in 1:k){
          ind<-as.character(i)
          if(length(which_tps[[ind]])>0){
            tmp[,,]<-aperm(X[,i,,],c(2,1,3))
            tmp[,,]<-dnorm(tmp[,which_tps[[ind]],],rep(tp_mus[[ind]],each=iter),rep(tp_sigs[[ind]],each=iter),log=T)
            out[,]<-out+as.vector(apply(tmp,c(1,3),sum))
          }
        }
      }else{
        X<-extract(stanfit,"X",permute=F,inc_warmup=T)
        X<-aperm(X,c(1,3,2))
        tmp<-X
        tmp[,,]<-dnorm(X[,dat$which_tp,],rep(dat$tp_mu,each=iter),rep(dat$tp_sig,each=iter),log=T)
        out[,]<-out+as.vector(apply(tmp,c(1,3),sum))
      }
    }else{
      misY<-extract(stanfit,"mis_Y",permute=F,inc_warmup=T)
      misY<-aperm(misY,c(1,3,2))
      tmp<-misY
      tmp[,,]<-dnorm(misY[,dat$which_tp,],rep(dat$tp_mu,each=iter),rep(dat$tp_sig,each=iter),log=T)
      out[,]<-out+as.vector(apply(tmp,c(1,3),sum))
    }
  }
  out
}

.pruning.alg<-function(trait.mat,fit,stanfit,dat,n,k,iter,chains){
  foo<-function(iter,chain,contra,inv.evocov){
    tmp.contra<-contra[,iter,chain]
    tmp.evocov<-inv.evocov[,,iter,chain]
    t(tmp.contra)%*%tmp.evocov%*%tmp.contra
  }
  if(k>1){
    chol.evocov<-extract(stanfit,"chol_Xcov",permute=F,inc_warmup=T)
    inv.evocov<-array(0,c(k,k,iter,chains))
    for(i in 1:k){
      for(j in 1:i){
        inv.evocov[j,i,,]<-as.vector(chol.evocov[,,paste('chol_Ycov[',j,',',i,']',sep='')])
      }
    }
    log.dets<-apply(evocov,c(3,4),function(ii) log(det(ii)))
    inv.evocov[,,,]<-unlist(lapply(asplit(inv.evocov,c(3,4)),function(ii) solve(t(ii)%*%ii)))
  }else{
    log.dets<-0
  }
  if(!fit$call$constrain.Rsig2|fit$call$trend){
    R<-extract(stanfit,"R",permute=F,inc_warmup=T)
    R<-aperm(R,c(1,3,2))
  }else{
    R<-extact(stanfit,"R_0",permute=F,inc_warmup=T)
    R<-aperm(array(R,dim=c(iter,chains,length(dat$real_e))),c(1,3,2))
  }
  X0<-extract(stanfit,"X0",permute=F,inc_warmup=T)
  X0<-aperm(X0,c(3,1,2))
  TT<-dat$prune_T
  SS<-array(0,dim=c(2*n-1,iter,chains))
  VV<-SS
  SS[dat$real_e,,]<-exp(aperm(R,c(2,1,3)))*TT[dat$real_e]
  XX<-array(NA,dim=c(k,2*n-1,iter,chains))
  XX[,dat$tip_e,,]<-aperm(trait.mat,c(2,1,3,4))
  tmp.contra<-array(NA,dim=c(iter,chains))
  LL<-array(NA,dim=c(length(dat$postorder),iter,chains))
  counter<-0
  contra<-array(NA,dim=c(k,iter,chains))
  sum.Vs<-array(NA,dim=c(iter,chains))
  for(i in dat$postorder){
    des.X<-XX[,dat$des_e[i,],,,drop=F]
    des.V<-VV[dat$des_e[i,],,,drop=F]+SS[dat$des_e[i,],,,drop=F]
    contra[,,]<-des.X[,2,,]-des.X[,1,,]
    if(k>1){
      for(j in 1:chains){
        tmp.contra[,j]<-unlist(lapply(1:iter,foo,chain=j,contra=contra,inv.evocov=inv.evocov))
      }
    }else{
      tmp.contra[,]<-contra^2
    }
    sum.Vs[,]<-des.V[1,,]+des.V[2,,]
    counter<-counter+1
    LL[counter,,]<- -0.5*(k*log(sum.Vs)+tmp.contra/sum.Vs)
    XX[,i,,]<-rep(des.V[2,,]/sum.Vs,each=k)*des.X[,1,,]+rep(des.V[1,,]/sum.Vs,each=k)*des.X[,2,,]
    VV[i,,]<-1/(1/des.V[1,,]+1/des.V[2,,])
  }
  contra[,,]<-as.vector(XX[,1,,])-X0
  if(k>1){
    for(j in 1:chains){
      tmp.contra[,j]<-unlist(lapply(1:iter,foo,chain=j,contra=contra,inv.evocov=inv.evocov))
    }
  }else{
    tmp.contra[,]<-contra^2
  }
  apply(LL,c(2,3),sum)-0.5*(n*log.dets+k*log(VV[1,,])+tmp.contra/VV[1,,])
}

.coerce.fit<-function(fit,stanfit,trans.const,dat){
  tip.names<-fit$call$tree$tip.label
  trait.names<-colnames(fit$call$trait.data)
  n<-length(tip.names);k<-length(trait.names);iter<-fit$sampler.control$iter;chains<-fit$sampler.control$chains
  X.param.names<-paste(rep(trait.names,n),'_',rep(tip.names,each=k),sep='')
  new.chains<-array(NA,c(iter,n*k+1,chains),
                    list(iterations=NULL,parameters=c(X.param.names,'R_0'),chains=paste('chains',1:chains)))
  if(fit$call$intra.var){
    out.X<-extract(stanfit,"X",permute=F,inc_warmup=T)*
      rep(sqrt(trans.const$X_sig2),each=iter*chains)+rep(trans.const$X_0,each=iter*chains)
    new.chains[,X.param.names,]<-aperm(out.X,c(1,3,2))
  }else if(length(dat$which_mis)>0){
    mis.Y<-extract(stanfit,"mis_Y",permute=F,inc_warmup=T)*
      rep(sqrt(trans.const$X_sig2),iter*chains*dat$k_mis)+rep(trans.const$X_0,iter*chains*dat$k_mis)
    tmp<-(dat$which_mis-1)*k+rep(1:k,dat$k_mis)
    new.chains[,X.param.names[tmp],]<-aperm(mis.Y,c(1,3,2))
  }
  new.chains[,'R_0',]<-extract(stanfit,"R0",permute=F,inc_warmup=T)
  incl.inds<-which(apply(new.chains[,,1],2,function(ii) !all(is.na(ii))))
  fit$chains<-.index.element(new.chains,incl.inds,2)
  fit
}
#need to have someway to specify graphing of sampler.params with trace.plot/prof.plot?