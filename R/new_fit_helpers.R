#works for the most part, except if new node has created new edge and position is specified to be past old
#edge...
#works in those cases now--simple sorting issue
#wrapper for bind.tip to allow adding multiple tips to a phylogeny at once
#' @export
multi.bind.tip<-function(tree,names,edge.lengths=NULL,nodes,positions=0){
  if(is.null(edge.lengths)){
    edge.lengths<-NA
  }
  args.list<-list(names=names,edge.lengths=edge.lengths,nodes=nodes,positions=positions)
  max.len<-max(sapply(args.list,length))
  for(i in 1:length(args.list)){
    args.list[[i]]<-rep(args.list[[i]],length.out=max.len)
  }
  nodes<-args.list$nodes
  positions<-args.list$positions
  names(nodes)<-ifelse(nodes<=length(tree$tip.label),nodes+tree$Nnode,nodes-length(tree$tip.label))
  args.list<-lapply(args.list,function(ii) ii[order(as.numeric(names(nodes)),-positions)])
  nodes<-nodes[order(as.numeric(names(nodes)))]
  int.nodes<-nodes>length(tree$tip.label)
  tips<-nodes<=length(tree$tip.label)
  for(i in 1:max.len){
    try.tree<-try(bind.tip(tree,args.list$names[i],
                           if(is.na(args.list$edge.lengths[i])) NULL else args.list$edge.lengths[i],
                           args.list$nodes[i],
                           args.list$positions[i]),silent=T)
    if(inherits(try.tree,'try-error')){
      warning("failed to bind tip '",args.list$names[i],"' to node ",nodes[i],' due to following error:\n',
              try.tree)
    }else{
      tree<-try.tree
      args.list$nodes[int.nodes]<-args.list$nodes[int.nodes]+if(args.list$positions[i]<=0) 1 else 2
      tmp<-which(tree$tip.label==args.list$names[i])
      tmp<-args.list$nodes>=tmp
      args.list$nodes[tips&tmp]<-args.list$nodes[tips&tmp]+1
    }
  }
  tree
}

#remove non-numeric/all NA rows, coerce to matrix with labelled rows/columns
#allow naming of prior components, and make sure to coerce them to have same number of columns as Y and
#such before exporting them from this function...
.format.trait.data<-function(Y,X0.prior.mu=NULL,X0.prior.sig=NULL,Xsig2.prior=NULL,Ysig2.prior=NULL,tree){
  #if any of the priors have names, attempt to reorder according to trait.data
  for(i in c('X0.prior.mu','X0.prior.sig','Xsig2.prior','Ysig2.prior')){
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
          for(i in c('X0.prior.mu','X0.prior.sig','Xsig2.prior','Ysig2.prior')){
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
        for(i in c('X0.prior.mu','X0.prior.sig','Xsig2.prior','Ysig2.prior')){
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
  for(i in c('X0.prior.mu','X0.prior.sig','Xsig2.prior','Ysig2.prior')){
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
    warning("corateBM currently doesn't support tips at root of tree: these tips were removed; please specify information regarding trait valuies at the root using the X0 prior",
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
      tp.sig<-tp.sig[sort(rownames(tp.sig)),]
    }
    if(length(missing.mu)>0){
      tp.mu<-do.call(rbind,c(list(tp.mu),rep(list(NA),length(missing.mu))))
      rownames(tp.mu)[-(length(missing.mu)-1):0+nrow(tp.mu)]<-missing.mu
      tp.mu<-tp.mu[sort(rownames(tp.mu)),]
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
                           X0.prior.mu=NULL,X0.prior.sig=NULL,Xsig2.prior=NULL,
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
  for(i in c('X0.prior.sig','Xsig2.prior','Ysig2.prior')){
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
       X0.prior.mu=X0.prior.mu,X0.prior.sig=X0.prior.sig,Xsig2.prior=Xsig2.prior,
       R0.prior.mu=R0.prior.mu,Rsig2.prior=Rsig2.prior,
       Ysig2.prior=Ysig2.prior,
       Rmu.prior.mu=Rmu.prior.mu,Rmu.prior.sig=Rmu.prior.sig)
}