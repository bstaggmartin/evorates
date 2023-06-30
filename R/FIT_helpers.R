#worth asking if it's worth it to create all named entries and sticking an NA in it if it's unspecified
#I don't think so--just do it as final step
.make.prior.list<-function(...){
  in.list<-list(...)
  names.table<-cbind(c('^(R[\\._]?0).*((mu)|(mean))','(R[\\._]?0).*((sig)|(sd))',
                       '^((intra)|(Y))[\\._]?((sig2)|(var)).*((sig)|(sd))','^((intra)|(Y))[\\._]?((sig2)|(var)).*((nu)|(df))',
                       '^R[\\._]?((sig2)|(var))',
                       '^((R[\\._]?mu)|(trend)).*((mu)|(mean))','^((R[\\._]?mu)|(trend)).*((sig)|(sd))'),
                     c('R0_prior_mu','R0_prior_sig','Ysig2_prior_sig','Ysig2_prior_df','Rsig2_prior_sig','Rmu_prior_mu','Rmu_prior_sig'))
  out.list<-list()
  for(i in 1:nrow(names.table)){
    tmp<-grepl(names.table[i,1],names(in.list))
    if(sum(tmp)>1){
      stop(names.table[i,2],
           ' matched with multiple arguments: please double-check input.')
    }
    out.list[[names.table[i,2]]]<-unname(unlist(in.list[tmp]))
  }
  #will have to make this more specific in the future with the generalization to multivariate traits
  #narrow this down to parameters that don't require multiple numbers
  conflicts.dups<-names(which(lengths(out.list)>1))
  if(length(conflicts.dups)>0){
    stop('Two numbers specified for prior parameter(s) ',
         paste(conflicts.dups,collapse=', '),
         ': please double-check input.')
  }
  #then check the multiple numbered entries have the right number of entires/are labelled/etc., etc.
  tmp.list<-out.list[grepl('sig$',names(out.list))]
  if(length(tmp.list)>0){
    conflicts.negs<-names(which(sapply(tmp.list,function(ii) any(ii<=0))))
    if(length(conflicts.negs)>0){
      out.list[conflicts.negs]<-lapply(out.list[conflicts.negs],function(ii) ifelse(ii<=0,NA,ii))
      warning('Negative numbers or 0 specified for prior scale parameter(s) ',
              paste(conflicts.negs,collapse=', '),
              ': all negative numbers set to defaults.',
              immediate.=TRUE)
    }
  }
  #convert Inf df to 0 (df <= 0 interpreted as Inf by Stan program)
  if(!is.null(out.list[["Ysig2_prior_df"]])){
    if(is.infinite(out.list[["Ysig2_prior_df"]])){
      out.list[["Ysig2_prior_df"]]<-0
    }
  }
  out.list
}

.coerce.trait.mat<-function(in.trait.mat){
  trait.mat<-in.trait.mat
  if(length(dim(trait.mat))<=1&is.numeric(trait.mat)){
    input.code<-1
    trait.mat<-as.matrix(trait.mat)
  }else if(is.data.frame(trait.mat)){
    input.code<-2
    #prefer tip label column, then rownames
    labels.col<-sapply(trait.mat,function(ii) is.character(ii)|is.factor(ii))
    if(sum(labels.col)>1){
      stop('Multiple columns of label-like data in ',
           deparse(substitute(in.trait.mat)),
           ': please provide only one column of character/factor data, which will be interpreted as tip labels.')
    }
    if(sum(labels.col)>0){
      new.trait.mat<-as.matrix(trait.mat[,!labels.col])
      rownames(new.trait.mat)<-as.character(trait.mat[,labels.col])
      trait.mat<-new.trait.mat
    }
  }else if(length(dim(trait.mat))==2&is.numeric(trait.mat)){
    input.code<-3
  }else{
    stop('Format of ',
         deparse(substitute(in.trait.mat)),
         'not recognized: please provide either a named vector of trait info, a rownamed matrix of trait info, or a data.frame with two columns consisting of tip labels and associated trait info.')
  }
  if(ncol(trait.mat)>1){
    stop('Multiple columns of numeric data in ',
         deparse(substitute(in.trait.mat)),
         ': please provide only one column of numeric data, which will be interpreted as trait information. EvoRate models do not yet support multivariate data.')
  }
  names.test<-rownames(trait.mat)
  if(is.null(names.test)){
    stop('No labels found in ',
         deparse(substitute(in.trait.mat)),
         ': please add a ',
         c('names attribute','column','rownames attribute')[input.code],
         ' providing tip labels to be matched with the tip.label element of the tree object.')
  }
  if(is.null(colnames(trait.mat))){
    colnames(trait.mat)<-'X1'
  }
  trait.mat
}

.coerce.tree<-function(tree,trait.data=NULL,trait.se=NULL){
  node.names<-as.character(1:tree$Nnode+length(tree$tip.label))
  node.matches<-unlist(sapply(paste('^',rownames(trait.data),'$',sep=''),
                              function(ii) grep(ii,node.names,value=TRUE)))
  if(length(node.matches)>0){
    tree<-multi.bind.tip(tree,paste0("t",node.matches),0,as.numeric(node.matches),0)
    inds<-rownames(trait.data)%in%node.matches
    rownames(trait.data)[inds]<-paste0("t",rownames(trait.data)[inds])
    inds<-rownames(trait.se)%in%node.matches
    rownames(trait.se)[inds]<-paste0("t",rownames(trait.se)[inds])
  }
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=1e-300)
  dist.mat<-cophenetic(tree)
  diag(dist.mat)<-NA
  if(any(dist.mat==0,na.rm=TRUE)){
    problem.indices<-which(dist.mat==0,arr.ind=TRUE)
    problem.indices<-t(apply(problem.indices,1,sort))
    problem.indices<-problem.indices[!duplicated(problem.indices),,drop=FALSE]
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
    warning('Detected group(s) of tips representing the same exact phylogenetic position(s): ',
            paste(sapply(collapses,paste,collapse=', '),collapse='; '),
            ': these tips were collapsed into single tips and any trait data was relabelled accordingly',
            immediate.=TRUE)
    new.tip.labels<-sapply(collapses,function(ii) paste(unique(ii),collapse='&'))
    for(i in 1:length(collapses)){
      if(!is.null(trait.data)){
        tmp.X<-grepl(paste(paste0('^',collapses[[i]],'$'),collapse='|'),rownames(trait.data))
        if(sum(tmp.X)>0){
          rownames(trait.data)[tmp.X]<-new.tip.labels[i]
        }
      }
      if(!is.null(trait.se)){
        tmp.SE<-grepl(paste(paste0('^',collapses[[i]],'$'),collapse='|'),rownames(trait.se))
        if(sum(tmp.SE)>0){
          rownames(trait.se)[tmp.SE]<-new.tip.labels[i]
        }
      }
    }
    tree$tip.label[label.pos]<-new.tip.labels
    tree<-drop.tip(tree,unlist(sapply(collapses,'[',-1)))
  }
  if(any(duplicated(tree$tip.label))){
    stop('Two tips share the same label: please relabel the tips of the tree such that each tip has a unique label')
  }
  list(tree=tree,trait.data=trait.data,trait.se=trait.se)
}

.coerce.trait.se<-function(trait.se,p_SE,has_intra,mis){
  #NAs = unfixed, estimated standard errors (coded -1)
  trait.se[is.na(trait.se)]<- -1
  #Infs = missing observations (coded -2)
  trait.se[is.infinite(trait.se)]<- -2
  #has_intra and mis for trait.se (keep duplicates for now, don't aggregate with p_SE yet)
  has_intra<-has_intra[rownames(trait.se)]
  mis<-mis[rownames(trait.se)]
  #Check for specification of fixed (or infinite) standard errors on tips with intra-tip observations
  tmp<-which(trait.se!=-1&has_intra)
  conflicts_intra<-unique(rownames(trait.se)[tmp])
  if(length(conflicts_intra)>0){
    trait.se[tmp]<- -1
    warning('Fixed standard errors are specified for tip(s) ',
            paste(conflicts_intra,collapse=', '),
            ', but tip(s) have intra-tip observations associated with them: defaulted to ignoring fixed standard errors. Specify a single observation for any tips you want to associate with a fixed standard error.',
            immediate.=TRUE)
  }
  #Check for specification of infinite standard errors on tips with observations
  tmp<-which(trait.se==-2&!mis)
  conflicts_mis<-unique(rownames(trait.se)[tmp])
  if(length(conflicts_mis)>0){
    trait.se[tmp]<- -1
    warning('Infinite standard errors are specified for tip(s) ',
            paste(conflicts_mis_X,collapse=', '),
            ', but tip(s) have observations associated with them: defaulted to using an unfixed standard error. Do not specify any observations for any tips you want to associate with an infinite standard error (i.e., missing data).',
            immediate.=TRUE)
  }
  #Check for specification of non-infinite (fixed or unfixed) standard errors on tips without observations
  tmp<-which(trait.se!=-2&mis)
  conflicts_mis_X<-unique(rownames(trait.se)[tmp])
  if(length(conflicts_mis_X)>0){
    trait.se[tmp]<- -2
    warning('Fixed standard errors are specified for tip(s) ',
            paste(conflicts_mis_X,collapse=', '),
            ', but tip(s) have no observations associated with them: defaulted to ignoring fixed standard errors. Specify a single observation for any tips you want to associate with a fixed standard error.',
            immediate.=TRUE)
  }
  trait.se<-trait.se[!duplicated(paste(rownames(trait.se),trait.se)),,drop=FALSE]
  #Check for specification of 2 different standard errors for same tip
  conflicts_dups<-rownames(trait.se)[duplicated(rownames(trait.se))]
  if(length(conflicts_dups)>0){
    stop('Multiple fixed standard errors are specified for tip(s)',
         paste(conflicts_intra,collapse=', '),
         ': please specify only one standard error for each tip.')
  }
  p_SE[rownames(trait.se)]<-trait.se
  p_SE
}

#might be better to account for standard error?
.get.trans.const<-function(tree,X,mis){
  trans.const<-list(hgt=NULL,X_sig2=rep(0,ncol(X)))
  trans.const$hgt<-max(node.depth.edgelength(tree))
  tree$edge.length<-tree$edge.length/trans.const$hgt
  for(i in 1:ncol(X)){
    tmp.X<-X[!mis,i]
    if(length(tmp.X)>1){
      trans.const$X_sig2[i]<-var(tmp.X)
    }
    if(trans.const$X_sig2[i]==0){
      warning('No information to scale trait data--every tip is missing and/or has the same mean trait value: defaulted to not scaling trait data but you should probably double-check your input!',
              immediate.=TRUE)
      trans.const$X_sig2[i]<-1
    }
  }
  trans.const
}

.prep.run<-function(input.evorates.obj,return.as.obj=TRUE,out.file=NULL,check.overwrite=TRUE,
                    constrain.Rsig2=NULL,trend=NULL,lik.power=NULL,prior.list=NULL,nchain=4){
  ##HANDLING FILE OUTPUTS##
  if(!return.as.obj|!is.null(out.file)){
    if(is.null(out.file)){
      message("no file output name specified, but model results aren't set to be returned in R working environment:\n picked file output name automatically based off current date and time")
      time<-Sys.time()
      out.file<-paste0('evorates_',gsub(' |-|:','_',time))
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
  dat<-input.evorates.obj$dat
  def.priors<-list('R0_prior_mu'=0,'R0_prior_sig'=10,
                   'Ysig2_prior_sig'=2,'Ysig2_prior_df'=1,
                   'Rsig2_prior_sig'=5,
                   'Rmu_prior_mu'=0,'Rmu_prior_sig'=10)
  #we'll need to tweak this for multivariate priors
  for(i in names(def.priors)){
    if(!is.null(prior.list[[i]])){
      if(!is.na(prior.list[[i]])){
        dat[[i]]<-prior.list[[i]]
      }
    }
    if(is.null(dat[[i]])){
      dat[[i]]<-def.priors[[i]]
    }else{
      if(is.na(dat[[i]])){
        dat[[i]]<-def.priors[[i]]
      }
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
              immediate.=TRUE)
      lik.power<-1
    }
    dat$lik_power<-lik.power
  }
  input.evorates.obj$dat<-dat
  
  ##CHOOSE MODEL/EXCLUDED PARAMETERS##
  #will become more complex/important with multivariate extension
  exclude.pars<-c('std_R0','tau_Ysig2','std_Ysig2','std_Ysig2_cauchy','std_Rsig2','std_Rmu','raw_R','SE','tmp_lik','tmp_mu')
  if(dat$n_mis_SE==0){
    exclude.pars<-c(exclude.pars,'Ysig2')
  }
  if(dat$constr_Rsig2){
    exclude.pars<-c(exclude.pars,'Rsig2')
  }
  if(dat$constr_Rmu){
    exclude.pars<-c(exclude.pars,'Rmu')
  }
  if(dat$constr_Rsig2&dat$constr_Rmu){
    exclude.pars<-c(exclude.pars,'R')
  }
  stanobj<-'univar_evorates_normpri_tpri'
  
  ##FORM OUTPUT FOR RUNNING MODEL##
  out<-list(return.as.obj=return.as.obj,stanobj=stanobj,exclude.pars=exclude.pars,out.file=out.file,
            input.evorates.obj=input.evorates.obj)
  if(!is.null(out.file)){
    out<-c(out,file.strings=list(file.strings),directory=directory)
  }
  out
}
