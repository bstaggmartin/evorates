.match.type<-function(type,choices=c('chains','quantiles','means','diagnostics')){
  type<-unlist(type[1],use.names=FALSE)
  if(is.numeric(type)){
    type<-choices[type][1]
  }
  matches<-pmatch(type,choices)
  if(is.na(matches)){
    stop("Parameter type should be one of ",paste0(choices,collapse=", "))
  }
  choices[matches]
}

#this can definitely be made simpler!
.check.edge.indices<-function(select,in.select,Nedges=NULL){
  if(length(select)!=0){
    old.select<-select
    if(!is.numeric(select)){
      if(is.character(select)){
        select<-gsub('R_','',select)
      }
      select<-suppressWarnings(as.numeric(select))
    }
    probs<-is.na(select)
    if(all(probs)){
      stop('format of ',in.select,' not recognized: the best practice is to use an integer vector to select edge indices')
    }
    if(any(probs)){
      tmp<-sum(probs)
      warning(paste0(old.select[probs][-tmp],collapse=', '),
              if(tmp==1) '' else ' and ',
              paste0(old.select[probs][tmp]),
              if(tmp==1) ' was ' else ' were ',
              'not recognized as ',
              if(tmp==1) ' a valid edge index ' else 'valid edge indices ',
              'and ignored: the best practice is to use an integer vector to select edge indices')
      select<-select[!probs]
    }
    all.neg<-all(select<=0)
    all.pos<-all(select>=0)
    if(all.neg|all.pos){
      if(all.pos){
        exceeds<-select>Nedges
        if(all(exceeds)){
          stop('all specified edge indices out of bounds (i.e., above the number of edges in tree)')
        }else if(any(exceeds)){
          select<-select[!exceeds]
          warning('some specified edge indices out of bounds (i.e., above the number of edges in tree): these indices were ignored')
        }
      }
      if(all.neg){
        select<-(1:Nedges)[select]
      }
    }else{
      stop('mixes of negative and postive integers to specify edge indices are unallowed')
    }
  }
  select[select!=0]
}

.make.edge.groups.list<-function(node.groups,edge.groups,tree,partial.match){
  if(is.null(node.groups)&is.null(edge.groups)){
    edge.groups<-list(1:Nedge(tree))
  }
  if(!is.list(node.groups)&length(node.groups)>0){
    node.groups<-list(node.groups)
  }
  node.groups<-lapply(node.groups,get.clade.edges,
                      tree=tree,
                      partial.match=partial.match)
  if(!is.list(edge.groups)&length(edge.groups)>0){
    edge.groups<-list(edge.groups)
  }
  edge.groups<-lapply(edge.groups,
                      function(ii) tryCatch(.check.edge.indices(ii,deparse(ii),Nedge(tree)),
                                            error=function(e) {warning(e);NULL}))
  edge.groups<-lapply(edge.groups,unique) #probably add a warning here...
  out<-c(node.groups,edge.groups)
  if(length(out)==0){
    stop("no valid groups of edge indices specified")
  }
  out
}