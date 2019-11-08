#' @export
thresh<-function(tree,sim,thresholds,n.sample=ncol(sim$nodes)){
  
  ##BASIC ERROR CHECKS AND FIXES##
  if(!inherits(tree,'phylo')){
    stop('tree should be object of class \"phylo\"')
  }
  if(!inherits(sim,'simBM')){
    stop('sim should be object of class \"simBM\"')
  }
  n.sample<-round(n.sample)
  if(n.sample<1){
    stop('n.sample is less than 1')
  }
  ####
  
  samp<-sample(1:ncol(sim$nodes),n.sample)
  
  ##CONVERTING NODE MATRIX TO EDGE MATRIX##
  time.vec<-sort(unique(round(node.depth.edgelength(tree),digits=10)))
  new.edges<-array(NA,dim=c(nrow(tree$edge),length(time.vec),length(samp)))
  n1<-tree$edge[,1];n2<-tree$edge[,2]
  t1<-round(node.depth.edgelength(tree)[n1],10);t2<-round(node.depth.edgelength(tree)[n2],10)
  new.edges[cbind(rep(rep(1:nrow(tree$edge),2),length(samp)),
                  rep(sapply(c(t1,t2),function(x) which(x==time.vec)),length(samp)),
                  rep(1:length(samp),each=nrow(tree$edge)*2))]<-
    sim$nodes[cbind(rep(c(n1,n2),length(samp)),
                    rep(samp,each=nrow(tree$edge)*2))]
  
  ##COLLATING NEW EDGE MATRIX WITH EXISTING ONE IN SIM OBJECT, IF IT EXISTS##
  if(!is.null(sim$edges)){
    new.edges<-abind::abind(new.edges,sim$edges[,,samp],along=2)
    time.vec<-c(time.vec,sim$ts)
    ord<-order(time.vec)
    new.edges<-new.edges[,ord,]
    time.vec<-time.vec[ord]
    dups<-which(duplicated(time.vec))
    new.edges[,dups-1,]<-ifelse(is.na(new.edges[,dups,]),new.edges[,dups-1,],new.edges[,dups,])
    new.edges<-new.edges[,-dups,]
    time.vec<-time.vec[-dups]
  }
  
  
  thresholds<-sort(thresholds)
  threshed<-sapply(thresholds,'<=',new.edges)
  states<-new.edges
  states[,,]<-rowSums(threshed)
  find.inter<-function(x,y,y.crit){
    return((y.crit-y[1])*(x[2]-x[1])/(y[2]-y[1])+x[1])
  }
  
  master.list<-vector(mode="list",length=length(samp))
  class(master.list)<-c("multiSimmap","multiPhylo")
  map.base<-list(edge=tree$edge,edge.length=tree$edge.length,Nnode=tree$Nnode,tip.label=tree$tip.label,
                 maps=vector(mode="list",length=nrow(tree$edge)),
                 mapped.edge=matrix(NA,nrow=nrow(tree$edge),ncol=length(thresholds)+1))
  rownames(map.base[["mapped.edge"]])<-paste(n1,n2,sep=",")
  colnames(map.base[["mapped.edge"]])<-0:length(thresholds)
  class(map.base)<-c("simmap","phylo")
  for(i in 1:length(samp)){
    master.list[[i]]<-map.base
    for(e in 1:nrow(tree$edge)){
      tmp.states<-states[e,,i]
      tmp.time.vec<-time.vec[which(!is.na(tmp.states))]
      tmp.new.edges<-new.edges[e,which(!is.na(tmp.states)),i]
      tmp.states<-tmp.states[!is.na(tmp.states)]
      crit.pts<-which(tmp.states[1:(length(tmp.states)-1)]!=tmp.states[2:length(tmp.states)])
      start.t<-round(node.depth.edgelength(tree)[tree$edge[e,1]],10)
      end.t<-round(node.depth.edgelength(tree)[tree$edge[e,2]],10)
      if(length(crit.pts)==0){
        tmp.vec<-end.t-start.t
        names(tmp.vec)<-tmp.states[which(tmp.time.vec==start.t)]
      }else{
        crossed.thresh<-mapply(seq,from=tmp.states[crit.pts],to=tmp.states[crit.pts+1],SIMPLIFY=F)
        crossed.thresh<-unlist(crossed.thresh)[-(sapply(crossed.thresh,which.min)+
                                                   c(0,cumsum(lengths(crossed.thresh)[-length(crossed.thresh)])))]
        new.crit.pts<-rep(crit.pts,abs(tmp.states[crit.pts]-tmp.states[crit.pts+1]))
        tmp.vec<-mapply(function(pt,thresh) find.inter(x=tmp.time.vec[c(pt,pt+1)],y=tmp.new.edges[c(pt,pt+1)],y.crit=thresh),
                        pt=new.crit.pts,thresh=thresholds[crossed.thresh])
        tmp.vec<-tmp.vec-start.t;tmp.vec<-c(tmp.vec,tree$edge.length[e]);tmp.vec<-c(tmp.vec[1],diff(tmp.vec))
        tmp.names<-mapply(seq,from=tmp.states[crit.pts],to=tmp.states[crit.pts+1],SIMPLIFY=F)
        #quirk with R -> -integer(0) returns integer(0)
        if(length(tmp.names)==1){
          names(tmp.vec)<-unlist(tmp.names)
        }else{
          names(tmp.vec)<-unlist(tmp.names)[-(cumsum(lengths(tmp.names)[-length(tmp.names)]))]
        }
      }
      master.list[[i]]$maps[[e]]<-tmp.vec
      master.list[[i]]$mapped.edge[e,]<-sapply(0:length(thresholds),function(x) sum(tmp.vec[names(tmp.vec)==x]))
    }
  }
  return(master.list)
}