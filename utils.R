#Draw from a bactrian distribution (see Yang et al. 2013 for details)
##works by drawing from a binomial distribution, then using that draw as an indicator to draw from either a "left" or "right" normal
##distribution
rbac<-function(n,x=0,sd=1,m=0.9){
  mu<-ifelse(rbinom(n,1,0.5)>0.5,x-m*sd,x+m*sd)
  rnorm(n,mu,sqrt(1-m^2)*sd)
}

#Extract edge-wise rates from the output of a relaxed.clock.BM mcmc
##helpful since the autocorrelated rate model output gives (log) node-wise rate parameters, rather than the edge-wise rates (which
##can be approximated by averaging across adjacent nodes)
get.edge.rates<-function(tree,par.mat){
  exp(sapply(1:ncol(par.mat),function(ii) {
    apply(cbind(par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,1]],par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,2]]),1,mean)
  }))
}

#Function to sample the indices of vector as "evenly" as possible by halving the vector into finer and finer increments
##I'm hoping this will make simulated BM bridges more noisy, a la true BM, but I still need to test this; note that this
##function is, unfortunately, much slower than sample, but I have yet to think of a more efficient way to do this...
sample.evenly<-function(vec){
  samp<-c(0,round(length(vec)/2),length(vec)+1)
  while(!all(1:length(vec)%in%samp)){
    space<-mean(sort(samp,decreasing=T)[-length(samp)]-sort(samp,decreasing=T)[-1])
    samp<-unique(c(samp,round(sort(samp)[-length(samp)]+space/2)))
  }
  vec[samp[-c(1,3)]]
}
##Update: did not seem to make the bridges more noisy...but keeping the function anyways; I commented out a line of code in the
##'SimBM' function so I can easily switch between random sampling and this type of sampling when simulating BM bridges

#Function to make a modify a vector of colors to have a certain level of transparency (or be darker/lighter)
##I find this is helpful in a number of visualization situations; use alph=NA to not modify original transparencies of colors
alter.cols<-function(x,alph=NA,mod.val=0){
  if(is.character(x)){
    cols<-col2rgb(x,alpha=T)/255
  }else{
    cols<-col2rgb(as.numeric(x),alpha=T)/255
  }
  cols[1:3,]<-cols[1:3,]+mod.val
  cols[1:3,]<-ifelse(cols[1:3,]<0,0,ifelse(cols[1:3,]>1,1,cols[1:3,]))
  alph<-ifelse(is.na(alph),cols[4,],alph)
  cols[4,]<-alph
  return(mapply(rgb,red=cols[1,],green=cols[2,],blue=cols[3,],alpha=cols[4,]))
}

#Function to produce an edge-wise color map for a given phylogenetic tree where monophyletic clades are assigned to a given set of
#colors based on MRCAs (note: if colors are assigned to nested clades, the color for the smaller clade will override that for the
#larger clade within the smaller clade, technically allowing paraphyletic clades to be assigned unqiue colors)
##you already have a working version of this integrated into the 'plot.simBM' method, but I think visualizing simulations in this
##package will go a lot smoother if you make coloring phylogeny edges a more modular process and make this a separate function;
##combined with evolve.colors below, as well as planned 'mask.clades/lineages' and 'jitter.colors' functions, users will have
##flexible, intuitive tools to coloring their simulations the way they see fit
##9/10 update: modified coloration system--moved this function from 'simBM' to here, made alter.cols handle more aspects of color
##management internally
##9/11 thoughts: modify the way this function (and its lineage counterpart) handles base colors such that it can take already
##made colmaps and recolor specified clades and lineages
color.clades<-function(tree,MRCAs,cols=palette()[2:(length(MRCAs)+1)],alph=NA,base.col=palette()[1],base.alph=NA){
  if(length(cols)<length(MRCAs)){
    cols<-rep(cols,length.out=length(MRCAs))
  }
  if(length(alph)<length(MRCAs)){
    alph<-rep(alph,length.out=length(MRCAs))
  }
  cols<-c(alter.cols(base.col,alph=base.alph),alter.cols(cols,alph=alph))
  clade.edges<-lapply(c(length(tree$tip.label)+1,MRCAs),function(nn) which(tree$edge[,2]%in%getDescendants(tree,nn)))
  inc.ord<-order(lengths(clade.edges),decreasing=T)
  clade.edges<-clade.edges[inc.ord]
  for(c in 1:length(clade.edges)){
    clade.edges[[c]]<-clade.edges[[c]][!(clade.edges[[c]]%in%unlist(clade.edges[-c]))]
  }
  clade.edges[inc.ord]<-clade.edges
  clade.score<-rep(1:length(clade.edges),lengths(clade.edges))
  edge.map<-cbind(clade.score,unlist(clade.edges))
  edge.map<-edge.map[order(edge.map[,2]),]
  colmap<-cols[edge.map[,1]]
  colmap
}

#Same function as above, but for lineages rather than clades
##In the case of overlapping lineages, the color is resolved by averaging rgb values
color.lineages<-function(tree,lins,cols=palette()[2:(length(lins)+1)],alph=NA,base.col=palette()[1],base.alph=NA){
  if(length(cols)<length(lins)){
    cols<-rep(cols,length.out=length(lins))
  }
  if(length(alph)<length(lins)){
    alph<-rep(alph,length.out=length(lins))
  }
  cols<-c(alter.cols(base.col,alph=base.alph),alter.cols(cols,alph=alph))
  if(is.character(lins)){
    lins<-match(lins,tree$tip.label)
  }
  get.lin.edges<-function(tree,lin){
    prev.len<-0
    pp<-which(tree$edge[,2]==lin)
    while((length(pp)-prev.len)>0){
      prev.len<-length(pp)
      pp<-c(pp,which(tree$edge[,2]==tree$edge[pp[length(pp)],1]))
    }
    pp
  }
  lin.edges<-lapply(lins,get.lin.edges,tree=tree)
  lin.scores<-rep(1:length(lin.edges),lengths(lin.edges))
  excl.edges<-(1:nrow(tree$edge))[!(1:nrow(tree$edge)%in%unlist(lin.edges))]
  lin.edges<-c(excl.edges,unlist(lin.edges))
  lin.scores<-c(rep(1,length(excl.edges)),lin.scores+1)
  lin.scores<-lin.scores[order(lin.edges)]
  tmp.cols<-t(col2rgb(cols[lin.scores],alpha=T))/255
  edge.wise.lin.scores<-split(tmp.cols,sort(lin.edges))
  colmap<-rgb(t(sapply(edge.wise.lin.scores,function(ii) apply(matrix(ii,ncol=4),2,mean))))
  colmap
}

#Function to produce an edge-wise color map for a given phylogenetic tree that exhibits phylogenetic signal
##Works by creating a 1 to 3-dimensional color space and simulating BM evolution within that space; to make sure lineages
##smoothly grade into one another, internal node values are estimated using ancestral state reconstruction, and edge-wise colors
##are averaged across adjacent nodes
evolve.colors<-function(tree,col.space.res=100,col.space.d=2,rate=0.1,root=T,
                        hlim=c(0,1),s=1,v=1,slim=c(0.4,1),vlim=c(0.4,1),circular.h=ifelse(all(hlim==c(0,1)),T,F),alph=1,
                        plot=F,...){
  if(col.space.d!=1){
    sv.vals<-switch(col.space.d-1,
                    cbind(seq(slim[1],slim[2],length.out=col.space.res),seq(vlim[1],vlim[2],length.out=col.space.res)),
                    expand.grid(seq(slim[1],slim[2],length.out=col.space.res),seq(vlim[1],vlim[2],length.out=col.space.res)))
    col.space<-array(rainbow(col.space.res,rep(sv.vals[,1],each=col.space.res),rep(sv.vals[,2],each=col.space.res),
                             hlim[1],hlim[2],alpha=alph),
                     dim=rep(col.space.res,col.space.d))
    if(is.vector(rate)){
      rate<-diag(rate,nrow=col.space.d,ncol=col.space.d)
    }
    node.cols<-matrix(rmvn(1,rep(0,col.space.d*length(tree$tip.label)),kronecker(rate,vcv(tree))),ncol=col.space.d)
    rownames(node.cols)<-tree$tip.label
    node.cols<-rbind(node.cols,anc.recon(node.cols,tree))
    edge.cols<-sapply(1:col.space.d,function(ii) apply(matrix(node.cols[as.vector(tree$edge),ii],ncol=2),1,mean))
    if(root){
      edge.cols<-rbind(edge.cols,node.cols[length(tree$tip.label)+1,])
    }
    if(circular.h){
      problem.indices<-which(edge.cols[,1]>hlim[2])
      while(!all(edge.cols[,1]<=hlim[2])){
        edge.cols[problem.indices,1]<-edge.cols[problem.indices,1]-(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols[,1]>hlim[2])
      }
      problem.indices<-which(edge.cols[,1]<hlim[1])
      while(!all(edge.cols[,1]>=hlim[1])){
        edge.cols[problem.indices,1]<-edge.cols[problem.indices,1]+(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols[,1]<hlim[1])
      }
      edge.cols[,1]<-round(edge.cols[,1]*(col.space.res-1)+1)
    }else{
      edge.cols[,1]<-round((edge.cols[,1]-min(edge.cols[,1]))/(max(edge.cols[,1])-min(edge.cols[,1]))*(col.space.res-1)+1)
    }
    for(i in 2:col.space.d){
      edge.cols[,i]<-round((edge.cols[,i]-min(edge.cols[,i]))/(max(edge.cols[,i])-min(edge.cols[,i]))*(col.space.res-1)+1)
    }
    edge.cols<-col.space[edge.cols]
    if(plot){
      plot(tree,edge.color=edge.cols,edge.width=4,...)
    }
    edge.cols
  }else{
    col.space<-rainbow(col.space.res,s,v,hlim[1],hlim[2])
    node.cols<-t(rmvn(1,rep(0,length(tree$tip.label)),rate*vcv(tree)))
    rownames(node.cols)<-tree$tip.label
    node.cols<-c(node.cols,anc.recon(node.cols,tree))
    edge.cols<-apply(matrix(node.cols[as.vector(tree$edge)],ncol=2),1,mean)
    if(root){
      edge.cols<-c(edge.cols,node.cols[length(tree$tip.label)+1])
    }
    if(circular.h){
      problem.indices<-which(edge.cols>hlim[2])
      while(!all(edge.cols<=hlim[2])){
        edge.cols[problem.indices]<-edge.cols[problem.indices]-(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols>hlim[2])
      }
      problem.indices<-which(edge.cols<hlim[1])
      while(!all(edge.cols>=hlim[1])){
        edge.cols[problem.indices]<-edge.cols[problem.indices]+(hlim[2]-hlim[1])
        problem.indices<-which(edge.cols<hlim[1])
      }
      edge.cols<-round(edge.cols*(col.space.res-1)+1)
    }else{
      edge.cols<-round((edge.cols-min(edge.cols))/(max(edge.cols)-min(edge.cols))*(col.space.res-1)+1)
    }
    edge.cols<-col.space[edge.cols]
    if(plot){
      plot(tree,edge.color=edge.cols,edge.width=4,...)
    }
    edge.cols
  }
}

#Function to modify colors of edges in specified clades (can be inverted)
##Helpful for modifying existing colmaps or masking everything except for certain focal clades
alter.clade<-function(tree,colmap,MRCAs,alph=NA,mod.val=0,col=NULL,invert=F){
  
}

#Function to modify colors of edges in specified lineages (can be inverted)
##Helpful for modifying existing colmaps or masking everything except for certain focal lineages
alter.lin<-function(tree,colmap,lins,alph=NA,mod.val=0,col=NULL,invert=F){
  
}
#You might be able to combine these two functions into an alter.colmap function?

#Function to "jitter" the RGB values of colors corresponding to specified clades or lineages
##Helpful if you want to be able to differentiate intra-clade/lineage edges
jitter.colors<-function(tree,colmap,MRCAs=NULL,lins=NULL,amount=0.1,red=T,green=T,blue=T,alpha=F){
  amount<-amount*c(red,green,blue,alpha)
  if(!is.null(lins)){
    get.lin.edges<-function(tree,lin){
      prev.len<-0
      pp<-which(tree$edge[,2]==lin)
      while((length(pp)-prev.len)>0){
        prev.len<-length(pp)
        pp<-c(pp,which(tree$edge[,2]==tree$edge[pp[length(pp)],1]))
      }
      pp
    }
    lin.edges<-lapply(lins,get.lin.edges,tree=tree)
  }else{
    lin.edges<-NULL
  }
  if(!is.null(MRCAs)){
    clade.edges<-lapply(MRCAs,function(nn) which(tree$edge[,2]%in%getDescendants(tree,nn)))
  }else{
    clade.edges<-NULL
  }
  edges<-unique(c(unlist(lin.edges),unlist(clade.edges)))
  new.cols<-t(col2rgb(colmap[edges],alpha=T)/255+rnorm(length(edges)*4,0,amount))
  new.cols[,]<-ifelse(as.vector(new.cols)<0,0,ifelse(as.vector(new.cols)>1,1,as.vector(new.cols)))
  colmap[edges]<-rgb(new.cols)
  colmap
}
