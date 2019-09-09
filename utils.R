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
##I find this is helpful in a number of visualization situations
alter.cols<-function(x,alph=0.5,mod.val=0,name=T){
  if(name){
    cols<-col2rgb(x,alpha=T)/255
  }else{
    cols<-col2rgb(as.numeric(sort(x)),alpha=T)/255
  }
  cols[1:3,]<-cols[1:3,]+mod.val
  cols[1:3,]<-ifelse(cols[1:3,]<0,0,ifelse(cols[1:3,]>1,1,cols[1:3,]))
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
color.clades<-function(tree,MRCAs,base.col=palette()[1],cols=palette()[-1],alph=1){
  
}

#Same function as above, but for lineages rather than clades
color.lineages<-function(tree,lins,base.col=palette()[1],cols=palette()[-1],alph=1){
  
}

#Function to produce an edge-wise color map for a given phylogenetic tree that exhibits phylogenetic signal
##Works by creating a 1 to 3-dimensional color space and simulating BM evolution within that space; to make sure lineages
##smoothly grade into one another, internal node values are estimated using ancestral state reconstruction, and edge-wise colors
##are averaged across adjacent nodes
evolve.colors<-function(tree,col.space.res=100,col.space.d=2,rate=0.1,root=T,
                        hlim=c(0,1),s=1,v=1,slim=c(0.4,1),vlim=c(0.4,1),circular.h=ifelse(all(hlim==c(0,1)),T,F),plot=F,...){
  if(col.space.d!=1){
    sv.vals<-switch(col.space.d-1,
                    cbind(seq(slim[1],slim[2],length.out=col.space.res),seq(vlim[1],vlim[2],length.out=col.space.res)),
                    expand.grid(seq(slim[1],slim[2],length.out=col.space.res),seq(vlim[1],vlim[2],length.out=col.space.res)))
    col.space<-array(rainbow(col.space.res,rep(sv.vals[,1],each=col.space.res),rep(sv.vals[,2],each=col.space.res),
                             hlim[1],hlim[2]),
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

#Function to modify colors of edges not in specified clades (can be inverted)
##Helpful for masking everything except for certain focal clades
alter.notclade<-function(tree,colmap,MRCAs,alph=0,mod.val=0,invert=F){
  
}

#Function to modify colors of edges not in specified lineages (can be inverted)
##Helpful for masking everything except for certain focal lineages
alter.notlin<-function(tree,colmap,lins,alph=0,mod.val=0,invert=F){
  
}

#Function to "jitter" the RGB values of colors corresponding to specified clades
##Helpful if you want to be able to differentiate intra-clade edges
jitter.colors<-function(tree,colmap,MRCAs,amount=0.1){
  
}