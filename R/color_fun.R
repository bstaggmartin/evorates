#Function to make a modify a vector of colors to have a certain level of transparency (or be darker/lighter)
##I find this is helpful in a number of visualization situations; use alph=NA to not modify original transparencies of colors
#' Modify the transparency and shade of colors
#'
#' This function takes a vector of colors (x) as input and returns a vector of colors of the same length. The function can be used to
#' either make the colors a certain transparency (alph) or make the colors lighter or darker (mod.val). alph and mod.val are recycled
#' or truncated to be the same length as x.
#'
#' @param x Vector of color values, which can either be a character vector of color names and/or hexadecimal codes, or a numeric
#' vector. In the latter case, the vecotr will be coerced to integers taken to correspond to colors currently defined in palette().
#' @param alpha Vector of alpha (transparency) values (numeric), with 0 corresponding to complete transparency and 1 to complete
#' opaqueness. Use NA to not modify transparency of original color vector. alph must consist of NA's and/or values between 0 and 1.
#' @param mod.val Vector of values (numeric) by which to brighten/darken colors, with positive and negative numbers corresponding to
#' brightening and darkening, respectively. In this context, rgb values range from 0 to 1, and any modified values exceeding
#' these boundaries will be rounded up or down, respectively.
#' @return a vector of color values in hexadecimal code (character)
#' @examples
#' colors<-c('red','green','blue')
#' alter.cols(colors,alph=c(0,1,0.5),mod.val=c(0,-0.5,-0.75))
#' #recycling behavior 
#' alter.cols(colors,alph=0.2,mod.val=-0.5)
#' 
#' @export
alter.cols<-function(x,alpha=NA,mod.val=0){
  if(is.character(x)){
    cols<-col2rgb(x,alpha=T)/255
  }else if(all(is.numeric(x)&x>0)){
    cols<-col2rgb(as.integer(x),alpha=T)/255
  }else{
    stop('x is not a recognizable vector of colors')
  }
  if(!all(is.na(alpha)|(alpha<=1&alpha>=0))){
    stop('alpha is outside permitted range (0-1 or NA)')
  }
  cols[1:3,]<-cols[1:3,]+mod.val
  cols[1:3,]<-ifelse(cols[1:3,]<0,0,ifelse(cols[1:3,]>1,1,cols[1:3,]))
  alpha<-ifelse(is.na(alpha),cols[4,],alpha)
  cols[4,]<-alpha
  mapply(rgb,red=cols[1,],green=cols[2,],blue=cols[3,],alpha=cols[4,])
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
##11/5 thoughts: if you make these functions capable of inversion it will render the 'alter functions' obsolete, which will
##simplify the coloring system to be more intuitive
###nah, I still like the alter.lin/clade idea, honestly--inversion doesn't make predictable sense in a scenario where you are
###multiple colors to multiple clades/lineages.
##11/5 update: completed 9/11 thoughts for color.clades, color.lineages, and alter.colmap
##basic system: generate colmaps using evolve.colors(), color.clades(), or color.lineages(), and further modify colmaps using
##color.clades()/color.lineages() (set base.cols equal to the colmap to be modified in this case), alter.colmap, and jitter.colors
#' @export
color.clades<-function(tree,MRCAs,cols=palette()[2:(length(MRCAs)+1)],alph=NA,base.cols=palette()[1],base.alph=NA){
  if(length(cols)<length(MRCAs)){
    cols<-rep(cols,length.out=length(MRCAs))
  }
  if(length(alph)<length(MRCAs)){
    alph<-rep(alph,length.out=length(MRCAs))
  }
  cols<-alter.cols(cols,alph=alph)
  clade.edges<-lapply(MRCAs,function(nn) which(tree$edge[,2]%in%phytools::getDescendants(tree,nn)))
  inc.ord<-order(lengths(clade.edges),decreasing=T)
  clade.edges<-clade.edges[inc.ord]
  for(c in 1:length(clade.edges)){
    clade.edges[[c]]<-clade.edges[[c]][!(clade.edges[[c]]%in%unlist(clade.edges[-c]))]
  }
  clade.edges[inc.ord]<-clade.edges
  clade.score<-rep(1:length(clade.edges),lengths(clade.edges))
  if(length(base.cols)<nrow(tree$edge)){
    base.cols<-rep(base.cols,length.out=nrow(tree$edge))
  }
  if(length(base.alph)<nrow(tree$edge)){
    base.alph<-rep(base.alph,length.out=nrow(tree$edge))
  }
  colmap<-alter.cols(base.cols,base.alph)
  colmap[unlist(clade.edges)]<-cols[clade.score]
  colmap
}

#Same function as above, but for lineages rather than clades
##In the case of overlapping lineages, the color is resolved by averaging rgb values
#' @export
color.lineages<-function(tree,lins,cols=palette()[2:(length(lins)+1)],alph=NA,base.cols=palette()[1],base.alph=NA){
  if(length(cols)<length(lins)){
    cols<-rep(cols,length.out=length(lins))
  }
  if(length(alph)<length(lins)){
    alph<-rep(alph,length.out=length(lins))
  }
  cols<-alter.cols(cols,alph=alph)
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
  tmp.cols<-t(col2rgb(cols[lin.scores],alpha=T))/255
  edge.wise.lin.scores<-split(tmp.cols,unlist(lin.edges))
  lin.scores<-as.numeric(names(edge.wise.lin.scores))
  cols<-rgb(t(sapply(edge.wise.lin.scores,function(ii) apply(matrix(ii,ncol=4),2,mean))))
  if(length(base.cols)<nrow(tree$edge)){
    base.cols<-rep(base.cols,length.out=nrow(tree$edge))
  }
  if(length(base.alph)<nrow(tree$edge)){
    base.alph<-rep(base.alph,length.out=nrow(tree$edge))
  }
  colmap<-alter.cols(base.cols,base.alph)
  colmap[lin.scores]<-cols
  colmap
}

#Function to produce an edge-wise color map for a given phylogenetic tree that exhibits phylogenetic signal
##Works by creating a 1 to 3-dimensional color space and simulating BM evolution within that space; to make sure lineages
##smoothly grade into one another, internal node values are estimated using ancestral state reconstruction, and edge-wise colors
##are averaged across adjacent nodes
#' @export
evolve.colors<-function(tree,col.space.res=100,col.space.d=2,rate=0.1,return.nodes=F,
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
    node.cols<-matrix(mvnfast::rmvn(1,rep(0,col.space.d*length(tree$tip.label)),kronecker(rate,vcv(tree))),ncol=col.space.d)
    rownames(node.cols)<-tree$tip.label
    node.cols<-rbind(node.cols,Rphylopars::anc.recon(node.cols,tree))
    edge.cols<-sapply(1:col.space.d,function(ii) apply(matrix(node.cols[as.vector(tree$edge),ii],ncol=2),1,mean))
    if(return.nodes){
      edge.cols<-rbind(edge.cols,node.cols)
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
    if(return.nodes){
      list(colmap=edge.cols[1:nrow(tree$edge)],node.cols=edge.cols[(nrow(tree$edge)+1):length(edge.cols)])
    }else{
      edge.cols
    }
  }else{
    col.space<-rainbow(col.space.res,s,v,hlim[1],hlim[2])
    node.cols<-t(mvnfast::rmvn(1,rep(0,length(tree$tip.label)),rate*vcv(tree)))
    rownames(node.cols)<-tree$tip.label
    node.cols<-c(node.cols,Rphylopars::anc.recon(node.cols,tree))
    edge.cols<-apply(matrix(node.cols[as.vector(tree$edge)],ncol=2),1,mean)
    if(return.nodes){
      edge.cols<-c(edge.cols,node.cols)
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
    if(return.nodes){
      list(colmap=edge.cols[1:nrow(tree$edge)],node.cols=edge.cols[(nrow(tree$edge)+1):length(edge.cols)])
    }else{
      edge.cols
    }
  }
}

#Function to modify colors of edges in specified lineages/clades (can be inverted)
##Helpful for modifying existing colmaps or masking everything except for certain focal lineages/clades
#' @export
alter.colmap<-function(tree,colmap,MRCAs=NULL,lins=NULL,alph=NA,mod.val=0,col=NULL,invert=F){
  if(!is.null(lins)){
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
  }else{
    lin.edges<-NULL
  }
  if(!is.null(MRCAs)){
    clade.edges<-lapply(MRCAs,function(nn) which(tree$edge[,2]%in%phytools::getDescendants(tree,nn)))
  }else{
    clade.edges<-NULL
  }
  edges<-unique(c(unlist(clade.edges),unlist(lin.edges)))
  if(invert){
    edges<-(1:nrow(tree$edge))[!(1:nrow(tree$edge)%in%edges)]
  }
  if(!is.null(col)){
    cols<-rep(col,length.out=length(edges))
  }else{
    cols<-colmap[edges]
  }
  colmap[edges]<-alter.cols(cols,alph=alph,mod.val=mod.val)
  colmap
}
#You might be able to combine these two functions into an alter.colmap function?
#done!

#Function to "jitter" the RGB values of colors corresponding to specified clades or lineages
##Helpful if you want to be able to differentiate intra-clade/lineage edges
#' @export
jitter.colors<-function(tree,colmap,MRCAs=NULL,lins=NULL,amount=0.1,red=T,green=T,blue=T,alpha=F){
  amount<-amount*c(red,green,blue,alpha)
  if(!is.null(lins)){
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
  }else{
    lin.edges<-NULL
  }
  if(!is.null(MRCAs)){
    clade.edges<-lapply(MRCAs,function(nn) which(tree$edge[,2]%in%phytools::getDescendants(tree,nn)))
  }else{
    clade.edges<-NULL
  }
  edges<-unique(c(unlist(lin.edges),unlist(clade.edges)))
  new.cols<-t(col2rgb(colmap[edges],alpha=T)/255+rnorm(length(edges)*4,0,amount))
  new.cols[,]<-ifelse(as.vector(new.cols)<0,0,ifelse(as.vector(new.cols)>1,1,as.vector(new.cols)))
  colmap[edges]<-rgb(new.cols)
  colmap
}

#' @export
get.node.colmap<-function(colmap,tree){
  edges<-tree$edge
  node.colmap<-matrix(NA,ncol=4,nrow=max(edges))
  for(n in 1:max(edges)){
    adj.edges<-c(which(edges[,1]==n),which(edges[,2]==n))
    adj.edge.cols<-t(col2rgb(colmap[adj.edges],alpha=T))/255
    if(is.vector(adj.edge.cols)){
      node.colmap[n,]<-adj.edge.cols
    }else{
      node.colmap[n,]<-apply(adj.edge.cols,2,mean)
    }
  }
  node.colmap<-mapply(rgb,red=node.colmap[,1],green=node.colmap[,2],blue=node.colmap[,3],alpha=node.colmap[,4])
  node.colmap
}
