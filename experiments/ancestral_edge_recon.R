library(phytools)
tree<-pbtree(n=50)
x<-fastBM(tree,internal=T)
x_trim<-x[1:length(tree$tip.label)]
library(contSimmap)
fastpgram(x_trim,tree)
library(Rphylopars)
recon<-anc.recon(x_trim,tree,vars=T)
recon$Yhat<-c(x_trim,as.vector(recon$Yhat));recon$vars<-c(rep(0,length(x_trim)),as.vector(recon$vars))
est.sig2<-mean(pic(x_trim,tree)^2)
edge<-list(expect=rep(NA,nrow(tree$edge)),var=rep(NA,nrow(tree$edge)))
for(i in 1:length(edge$expect)){
  edge$expect[i]<-mean(recon$Yhat[tree$edge[i,]])
  edge$var[i]<-sum(recon$var[tree$edge[i,]]+min(recon$var[tree$edge[i,]]))/4+tree$edge.length[i]*est.sig2/12
  #min(recon$var[tree$edge[i,]]) for covariance
}
tmp.x<-node.depth.edgelength(tree)[tree$edge[,1]]+tree$edge.length/2
points(edge$expect~tmp.x,pch=16)
segments(y0=edge$expect-sqrt(edge$var)*1.96,y1=edge$expect+sqrt(edge$var)*1.96,x0=tmp.x)
segments(y0=recon$Yhat-sqrt(recon$vars)*1.96,y1=recon$Yhat+sqrt(recon$vars)*1.96,x0=node.depth.edgelength(tree))
#randomization
#more complicated than I initially anticipated--sister edges must be jointly drawn from a multivariate normal since they have a non-zero
#covariance relative to the average of the ancestral edge
rand.edge<-edge
#eliminate edge variance associated with variance(+covariance) of the ancestral node value
rand.edge$var[-1]<-rand.edge$var[-1]-recon$var[tree$edge[-1,1]]/4-apply(matrix(recon$var[tree$edge[-1,]],ncol=2),1,min)/2
#add variance associated with node now that edge mean has been fixed
tmp<-sapply(2:nrow(tree$edge),function(ii) which(tree$edge[,2]==tree$edge[ii,1]));tmp[lengths(tmp)==0]<-1;tmp<-unlist(tmp)
tmp2<-node.depth.edgelength(tree)[tree$edge[-1,1]]
tmp3<-node.depth.edgelength(tree)[tree$edge[-1,2]]

#here's the formula that has to be fixed...
tmp.var<-est.sig2*(tree$edge.length[tmp]/2)*(tmp3-tmp2)/(tmp3-tmp2+tree$edge.length[tmp]/2)+rand.edge$var[-1]


rand.edge$var[-1]<-rand.edge$var[-1]+tmp.var/4+
  apply(matrix(c(recon$var[tree$edge[-1,2]],tmp.var),ncol=2),1,min)/2
#this doesn't seem to work actually; problem is that the node values themselves have artificially low variances since their descendant
#node values are known--you may have to subtract something to do with the covariance of ancestral-descendant node values...
rand.edge$expect<-rep(NA,nrow(tree$edge))
library(mvnfast)
for(i in 1:nrow(tree$edge)){
  if(i==1){
    rand.edge$expect[i]<-rnorm(1,edge$expect[i],sqrt(rand.edge$var[i]))
    next
  }
  if(!is.na(rand.edge$expect[i])){
    next
  }
  a.n<-tree$edge[i,1]
  d.n<-tree$edge[i,2]
  if(any(tree$edge[,2]==n)){
    p<-which(tree$edge[,2]==a.n)
    if(length(p)==0){
      p<-1
    }
    d<-which(tree$edge[,1]==d.n)
    if(length(d)==0){
      rand.edge$expect[i]<-edge$expect[p]*2/(tree$edge.length[p]+tree$edge.length[i])/
        (2/(tree$edge.length[p]+tree$edge.length[i])+2/tree$edge.length[i])+
        recon$Yhat[d.n]*2/tree$edge.length[i]/
        (2/(tree$edge.length[p]+tree$edge.length[i])+2/tree$edge.length[i])
    }else{
      rand.edge$expect[i]<-sum(edge$expect[c(p,d)]*2/(tree$edge.length[c(p,d)]+tree$edge.length[i])/
                                 sum(2/(tree$edge.length[c(p,d)]+tree$edge.length[i])))
    }
    s<-which(tree$edge[,1]==a.n);s<-s[s!=i]
    if(s==1){
      s<-integer(0)
    }
    if(length(s)>0){
      s.d.n<-tree$edge[s,2]
      s.d<-which(tree$edge[,1]==s.d.n)
      if(length(s.d)==0){
        rand.edge$expect[s]<-edge$expect[p]*2/(tree$edge.length[p]+tree$edge.length[s])/
          (2/(tree$edge.length[p]+tree$edge.length[s])+2/tree$edge.length[s])+
          recon$Yhat[s.d.n]*2/tree$edge.length[s]/
          (2/(tree$edge.length[p]+tree$edge.length[s])+2/tree$edge.length[s])
      }else{
        rand.edge$expect[s]<-sum(edge$expect[c(p,s.d)]*2/(tree$edge.length[c(p,s.d)]+tree$edge.length[s])/
                                   sum(2/(tree$edge.length[c(p,s.d)]+tree$edge.length[s])))
      }
      rand.edge$expect[c(i,s)]<-rmvn(1,rand.edge$expect[c(i,s)],matrix(c(rand.edge$var[i],0,0,
                                                                         rand.edge$var[s]),ncol=2))
    }else{
      rand.edge$expect[i]<-rnorm(1,edge$expect[i],sqrt(rand.edge$var[i]))
    }
  }
  
}
#plot(rand.edge$expect~apply(matrix(node.depth.edgelength(tree)[as.vector(tree$edge)],ncol=2),1,mean))
edge.t.ranges<-matrix(node.depth.edgelength(tree)[as.vector(tree$edge)],ncol=2)
interp.nod.vals<-rep(NA,max(tree$edge))
for(n in 1:length(interp.nod.vals)){
  if(n<=length(tree$tip.label)){
    interp.nod.vals[n]<-recon$Yhat[n]
    next
  }
  involved.edges<-c(which(tree$edge[,2]==n),which(tree$edge[,1]==n))
  weights<-(2/tree$edge.length[involved.edges])/sum(2/tree$edge.length[involved.edges])
  interp.nod.vals[n]<-sum(weights*rand.edge$expect[involved.edges])
}
plot(interp.nod.vals~node.depth.edgelength(tree),col='white')
segments(x0=c(edge.t.ranges[,1],edge.t.ranges[,1]+tree$edge.length/2),x1=c(edge.t.ranges[,1]+tree$edge.length/2,edge.t.ranges[,2]),
         y0=c(interp.nod.vals[tree$edge[,1]],rand.edge$expect),y1=c(rand.edge$expect,interp.nod.vals[tree$edge[,2]]))
#this is unwieldly and unnecessarily complicated--I would just do randomization at nodes followed by random interpolations...
#part where it breaks down is the covariance between sister taxa --> your calculations of edge means variances conditioned on previous
#edge mean being fixed isn't correct yet; the fromula yields distinct variance for the node in question depending on whether you're
#focusing on one sister branch or the other. The node variance, however, must be conditioned on BOTH for a correct calculation...
#yeah, I think the above is technically possible, but unnecessarily complicated...

#okay, honestly--it's a reasonably good approximation! Still, covariances aren't calculated exactly correctly, and it's a bit noisier
#than 'natural'  BM...

#let's try the 'node-first' approach
rand.recon<-recon
root<-length(tree$tip.label)+1
mu<-c(x_trim,rep(NA,tree$Nnode))
p<-rep(NA,max(tree$edge))
for(ee in nrow(tree$edge):0){
  if(ee==0){
    nn<-root
    dd<-tree$edge[which(tree$edge[,1]==nn),2]
    p.a<-sum(p[dd],na.rm=T)
    mu[nn]<-sum(mu[dd]*p[dd]/p.a,na.rm=T)
    t<-0
    p[nn]<-p.a/(1+t*p.a)
    break
  }
  nn<-tree$edge[ee,2]
  if(length(which(tree$edge[,1]==nn))==0){
    if(is.na(mu[nn])){
      next
    }else{
      t<-tree$edge.length[ee]
      p[nn]<-1/t
    }
  }
  if(length(which(tree$edge[,1]==nn))>0){
    dd<-tree$edge[which(tree$edge[,1]==nn),2]
    p.a<-sum(p[dd],na.rm=T)
    mu[nn]<-sum(mu[dd]*p[dd]/p.a,na.rm=T)
    t<-tree$edge.length[ee]
    p[nn]<-p.a/(1+t*p.a)
  }
}
p[1:length(tree$tip.label)]<-0
a.nodes<-sapply((1:max(tree$edge))[-(0:length(tree$tip.label)+1)],function(ii) tree$edge[which(tree$edge[,2]==ii),1])
a.edges<-sapply((1:max(tree$edge))[-(0:length(tree$tip.label)+1)],function(ii) which(tree$edge[,2]==ii))
rand.recon$vars[-(0:length(tree$tip.label)+1)]<-est.sig2*1/(est.sig2/recon$vars[-(0:length(tree$tip.label)+1)]-
                                                              (est.sig2/recon$vars[a.nodes]-p[-(0:length(tree$tip.label)+1)])/
                                                              (1+tree$edge.length[a.edges]*
                                                                 (est.sig2/recon$vars[a.nodes]-p[-(0:length(tree$tip.label)+1)]))+
                                                              1/tree$edge.length[a.edges])
nsim<-100
rand.recon$Yhat<-matrix(rep(mu,nsim),ncol=nsim)
rand.recon$Yhat[root,]<-rnorm(nsim,recon$Yhat[root],sqrt(rand.recon$vars[root]))
for(ee in 1:nrow(tree$edge)){
  dn<-tree$edge[ee,2];an<-tree$edge[ee,1]
  if(dn<root){
    next
  }
  rand.recon$Yhat[dn,]<-rand.recon$Yhat[dn,]*p[dn]*tree$edge.length[ee]+rand.recon$Yhat[an,]-rand.recon$Yhat[an,]*p[dn]*tree$edge.length[ee]
  rand.recon$Yhat[dn,]<-rnorm(nsim,rand.recon$Yhat[dn,],sqrt(rand.recon$vars[dn]))
}
fastpgram(rand.recon$Yhat[,sample(nsim,1)],tree)
fastpgram(rand.recon$Yhat[,sample(nsim,1)],tree,add=T)
rand.edge<-edge
rand.edge$var<-est.sig2*tree$edge.length/12
rand.edge$expect<-apply(rand.recon$Yhat,2,function(ii) sapply(1:nrow(tree$edge),function(jj) mean(ii[tree$edge[jj,]])))
rand.edge$expect<-rand.edge$expect+rnorm(nsim*nrow(tree$edge),0,sqrt(rand.edge$var))
interp.nod.vals<-matrix(NA,nrow=max(tree$edge),ncol=nsim)
for(n in 1:nrow(interp.nod.vals)){
  if(n<=length(tree$tip.label)){
    interp.nod.vals[n,]<-recon$Yhat[n]
    next
  }
  involved.edges<-c(which(tree$edge[,2]==n),which(tree$edge[,1]==n))
  weights<-(2/tree$edge.length[involved.edges])/sum(2/tree$edge.length[involved.edges])
  interp.nod.vals[n,]<-apply(rand.edge$expect[involved.edges,],2,function(ii) sum(weights*ii))
}

plot(0,ylim=range(interp.nod.vals),xlim=range(node.depth.edgelength(tree)),col='white')
samp<-sample(nsim,1)
segments(x0=c(edge.t.ranges[,1],edge.t.ranges[,1]+tree$edge.length/2),x1=c(edge.t.ranges[,1]+tree$edge.length/2,edge.t.ranges[,2]),
         y0=c(interp.nod.vals[tree$edge[,1],samp],rand.edge$expect[,samp]),y1=c(rand.edge$expect[,samp],interp.nod.vals[tree$edge[,2],samp]))
#works well--deriving a direct--ML edge recon to randomized edge recon algorithm isn't necessary for now, maybe helpful for speeding
#up future calculations...? But more complicated than I initially thought.

colramp<-colorRampPalette(c('blue','cyan','gray','orange','red'))(100)
cols<-matrix(
  colramp[round((rand.edge$expect-min(rand.edge$expect))/(max(rand.edge$expect)-min(rand.edge$expect))*(length(colramp)-1)+1)],
  ncol=nsim)
plot(tree,edge.color=cols[,sample(nsim,1)],edge.width=4)

#variance seems slightly off--a bit underestimated. Odd, but probably a simple arithmetic error somewhere earlier...
#tested with original simBM function (corrected to sqrt the variance) and the function works better
plot(recon$vars~apply(rand.recon$Yhat,1,var))
plot(edge$var~apply(rand.edge$expect,1,var))
#got it, simple algebra error (to go from MLE var to p --> sigma/MLE var, NOT 1/MLE var/sigma)