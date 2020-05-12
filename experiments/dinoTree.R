DinoTrees<-read.tree('~/../Desktop/Dinosaur trees cal3 Benson 2017.tre')
#cal3 doesn't really work correctly with sim since the inter-node distances are so short as to be interpreted as polytomies? That seems
#to be the case anyways--should explore this more
#figured out the problem--just look at tree$edge.length--some edge lengths are 0!!! the way I see it, there are two answers to this
#problem: 1) make sure the function works with polytomies (probably the best), 2) fudge edge.lengths to make sure they are all of an
#infinitesmally small length (kinda hacky, but might be appropriate in some instances)
#turns out you may have to do both in this case--simBM DOES work with multifurcating trees, but this doesn't get rid of all 0 edge
#lengths
DinoData<-read.delim('~/../Desktop/Dataset S1.txt')
DinoSize<-setNames(DinoData$Mass,DinoData$Name_in_tree)
tree<-DinoTrees[[1]]
x<-log(DinoSize[DinoTree$tip.label])
tree<-di2multi(tree)
tree$edge.length[tree$edge.length==0]<-1e-6
sim<-test(tree,x,n.sim=10)
plot(sim,DinoTree3,alph=0.1,n.sample=100)

#Okay, here are some more controlled tests--it works fine with polytomies; that doesn't seem to be the issue
tree<-pbtree(n=50)
zero.edges<-sample(1:nrow(tree$edge),10)
tree$edge.length[zero.edges]<-0
tree<-force.ultrametric(tree,method='extend')
plot(tree)
tree<-di2multi(tree)
x<-fastBM(tree)
sim<-simBM(tree,x,n.sim=10)
plot(sim,tree)
#It has problems with 0 edge lengths, since some of the calculations inevitably take the form Inf/Inf, which gives NaN
#For now I added an additional argument to simBM to increase 0 edge lengths to arbitrarily small numbers
tree<-pbtree(n=50)
zero.edges<-sample(1:nrow(tree$edge),10)
tree$edge.length[zero.edges]<-0
x<-fastBM(tree)
fastpgram(x,tree)
sim<-simBM(tree,x,n.sim=10)

np<-100
nsim<-100
ps<-matrix(rnorm(np*nsim,100,1),ncol=np)
ps[,1]<-0
BM<-t(apply(ps,1,cumsum))
plot(BM[1,],type='l',ylim=c(min(BM),max(BM)))
for(i in 2:nsim){
  lines(BM[i,])
}
interpolation<-t(apply(BM,1,function(ii) (ii[np]-ii[1])/(np-1)*(1:np)+ii[1]))
plot(interpolation[1,],type='l',ylim=c(min(interpolation),max(interpolation)))
for(i in 2:nsim){
  lines(interpolation[i,])
}
plot(BM[1,]-interpolation[1,],type='l',ylim=c(min(BM-interpolation),max(BM-interpolation)))
for(i in 2:nsim){
  lines(BM[i,]-interpolation[i,])
}

ps_0<-matrix(rnorm(np*nsim,0,1),ncol=np)
ps_0<-0
BM_0<-t(apply(ps_0,1,cumsum))
plot(BM_0[1,],type='l',ylim=c(min(BM_0),max(BM_0)))
for(i in 2:nsim){
  lines(BM_0[i,])
}
interpolation_0<-t(apply(BM_0,1,function(ii) (ii[np]-ii[1])/(np-1)*(1:np)+ii[1]))
plot(interpolation_0[1,],type='l',ylim=c(min(interpolation_0),max(interpolation_0)))
for(i in 2:nsim){
  lines(interpolation_0[i,])
}
plot(BM_0[1,]-interpolation_0[1,],type='l',ylim=c(min(BM_0-interpolation_0),max(BM_0-interpolation_0)))
for(i in 2:nsim){
  lines(BM_0[i,]-interpolation_0[i,])
}

plot(BM[1,]-interpolation[1,],type='l',ylim=c(min(BM-interpolation),max(BM-interpolation)))
for(i in 2:nsim){
  lines(BM[i,]-interpolation[i,])
  lines(BM_0[i,]-interpolation_0[i,],col='red')
}
