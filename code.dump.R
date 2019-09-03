###SEEING IF REDUCED MVNORMAL WORKS TO SIMPLIFY LR RATIO CALC AND SPEED UP MCMC###
get.adjacent<-function(n,tree){
  P.edges<-sapply(n,function(nn) which(tree$edge[,2]==nn))
  P.edges[lengths(P.edges)<1]<-NA;P.edges<-unlist(P.edges)
  P.nodes<-tree$edge[P.edges,1]
  D.edges<-sapply(n,function(nn) which(tree$edge[,1]==nn))
  D.edges[lengths(D.edges)<1]<-list(c(NA,NA));D.edges<-matrix(unlist(D.edges),ncol=2,byrow=T)
  D.nodes<-matrix(tree$edge[as.vector(D.edges),2],ncol=2)
  out.nodes<-cbind(n,P.nodes,D.nodes);colnames(out.nodes)<-c("N","P","D1","D2")
  out.edges<-cbind(P.edges,D.edges);colnames(out.edges)<-c("P","D1","D2");rownames(out.edges)<-n
  return(list(nodes=out.nodes,edges=out.edges))
}
full0<-dmvnorm(corateBM$ln.node.rates[-(length(tree$tip.label)+1)],
               mean=corateBM$params$r0+corateBM$params$rtrend*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
               sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(corateBM$params$rvar,dim(C)[3])),log=T)
node.rates<-corateBM$ln.node.rates;node.rates[c(25,35)]<-node.rates[c(25,35)]+rnorm(1,0,1)
full1<-dmvnorm(node.rates[-(length(tree$tip.label)+1)],
               mean=corateBM$params$r0+corateBM$params$rtrend*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
               sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(corateBM$params$rvar,dim(C)[3])),log=T)
test<-as.vector(get.adjacent(c(25,35),tree)$nodes)
start<-proc.time()
red0<-dmvnorm(corateBM$ln.node.rates[test],
              mean=corateBM$ln.node.rates[21]+corateBM$params$rtrend*node.depth.edgelength(tree)[test],
              sigma=trans.C(C[test,test,],rep(corateBM$params$rvar,dim(C)[3])),log=T)
red1<-dmvnorm(node.rates[test],
              mean=node.rates[21]+corateBM$params$rtrend*node.depth.edgelength(tree)[test],
              sigma=trans.C(C[test,test,],rep(corateBM$params$rvar,dim(C)[3])),log=T)
proc.time()-start
##It does!
##it seems that you can crop out edges in addition to nodes when doing reduced calculations, but you have to cut out the "new
##root nodes"...
##Function for identifying new roots given adjacent node/edge lookup table and subsequently identifying the subtrees associated
##w/ each new root
test<-get.adjacent(1:(tree$Nnode*2+1),tree)
focal.nodes<-c(25,35)
nodes<-as.vector(test$nodes[focal.nodes,]);nodes<-nodes[!(is.na(nodes))]
new.roots<-unique(nodes[which(!(test$nodes[nodes,2]%in%nodes))])
subtrees<-list(NULL)
subtrees.edges<-list(NULL)
for(i in 1:length(new.roots)){
  subtrees[[i]]<-new.roots[i]
  subtrees.edges[[i]]<-vector()
  no.D<-1
  while(no.D>0){
    D<-test$nodes[subtrees[[i]][(length(subtrees[[i]])-no.D+1):length(subtrees[[i]])],3:4]
    D.edges<-test$edges[subtrees[[i]][(length(subtrees[[i]])-no.D+1):length(subtrees[[i]])],2:3]
    D.edges<-D.edges[D%in%nodes];D<-D[D%in%nodes]
    no.D<-length(D)
    if(no.D>0){
      subtrees[[i]]<-c(subtrees[[i]],D)
      subtrees.edges[[i]]<-c(subtrees.edges[[i]],D.edges)
    }
  }
}
##using subtrees to speed up calculations?
start<-proc.time()
red0<-red1<-0
for(i in 1:length(subtrees)){
  red0<-sum(red0,dmvnorm(corateBM$ln.node.rates[subtrees[[i]][-1]],
                         mean=corateBM$ln.node.rates[subtrees[[i]][1]]+
                           corateBM$params$rtrend*
                           (node.depth.edgelength(tree)[subtrees[[i]][-1]]-node.depth.edgelength(tree)[subtrees[[i]][1]]),
                         sigma=trans.C(C[subtrees[[i]][-1],subtrees[[i]][-1],subtrees.edges[[i]]],
                                       rep(corateBM$params$rvar,length(subtrees.edges[[i]]))),
                         log=T))
  red1<-sum(red1,dmvnorm(node.rates[subtrees[[i]][-1]],
                         mean=corateBM$ln.node.rates[subtrees[[i]][1]]+
                           corateBM$params$rtrend*
                           (node.depth.edgelength(tree)[subtrees[[i]][-1]]-node.depth.edgelength(tree)[subtrees[[i]][1]]),
                         sigma=trans.C(C[subtrees[[i]][-1],subtrees[[i]][-1],subtrees.edges[[i]]],
                                       rep(corateBM$params$rvar,length(subtrees.edges[[i]]))),
                         log=T))
}
proc.time()-start
##doesn't produce any benefit, unfortunately (probably because it has to invert multiple reduced matrices rather than a single
##one)
##keeping the subtree ID functions-->could prove helpful in other calculations
##Also keep in mind the above approach needs to be tweaked to accomodate cases where the root or a tip node is the focal node
##(in the latter case, the C array is simplified to a matrix because it only consists of one edge, and thus the trans.C 
##function breaks--altering the trans.C function could easily allow it to accomodate matrices, though)
##seeing if a reduced mvnorm works (without subtrees simplification) works in the case of updating r0 parameter...
full0<-dmvnorm(corateBM$ln.node.rates[-(length(tree$tip.label)+1)],
               mean=corateBM$ln.node.rates[21]+corateBM$params$rtrend*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
               sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(corateBM$params$rvar,dim(C)[3])),log=T)
node.rates<-corateBM$ln.node.rates;node.rates[21]<-node.rates[21]+rnorm(1,0,1)
full1<-dmvnorm(node.rates[-(length(tree$tip.label)+1)],
               mean=node.rates[21]+corateBM$params$rtrend*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
               sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(corateBM$params$rvar,dim(C)[3])),log=T)
test<-c(as.vector(get.adjacent.nodes(21,tree)))[-1]
test<-test[!is.na(test)]
red0<-dmvnorm(corateBM$ln.node.rates[test],
              mean=corateBM$ln.node.rates[21]+corateBM$params$rtrend*node.depth.edgelength(tree)[test],
              sigma=trans.C(C[test,test,],rep(corateBM$params$rvar,dim(C)[3])),log=T)
red1<-dmvnorm(node.rates[test],
              mean=node.rates[21]+corateBM$params$rtrend*node.depth.edgelength(tree)[test],
              sigma=trans.C(C[test,test,],rep(corateBM$params$rvar,dim(C)[3])),log=T)
##it does!!!
##seeing if a reduced mvnorm works for trait values
full0<-dmvnorm(corateBM$traits,
               mean=corateBM$params$x0+corateBM$params$xtrend*node.depth.edgelength(tree)[1:length(tree$tip.label)],
               sigma=trans.C(C[1:length(tree$tip.label),1:length(tree$tip.label),],corateBM$edge.rates),log=T)
edge.rates<-corateBM$edge.rates;edge.rates[5]<-edge.rates[5]*exp(rnorm(1,0,1))
full1<-dmvnorm(corateBM$traits,
               mean=corateBM$params$x0+corateBM$params$xtrend*node.depth.edgelength(tree)[1:length(tree$tip.label)],
               sigma=trans.C(C[1:length(tree$tip.label),1:length(tree$tip.label),],edge.rates),log=T)
test<-getDescendants(tree,23)
test<-test[test<=length(tree$tip.label)]
red0<-dmvnorm(corateBM$traits[test],
              mean=corateBM$params$x0+corateBM$params$xtrend*node.depth.edgelength(tree)[test],
              sigma=trans.C(C[test,test,],corateBM$edge.rates),log=T)
red1<-dmvnorm(corateBM$traits[test],
              mean=corateBM$params$x0+corateBM$params$xtrend*node.depth.edgelength(tree)[test],
              sigma=trans.C(C[test,test,],edge.rates),log=T)
red1-red0;full1-full0
#only apporximately, and the approximation gets better and better as you go towards roots--like to do with the fact that internal
#nodes are not being directly estimated...
###TAKE HOMES###
##Use C arrays with dimensions 1 and 2 reduced to the appropriate nodes
##No need to identify appropriate edges and simplify to subtrees, calculating likelihoods for each subtree--while a neat trick,
##it doesn't speed up calculations meaningfully; such a thing could become relevant to many-edged "supertrees", but, in this case,
##R doesn't even have the memory allocation to create the appropriate C arrays anyways...stress tests seem to indicate that R
##would have trouble with a tree ~500 tips...

###MUSINGS ON SIMBM AND "OVER-OPTIMIZATION" OF ANCESTRAL STATES###
##So I figured out why simBM often returns ancestral state estimates wayyy different from those the phylogeny was simulatedunder:
tree<-pbtree(n=25)
ntax<-length(tree$tip.label);nnode<-tree$Nnode*2+1
x<-fastBM(tree,internal=T,sig2=500)
n.sim=1000
test<-simBM(tree,x[1:length(tree$tip.label)],n.sim=n.sim)
plot(test,tree,n.sample=20,alph=0.1)
segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
         y0=x[tree$edge[,1]],y1=x[tree$edge[,2]])
plot(as.vector(test$nodes[(ntax+1):nnode,])~rep(x[(ntax+1):nnode],n.sim),pch=16,col=rgb(0,0,0,0.1))
abline(0,1)
library(mvnfast)
C<-array(0,dim=c(rep(tree$Nnode*2+1,2),nrow(tree$edge)))
for(e in 1:nrow(tree$edge)){
  D<-unique(c(tree$edge[e,2],getDescendants(tree,tree$edge[e,2])))
  C[as.matrix(cbind(expand.grid(D,D),rep(e,length(D)^2)))]<-tree$edge.length[e]
}
C<-apply(C,c(1,2),sum)
C<-C[-(ntax+1),-(ntax+1)]
loglik<-dmvn(x[-(ntax+1)],rep(0,nnode-1),C,log=T)
sim.logliks<-dmvn(t(test$nodes[-(ntax+1),]),rep(0,nnode-1),C,log=T)
hist(sim.logliks,xlim=c(min(sim.logliks,loglik),max(sim.logliks,loglik)))
abline(v=loglik)
##The generating model is almost never anywhere near an optimal likelihood value; really, all you have is one, highly
##multidimensional data point!
##Perhaps you could try to adjust for that? Adjust the error around each node, ajust for the fact that some nodes have fewer
##descendants and thus we can call their ancestral state with less certainty?
plot(as.vector(test$nodes[(ntax+1):nnode,])~jitter(rep(node.depth.edgelength(tree)[(ntax+1):nnode],n.sim),100),pch=16,col=rgb(0,0,0,0.05))
segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
          y0=x[tree$edge[,1]],y1=x[tree$edge[,2]])
points(x[(ntax+1):nnode]~node.depth.edgelength(tree)[(ntax+1):nnode],pch=16,col="red")

###IS IT APPROPRIATE TO AVERAGE NODES TO APPROXIMATE AVERAGE EDGE VALUES?###
##quick test regarding the appropriateness of averaging nodes to approximate average edge values...
BM<-matrix(NA,nrow=1000,ncol=1000)
for(i in 1:1000){
  BM[i,]<-cumsum(rnorm(1000))
}
plot(apply(BM,1,mean)~apply(BM[,c(1,1000)],1,mean))
summary(lm(apply(BM,1,mean)~apply(BM[,c(1,1000)],1,mean)))
##seems the expected value of a brownian bridge is probably the mean of the start and the end point, but there's considerable
##variance...is there a analytical solution for this? R^2 is only ~0.75...
