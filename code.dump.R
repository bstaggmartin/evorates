

get.edge.rates<-function(tree,par.mat){
  return(exp(sapply(1:ncol(par.mat),function(ii) apply(cbind(par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,1]],
                                                             par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,2]]),
                                                       1,mean))))
}


start<-proc.time()
par.mat<-relaxed.clock.BM(tree,corateBM$traits,n.iter=1e5,thin=25,tune.period=1e4,report.every=50,win.update=100,block.size=39)
elapsed<-proc.time()-start

thin<-par.mat[,seq(1200,ncol(par.mat),2)]
edge.rates<-get.edge.rates(tree,thin)
pdf("promising initial results2.pdf",height=10,width=10)
par(mfrow=c(2,2))
plot(log(apply(edge.rates,1,median))~log(corateBM$edge.rates),pch=19,xlab="true ln(rate)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=rgb(1,0,0,0.5))
segments(x0=log(corateBM$edge.rates),x1=log(corateBM$edge.rates),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
abline(0,1,lwd=4)
plot(corateBM,tree,norm.lb=-4,norm.ub=4,
     cols=colorRampPalette(c("skyblue2","cyan2","chartreuse2","goldenrod","red"))(20),log=T,edge.width=5,main="true ln(rates)")
frame()
med.edge.rates<-apply(edge.rates,1,median)
plot.corateBM(med.edge.rates,tree,norm.lb=-4,norm.ub=4,
              cols=colorRampPalette(c("skyblue2","cyan2","chartreuse2","goldenrod","red"))(20),log=T,edge.width=5,
              main="median ln(rates) from posterior sample")
dev.off()

#seeing if a reduced mvnorm works
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
#using subtrees to speed up calculations?
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
#doesn't produce any benefit, unfortunately (probably because it has to invert multiple reduced matrices rather than a single one)
#keeping the subtree ID functions-->could prove helpful in other calculations
###Also keep in mind the above approach needs to be tweaked to accomodate cases where the root or a tip node is the focal node
###(in the latter case, the C array is simplified to a matrix because it only consists of one edge, and thus the trans.C 
##function breaks--altering the trans.C function could easily allow it to accomodate matrices, though)

#seeing if a reduced mvnorm works (without subtrees simplification) works in the root case...
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
#it does!!!

#seeing if a reduced mvnorm works for trait values
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


#Use C arrays with dimensions 1 and 2 reduced to the appropriate nodes
#No need to identify appropriate edges and simplify to subtrees, calculating likelihoods for each subtree--while a neat trick,
#it doesn't speed up calculations meaningfully; such a thing could become relevant to many-edged "supertrees", but, in this case,
#R doesn't even have the memory allocation to create the appropriate C arrays...stress tests seem to indicate that R would have
#trouble with a tree ~500 tips...