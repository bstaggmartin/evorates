###Michael's histogram visuialization idea###
#Have color-coded histograms that track trait value distributions over time for given lineages
##Maybe even make the colors evolve according to brownian motion to improve visualization?
color.scale<-sapply(seq(0.4,1,length.out=100),function(ii) rainbow(100,s=ii,v=ii))
plot(expand.grid(1:100,1:100),pch=16,col=as.vector(color.scale))
library(phytools)
tree<-pbtree(n=25)
#colx<-cbind(runif(length(tree$tip.label),1,100),runif(length(tree$tip.label),1,100))
colx<-fastBM(tree,a=0,sig2=10,nsim=2)
colx<-rbind(colx,cbind(fastAnc(tree,colx[,1]),fastAnc(tree,colx[,2])))
edge.avg<-cbind(apply(cbind(colx[tree$edge[,1],1],colx[tree$edge[,2],1]),1,mean),
                apply(cbind(colx[tree$edge[,1],2],colx[tree$edge[,2],2]),1,mean))
std.edge.avg<-apply(edge.avg,2,function(xx) round((xx-min(xx))/(max(xx)-min(xx))*99+1))
plot(tree,edge.width=4,edge.color=color.scale[cbind(std.edge.avg[,1],std.edge.avg[,2])])
plot(std.edge.avg[,1]~std.edge.avg[,2],col=color.scale[cbind(std.edge.avg[,1],std.edge.avg[,2])],pch=16)
colmap<-color.scale[cbind(std.edge.avg[,1],std.edge.avg[,2])]
base.edges<-which(tree$edge[,1]==(length(tree$tip.label)+1))
colmap<-c(colmap,color.scale[cbind(round(mean(std.edge.avg[base.edges,1])),round(mean(std.edge.avg[base.edges,2])))])
x<-fastBM(tree)
source("simBM.R");source("utils.R");library(truncnorm);library(abind)
test<-simBM(tree,x)



plot.trait.profile<-function(sim,tree,t,colmap,alph=0.2){
  plot(0,xlim=c(min(simBM$edges,na.rm=T),max(simBM$edges,na.rm=T)),ylim=c(0,3),col='white',xlab="trait val",ylab="prob density")
  if(t==1){
    polygon(density(simBM$nodes[length(tree$tip.label)+1,]),col=alter.cols(colmap[length(colmap)],alph,name=T),border=NA)
  }else if(t==length(simBM$ts)){
    segments(y0=rep(0,length(tree$tip.label)),y1=rep(1e6,length(tree$tip.label)),
             x0=simBM$nodes[1:length(tree$tip.label),1],x1=simBM$nodes[1:length(tree$tip.label),1],
             col=alter.cols(colmap[sapply(1:length(tree$tip.label),function(ii) which(tree$edge[,2]==ii))],alph,name=T),lwd=4)
  }else{
    for(i in 1:dim(simBM$edges)[1]){
      if(!all(is.na(simBM$edges[i,t,]))){
        polygon(density(simBM$edges[i,t,]),col=alter.cols(colmap[i],alph,name=T),border=NA)
      }
    }
  }
}
plot.trait.profile(test,tree,99,colmap)


library(magick)
setwd("render_test")
file.remove(list.files(pattern=".png"))
png(paste("profile","%02d.png",sep=""),width=600,height=600)
par(mfrow=c(2,1))
for(t in 1:dim(test$edges)[2]){
  par(mar=c(0.2,5,1,3))
  #plot(tree,edge.color=colmap,edge.width=4)
  plot(0,xlim=c(min(node.depth.edgelength(tree)),max(node.depth.edgelength(tree))),
       ylim=c(min(test$edges,na.rm=T),max(test$edges,na.rm=T)),
       axes=F,xlab='',ylab='',bty='n',col=rgb(0,0,0,0))
  segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
           y0=apply(test$nodes[tree$edge[,1],],1,mean),y1=apply(test$nodes[tree$edge[,2],],1,mean),
           col=alter.cols(colmap,name=T),lwd=4)
  abline(v=test$ts[t],lwd=3,col=rgb(0,0,0,0.5))
  par(mar=c(5,5,0.2,3))
  plot.trait.profile(test,tree,t,colmap)
}
dev.off()
profiles<-paste("profile",c(paste("0",1:9,sep=""),10:dim(test$edges)[2]),".png",sep="")
images<-list(NULL)
for(i in 1:dim(test$edges)[2]){
  images[[i]]<-image_read(profiles[i])
}
frames<-image_join(c(images,rep(images[[length(images)]],12),rev(images),rep(images[[1]],12)))
animation<-image_animate(frames,fps=10)
image_write(animation,paste("render",".gif",sep=""))
file.remove(list.files(pattern=".png"))






#look into edge collating function--turned out to produce some odd behaviors when you played around with it here
#also, simBM returns the edge matrix  with first and last column containing all NAs--makes sense given the sequence of time points
#under which you simulate value, but also those columns are kinda pointless...probably partially accounts for unexpected behavior
#of edge collater
samp<-1:1000
time.vec<-sort(unique(round(node.depth.edgelength(tree),digits=6)))
new.edges<-array(NA,dim=c(nrow(tree$edge),length(time.vec),length(samp)))
n1<-tree$edge[,1];n2<-tree$edge[,2]
t1<-round(node.depth.edgelength(tree)[n1],6);t2<-round(node.depth.edgelength(tree)[n2],6)
new.edges[cbind(rep(rep(1:nrow(tree$edge),2),length(samp)),
                rep(sapply(c(t1,t2),function(x) which(x==time.vec)),length(samp)),
                rep(1:length(samp),each=nrow(tree$edge)*2))]<-
  test$nodes[cbind(rep(c(n1,n2),length(samp)),
                   rep(samp,each=nrow(tree$edge)*2))]
new.edges<-abind(new.edges,test$edges[,,samp],along=2)
time.vec<-c(time.vec,test$ts)
ord<-order(time.vec)
new.edges<-new.edges[,ord,]
time.vec<-time.vec[ord]
dups<-which(duplicated(time.vec))
new.edges[,dups-1,]<-ifelse(is.na(new.edges[,dups,]),new.edges[,dups-1,],new.edges[,dups,])
new.edges<-new.edges[,-dups,]
time.vec<-time.vec[-dups]
na.entries<-sapply(1:dim(new.edges)[2],function(ii) all(is.na(new.edges[,ii,])))
new.edges<-new.edges[,-na.entries,]
time.vec<-time.vec[-na.entries]
new.edges
time.vec