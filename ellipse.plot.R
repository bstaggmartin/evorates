ellipse.plot<-function(fit,traits=colnames(fit$call$trait.data),
                       element=c('chains','quantiles','means','MAPs'),
                       type=c('evocov','intracov'),
                       select.extra=NULL,res=25,sample=500,
                       kernel.res=100,density.res=100,plot.type='lines',add=F,level=0.95){
  call<-as.list(match.call(ellipse.plot))
  #handle matching element before call...
  cov.mats<-do.call(get.cov.mat,c(call[names(call)%in%names(formals(get.cov.mat))],simplify=F))
  if(plot.type=='lines'){
    angles<-seq(0,2*pi,length.out=res)
    angles<-rbind(cos(angles),sin(angles))
  }
  for(i in 1:dim(cov.mats)[1]){
    for(j in 1:i){
      if(j==i){
        next
      }
      if(dim(cov.mats)[3]>sample){
        tmp<-.index.element(cov.mats,list(c(i,j),c(i,j),sample(dim(cov.mats)[3],sample)),c(1,2,3))
      }else{
        tmp<-.index.element(cov.mats,c(i,j),c(1,2))
      }
      eigen.decomp<-array(apply(tmp,3:4,eigen),dim(tmp)[3:4],dimnames(tmp)[3:4])
      coords<-array(list(matrix(NA,res,2)),dim(eigen.decomp))
      for(k in 1:ncol(eigen.decomp)){
        for(l in 1:nrow(eigen.decomp)){
          if(any(eigen.decomp[[l,k]]$values<0)){
            warning('oh no')
            next
          }
          if(plot.type!='lines'){
            angles<-sort(runif(res,0,2*pi))
            angles<-rbind(cos(angles),sin(angles))
          }
          coords[[l,k]]<-t(eigen.decomp[[l,k]]$vectors%*%(angles*sqrt(eigen.decomp[[l,k]]$values)))
        }
      }
      xlim<-range(unlist(coords)[rep(c(T,F),each=res)],na.rm=T)
      ylim<-range(unlist(coords)[rep(c(F,T),each=res)],na.rm=T)
      if(!add){
        plot(0,col='white',
             xlim=xlim,
             ylim=ylim,
             xlab=dimnames(cov.mats)[[1]][j],ylab=dimnames(cov.mats)[[1]][i])
      }
      for(k in 1:ncol(eigen.decomp)){
        if(plot.type=='lines'){
          for(l in 1:nrow(eigen.decomp)){
            lines(coords[[l,k]],lty=k,col=k)
          }
        }else if(plot.type=='points'){
          points(unlist(coords[,k])[rep(c(F,T),each=res)]~unlist(coords[,k])[rep(c(T,F),each=res)],col=k)
        }else if(plot.type=='density'){
          tmp.dens<-MASS::kde2d(unlist(coords[,k])[rep(c(T,F),each=res)],
                                unlist(coords[,k])[rep(c(F,T),each=res)],
                                n=kernel.res,
                                lims=c(xlim+c(-1,1)*0.1*diff(xlim),ylim+c(-1,1)*0.1*diff(ylim)))
          image(tmp.dens,col=contSimmap::alter.cols(rep(k,length(level)),alpha=seq(0.2,0.4,length.out=length(level))),
                  #colorRampPalette(contSimmap::alter.cols(rep(k,2),alpha=c(0,0.2)),alpha=T)(1),
                breaks=c(dnorm(qnorm((1-level)/2))/dnorm(0)*max(tmp.dens$z),max(tmp.dens$z)),add=T)
        }
      }
    }
  }
}
#for most applications, res of 40 and sample of 500 seems to do the trick...
#disable diagnostics and (maybe) quantiles
#add support for line, points, or density plotting...
tree<-phytools::pbtree(n=25)
X<-phytools::fastBM(tree,nsim=3)
fit<-contSimmap::fit.corateBM(tree,X,chains=2,cores=2)
ellipse.plot(fit,sample=500,res=25)

plot(0,col='white',xlim=c(-2,2),ylim=c(-2,2))
image(hmm,col = colorRampPalette(c(rgb(0,0,0,0),palette()[k]),alpha=T)(100))

test<-contourLines(hmm)
levs<-sort(unique(sapply(test,'[[','level')))
cols<-colorRampPalette(c(rgb(0,0,0,0),palette()[k]),alpha=T)(length(levs))
names(cols)<-levs
for(i in 1:length(test)){
  polygon(test[[i]]$x,test[[i]]$y,col=cols[as.character(test[[i]]$level)],border=NA)
}

res<-100
zlim<-range(hmm$z)
colramp<-colorRampPalette(c(rgb(0,0,0,0),palette()[k]),alpha=T)(res)
inds<-round((hmm$z-zlim[1])/(zlim[2]-zlim[1])*(res-1))+1
inds[inds<1]<-1;inds[inds>res]<-res
colvec<-colramp[inds]
off.x<-diff(hmm$x)[1]/2
off.y<-diff(hmm$y)[1]/2
for(i in hmm$x){
  for(j in hmm$y){
    polygon(c(i-off.x,i+off.x,i+off.x,i-off.x),c(j-off.y,j-off.y,j+off.y,j+off.y),col=colvec[(i-1)*100+j],border=NA)
  }
}

level<-0.95
dnorm(qnorm((1-level)/2))/dnorm(0)
tmp<-tmp.dens$z/max(tmp.dens$z)
tmp[,]<-ifelse(tmp>dnorm(qnorm((1-level)/2))/dnorm(0),1,0)
image(tmp,)
plot(0,col='white',xlim=xlim,ylim=ylim)
testy<-contourLines(list(x=tmp.dens$x,y=tmp.dens$y,z=tmp),levels=dnorm(qnorm((1-level)/2))/dnorm(0))
polygon(testy[[2]]$x,testy[[2]]$y)

xx<-mean(unlist(sapply(testy,'[','x')))
yy<-mean(unlist(sapply(testy,'[','y')))
r<-NULL
theta<-NULL
for(i in 1:length(testy[[1]]$x)){
  dy<-testy[[1]]$y[i]-yy
  dx<-testy[[1]]$x[i]-xx
  r[i]<-sqrt((dy)^2+(dx)^2)
  theta[i]<-atan(dy/dx)
  if(dx<0){
    theta[i]<-theta[i]+pi
  }
}
r2<-NULL
theta2<-NULL
for(i in 1:length(testy[[2]]$x)){
  dy<-testy[[2]]$y[i]-yy
  dx<-testy[[2]]$x[i]-xx
  r2[i]<-sqrt((dy)^2+(dx)^2)
  theta2[i]<-atan(dy/dx)
  if(dx<0){
    theta2[i]<-theta2[i]+pi
  }
}
yyy<-r*sin(theta)
xxx<-r*cos(theta)
plot(yyy~xxx,type='l')

plot(r~theta)
plot(r2~theta2)

r<-r[order(theta)]
theta<-sort(theta)
r2<-r2[order(theta2)]
theta2<-sort(theta2)
angs<-seq(min(c(theta,theta2)),max(c(theta,theta2)),length.out=22)[-c(1,22)]
rs<-matrix(NA,20,2)
for(i in 1:length(angs)){
  ang11<-min(which(theta>angs[i]))
  ang12<-max(which(theta<angs[i]))
  ang21<-min(which(theta2>angs[i]))
  ang22<-max(which(theta2<angs[i]))
  rs[i,]<-c((r[ang12]-r[ang11])/(theta[ang12]-theta[ang11])*(angs[i]-theta[ang11])+r[ang11],
            (r2[ang22]-r2[ang21])/(theta2[ang22]-theta2[ang21])*(angs[i]-theta2[ang21])+r2[ang21])
}
plot(rs[,1]~angs)
plot(rs[,2]~angs)
yyy<-rs[,1]*sin(angs)
xxx<-rs[,1]*cos(angs)
plot(yyy~xxx)
#function returns NAs--could actually use this...
#basically a polar sweep line algorithm

apply(rs,1,function(ii) ii[1]>ii[2])
#if all NAs--no intersections
#if any TRUE, none FALSE--shape 1 encompasses shape 2
#if any FALSE, none TRUE--shape 2 encompasses shape 1
new.testy<-testy
tmp.y1<-abs(new.testy[[1]]$y-yy)
out.slice1<-which(tmp.y1==min(tmp.y1[new.testy[[1]]$x-xx<0])&new.testy[[1]]$x-xx<0)
tmp.y2<-abs(new.testy[[2]]$y-new.testy[[1]]$y[out.slice1])
in.slice1<-which(tmp.y2==min(tmp.y2[new.testy[[2]]$x-xx<0])&new.testy[[2]]$x-xx<0)
in.slice2<-which(tmp.y2==min(tmp.y2[new.testy[[2]]$x-xx>0])&new.testy[[2]]$x-xx>0)
out.slice2<-which(tmp.y1==min(tmp.y1[new.testy[[1]]$x-xx>0])&new.testy[[1]]$x-xx>0)
vec1<-out.slice2:out.slice1
vec2<-in.slice1:in.slice2
new.testy[[1]]$y<-c(testy[[1]]$y[vec1],testy[[2]]$y[vec2])
new.testy[[1]]$x<-c(testy[[1]]$x[vec1],testy[[2]]$x[vec2])
new.testy[[2]]$y<-c(testy[[1]]$y[-vec1[-c(1,length(vec1))]],testy[[2]]$y[-vec2[-c(1,length(vec2))]])
new.testy[[2]]$x<-c(testy[[1]]$x[-vec1[-c(1,length(vec1))]],testy[[2]]$x[-vec2[-c(1,length(vec2))]])


#mix of FALSE and TRUE--shape 1 and 2 intersect one another... (shouldn't happen with contour lines)
#should vectorize  to make useful...
#can use this for general-purpose drawing polygons for levels--might be faster than hi-res images (will look
#nicer, in any case)

#strive for this--but hold for now--would need to generalize to cases with multiple holes, probably...

#quantile options--either select quantiles by covariance or variance?



#figured it out!!!
plot(0,col='white',xlim=c(-2,2),ylim=c(-2,2))
tmp.coords<-rbind(cbind(testy[[1]]$x,testy[[1]]$y),cbind(testy[[2]]$x,testy[[2]]$y))
tmp.coords<-rbind(cbind(testy[[1]]$x,testy[[1]]$y),cbind(c(-0.1,0.1,0.1,-0.1,-0.1)-0.5,c(0.9,0.9,1.1,1.1,0.9)),cbind(c(-0.1,0.1,0.1,-0.1,-0.1),c(0.9,0.9,1.1,1.1,0.9)),cbind(-0.6,0.9),cbind(testy[[1]]$x[length(testy[[1]]$x)],testy[[1]]$y[length(testy[[1]]$x)]))
#outer closed polygon, all inner closed polygons..., end points of all polygons except last in reverse sequence
polygon(tmp.coords,fillOddEven=T,col='black')
#find larger polygon

coords<-cbind(c(-0.1,0.1,0.1,-0.1,-0.1),c(0.9,0.9,1.1,1.1,0.9))
coords<-cbind(testy[[2]]$x,testy[[2]]$y)
abs(0.5*sum(coords[-nrow(coords),1]*coords[-1,2]-coords[-nrow(coords),2]*coords[-1,1]))
get.area<-function(ls){
  out<-rep(NA,length(ls))
  for(i in 1:length(ls)){
    len<-length(ls[[i]]$x)
    out[i]<-abs(0.5*sum(ls[[i]]$x[-len]*ls[[i]]$y[-1]-ls[[i]]$y[-len]*ls[[i]]$x[-1]))
  }
  out
}
testy
get.area(testy)
make.coords.list<-function(ls){
  xx<-sapply(ls,'[','x');yy<-sapply(ls,'[','y')
  out<-cbind(unlist(xx),unlist(yy))
  out<-rbind(out,out[cumsum(lengths(xx))[length(ls):1],])
  out
}

fit<-readRDS('large_fit')

testy<-contourLines(tmp.dens,levels=dnorm(qnorm((1-c(0.99,0.95,0.9,0.8))/2))/dnorm(0)*max(tmp.dens$z))
new.testy<-split(testy,sapply(testy,'[[','level'))
new.testy2<-lapply(new.testy,make.coords.list)
colvec<-contSimmap::alter.cols(rep(1,length(new.testy)+1),alpha=seq(0,0.4,length.out=length(new.testy)+1))[-1]
plot(0,col='white',xlim=c(-1.5,1.5),ylim=c(-1.6,1.6))
for(i in 1:length(new.testy2)){
  polygon(new.testy2[[i]],col=colvec[i],border=NA)
}
for(i in 1:length(new.testy)){
  for(j in 1:length(new.testy[[i]])){
    lines(new.testy[[i]][[j]])
  }
}
#perfect--contour lines already organizes by area
#got it all working--should be generalizable to any group of non-intersecting polygons...overlap doesn't seem
#to create much visual issue
#about 2x faster than image() for a limited number of polygons...and much smoother visuals
#could improve by mixing each level wih next level up to cut out holes in lower polygons...
get.area<-function(ls){
  out<-rep(NA,length(ls))
  for(i in 1:length(ls)){
    len<-length(ls[[i]]$x)
    out[i]<-abs(0.5*sum(ls[[i]]$x[-len]*ls[[i]]$y[-1]-ls[[i]]$y[-len]*ls[[i]]$x[-1]))
  }
  out
}
expand.range<-function(x,fac){
  ran<-range(x)
  ran+c(-1,1)*fac*abs(diff(ran))
}
make.coords.mat<-function(ls){
  xx<-sapply(ls,'[','x');yy<-sapply(ls,'[','y')
  out<-cbind(unlist(xx),unlist(yy))
  out<-rbind(out,out[cumsum(lengths(xx))[length(ls):1],])
  out
}
interpolate<-function(x,length.out){
  xx<-seq(1,length(x),length.out=length.out)
  in.inds<-xx%in%(1:length(xx))
  out<-rep(NA,length.out)
  out[in.inds]<-x[xx[in.inds]]
  xx<-xx[!in.inds]
  out[!in.inds]<-(x[ceiling(xx)]-x[floor(xx)])/(ceiling(xx)-floor(xx))*(xx-floor(xx))+x[floor(xx)]
  out
}
hist2d<-function(x,y,breaks=NULL,nbreaks=NULL,res=25,as.dens=F,range.expand.fac=0.1,add=F,...){
  dens<-MASS::kde2d(x,y,n=res,lims=c(expand.range(x,range.expand.fac),expand.range(y,range.expand.fac)))
  if(is.null(nbreaks)){
    nbreaks<-nclass.Sturges(dens$z)
  }
  if(is.null(breaks)){
    breaks<-pretty(dens$z,nbreaks)
  }else{
    breaks<-sort(unique(breaks),T)
    breaks<-dnorm(qnorm((1-breaks)/2))/dnorm(0)*max(dens$z)
    if(breaks[1]!=-1){
      breaks<-c(-1,breaks)
    }else if(breaks[1]==0){
      breaks[1]<- -1
    }
    if(breaks[length(breaks)]!=max(dens$z)+1){
      breaks<-c(breaks,max(dens$z)+1)
    }else if(breaks[length(breaks)]==max(dens$z)){
      breaks[length(breaks)]<-max(dens$z)+1
    }
  }
  col.list<-list()
  for(i in c('col','border')){
    if(is.null(list(...)[[i]])){
      tmp<-NA
    }else{
      tmp<-list(...)[[i]]
    }
    if(all(is.na(tmp))){
      col.list[[i]]<-rep(NA,length(breaks)-1)
    }else{
      col.list[[i]]<-c(rgb(0,0,0,0),colorRampPalette(tmp,alpha=T)(length(breaks)-2))
    }
  }
  if(!add){
    do.call(plot,
            c(x=list(x),
              y=list(y),
              col='white',
              type='p',
              pch=1,
              xlab='x',
              ylab='y',
              list(...)[!(names(list(...))%in%c('col','type','pch'))]))
  }
  if(as.dens){
    poly.list<-list()
    for(i in c('density','angle')){
      if(is.null(list(...)[[i]])){
        tmp<-NA
      }else{
        tmp<-list(...)[[i]]
      }
      if(all(is.na(tmp))){
        poly.list[[i]]<-rep(NA,length(breaks)-1)
      }else{
        poly.list[[i]]<-c(0,interpolate(tmp,length(breaks)-2))
      }
    }
    poly.list[['angle']]<-ifelse(is.na(poly.list[['angle']]),45,poly.list[['angle']])
    cont.lines<-contourLines(dens,levels=breaks)
    levs<-sapply(cont.lines,'[[','level')
    cont.areas<-get.area(cont.lines)
    cont.areas<-split(cont.areas,levs)
    cont.lines<-split(cont.lines,levs)
    cont.poly<-cont.lines
    for(i in 1:(length(cont.poly)-1)){
      tmp.areas<-c(cont.areas[[i]],cont.areas[[i+1]][cont.areas[[i+1]]<max(cont.areas[[i]])])
      tmp.areas.ord<-order(tmp.areas,decreasing=T)
      cont.poly[[i]]<-c(cont.poly[[i]],cont.poly[[i+1]][cont.areas[[i+1]]<max(cont.areas[[i]])])[tmp.areas.ord]
    }
    cont.poly<-lapply(cont.poly,make.coords.mat)
    if(!all(is.na(col.list['col']))){
      for(i in 1:length(cont.poly)){
        do.call(polygon,
                c(x=list(cont.poly[[i]][,1]),
                  y=list(cont.poly[[i]][,2]),
                  col=col.list[['col']][-1][i],
                  density=poly.list[['density']][-1][i],
                  angle=poly.list[['angle']][-1][i],
                  border=NA,
                  fillOddEven=T,
                  list(...)[!(names(list(...))%in%c('col','density','angle','border','fillOddEven'))]))
      }
    }
    for(i in 1:length(cont.lines)){
      for(j in 1:length(cont.lines[[i]])){
        do.call(lines,
                c(x=list(cont.lines[[i]][[j]]$x),
                  y=list(cont.lines[[i]][[j]]$y),
                  col=col.list[['border']][-1][i],
                  list(...)[!(names(list(...))%in%c('col','border'))]))
      }
    }
  }else{
    do.call(image,
            c(x=list(dens$x),
              y=list(dens$y),
              z=list(dens$z),
              col=list(col.list[['col']]),
              breaks=breaks,
              add=T,
              list(...)[!(names(list(...))%in%c('z','col','breaks','add'))]))
  }
}

#fuck--.filled.contour works x.x
#we'll see--your function is slower, but more flexible, I think...

#it works, just a few additional things:
#cut out polygon of the next level by sorting them into poly list BY area, excluding any polygons larger than
#largest area of current level
#work on cases with cutoff--causes weird fill artifacts...

test.fun<-function(...){
  col.list<-list()
  for(i in c('col','border')){
    if(!hasArg(i)){
      tmp<-NA
    }else{
      tmp<-list(...)[[i]]
      # if(is.null(get(i))){
      #   tmp<-NA
      # }else{
      #   tmp<-get(i)
      # }
    }
    # if(all(is.na(tmp))){
    #   col.list[[i]]<-rep(NA,length(breaks)-1)
    # }else{
    #   col.list[[i]]<-c(rgb(0,0,0,0),colorRampPalette(tmp,alpha=T)(length(breaks)-2))
    # }
  }
  tmp
  # col.list
}


hist2d<-function(x,y,breaks=NULL,nbreaks=NULL,res=25,as.dens=F,range.expand.fac=0.1,add=F,...){
  geom.formals<-c('density','angle','type','lty','lwd','pch','lend','lmitre')
  dens<-MASS::kde2d(x,y,n=res,lims=c(expand.range(x,range.expand.fac),expand.range(y,range.expand.fac)))
  if(is.null(nbreaks)){
    nbreaks<-nclass.Sturges(dens$z)
  }
  if(is.null(breaks)){
    breaks<-pretty(dens$z,nbreaks)
  }else{
    breaks<-sort(unique(breaks),T)
    breaks<-dnorm(qnorm((1-breaks)/2))/dnorm(0)*max(dens$z)
    if(breaks[1]!=-1){
      breaks<-c(-1,breaks)
    }else if(breaks[1]==0){
      breaks[1]<- -1
    }
    if(breaks[length(breaks)]!=max(dens$z)+1){
      breaks<-c(breaks,max(dens$z)+1)
    }else if(breaks[length(breaks)]==max(dens$z)){
      breaks[length(breaks)]<-max(dens$z)+1
    }
  }
  col.list<-list()
  for(i in c('col','border')){
    if(is.null(list(...)[[i]])){
      tmp<-NA
    }else{
      tmp<-list(...)[[i]]
    }
    if(all(is.na(tmp))){
      col.list[[i]]<-rep(NA,length(breaks)-1)
    }else{
      col.list[[i]]<-c(rgb(0,0,0,0),colorRampPalette(tmp,alpha=T)(length(breaks)-2))
    }
  }
  if(!add){
    do.call(plot,
            c(x=list(x),
              y=list(y),
              col='white',
              type='p',
              pch=1,
              xlab='x',
              ylab='y',
              list(...)[!(names(list(...))%in%c('col','type','pch'))]))
  }
  if(as.dens){
    args.list<-list(...)[names(list(...))%in%geom.formals]
    args.list<-lapply(args.list,function(ii) c(0,interpolate(ii,length(breaks)-2)))
    args.list
    cont.lines<-contourLines(dens,levels=breaks)
    levs<-sapply(cont.lines,'[[','level')
    cont.areas<-get.area(cont.lines)
    cont.areas<-split(cont.areas,levs)
    cont.lines<-split(cont.lines,levs)
    cont.poly<-cont.lines
    args.list
    if(length(cont.poly)>1){
      for(i in 1:(length(cont.poly)-1)){
        tmp.areas<-c(cont.areas[[i]],cont.areas[[i+1]][cont.areas[[i+1]]<max(cont.areas[[i]])])
        tmp.areas.ord<-order(tmp.areas,decreasing=T)
        cont.poly[[i]]<-c(cont.poly[[i]],cont.poly[[i+1]][cont.areas[[i+1]]<max(cont.areas[[i]])])[tmp.areas.ord]
      }
    }
    cont.poly<-lapply(cont.poly,make.coords.mat)
    if(!all(is.na(col.list['col']))){
      for(i in 1:length(cont.poly)){
        do.call(polygon,
                c(x=list(cont.poly[[i]][,1]),
                  y=list(cont.poly[[i]][,2]),
                  col=col.list[['col']][-1][i],
                  border=NA,
                  fillOddEven=T,
                  lapply(args.list[names(args.list)%in%geom.formals],function(ii) ii[-1][i])))
      }
    }
    for(i in 1:length(cont.lines)){
      for(j in 1:length(cont.lines[[i]])){
        do.call(lines,
                c(x=list(cont.lines[[i]][[j]]$x),
                  y=list(cont.lines[[i]][[j]]$y),
                  col=col.list[['border']][-1][i],
                  lapply(args.list[names(args.list)%in%geom.formals],function(ii) ii[-1][i])))
      }
    }
  }else{
    do.call(image,
            c(x=list(dens$x),
              y=list(dens$y),
              z=list(dens$z),
              col=list(col.list[['col']]),
              breaks=breaks,
              add=T,
              list(...)[!(names(list(...))%in%c('z','col','breaks','add'))]))
  }
}
#for some reason, polygon won't accept lwd of 0...
#fixed
#now just need to clean up argument matching and bandwidth specification...
#need to work on breaks--something still not quite right here...