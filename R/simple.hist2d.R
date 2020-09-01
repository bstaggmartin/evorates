simple.hist2D<-function(x,y,p=0.5,...){
  
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