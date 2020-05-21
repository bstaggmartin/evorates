trace.plot<-function(fit,select='.',separate.R=F,...){#,
                     # args.R0=list(),args.Rsig2=list(),args.X0=list(),args.bg.rate=list(),args.R=list(),
                     # separate.X=F,args.Xsig2=list(),args.X=list()){
  # args.master<-list(args.R0=args.R0,args.Rsig2=args.Rsig2,args.X0=args.X0,args.bg.rate=args.bg.rate,args.R=args.R,
  #                   args.Xsig2=args.Xsig2,args.X=args.X)
  # for(i in 1:length(args.master)){
  #   tmp<-list(...)
  #   tmp<-tmp[!(names(tmp)%in%names(args.master[i]))]
  #   args.master[i]<-c(args.master[i],tmp)
  # }
  #actually, restructure the 'args' bit --> first select the parameters, than iterate through any  possible args, check if the
  #function has those arguments, and construct the args.master list...
  fit$chains<-fit%chains%select
  chains<-check.n.proc(fit,'chains')
  
  #do it here
  
  for(i in dimnames(chains)[[2]]){
    if(!separate.R){
      if(R.flag){
        next
      }
      if(grepl('^R\\[\\d+\\]$',i)){
        matplot(fit%chains%'R[',xlab='iterations',ylab='R',type='l',args.R)
      }
    }
    if(!separate.X){
      
    }
  }
    
  
  if(is.null(select)){
    plot(fit%chains%'R0',type='l',xlab='iterations',ylab='R0',type='l',args.R0)
    plot(fit%chains%'Rsig2',xlab='iterations',ylab='Rsig2',type='l',args.Rsig2)
    plot(fit%chains%'X0',xlab='iterations',ylab='X0',type='l',args.X0)
    plot(fit%chains%'bg.rate',xlab='iterations',ylab='bg.rate',type='l',args.bg.rate)
    if(separate.R){
      for(i in blah){
        
      }
    }else{
      
    }
  }else{
    
  }
}
