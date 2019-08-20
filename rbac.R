rbac<-function(n,x=0,sd=1,m=0.9){
  mu<-ifelse(rbinom(n,1,0.5)>0.5,x-m*sd,x+m*sd)
  rnorm(n,mu,sqrt(1-m^2)*sd)
}