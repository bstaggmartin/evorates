#'@title Draw from a bactrian distribution (see Yang et al. 2013 for details)
#'@name rbac
#'@author Bruce
#'
#'@description works by drawing from a binomial distribution, then using that draw as an indicator to draw from either a "left" or "right" normal
#'distribution
#'
#'@param n number of observations 
#'@param x 
#'@param sd standard deviation
#'@param m 

rbac<-function(n,x=0,sd=1,m=0.9){
  mu<-ifelse(rbinom(n,1,0.5)>0.5,x-m*sd,x+m*sd)
  rnorm(n,mu,sqrt(1-m^2)*sd)
}