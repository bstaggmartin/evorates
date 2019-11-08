#'@title Function to sample the indices of vector as "evenly" as possible by halving the vector into finer and finer increments
#'@name sample.evenly
#'@author Bruce
#'
#'@description I'm hoping this will make simulated BM bridges more noisy, a la true BM, but I still need to test this; note that this
#' function is, unfortunately, much slower than sample, but I have yet to think of a more efficient way to do this...
#' 
#'@param vec numeric vector 
#'@return vector of evenly sampled indices

sample.evenly<-function(vec){
  samp<-c(0,round(length(vec)/2),length(vec)+1)
  while(!all(1:length(vec)%in%samp)){
    space<-mean(sort(samp,decreasing=T)[-length(samp)]-sort(samp,decreasing=T)[-1])
    samp<-unique(c(samp,round(sort(samp)[-length(samp)]+space/2)))
  }
  vec[samp[-c(1,3)]]
}