#Function to make a modify a vector of colors to have a certain level of transparency (or be darker/lighter)
##I find this is helpful in a number of visualization situations; use alph=NA to not modify original transparencies of colors
#' Modify the transparency and shade of colors
#'
#' This function takes a vector of colors (x) as input and returns a vector of colors of the same length. The function can be used to
#' either make the colors a certain transparency (alph) or make the colors lighter or darker (mod.val). alph and mod.val are recycled
#' or truncated to be the same length as x.
#'
#' @param x Vector of color values, which can either be a character vector of color names and/or hexadecimal codes, or a numeric
#' vector. In the latter case, the vecotr will be coerced to integers taken to correspond to colors currently defined in palette().
#' @param alpha Vector of alpha (transparency) values (numeric), with 0 corresponding to complete transparency and 1 to complete
#' opaqueness. Use NA to not modify transparency of original color vector. alph must consist of NA's and/or values between 0 and 1.
#' @param mod.val Vector of values (numeric) by which to brighten/darken colors, with positive and negative numbers corresponding to
#' brightening and darkening, respectively. In this context, rgb values range from 0 to 1, and any modified values exceeding
#' these boundaries will be rounded up or down, respectively.
#' @return a vector of color values in hexadecimal code (character)
#' @examples
#' colors<-c('red','green','blue')
#' alter.cols(colors,alph=c(0,1,0.5),mod.val=c(0,-0.5,-0.75))
#' #recycling behavior 
#' alter.cols(colors,alph=0.2,mod.val=-0.5)
#' 
#' @export
alter.cols<-function(x,alpha=NA,mod.val=0){
  if(is.character(x)){
    cols<-col2rgb(x,alpha=T)/255
  }else if(all(is.numeric(x)&x>0)){
    cols<-col2rgb(as.integer(x),alpha=T)/255
  }else{
    stop('x is not a recognizable vector of colors')
  }
  if(!all(is.na(alpha)|(alpha<=1&alpha>=0))){
    stop('alpha is outside permitted range (0-1 or NA)')
  }
  cols[1:3,]<-cols[1:3,]+mod.val
  cols[1:3,]<-ifelse(cols[1:3,]<0,0,ifelse(cols[1:3,]>1,1,cols[1:3,]))
  alpha<-rep(alpha,length.out=length(x))
  alpha<-ifelse(is.na(alpha),cols[4,],alpha)
  cols[4,]<-alpha
  mapply(rgb,red=cols[1,],green=cols[2,],blue=cols[3,],alpha=cols[4,])
}