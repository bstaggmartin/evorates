#Function to modify a vector of colors via mixing and/or altering transparency
##I find this is helpful in a number of visualization situations; use alph=NA to not modify original transparencies of colors
#' Modify vectors of colors
#'
#' This function alters a vector of colors by modifying their transparency and/or mixing them with other colors.
#'
#' @param x A vector of colors specified either by name, hexadecimal code, or numbers (corresponding to colors defined in
#' \code{palette()}). RGB values can be coerced to hexadecimal codes via the \code{rgb()}.
#' @param alpha A numeric vector of alpha (transparency) values, with 0 corresponding to complete transparency and 1 to complete
#' opaqueness. Use NA to not modify transparency of colors in \code{x}.
#' @param mix Another vector of colors to mix with the colors in \code{x}. Use NA to prevent any mixing.
#' @param wgt A numeric vector of weights to control how colors in \code{mix} are combined with colors in \code{x}, with 0 
#' corresponding to no mixing and 1 to complete replacement.
#' @return a character vector of color values in hexadecimal code.
#' @examples
#' colors<-c('red','green','blue')
#' alter.cols(colors,alpha=c(0,1,0.5),mix=c('white','black','gray'),wgt=c(0,0.5,1))
#' #recycling behavior--recycled to maximum possible length
#' alter.cols(colors,alpha=0.2,mod.val=c('black','white'),wgt=c(0.1,0.2,0.8,0.9))
#' 
#' @export
alter.cols<-function(x,alpha=NA,mix=NA,wgt=0.5){
  lens<-lengths(list(x,alpha,mix,wgt))
  #prepare mix
  if(any(!is.na(mix))){
    mix.len<-max(lens[-c(1,2)])
    out.len<-max(lens)
    tmp<-rep(mix,length.out=mix.len)
    wgt[wgt<0]<-0
    wgt[wgt>1]<-1
    wgt<-rep(wgt,length.out=mix.len)
    wgt[is.na(tmp)]<-0
    mix<-col2rgb(mix,alpha=TRUE)[,rep(seq_len(lens[3]),length.out=mix.len),drop=FALSE]
    mix<-sweep(mix,2,wgt,'*')
    mix<-mix[,rep(seq_len(mix.len),length.out=out.len),drop=FALSE]
  }else{
    out.len<-max(lens[c(1,2)])
    mix<-NULL
  }
  #prepare x
  x<-col2rgb(x,alpha=TRUE)[,rep(seq_len(lens[1]),length.out=out.len),drop=FALSE]
  if(!is.null(mix)){
    x<-sweep(x,2,1-wgt,'*')
    x<-x+mix
  }
  #prepare alpha
  if(any(!is.na(alpha))){
    pre<-x[4,,drop=FALSE]
    alpha[alpha<0]<-0
    alpha[alpha>1]<-1
    x[4,]<-255*alpha
    nas<-is.na(x[4,,drop=FALSE])
    x[4,nas]<-pre[nas]
  }
  do.call(rgb,c(asplit(x,1),maxColorValue=255))
}