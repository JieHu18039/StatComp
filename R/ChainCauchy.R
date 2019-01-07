#' @title The MC chain of Cauchy distribution
#' @description Give a MC chain of Cauchy distribution
#' with location parameter \code{location} and scale parameter \code{scale}.
#' @param n the length of the chain.
#' @param sigma The standard deviation of the proposed distribution.
#' @param x0	The initial value of the chain.
#' @param location,scale  location and scale parameters.
#' @details If \code{location} or \code{scale} are not specified, they assume the
#' default values of 0 and 1 respectively.
#'
#' The Cauchy distribution with location \code{l} and scale \code{s} has density
#' \code{f(x)=1/(pi s(1+((x-l)/s)^2))} for all \code{x}.
#' @return \code{x}  return the chain.
#' \code{k}  Returns the number of rejections.
#' @examples
#' \dontrun{
#' n<-10000
#' cauchy<-cauchy.chain(n,sigma=1,x0=0,location=0,scale=1)
#' cauchy$k/n
#' hist(cauchy$x[1001:n],freq=FALSE,main="cauchy密度图",breaks=60)
#' curve(dcauchy(x,0,1),add=TRUE)
#' }
#' @export
cauchy.chain<-function(n,sigma,x0,location,scale){
  x<-NULL
  x[1]<-x0
  e=runif(n)
  k<-0
  for(i in 2:n){
    y<-rnorm(1,x[i-1],sigma)
    if(e[i]<=(dcauchy(y,location,scale)/dcauchy(x[i-1],location,scale)))
    {x[i]<-y}
    else
    {x[i]<-x[i-1]
    k<-k+1}
  }
  return(list(x=x,k=k))
}
