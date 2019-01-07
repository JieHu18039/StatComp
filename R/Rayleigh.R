#' @title The Reyleigh distribution
#' @aliases dRayleigh pRayleigh qRayleigh rRayleigh
#' @usage dRayleigh(x,sigma)
#' pRayleigh(q,sigma,lower.tail=TRUE)
#' qRayleigh(p,sigma,lower.tail=TRUE)
#' rRayleigh(n,sigma)
#' @description Density, distribution function, quantile function and
#'  random generation for the Reyleigh distribution.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n the number of samples.
#' @param sigma the scale parameter of Reyleigh distribution.
#' @param lower.tail	 logical; if TRUE (default), probabilities are \code{P[X<=x]}, otherwise, \code{P[X > x]}.
#' @return dt gives the density, pt gives the distribution function, qt gives
#'  the quantile function, and rt generates random deviates.
#' @examples
#' \dontrun{
#' n<-10000
#' y<-rRayleigh(n,sigma=2)
#' hist(y,breaks=seq(0,ceiling(max(y)),0.1),freq = FALSE,
#' main = "Histogram of Reyleigh",xlab = " ")
#' x<-seq(0,8.5,0.05)
#' lines(x,dRayleigh(x,sigma=2),col="red")
#' qRayleigh(0.05,sigma=2)
#' pRayleigh(0.6405828,sigma = 2)
#' }
#' @export
dRayleigh<-function(x,sigma){
  x<-x*I(x>=0)
  x*exp(-x^2/(2*sigma^2))/(sigma^2)
}
#' @export
pRayleigh<-function(q,sigma,lower.tail=TRUE){
  q<-q*I(q>=0)
  if(lower.tail) 1-exp(-q^2/(2*sigma^2))
  else exp(-q^2/(2*sigma^2))
}

#' @export
qRayleigh<-function(p,sigma,lower.tail=TRUE){
  if(any(p<=0) || any(p>=1))
    stop("all 'p' must between 0 with 1")
  if(lower.tail) sigma*sqrt(-2*log(1-p))
  else sigma*sqrt(-2*log(p))
}

#' @export
rRayleigh<-function(n,sigma){
 sigma*sqrt(-2*log(runif(n)))
}

