#' @title Fibonacci sequence
#' @description Compute item \code{n} of the Fibonacci sequence using R
#' @param n item \code{n}
#' @return The value of item \code{n} of the Fibonacci sequence
#' @examples
#' \dontrun{
#' fibR(29)
#' }
#' @export
fibR<- function(n){
  if(n==1||n==2) return(1)
  return(fibR(n-1)+fibR(n-2))
}
