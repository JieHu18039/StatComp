#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{system.time} to compare the performance of C functions (\code{fibC} and \code{fibR}).
#' @importFrom stats rnorm runif dcauchy
#' @importFrom Rcpp evalCpp
#' @useDynLib StatComp18039
#' @examples
#' \dontrun{
#' system.time(fibR(35))
#' system.time(fibC(35))
#' }
NULL
