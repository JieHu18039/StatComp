#include <Rcpp.h>
using namespace Rcpp;

//' @title Fibonacci sequence
//' @description Compute item \code{n} of the Fibonacci sequence using c++
//' @param n  item \code{n}
//' @return The value of item \code{n} of the Fibonacci sequence
//' @examples
//' \dontrun{
//' fibC(29)
//' }
//' @export
// [[Rcpp::export]]
int fibC(int n)
{
  if(n==1||n==2) return 1;
  return fibC(n-1)+fibC(n-2);
}

