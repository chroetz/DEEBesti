#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int whichMinDist(NumericMatrix x, NumericVector y) {
  int nrow = x.nrow(), ncol = x.ncol();
  double minDist = INFINITY;
  double minIdx;
  double dst;
  double v;
  for (int i = 0; i < nrow; ++i) {
    dst = 0;
    for (int j = 0; j < ncol; ++j) {
      v = x(i, j) - y(j);
      dst += v*v;
    }
    if (dst < minDist) {
      minDist = dst;
      minIdx = i;
    }
  }
  return minIdx+1;
}

// [[Rcpp::export]]
NumericVector distToVec(NumericMatrix x, NumericVector y) {
  int n = x.nrow(), d = x.ncol();
  NumericVector out(n);
  double dst;
  double v;
  for (int i = 0; i < n; ++i) {
    dst = 0;
    for (int j = 0; j < d; ++j) {
      v = x(i, j) - y(j);
      dst += v*v;
    }
    out[i] = sqrt(dst);
  }
  return out;
}
