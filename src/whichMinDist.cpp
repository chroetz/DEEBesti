#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int whichMinDist(NumericMatrix x, NumericVector y) {
  int nrow = x.nrow(), ncol = x.ncol();
  int minI = 0;
  double minV = 0;
  for (int j = 0; j < ncol; ++j) {
    double v = x(0, j)-y(j);
    minV += v*v;
  }
  double z = 0;
  for (int i = 1; i < nrow; ++i) {
    z = 0;
    for (int j = 0; j < ncol; ++j) {
      double v = x(i, j)-y(j);
      z += v*v;
    }
    if (z < minV) {
      minV = z;
      minI = i;
    }
  }
  return minI+1;
}

// [[Rcpp::export]]
NumericVector distToVec(NumericMatrix x, NumericVector y) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  for (int i = 1; i < nrow; ++i) {
    double z = 0;
    for (int j = 0; j < ncol; ++j) {
      double v = x(i, j)-y(j);
      z += v*v;
    }
    out[i] = sqrt(z);
  }
  return out;
}
