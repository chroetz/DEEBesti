#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix expKernelMatrix(NumericMatrix state, double bandwidth, double regulation) {
  int n = state.nrow();
  int d = state.ncol();
  NumericMatrix out(n,n);
  double bwSqr = bandwidth*bandwidth;
  double dst, v;

  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      dst = 0;
      for (int k = 0; k < d; ++k) {
        v = state(i,k)-state(j,k);
        dst += v*v;
      }
      dst = exp(-0.5*dst/bwSqr);
      out(i, j) = dst;
      out(j, i) = dst;
    }
  }
  out.fill_diag(1 + regulation);

  return out;
}

// [[Rcpp::export]]
NumericVector expKernelVectorFromDistSqr(NumericVector distSqr, double bandwidth) {
  int n = distSqr.length();
  NumericVector out(n);
  double bwSqr = bandwidth*bandwidth;
  for (int i = 0; i < n; ++i) {
    out(i) = exp(-0.5*distSqr(i)/bwSqr);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix expKernelMatrix1D(NumericVector x, double bandwidth, double regulation) {
  int n = x.length();
  NumericMatrix out(n,n);
  double bwSqr = bandwidth*bandwidth;
  double dst, v;

  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      v = x(i)-x(j);
      dst = exp(-0.5*v*v/bwSqr);
      out(i, j) = dst;
      out(j, i) = dst;
    }
  }
  out.fill_diag(1 + regulation);

  return out;
}

// [[Rcpp::export]]
NumericMatrix expKernelVectors1D(NumericVector x, NumericVector xout, double bandwidth) {
  int n = x.length();
  int m = xout.length();
  NumericMatrix out(n, m);
  double bwSqr = bandwidth*bandwidth;
  double v;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      v = x(i)-xout(j);
      out(i, j) = exp(-0.5*v*v/bwSqr);
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix expKernelDerivVectors1D(NumericVector x, NumericVector xout, double bandwidth) {
  int n = x.length();
  int m = xout.length();
  NumericMatrix out(n, m);
  double bwSqr = bandwidth*bandwidth;
  double v;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      v = x(i)-xout(j);
      out(i, j) = v/bwSqr*exp(-0.5*v*v/bwSqr);
    }
  }

  return out;
}
