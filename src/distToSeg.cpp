#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double distSqrToSeg(NumericVector u, NumericVector v0, NumericVector v1) {
  int d = u.length();
  double* del = new double[d];
  double dist;
  double lenSqr;
  double ds;
  double prod;
  double t;

  lenSqr = 0;
  for (int j = 0; j < d; ++j) {
    ds = v1[j] - v0[j];
    del[j] = ds;
    lenSqr += ds*ds;
  }
  if (lenSqr < 1e-15) {
    dist = 0;
    for (int j = 0; j < d; ++j) {
      ds = u[j] - v0[j];
      dist += ds*ds;
    }
    return dist;
  }

  prod = 0;
  for (int j = 0; j < d; ++j) prod += (u[j] - v0[j])*del[j];

  t = prod / lenSqr;
  if (t > 1) t = 1;
  if (t < 0) t = 0;

  dist = 0;
  for (int j = 0; j < d; ++j) {
    ds = u[j] - (v0[j] + t * del[j]);
    dist += ds*ds;
  }

  delete[] del;

  return dist;
}

// [[Rcpp::export]]
NumericVector distSqrToPwLin(NumericMatrix path, NumericVector query) {
  int d = query.length();
  int n = path.nrow()-1;

  double dist;
  double lenSqr;
  double ds;
  double prod;
  double t;
  double* del = new double[d];

  NumericVector dsts(n);

  for (int i = 0; i < n; ++i) {

    lenSqr = 0;
    for (int j = 0; j < d; ++j) {
      ds = path(i+1,j) - path(i,j);
      del[j] = ds;
      lenSqr += ds * ds;
    }
    if (lenSqr < 1e-15) {
      dist = 0;
      for (int j = 0; j < d; ++j) {
        ds = query[j] - path(i,j);
        dist += ds*ds;
      }
      dsts[i] = dist;
      continue;
    }

    prod = 0;
    for (int j = 0; j < d; ++j) prod += (query[j] - path(i,j))*del[j];

    t = prod / lenSqr;
    if (t > 1) t = 1;
    if (t < 0) t = 0;

    dist = 0;
    for (int j = 0; j < d; ++j) {
      ds = query[j] - (path(i,j) + t * del[j]);
      dist += ds*ds;
    }
    dsts[i] = dist;

  }

  delete[] del;

  return dsts;
}

// expects points of the same trajectory to be next to each other, i.e.,
// when id changes from one index to the next, a completely new trajectory starts
// [[Rcpp::export]]
double whichMinDistToPwLin(NumericMatrix path, NumericVector id, NumericVector query) {
  int d = query.length();
  int n = path.nrow()-1;

  double dist;
  double lenSqr;
  double ds;
  double prod;
  double t;
  double* p = new double[d];
  double* del = new double[d];

  double minDist = INFINITY;
  double minIdx;

  for (int i = 0; i < n; ++i) {

    if (id[i] != id[i+1]) continue;

    lenSqr = 0;
    for (int j = 0; j < d; ++j) {
      ds = path(i+1,j) - path(i,j);
      del[j] = ds;
      lenSqr += ds * ds;
    }
    if (lenSqr < 1e-15) {
      dist = 0;
      for (int j = 0; j < d; ++j) {
        ds = query[j] - path(i,j);
        dist += ds*ds;
      }
      if (dist < minDist) {
        minDist = dist;
        minIdx = i + 0.5;
      }
      continue;
    }

    prod = 0;
    for (int j = 0; j < d; ++j) prod += (query[j] - path(i,j))*del[j];

    t = prod / lenSqr;
    if (t > 1) t = 1;
    if (t < 0) t = 0;

    for (int j = 0; j < d; ++j) p[j] = path(i,j) + t * del[j];

    dist = 0;
    for (int j = 0; j < d; ++j) {
      ds = query[j] - p[j];
      dist += ds*ds;
    }

    if (dist < minDist) {
      minDist = dist;
      minIdx = i+t;
    }

  }

  delete[] del;
  delete[] p;

  return minIdx + 1;
}

