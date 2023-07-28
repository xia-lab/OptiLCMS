#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector DescendZero(NumericVector yvals, int numin, int istart) {
  IntegerVector res(2);
  int ilower = 0, iupper=0;
  int i;
  for (i = istart; i >= 0; i--) {
    if (yvals[i] < 0) {
      break;
    }
  }
  ilower = i + 1;
  
  for (i = istart; i < numin; i++) {
    if (yvals[i] < 0) {
      break;
    }
  }
  iupper = i - 1;
  
  res[0] = ilower;
  res[1] = iupper;
  //std::cout << "ilower=>" << ilower << std::endl;
  //std::cout << "iupper=>" << iupper << std::endl;
  return res;
}


// [[Rcpp::export]]
NumericVector ColMax(const NumericVector& inval, const int& n, const int& dn) {
  
  NumericVector out(dn);
  
  for (int i = 0; i < dn; i++) {
    //std::cout << "running into here ==> " << i << std::endl;
    out[i] = inval[n*i];
    for (int j = 1; j < n; j++) {
      if (inval[n*i + j] > out[i]) {
        out[i] = inval[n*i + j];
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
IntegerVector DescendMin(NumericVector yvals, int numin, int istart) {
  IntegerVector res(2);
  int i;
  
  for (i = istart; i > 0; i--) {
    if (yvals[i-1] >= yvals[i]) {
      break;
    }
  }
  res[0] = i;
  
  for (i = istart; i < numin-1; i++) {
    if (yvals[i+1] >= yvals[i]) {
      break;
    }
  }
  res[1] = i;
  return res;
}


// [[Rcpp::export]]
IntegerVector WhichColMax(const NumericVector& inval, const int& n, const int& dn) {
  IntegerVector out(dn);
  for (int i = 0; i < dn; i++) {
    out[i] = n*i;
    for (int j = 1; j < n; j++) {
      if (inval[n*i + j] > inval[out[i]]) {
        out[i] = n*i + j;
      }
    }
  }
  for (int i = 0; i < dn; i++) {
    out[i]++;
  }
  return(out);
}


// [[Rcpp::export]]
IntegerVector FindEqualGreaterM(const NumericVector& inval, const int& size,
                                       const NumericVector& values, const int& valsize) {
  
  IntegerVector index(valsize);
  int idx = 0;
  
  for (int i = 0; i < valsize; i++) {
    while (idx < size && inval[idx] < values[i]) {
      idx++;
    }
    index[i] = idx;
  }
  return index;
}


// [[Rcpp::export]]
IntegerVector RectUnique(const NumericVector& m, const IntegerVector& order,
                                const int& nrow, const int& ncol,
                                const double& xdiff, const double& ydiff) {
  
  IntegerVector keep(nrow);
  
  int x1 = 0, x2 = nrow, y1 = nrow*2, y2 = nrow*3;
  
  for (int i = 0; i < nrow; i++) {
    int io = order[i];
    keep[io] = 1;
    for (int j = 0; j < i; j++) {
      int jo = order[j];
      if (keep[jo] &&
          !(m[x1+io] - m[x2+jo] > xdiff || m[x1+jo] - m[x2+io] > xdiff ||
          m[y1+io] - m[y2+jo] > ydiff || m[y1+jo] - m[y2+io] > ydiff)) {
        keep[io] = 0;
        break;
      }
    }
  }
  return keep;
}


// [[Rcpp::export]]
int continuousPtsAboveThreshold(NumericVector x, int istart, int numin, double threshold, int num) {
  
  int cnt = 0, n;
  
  for (int i = istart; i < numin; i++) {
    if (x[i] > threshold) cnt++;
    else cnt = 0;
    if (cnt >= num) {
      n = cnt;
      return n;
    }
  }
  return 0;
}



// [[Rcpp::export]]
IntegerVector continuousPtsAboveThresholdIdx(NumericVector x, int istart, int numin, double threshold, int num) {
  
  int cnt = 0;
  int stidx = 0;
  int enidx = 0;
  IntegerVector n(numin);
  
  for (int i = istart; i < numin; i++) {
    
    if (x[i] > threshold) {
      cnt++;
      if (cnt == 1) stidx=i;
      else enidx=i;
    }  else cnt = 0;
    
    if ((cnt==0 || i == numin - 1) && (enidx - stidx + 1) >= num) {
      for (int j = stidx; j<= enidx; j++) {
        n[j] = 1;
      }
      stidx = 0;
      enidx = 0;
    }
  }
  return n;
  
}


