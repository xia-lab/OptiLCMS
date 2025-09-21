#include <Rcpp.h>
#include <vector>
#include <cfloat>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector R_mzClust_hclust_rcpp(NumericVector x, int num, NumericVector d, double eppm, double eabs) {
  int n = num;
  int m = n + 2;
  
  // prepare d index offsets (each row i has length n-i-1)
  std::vector<int> d_offset(n, 0);
  int pos = 0;
  for (int i = 0; i < n - 1; ++i) {
    d_offset[i] = pos;
    pos += (n - i - 1);
  }
  
  std::vector<double> means(n);
  double mrange[2];
  int overlimit = 0;
  
  std::array<int,2> minclust = {0,0};
  double mindst = DBL_MAX;
  
  // loop vars
  int i, j, z;
  
  std::vector<int> clust(n * m, 0);
  
  // init means
  for (i = 0; i < n; i++) means[i] = x[i];
  
  // init cluster matrix: first element is length (1), second is index
  for (z = 0; z < n; z++) {
    clust[z * m] = 1;
    clust[z * m + 1] = z;
  }
  
  for (z = 0; z < n - 1; z++) {
    // find minimal distance among active clusters
    for (i = 0; i < n; i++) {
      if (clust[i * m] <= 0) continue;
      for (j = i + 1; j < n; j++) {
        if (clust[j * m] <= 0) continue;
        double dij = d[d_offset[i] + (j - i - 1)];
        if (dij < mindst) {
          mindst = dij;
          minclust[0] = i;
          minclust[1] = j;
        }
      }
    }
    
    if (mindst == DBL_MAX) break; // no clusters left
    
    // calculate means for merged cluster (store in minclust[0])
    means[minclust[0]] = (means[minclust[0]] * clust[minclust[0] * m]
                            + means[minclust[1]] * clust[minclust[1] * m])
      / (double)(clust[minclust[0] * m] + clust[minclust[1] * m]);
      
      mrange[0] = means[minclust[0]] - means[minclust[0]] * eppm - eabs;
      mrange[1] = means[minclust[0]] + means[minclust[0]] * eppm + eabs;
      
      // find outliers in cluster 0
      for (i = 1; i <= clust[minclust[0] * m]; i++) {
        int idx = clust[minclust[0] * m + i];
        if (x[idx] < mrange[0] || x[idx] > mrange[1]) {
          overlimit = 1;
          break;
        }
      }
      // if still ok, check cluster 1
      if (!overlimit) {
        for (i = 1; i <= clust[minclust[1] * m]; i++) {
          int idx = clust[minclust[1] * m + i];
          if (x[idx] < mrange[0] || x[idx] > mrange[1]) {
            overlimit = 1;
            break;
          }
        }
      }
      
      if (overlimit) {
        clust[minclust[0] * m] *= -1;
        clust[minclust[1] * m] *= -1;
        overlimit = 0;
      } else {
        // merge: append elements of cluster1 to cluster0
        int len0 = clust[minclust[0] * m];
        int len1 = clust[minclust[1] * m];
        for (i = 1; i <= len1; i++) {
          clust[minclust[0] * m + len0 + i] = clust[minclust[1] * m + i];
        }
        clust[minclust[0] * m] += clust[minclust[1] * m];
        clust[minclust[1] * m] = 0;
        
        // update distance: max linkage between new cluster and others
        for (i = 0; i < n; i++) {
          if (clust[i * m] <= 0 || i == minclust[0]) continue;
          
          // helper to get dij between a and b where a < b
          auto get_d = [&](int a, int b)->double {
            if (a >= b) return DBL_MAX;
            return d[d_offset[a] + (b - a - 1)];
          };
          auto set_d = [&](int a, int b, double val){
            if (a >= b) return;
            d[d_offset[a] + (b - a - 1)] = val;
          };
          
          if (i < minclust[0]) {
            if (i < minclust[1]) {
              double v1 = get_d(i, minclust[0]);
              double v2 = get_d(i, minclust[1]);
              if (v1 < v2) set_d(i, minclust[0], v2);
            } else {
              double v1 = get_d(i, minclust[0]);
              double v2 = get_d(minclust[1], i);
              if (v1 < v2) set_d(i, minclust[0], v2);
            }
          } else {
            if (i < minclust[1]) {
              double v1 = get_d(minclust[0], i);
              double v2 = get_d(i, minclust[1]);
              if (v1 < v2) set_d(minclust[0], i, v2);
            } else {
              double v1 = get_d(minclust[0], i);
              double v2 = get_d(minclust[1], i);
              if (v1 < v2) set_d(minclust[0], i, v2);
            }
          }
        }
      }
      
      mindst = DBL_MAX;
  }
  
  // create output g
  IntegerVector g(n, NA_INTEGER);
  int gnum = 1;
  for (i = 0; i < n; i++) {
    int tmp = std::abs(clust[i * m]);
    if (tmp == 0) continue;
    for (j = 1; j <= tmp; j++) {
      int idx = clust[i * m + j];
      g[idx] = gnum;
    }
    gnum++;
  }
  
  return g;
}
