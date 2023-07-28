#ifndef UTILITY_H
#define UTILITY_H

#include <RcppArmadillo.h>
#include <iostream>
#include <math.h>
#include "lowess.h"

using namespace Rcpp;
using namespace std;

NumericVector SmoothLoess(NumericMatrix &eic, double span);

vector<double> lowessCpp(IntegerVector x, NumericVector y, double spanVal);

bool checkContinuousPtsAboveThr(NumericVector v, int iStart, double num, double thr, int nSkipMax);

NumericVector getContinuousPtsAboveThrIdx(NumericVector v, int iStart, int num, double thr, int nSkipMax);

IntegerVector whichTrue(LogicalVector vecValues);

int whichTrue1(LogicalVector vecValues);

IntegerVector GetRoi(NumericVector is_roi, int idx_apex_eic);

double EstimateChromNoise(NumericVector &x, double trim, int min_pts);

NumericVector GetLocalNoiseEstimate(NumericVector d, IntegerVector idx_fr_roi, double noiserange_min, double noiserange_max, 
                                    int Nscantime, double threshold, int num);

NumericVector CalculateBL(NumericVector d, IntegerVector drange, double threshold, int num, int n_skip_max, double noiserange_min);

IntegerVector FindLocalMax(NumericVector x, int m, double v);

IntegerVector FindLocalMin(NumericVector x, int m);

NumericMatrix mergeEIC(NumericMatrix x, NumericMatrix y);

double Gauss(int x, int h, int mu, int sigma);

double cor_fast (NumericVector x, NumericVector y);

double GetDistantP(NumericMatrix peak1, NumericMatrix peak2);

// double funOptimc(NumericVector x, NumericVector M, NumericMatrix S);
// 
// NumericVector optim_real(NumericMatrix mpk_mtx, NumericVector eic);
// 
// NumericVector optim_ultra(NumericMatrix mpk_mtx, NumericVector eic);
// 

// List extractEIC(List specExp, NumericMatrix mzRange, NumericVector mz);

#endif