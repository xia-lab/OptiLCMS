#ifndef OPTIM_H
#define OPTIM_H

#include "optim_src.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace roptim;

double optim_ultra(NumericMatrix mpkmtx, NumericVector vec_eic, int main_idx);


#endif