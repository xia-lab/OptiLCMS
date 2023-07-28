#ifndef DDAUTILITY_H
#define DDAUTILITY_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

List MSCentroidsGrouping(NumericVector allMZs);

NumericMatrix row_erase (NumericMatrix& x, IntegerVector& rowID);

NumericMatrix cosineSimilarity(NumericMatrix Xr);

float dot(NumericVector a, NumericVector b, bool norm=true);

double entropy(NumericMatrix peaks_a, NumericMatrix peaks_b);

NumericMatrix ms2peak_parse(string text);

double spectrumSimilarity(NumericMatrix mtx1, NumericMatrix mtx2, double ppm_ms2);

double entropySimilarity(NumericMatrix mtx1, NumericMatrix mtx2, double ppm_ms2);
double get_mass_sodium();

double get_mass_potassium();

vector<int> parse_formula(std::string formula_txt);

double neutral_loss_similarity(NumericVector exp_mzs, NumericVector exp_ints, 
                               NumericVector ref_mzs, NumericVector ref_ints,
                               double ppm2);

NumericVector unique_num(NumericVector x);

#endif