// Li, Y., Kind, T., Folz, J. et al. Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. Nat Methods 18, 1524–1531 (2021). https://doi.org/10.1038/s41592-021-01331-z


#define SPEC_TYPE
typedef double float_spec;

#include <Rcpp.h>
#include "CleanSpectrum.h"
#include "SpectralEntropy.h"

using namespace Rcpp;


#ifndef SPEC_TYPEX
#define SPEC_TYPEX


Rcpp::NumericVector convert_matrix_to_vector(const Rcpp::NumericMatrix peaks);

Rcpp::NumericMatrix convert_vector_to_matrix(const Rcpp::NumericVector peaks, int nrow);


double r_calculate_unweighted_entropy_similarity(const Rcpp::NumericMatrix peaks_a,
                                                 const Rcpp::NumericMatrix peaks_b,
                                                 float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
                                                 bool clean_spectra,
                                                 float min_mz, float max_mz,
                                                 float noise_threshold,
                                                 int max_peak_num);


double r_calculate_entropy_similarity(const Rcpp::NumericMatrix peaks_a,
                                      const Rcpp::NumericMatrix peaks_b,
                                      float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
                                      bool clean_spectra,
                                      float min_mz, float max_mz,
                                      float noise_threshold,
                                      int max_peak_num);

#endif