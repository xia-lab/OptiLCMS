#ifndef DSPEC_H
#define DSPEC_H

#include "DetectPeaks.h"

using namespace Rcpp;
using namespace std;

List DecoSpectra(int idx_pg, 
                 List spectra_eics,
                 NumericVector peak_ms1,
                 int num_scantime,
                 int idx_apex_eic,
                 NumericVector info_pk_ms1,
                 double peakwidth_min,
                 double snthr,
                 bool is_dec_smoothed);

// 
// List DecoSpectra(int idx_pg, 
//                  String nm_smp,
//                  List spectra_eics,
//                  NumericMatrix peak_ms1,
//                  int num_scantime,
//                  int idx_apex_eic,
//                  NumericVector info_pk_ms1,
//                  int peakwidth_min,
//                  int snthr,
//                  bool isFWHM,
//                  bool is_dec_all,
//                  bool is_dec_smoothed);

#endif