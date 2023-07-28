#ifndef DPEAKS_H
#define DPEAKS_H

#include "utilities.h"

using namespace Rcpp;
using namespace std;

List DetectPeaks(NumericMatrix eic, double peakwidth_min, int num_scantime, int idx_apex_eic,
                 double snthr, bool is_smooth, int n_skip_max);

#endif