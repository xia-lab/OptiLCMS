//
//  r_utilities.cc
//  pense
//
//  Created by David Kepplinger on 2019-05-12.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "r_utilities.hpp"

#include <cmath>

namespace pense {
namespace r_interface {
//! Approximate value matching.
//!
//! Returns a vector of 1-based positions of the (first) matches of `x` in `table`.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @return a vector the same lenght of `x` with integers giving the position in `table` of the first match
//!         if there is a match, or `NA_integer_` otherwise.
SEXP ApproximateMatch(SEXP r_x, SEXP r_table, SEXP r_eps) noexcept {
  const R_xlen_t len_x = Rf_xlength(r_x);
  const int len_table = Rf_length(r_table);
  SEXP r_matches = PROTECT(Rf_allocVector(INTSXP, len_x));
  int* matches = INTEGER(r_matches);
  double const * x = REAL(r_x);
  double const * table = REAL(r_table);
  const double eps = *REAL(r_eps);

  for (R_xlen_t i = 0; i < len_x; ++i) {
    matches[i] = NA_INTEGER;
    for (int j = 0; j < len_table; ++j) {
      if (std::abs(x[i] - table[j]) < eps) {
        matches[i] = j + 1;
        break;
      }
    }
  }

  UNPROTECT(1);
  return r_matches;
}

}  // namespace r_interface
}  // namespace pense
