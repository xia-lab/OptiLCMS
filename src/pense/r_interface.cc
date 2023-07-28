//
//  r_interface.cc
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifdef HAVE_RCPP
#include <R_ext/Rdynload.h>

#include "rcpp_integration.hpp"
#include "r_en_regression.hpp"
#include "r_pense_regression.hpp"
#include "r_mesten_regression.hpp"
#include "r_robust_utils.hpp"
#include "r_enpy.hpp"
#include "r_utilities.hpp"

extern "C" SEXP run_testthat_tests() noexcept;

//! R initialzing function (must be in the global namespace).
extern "C" void R_init_pense(DllInfo *dll) noexcept;

using namespace pense::r_interface;

namespace {
//! Exported methods

}  // namespace



#endif  // HAVE_RCPP
