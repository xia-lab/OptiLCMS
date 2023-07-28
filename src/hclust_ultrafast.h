#ifndef ultrafastcluste_H
#define ultrafastcluste_H

#include <vector>
#include <cmath> // for std::pow, std::sqrt
#include <cstddef> // for std::ptrdiff_t
#include <limits> // for std::numeric_limits<...>::infinity()
#include <algorithm> // for std::fill_n
#include <stdexcept> // for std::runtime_error
#include <string> // for std::string
#include <cfloat> // also for DBL_MAX, DBL_MIN

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <RcppArmadillo.h>

using namespace Rcpp;


// Older versions of Microsoft Visual Studio do not have the fenv header.
#ifdef _MSC_VER
#if (_MSC_VER == 1500 || _MSC_VER == 1600)
#define NO_INCLUDE_FENV
#endif
#endif
// NaN detection via fenv might not work on systems with software
// floating-point emulation (bug report for Debian armel).
#ifdef __SOFTFP__
#define NO_INCLUDE_FENV
#endif
#ifdef NO_INCLUDE_FENV
#pragma message("Do not use fenv header.")
#else
#include <fenv.h>
#endif


#ifndef DBL_MANT_DIG
#error The constant DBL_MANT_DIG could not be defined.
#endif
#define T_FLOAT_MANT_DIG DBL_MANT_DIG

#ifndef LONG_MAX
#include <climits>
#endif
#ifndef LONG_MAX
#error The constant LONG_MAX could not be defined.
#endif
#ifndef INT_MAX
#error The constant INT_MAX could not be defined.
#endif

#ifndef INT32_MAX
#ifdef _MSC_VER
#if _MSC_VER >= 1600
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#else
typedef __int32 int_fast32_t;
typedef __int64 int64_t;
#endif
#else
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#endif
#endif

#define FILL_N std::fill_n
#ifdef _MSC_VER
#if _MSC_VER < 1600
#undef FILL_N
#define FILL_N stdext::unchecked_fill_n
#endif
#endif

// Suppress warnings about (potentially) uninitialized variables.
#ifdef _MSC_VER
#pragma warning (disable:4700)
#endif

#ifndef HAVE_DIAGNOSTIC
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 6))
#define HAVE_DIAGNOSTIC 1
#endif
#endif

#ifndef HAVE_VISIBILITY
#if __GNUC__ >= 4
#define HAVE_VISIBILITY 1
#endif
#endif

/* Since the public interface is given by the Python respectively R interface,
 * we do not want other symbols than the interface initalization routines to be
 * visible in the shared object file. The "visibility" switch is a GCC concept.
 * Hiding symbols keeps the relocation table small and decreases startup time.
 * See http://gcc.gnu.org/wiki/Visibility
 */
#if HAVE_VISIBILITY
#pragma GCC visibility push(hidden)
#endif

typedef int_fast32_t t_index;
#ifndef INT32_MAX
#define MAX_INDEX 0x7fffffffL
#else
#define MAX_INDEX INT32_MAX
#endif
#if (LONG_MAX < MAX_INDEX)
#error The integer format "t_index" must not have a greater range than "long int".
#endif
#if (INT_MAX > MAX_INDEX)
#error The integer format "int" must not have a greater range than "t_index".
#endif
typedef double t_float;

NumericVector auto_hclust(NumericVector x0);

NumericVector auto_hclust_median(NumericVector x0);
  
NumericVector ultra_hclust(NumericVector x0, int n_clusts);

NumericVector matrix_hclust(NumericMatrix data_mtx);

#endif
