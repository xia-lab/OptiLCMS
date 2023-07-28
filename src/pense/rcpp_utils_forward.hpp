//
//  rcpp_utils_forward.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//
#ifndef RCPP_UTILS_FORWARD_HPP_
#define RCPP_UTILS_FORWARD_HPP_

#include <memory>

#include "nsoptim_forward.hpp"
#include "alias.hpp"

namespace pense {
//! PY Result Structure
//! Contains a list of initial estimates and the associated metrics.
template<typename T> struct PyResult;
}  // namespace pense

namespace Rcpp {
//! Wrap a pense::alias::FwdList (aka std::forward_list) into an R list.
//!
//! @param list forward list.
//! @return an R list.
template<typename T> SEXP wrap(const pense::alias::FwdList<T>& list);

//! Wrap an nsoptim::Metrics object into an R list.
//!
//! @param metrics the metrics object.
//! @return an R list.
template<> SEXP wrap(const nsoptim::Metrics& metrics);

//! Wrap a PyResult into an R list.
//!
//! @param py_result PyResult structure.
//! @return an R list.
template<typename T> SEXP wrap(const pense::PyResult<T>& py_result);

namespace traits {
//! Create an exporter for any `std::forward_list` where the elements are supported by Rcpp::as.
//! This does not use the same functionality as in `Rcpp/internal/Exporter.h` because forward_lists are better
//! created sequentially.
template <typename T>
class Exporter< std::forward_list<T> >;

//! Converter for an R-list to an EN penalty.
template<> class Exporter< nsoptim::EnPenalty >;
//! Converter for an R-list to a LASSO penalty.
template<> class Exporter< nsoptim::LassoPenalty >;
//! Converter for an R-list to a Ridge penalty.
template<> class Exporter< nsoptim::RidgePenalty >;

}  // namespace traits
}  // namespace Rcpp

#endif  // RCPP_UTILS_FORWARD_HPP_
