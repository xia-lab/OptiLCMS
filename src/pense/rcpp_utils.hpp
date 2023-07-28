//
//  rcpp_utils.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef RCPP_UTILS_HPP_
#define RCPP_UTILS_HPP_

#include <string>
#include <memory>
#include <type_traits>

#include "nsoptim.hpp"
#include "constants.hpp"
#include "enpy_types.hpp"

namespace pense {
//! Get an item from the list or use the fallback if the item does not exist.
//!
//! @param list the R list object to extract the data from.
//! @param name the name of the object to extract.
//! @param fallback the fallback value if the list does not contain an element with name `name`.
//! @return either the element from the list with name `name`, or `fallback` if the item does not exist.
template<typename T>
T GetFallback(const Rcpp::List& list, const std::string& name, const T fallback) noexcept {
  try {
    // Check if the element exists to avoid unnecessary exceptions.
    // An unsupported cast to `T` still triggers an exception, but this shouldn't happen very often!
    if (list.containsElementNamed(name.c_str())) {
      return Rcpp::as<T>(list[name]);
    }
  } catch (...) {}
  return fallback;
}

//! enum-specific overload
template<>
inline pense::EnAlgorithm GetFallback<pense::EnAlgorithm>(const Rcpp::List& list, const std::string& name,
                                                          const pense::EnAlgorithm fallback) noexcept {
  try {
    // Check if the element exists to avoid unnecessary exceptions.
    // An unsupported cast to `T` still triggers an exception, but this shouldn't happen very often!
    if (list.containsElementNamed(name.c_str())) {
      return static_cast<pense::EnAlgorithm>(Rcpp::as<int>(list[name]));
    }
  } catch (...) {}
  return fallback;
}
//! enum-specific overload
template<>
inline pense::PenseAlgorithm GetFallback<pense::PenseAlgorithm>(const Rcpp::List& list, const std::string& name,
                                                                const pense::PenseAlgorithm fallback) noexcept {
  try {
    // Check if the element exists to avoid unnecessary exceptions.
    // An unsupported cast to `T` still triggers an exception, but this shouldn't happen very often!
    if (list.containsElementNamed(name.c_str())) {
      return static_cast<pense::PenseAlgorithm>(Rcpp::as<int>(list[name]));
    }
  } catch (...) {}
  return fallback;
}

//! enum-specific overload
template<>
inline pense::RhoFunctionType GetFallback<pense::RhoFunctionType>(const Rcpp::List& list, const std::string& name,
                                                                  const pense::RhoFunctionType fallback) noexcept {
  try {
    // Check if the element exists to avoid unnecessary exceptions.
    // An unsupported cast to `T` still triggers an exception, but this shouldn't happen very often!
    if (list.containsElementNamed(name.c_str())) {
      return static_cast<pense::RhoFunctionType>(Rcpp::as<int>(list[name]));
    }
  } catch (...) {}
  return fallback;
}

//! enum-specific overload
template<>
inline nsoptim::MMConfiguration::TighteningType GetFallback<nsoptim::MMConfiguration::TighteningType>(
  const Rcpp::List& list, const std::string& name, const nsoptim::MMConfiguration::TighteningType fallback) noexcept {
  try {
    // Check if the element exists to avoid unnecessary exceptions.
    // An unsupported cast to `T` still triggers an exception, but this shouldn't happen very often!
    if (list.containsElementNamed(name.c_str())) {
      return static_cast<nsoptim::MMConfiguration::TighteningType>(Rcpp::as<int>(list[name]));
    }
  } catch (...) {}
  return fallback;
}

//! Wrap an Optimum for any EN-type penalty function into an R list.
//!
//! @param optimium the Optimum object.
//! @return the optimum as Rcpp::List.
template <typename T>
Rcpp::List WrapOptimum(const T& optimum) {
  using Rcpp::Named;
  return Rcpp::List::create(Named("alpha") = optimum.penalty.alpha(),
                            Named("lambda") = optimum.penalty.lambda(),
                            Named("objf_value") = optimum.objf_value,
                            Named("statuscode") = static_cast<int>(optimum.status),
                            Named("status") = optimum.message,
                            Named("intercept") = optimum.coefs.intercept,
                            Named("beta") = optimum.coefs.beta);
}

//! Wrap a list of optima for any EN-type penalty function into an R list.
//! @param optima list of Optimum objects.
//! @return the optimum as Rcpp::List.
template <typename T>
Rcpp::List WrapOptima(const std::forward_list<T>& optima) {
  Rcpp::List output_list;
  for (auto&& optimum : optima) {
    output_list.push_back(WrapOptimum(optimum));
  }
  return output_list;
}
}  // namespace pense

namespace Rcpp {
//! Wrap a alias::FwdList (aka std::forward_list) into an R list.
//!
//! @param list forward list.
//! @return an R list.
template<typename T> SEXP wrap(const pense::alias::FwdList<T>& list) {
  List r_list;
  for (auto&& element : list) {
    r_list.push_back(wrap(element));
  }
  return r_list;
}

//! Wrap a PyResult into an R list.
//!
//! @param py_result PyResult structure.
//! @return an R list.
template<typename T>
SEXP wrap(const pense::PyResult<T>& py_result) {
  using Rcpp::Named;
  return Rcpp::List::create(Named("metrics") = py_result.metrics,
                            Named("estimates") = pense::WrapOptima(py_result.initial_estimates));
}

//! Wrap an nsoptim::Metrics object into an R list.
//!
//! @param metrics the metrics object.
//! @return an R list.
template<>
inline SEXP wrap(const nsoptim::Metrics& metrics) {
  List r_list;
  List sub_metrics;
  r_list["name"] = metrics.name();

  for (auto&& metric : metrics.DoubleMetrics()) {
    r_list[metric.name] = metric.value;
  }
  for (auto&& metric : metrics.IntegerMetrics()) {
    r_list[metric.name] = metric.value;
  }
  for (auto&& metric : metrics.StringMetrics()) {
    r_list[metric.name] = metric.value;
  }
  for (auto&& sub_metric : metrics.SubMetrics()) {
    sub_metrics.push_back(sub_metric);
  }
  if (sub_metrics.size() > 0) {
    r_list["sub_metrics"] = sub_metrics;
  }

  return wrap(r_list);
}

namespace traits {
//! Converter for an R-list to an EN penalty.
template<> class Exporter< nsoptim::EnPenalty > {
 public:
  explicit Exporter(SEXP r_obj) noexcept : r_obj_(r_obj) {}

  nsoptim::EnPenalty get() const {
    // We only create the penalty once we actually need it.
    List list(r_obj_);
    return nsoptim::EnPenalty(as<double>(list["alpha"]), as<double>(list["lambda"]));
  }

 private:
  SEXP r_obj_;
};

//! Converter for an R-list to a LASSO penalty.
template<> class Exporter< nsoptim::LassoPenalty > {
 public:
  explicit Exporter(SEXP r_obj) noexcept : r_obj_(r_obj) {}

  nsoptim::LassoPenalty get() const {
    List list(r_obj_);
    return nsoptim::LassoPenalty(as<double>(list["lambda"]));
  }

 private:
  SEXP r_obj_;
};

//! Converter for an R-list to a Ridge penalty.
template<> class Exporter< nsoptim::RidgePenalty > {
 public:
  explicit Exporter(SEXP r_obj) noexcept : r_obj_(r_obj) {}

  nsoptim::RidgePenalty get() const {
    List list(r_obj_);
    return nsoptim::RidgePenalty(as<double>(list["lambda"]));
  }

 private:
  SEXP r_obj_;
};

//! Create an exporter for any `std::forward_list` where the elements are supported by Rcpp::as.
//! This does not use the same functionality as in `Rcpp/internal/Exporter.h` because forward_lists are better
//! created sequentially.
template <typename T>
class Exporter< std::forward_list<T> > {
 public:
  explicit Exporter(SEXP r_obj) noexcept : r_obj_(r_obj) {}

  std::forward_list<T> get() const {
    std::forward_list<T> ret_list;
    auto ret_list_it = ret_list.before_begin();

    for (auto&& list_item : Rcpp::List(r_obj_)) {
      ret_list_it = ret_list.insert_after(ret_list_it, as<T>(list_item));
    }

    return ret_list;
  }
 private:
  SEXP r_obj_;
};

}  // namespace traits
}  // namespace Rcpp

#endif  // RCPP_UTILS_HPP_
