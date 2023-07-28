//
//  enpy_types.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef ENPY_TYPES_HPP_
#define ENPY_TYPES_HPP_

#include "nsoptim.hpp"
#include "alias.hpp"

namespace pense {
//! PY Result Structure
//! Contains a list of initial estimates and the associated metrics.
template<typename T>
struct PyResult {
  PyResult() noexcept : metrics("enpy_initest") {}
  explicit PyResult(nsoptim::Metrics&& _metrics) noexcept : metrics(std::move(_metrics)) {}

  nsoptim::Metrics metrics;
  alias::Optima<T> initial_estimates;
};

}  // namespace pense

#endif  // ENPY_TYPES_HPP_
