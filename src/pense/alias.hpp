//
//  alias.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef ALIAS_HPP_
#define ALIAS_HPP_

#include <forward_list>
#include <memory>
#include "nsoptim_forward.hpp"

namespace pense {
namespace alias {
//! Alias for pointers to regression data.
using RegressionDataPtr = std::shared_ptr<nsoptim::PredictorResponseData>;
using ConstRegressionDataPtr = std::shared_ptr<const nsoptim::PredictorResponseData>;

//! Alias for std::forward_list used throughout the codebase
template<typename T>
using FwdList = std::forward_list<T>;

//! Alias for a list of optima.
template<typename Optimizer>
using Optima = FwdList<typename Optimizer::Optimum>;

}  // namespace alias
}  // namespace pense

#endif  // ALIAS_HPP_
