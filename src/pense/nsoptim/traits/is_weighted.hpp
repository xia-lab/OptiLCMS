//
//  is_weighted.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_WEIGHTED_HPP_
#define NSOPTIM_TRAITS_IS_WEIGHTED_HPP_

#include <utility>
#include "../armadillo.hpp"
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename T>
static auto test_is_weighted(int) -> sfinae_method_type<decltype(std::declval<T>().weights()), arma::vec>;

template<typename T>
static auto test_is_weighted(double) -> std::false_type;
}  // namespace internal

//! Type trait for weighted loss functions.
//! Weighted loss functions must have a member function `weights` to access the weights.
template<typename T>
struct is_weighted : decltype(internal::test_is_weighted<T>(0)) {};
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_IS_WEIGHTED_HPP_
