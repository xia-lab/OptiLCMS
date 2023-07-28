//
//  is_adaptive.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_ADAPTIVE_HPP_
#define NSOPTIM_TRAITS_IS_ADAPTIVE_HPP_

#include <utility>
#include "../armadillo.hpp"
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename T>
static auto test_has_loadings(int) -> sfinae_method_type<decltype(std::declval<T>().loadings()), arma::vec>;

template<typename T>
static auto test_has_loadings(double) -> std::false_type;
}  // namespace internal

//! Type trait for adaptive penalty functions.
//! Adaptive penalty functions have a member function `loadings` to access the penalty loadings.
template<typename T>
struct is_adaptive : decltype(internal::test_has_loadings<T>(0)) {};
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_IS_ADAPTIVE_HPP_
