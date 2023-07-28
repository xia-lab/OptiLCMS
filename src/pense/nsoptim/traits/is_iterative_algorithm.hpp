//
//  is_iterative.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_ITERATIVE_ALGORITHM_HPP_
#define NSOPTIM_TRAITS_IS_ITERATIVE_ALGORITHM_HPP_

#include <utility>
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename T>
static auto test_is_iterative(int) -> sfinae_method_any<decltype(std::declval<T>().convergence_tolerance(1.0))>;

template<typename T>
static auto test_is_iterative(double) -> std::false_type;
}  // namespace internal

//! Type trait for iterative algorithms.
//! Iterative algorithms support changes to the convergence threshold and calling Optimize with a maximum number
//! of iterations.
template<typename T>
struct is_iterative_algorithm : decltype(internal::test_is_iterative<T>(0)) {};
}  // namespace traits
}  // namespace nsoptim
#endif  // NSOPTIM_TRAITS_IS_ITERATIVE_ALGORITHM_HPP_
