//
//  can_optimize.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_CAN_OPTIMIZE_HPP_
#define NSOPTIM_TRAITS_CAN_OPTIMIZE_HPP_

#include <utility>
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename T, typename U>
static auto test_can_optimize_empty(double) -> std::false_type;

template<typename T, typename U>
static auto test_can_optimize_empty(int) -> sfinae_method_type<decltype(std::declval<T>().Optimize()),
                                                               typename T::Optimum>;


template<typename T, typename U>
static auto test_can_optimize_start(double) -> std::false_type;

template<typename T, typename U>
static auto test_can_optimize_start(int) -> sfinae_method_type<
  decltype(std::declval<T>().Optimize(std::declval<U>())), typename T::Optimum>;
}  // namespace internal

//! Type trait for optimizer.
//! Checks whether the optimizer T can optimize for coefficients U
template<typename T, typename U>
struct can_optimize : internal::tf_switch<decltype(internal::test_can_optimize_empty<T, U>(0))::value &&
                                          decltype(internal::test_can_optimize_start<T, U>(0))::value> {};
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_CAN_OPTIMIZE_HPP_
