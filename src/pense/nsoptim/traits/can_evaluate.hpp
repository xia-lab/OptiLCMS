//
//  can_evaluate.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_CAN_EVALUATE_FUNCTION_HPP_
#define NSOPTIM_TRAITS_CAN_EVALUATE_FUNCTION_HPP_

#include <type_traits>
#include <utility>
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
//! Test if the type T can be evaluated with the coefficient type U
template<typename, typename>
static auto test_can_evaluate(double) -> std::false_type;

//! Test if the type T can be evaluated with the coefficient type U
template<typename T, typename U>
static auto test_can_evaluate(int) -> sfinae_method_type<decltype(std::declval<T>()(std::declval<U>())), double>;

}  // namespace internal

//! Type trait if the type T supports evaluation of the coefficient type `U`.
template<typename T, typename U>
struct can_evaluate : decltype(internal::test_can_evaluate<T, U>(0)) {};

}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_CAN_EVALUATE_HPP_
