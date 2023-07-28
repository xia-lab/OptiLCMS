//
//  is_loss_function.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_LOSS_FUNCTION_HPP_
#define NSOPTIM_TRAITS_IS_LOSS_FUNCTION_HPP_

#include <type_traits>
#include <utility>

#include "sfinae_types.hpp"
#include "can_evaluate.hpp"
#include "../objective/loss.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
//! Test if the loss function T uses data type U
template<typename, typename>
static auto test_loss_supports_data(double) -> std::false_type;

//! Test if the loss function T uses data type U
template<typename T, typename U>
static auto test_loss_supports_data(int) -> sfinae_method_type<decltype(std::declval<T>().data()), U>;

//! Test if the loss function T can create a "zero" coefficient object of type U.
template<typename, typename>
static auto test_loss_supports_zero_coefs(double) -> std::false_type;

//! Test if the loss function T can create a "zero" coefficient object of type U.
template<typename T, typename U>
static auto test_loss_supports_zero_coefs(int) -> sfinae_method_type<
    decltype(std::declval<T>().template ZeroCoefficients<U>()), U>;
}  // namespace internal

//! Type trait for a loss functions.
//! Tests if the given type implements the LossFunction interface.
template<typename T, typename U = typename T::DataType>
struct has_data_member : decltype(internal::test_loss_supports_data<T, U>(0))::type {};

//! Type trait if a loss function supports evaluation of the coefficient type `U`.
template<typename T, typename U>
struct loss_supports_evaluation : internal::tf_switch<
  decltype(internal::test_loss_supports_zero_coefs<T, U>(0))::value &&
  can_evaluate<T, U>::value> {};

//! Type trait if a type implements the LossFunction interface.
template<typename T> struct is_loss_function : internal::tf_switch<
  std::is_copy_constructible<T>::value && std::is_base_of<LossFunction<typename T::DataType>, T>::value> {};

}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_IS_LOSS_FUNCTION_HPP_
