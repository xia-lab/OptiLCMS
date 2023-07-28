//
//  is_differentiable.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_DIFFERENTIABLE_HPP_
#define NSOPTIM_TRAITS_IS_DIFFERENTIABLE_HPP_

#include <utility>
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename, typename>
static auto test_is_differentiable(double) -> std::false_type;

template<typename T, typename U>
static auto test_is_differentiable(int) -> sfinae_method_type<decltype(std::declval<T>().Gradient(std::declval<U>())),
                                                              typename T::template GradientType<U>>;

}  // namespace internal

//! Type trait for differentiable loss & penalty functions.
//! A differentiable loss/penalty function supports computing of the gradient.
template<typename T, typename U>
struct is_differentiable : decltype(internal::test_is_differentiable<T, U>(0)) {};
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_IS_DIFFERENTIABLE_HPP_
