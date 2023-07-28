//
//  has_convex_surrogate.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_HAS_CONVEX_SURROGATE_HPP_
#define NSOPTIM_TRAITS_HAS_CONVEX_SURROGATE_HPP_

#include <utility>
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {

template<typename, typename>
static auto test_has_convex_surrogate(double) -> std::false_type;

template<typename T, typename U>
static auto test_has_convex_surrogate(int) -> sfinae_method_type<
  decltype(std::declval<T>().GetConvexSurrogate(std::declval<U>())), typename T::ConvexSurrogateType>;

}  // namespace internal

//! Type trait for loss & penalty functions that have a convex surrogate.
//! A loss/penalty function which has a convex surrogate must have a member function `ConvexSurrogate`.
template<typename T, typename U>
struct has_convex_surrogate : decltype(internal::test_has_convex_surrogate<T, U>(0)) {};
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_HAS_CONVEX_SURROGATE_HPP_
