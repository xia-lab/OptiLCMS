//
//  is_en_penalty.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_EN_PENALTY_HPP_
#define NSOPTIM_TRAITS_IS_EN_PENALTY_HPP_

#include <utility>
#include "sfinae_types.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename T>
static auto test_is_en_penalty(sfinae_type_wrapper<typename T::is_en_penalty_tag>*) -> std::true_type;

template<typename>
static auto test_is_en_penalty(...) -> std::false_type;

}  // namespace internal

//! Type trait to identify a penalty function as "elastic net"-like.
template<typename T>
struct is_en_penalty : decltype(internal::test_is_en_penalty<T>(0)) {};
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_IS_EN_PENALTY_HPP_
