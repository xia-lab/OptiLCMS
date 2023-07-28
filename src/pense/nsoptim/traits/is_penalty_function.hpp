//
//  is_penalty_function.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_PENALTY_FUNCTION_HPP_
#define NSOPTIM_TRAITS_IS_PENALTY_FUNCTION_HPP_

#include <type_traits>
#include <utility>

#include "sfinae_types.hpp"
#include "can_evaluate.hpp"
#include "../objective/penalty.hpp"

namespace nsoptim {
namespace traits {
//! Type trait if a type implements the LossFunction interface.
// template<typename> struct is_penalty_function : std::false_type {};
template<typename T> struct is_penalty_function : internal::tf_switch<std::is_base_of<PenaltyFunction, T>::value &&
                                                                      std::is_copy_constructible<T>::value> {};
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_IS_PENALTY_FUNCTION_HPP_
