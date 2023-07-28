//
//  is_ls_regression_loss.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_IS_LS_REGRESSION_LOSS_HPP_
#define NSOPTIM_TRAITS_IS_LS_REGRESSION_LOSS_HPP_

#include <utility>
#include "sfinae_types.hpp"
#include "is_loss_function.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename T>
static auto test_is_ls_regression_loss(sfinae_type_wrapper<typename T::is_ls_regression_loss_tag>*) -> std::true_type;

template<typename>
static auto test_is_ls_regression_loss(...) -> std::false_type;

}  // namespace internal

//! Type trait to identify a penalty function as "elastic net"-like.
// template<typename T>
// struct is_ls_regression_loss : decltype(internal::test_is_ls_regression_loss<T>(0)) {};

template<typename T>
struct is_ls_regression_loss : internal::tf_switch<decltype(internal::test_is_ls_regression_loss<T>(0))::value &&
                                                            is_loss_function<T>::value>::type {};

}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_IS_LS_REGRESSION_LOSS_HPP_
