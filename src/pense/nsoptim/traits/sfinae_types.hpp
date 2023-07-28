//
//  sfinae_types.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-01-26.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_TRAITS_SFINAE_TYPES_HPP_
#define NSOPTIM_TRAITS_SFINAE_TYPES_HPP_

#include <type_traits>

#include "../armadillo.hpp"

namespace nsoptim {
namespace traits {
namespace internal {
template<typename T, typename U,
         typename = typename std::enable_if<std::is_same<typename std::decay<T>::type, U>::value, void>::type>
struct sfinae_method_type : std::true_type {};

template<typename T, template<typename> class Tag>
struct sfinae_method_tagged : Tag<typename std::decay<T>::type>::type {};

template<typename T, typename... Ts>
struct sfinae_method_any : std::true_type {};

template<typename T> struct sfinae_type_wrapper {};

template<bool IsTrue>
using tf_switch = typename std::conditional<IsTrue, std::true_type, std::false_type>::type;

}  // namespace internal
}  // namespace traits
}  // namespace nsoptim

#endif  // NSOPTIM_TRAITS_SFINAE_TYPES_HPP_
