//
//  convex.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OBJECTIVE_CONVEX_HPP_
#define NSOPTIM_OBJECTIVE_CONVEX_HPP_

namespace nsoptim {

//! CRTP helper class for convex functions which returns the object itself as convex surrogate.
template <typename Function>
class ConvexFunction {
 public:
  using ConvexSurrogateType = Function;

  template<typename T>
  Function& GetConvexSurrogate(const T&) {
    return static_cast<Function&>(*this);
  }
};

}  // namespace nsoptim

#endif  // NSOPTIM_OBJECTIVE_CONVEX_HPP_
