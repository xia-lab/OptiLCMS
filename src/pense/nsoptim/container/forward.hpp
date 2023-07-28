//
//  forward.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_CONTAINER_FORWARD_HPP_
#define NSOPTIM_CONTAINER_FORWARD_HPP_

#include "../config.hpp"
#include "../armadillo_forward.hpp"

namespace nsoptim {
//! Full definition at nsoptim/container/regression_coefficients.hpp
template <class T> class RegressionCoefficients;

//! Full definition at nsoptim/container/regression_coefficients.hpp
template <class T> class RegressionCoefficients;

//! Full definition at nsoptim/container/data.hpp
class PredictorResponseData;

namespace _metrics_internal {
//! Full definition at nsoptim/container/metrics.hpp
template<int> class Metrics;
}  // namespace _metrics_internal

//! Export the correct Metrics collection based on the configuration.
using Metrics = _metrics_internal::Metrics<NSOPTIM_METRICS_LEVEL>;
}  // namespace nsoptim

#endif  // NSOPTIM_CONTAINER_FORWARD_HPP_
