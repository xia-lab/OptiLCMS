//
//  armadillo.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_CONFIG_HPP_
#define NSOPTIM_CONFIG_HPP_

#define NSOPTIM_METRICS_LEVEL 1

#ifdef NSOPTIM_METRICS_DETAILED
# undef NSOPTIM_METRICS_LEVEL
# define NSOPTIM_METRICS_LEVEL 2
#endif

#ifdef NSOPTIM_METRICS_ENABLED
# undef NSOPTIM_METRICS_LEVEL
# define NSOPTIM_METRICS_LEVEL 1
#endif

#ifdef NSOPTIM_METRICS_DISABLED
# undef NSOPTIM_METRICS_LEVEL
# define NSOPTIM_METRICS_LEVEL 0
#endif

#endif  // NSOPTIM_CONFIG_HPP_
