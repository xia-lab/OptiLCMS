//
//  regression_coefficients.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_CONTAINER_REGRESSION_COEFFICIENTS_HPP_
#define NSOPTIM_CONTAINER_REGRESSION_COEFFICIENTS_HPP_

#include <memory>
#include <type_traits>

#include "../armadillo.hpp"

namespace nsoptim {

//! Regression coefficients in the form of a scalar intercept and a (sparse) vector of slope coefficients.
//! The slope coefficients must be either of type `arma::vec` or `arma::sp_vec`.
template <class T>
class RegressionCoefficients {
  static_assert(std::is_same<T, arma::vec>::value || std::is_same<T, arma::sp_vec>::value,
                "T must be a (sparse) vector.");

 public:
  //! The dimensions of regression coefficients is the number if *slope* parameters.
  //! The intercept is not included in this number!
  using Dimensions = int;

  //! The type of the regression slope vector.
  using SlopeCoefficient = T;

  //! Initialize an empty coefficient vector and an intercept of 0.
  RegressionCoefficients() noexcept : intercept(0), beta() {}

  //! Initialize the coefficients as 0-vector
  //!
  //! @param n_pred the number of predictors, i.e., number of slope coefficients.
  explicit RegressionCoefficients(const arma::uword n_pred) noexcept
    : RegressionCoefficients(n_pred, typename std::is_same<T, arma::sp_vec>::type{}) {}

  //! Copy the slope coefficients from the given vector. The intercept will be set to 0.
  //!
  //! @param beta the value for the slope coefficients.
  explicit RegressionCoefficients(const SlopeCoefficient& beta) noexcept : intercept(0), beta(beta) {}

  //! Copy the coefficients from the given intercept and slope.
  //!
  //! @param intercept the value for the intercept.
  //! @param beta the value for the slope coefficients.
  RegressionCoefficients(const double intercept, const SlopeCoefficient& beta) noexcept
      : intercept(intercept), beta(beta) {}

  //! Copy constructor
  //!
  //! @param coef the RegressionCoefficients object to copy.
  RegressionCoefficients(const RegressionCoefficients& coef) = default;

  //! Move constructor.
  //!
  //! @param coef the RegressionCoefficients object to move.
  RegressionCoefficients(RegressionCoefficients&& coef) = default;

  //! Copy assignment
  //!
  //! @param coef the RegressionCoefficients object to copy.
  RegressionCoefficients& operator=(const RegressionCoefficients& other) = default;

  //! Move assignment
  //!
  //! @param coef the RegressionCoefficients object to move.
  RegressionCoefficients& operator=(RegressionCoefficients&& other) = default;

  //! Reset the coefficients to an empty vector and an intercept of 0.
  inline void Reset() noexcept {
    intercept = 0;
    beta.reset();
  }

  //! Print the coefficients to the given output stream
  //!
  //! @param user_stream the stream to print the coefficients to.
  inline void Print(std::ostream &user_stream) const {
    user_stream << "Intercept: " << intercept << "\n";
    beta.t().print(user_stream, "Beta:");
  }

  //! Print the coefficients to stdout
  inline void Print() const {
    std::cout << "Intercept: " << intercept << "\n";
    beta.t().print("Beta:");
  }

  double intercept;       //< The intercept coefficient
  SlopeCoefficient beta;  //< The slope coefficients

 private:
  // Initialize the regression coefficients with all zeros if the slope is a dense vector.
  RegressionCoefficients(const arma::uword n_pred, std::false_type) noexcept
    : intercept(0), beta(n_pred, arma::fill::zeros) {}
  // Initialize the regression coefficients with all zeros if the slope is a sparse vector.
  RegressionCoefficients(const arma::uword n_pred, std::true_type) noexcept
    : intercept(0), beta(n_pred) {}
};
}  // namespace nsoptim

#endif  // NSOPTIM_CONTAINER_REGRESSION_COEFFICIENTS_HPP_
