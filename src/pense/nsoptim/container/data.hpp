//
//  data.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_CONTAINER_DATA_HPP_
#define NSOPTIM_CONTAINER_DATA_HPP_

#include <iostream>

#include "../armadillo.hpp"
#include "../utilities.hpp"

namespace nsoptim {

//! Simple structure holding a matrix with predictor data *x* (with `n_obs` rows and `n_pred` colums)
//! and the response vector *y* (with `n_obs` rows)
class PredictorResponseData {
 public:
  //! Initialize the predictor-response data with empty x and y
  PredictorResponseData() noexcept : x_(), y_(), n_obs_(0), n_pred_(0) {}

  //! Initialize predictor-response data with the given x and y.
  //! @note the given data will be copied!
  //!
  //! @param other_x predictor matrix to copy.
  //! @param other_y response vector to copy.
  PredictorResponseData(const arma::mat& other_x, const arma::vec& other_y) noexcept
    : x_(other_x), y_(other_y), n_obs_(other_x.n_rows), n_pred_(other_x.n_cols) {}

  //! Initialize predictor-response data with the given x and y.
  //! @note the given data will be moved to this container!
  //!
  //! @param other_x predictor matrix to move.
  //! @param other_y response vector to move.
  PredictorResponseData(arma::mat&& other_x, arma::vec&& other_y) noexcept
    : x_(std::move(other_x)), y_(std::move(other_y)), n_obs_(x_.n_rows), n_pred_(x_.n_cols) {}

  //! Copy the given predictor-response data, but pointing to the same underlying data!
  //!
  //! @param other predictor-response data to copy.
  PredictorResponseData(const PredictorResponseData& other) = default;

  PredictorResponseData& operator=(const PredictorResponseData& other) = default;

  //! Get a data set with the observations at the requested indices.
  //!
  //! @param indices the indicies of the observations to get.
  //! @return the subset of the data with the requested observations.
  PredictorResponseData Observations(const arma::uvec& indices) const {
    return PredictorResponseData(x_.rows(indices), y_.rows(indices));
  }

  //! Get a data set with the given observation removed.
  //!
  //! @param index the index of the observation to remove.
  //! @return the subset of the data with the observation removed.
  PredictorResponseData RemoveObservation(const arma::uword index) const {
    return PredictorResponseData(arma::join_vert(x_.head_rows(index), x_.tail_rows(n_obs_ - index - 1)),
                                 arma::join_vert(y_.head(index), y_.tail(n_obs_ - index - 1)));
  }

  //! Get a data set with the first `n_obs` observations of the data.
  //!
  //! @param n_obs number of observations to extract.
  //! @return the subset of the data with the requested observations.
  PredictorResponseData HeadRows(const arma::uword n_obs) const {
    return PredictorResponseData(x_.head_rows(n_obs), y_.head_rows(n_obs));
  }

  //! Get a data set with the last `n_obs` observations of the data.
  //!
  //! @param n_obs number of observations to extract.
  //! @return the subset of the data with the requested rows.
  PredictorResponseData TailRows(const arma::uword n_obs) const {
    return PredictorResponseData(x_.tail_rows(n_obs), y_.tail_rows(n_obs));
  }

  //! Compare two data containers based on their object ID.
  //!
  //! This does not compare the actual *data*, it only compares their ID. The ID is only equal if the objects
  //! are the same object or exact copies of each other.
  //!
  //! @param other the other data container.
  //! @return true if both container have different ids (i.e., they are different data objects!)
  bool operator!=(const PredictorResponseData& other) const noexcept {
    return id_ != other.id_;
  }

  //! Compare two data containers based on their object ID.
  //!
  //! This does not compare the actual *data*, it only compares their ID. The ID is only equal if the objects
  //! are the same object or exact copies of each other.
  //!
  //! @param other the other data container.
  //! @return true if both container have the same id (i.e., they are the same data object!)
  bool operator==(const PredictorResponseData& other) const noexcept {
    return id_ == other.id_;
  }

  //! Get the ID of the data container.
  //!
  //! @return a program-unique ID for the data container.
  ObjectId id() const noexcept {
    return id_;
  }

  //! Get a constant reference to the predictor matrix.
  //! Only valid as long as the PredictorResponseData object is in scope.
  //!
  //! @return constant reference to the predictor matrix
  const arma::mat& cx() const noexcept {
    return x_;
  }

  //! Get a constant reference to the response vector.
  //! Only valid as long as the PredictorResponseData object is in scope.
  //!
  //! @return constant reference to the response vector
  const arma::vec& cy() const noexcept {
    return y_;
  }

  //! Get non-const references to the data
  //! Get a reference to the predictor matrix.
  //! Only valid as long as the PredictorResponseData object is in scope.
  //!
  //! @return reference to the predictor matrix
  arma::mat& x() noexcept {
    return x_;
  }

  //! Get a reference to the response vector.
  //! Only valid as long as the PredictorResponseData object is in scope.
  //!
  //! @return reference to the response vector
  arma::vec& y() noexcept {
    return y_;
  }

  //! Get the number of observations in this data set.
  arma::uword n_obs() const noexcept {
    return n_obs_;
  }
  //! Get the number of observations in this data set.
  arma::uword n_pred() const noexcept {
    return n_pred_;
  }

 private:
  ObjectId id_;
  arma::mat x_;
  arma::vec y_;
  arma::uword n_obs_;   //< The number of observations in the data.
  arma::uword n_pred_;  //< The number of variables in the data.
};

}  // namespace nsoptim

#endif  // NSOPTIM_CONTAINER_DATA_HPP_
