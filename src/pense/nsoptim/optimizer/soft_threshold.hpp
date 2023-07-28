//
//  soft_threshold.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_SOFT_THRESHOLD_HPP_
#define NSOPTIM_OPTIMIZER_SOFT_THRESHOLD_HPP_

#include <algorithm>

#include "../armadillo.hpp"

namespace nsoptim {
//! Soft-thresholding sign(z) * max(0, |z| - gamma) for a scalar.
//!
//! @param z the value to soft-threshold.
//! @param gamma the threshold.
//! @return the soft-thresholded value.
inline double SoftThreshold(const double z, const double gamma) noexcept {
  if (std::fabs(z) <= gamma) {
    return 0.;
  } else if (z < 0) {
    return z + gamma;
  }
  return z - gamma;
}

//! Soft-thresholding sign(z) * max(0, |z| - gamma) for a vector of values using the same threshold.
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the threshold.
inline void SoftThreshold(const double gamma, arma::vec* z) noexcept;

//! Soft-thresholding sign(z) * max(0, |z| - gamma) for a vector of values using different thresholds.
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the thresholds.
inline void SoftThreshold(const arma::vec& gamma, arma::vec* z) noexcept;

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a dense vector `z1`.
//!
//! @param z1 a vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma the threshold.
//! @return the vector of thresholded values.
inline arma::vec SoftThreshold(arma::vec z1, const double c, const arma::vec& z2,
                               const double gamma) noexcept;

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1`.
//!
//! @param z1 a sparse vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma the threshold.
//! @return the sparse vector of thresholded values.
inline arma::sp_vec SoftThreshold(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                  const double gamma) noexcept;

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma_j) for a dense vector `z1`.
//!
//! @param z1 a sparse vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma a vector of different threshold values.
//! @return the sparse vector of thresholded values.
inline arma::vec SoftThreshold(arma::vec z1, const double c, const arma::vec& z2,
                               const arma::vec& gamma) noexcept;

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma_j) for a sparse vector `z1`.
//!
//! @param z1 a sparse vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma a vector of different threshold values.
//! @return the sparse vector of thresholded values.
inline arma::sp_vec SoftThreshold(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                  const arma::vec& gamma) noexcept;

namespace soft_threshold {
constexpr float kVecSoftthreshInplaceSparse = 1.5;   //< Vectors with more than 2/3 non-zero elements are always
                                                     //< considered dense.

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1` that is *actually*
//! sparse.
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the thresholds.
inline arma::sp_vec SoftThresholdSparse(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                        const double gamma) noexcept;

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1` that is
//! *actually dense*
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the thresholds.
inline arma::sp_vec SoftThresholdDense(const arma::sp_vec& sz1, const double c, const arma::vec& z2,
                                       const double gamma) noexcept;

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1` that is *actually*
//! sparse.
//!
//! @param gamma the thresholds.
inline arma::sp_vec SoftThresholdSparse(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                        const arma::vec& gamma) noexcept;

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma_j) for a sparse vector `z1` that is
//! *actually dense*
//!
//! @param gamma the thresholds.
inline arma::sp_vec SoftThresholdDense(const arma::sp_vec& sz1, const double c, const arma::vec& z2,
                                       const arma::vec& gamma) noexcept;

}  // namespace soft_threshold

//! Soft-thresholding sign(z) * max(0, |z| - gamma) for a vector of values using the same threshold.
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the threshold.
inline void SoftThreshold(const double gamma, arma::vec* z) noexcept {
  for (auto el_it = z->begin(), z_end = z->end(); el_it != z_end; ++el_it) {
    (*el_it) = SoftThreshold(*el_it, gamma);
  }
}

//! Soft-thresholding sign(z) * max(0, |z| - gamma) for a vector of values using different thresholds.
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the thresholds.
inline void SoftThreshold(const arma::vec& gamma, arma::vec* z) noexcept {
  auto gamma_it = gamma.begin();
  for (auto el_it = z->begin(), z_end = z->end(); el_it != z_end; ++el_it, ++gamma_it) {
    (*el_it) = SoftThreshold(*el_it, *gamma_it);
  }
}

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a dense vector `z1`.
//!
//! @param z1 a vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma the threshold.
//! @return the vector of thresholded values.
inline arma::vec SoftThreshold(arma::vec z1, const double c, const arma::vec& z2,
                               const double gamma) noexcept {
  auto z2_it = z2.begin();

  const auto z1_end = z1.end();
  for (auto z1_it = z1.begin(); z1_it != z1_end; ++z1_it, ++z2_it) {
    *z1_it += c * (*z2_it);
    if (*z1_it > gamma) {
      *z1_it -= gamma;
    } else if (*z1_it < -gamma) {
      *z1_it += gamma;
    } else {
      *z1_it = 0;
    }
  }
  return z1;
}

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1`.
//!
//! @param z1 a sparse vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma the threshold.
//! @return the sparse vector of thresholded values.
inline arma::sp_vec SoftThreshold(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                  const double gamma) noexcept {
  if (z1.n_elem > soft_threshold::kVecSoftthreshInplaceSparse * z1.n_nonzero) {
    // This branch is better for large sparse vectors
    return soft_threshold::SoftThresholdSparse(z1, c, z2, gamma);
  } else {
    // This branch is better for all other vectors
    return soft_threshold::SoftThresholdDense(z1, c, z2, gamma);
  }
}

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma_j) for a dense vector `z1`.
//!
//! @param z1 a vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma a vector of different threshold values.
//! @return the vector of thresholded values.
inline arma::vec SoftThreshold(arma::vec z1, const double c, const arma::vec& z2,
                               const arma::vec& gamma) noexcept {
  auto z2_it = z2.begin();
  auto gamma_it = gamma.begin();

  const auto z1_end = z1.end();
  for (auto z1_it = z1.begin(); z1_it != z1_end; ++z1_it, ++z2_it, ++gamma_it) {
    *z1_it += c * (*z2_it);
    if (*z1_it > *gamma_it) {
      *z1_it -= *gamma_it;
    } else if (*z1_it < -*gamma_it) {
      *z1_it += *gamma_it;
    } else {
      *z1_it = 0;
    }
  }
  return z1;
}

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma_j) for a sparse vector `z1`.
//!
//! @param z1 a sparse vector.
//! @param c a constant multiplicative value.
//! @param z2 a vector.
//! @param gamma a vector of different threshold values.
//! @return the sparse vector of thresholded values.
inline arma::sp_vec SoftThreshold(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                  const arma::vec& gamma) noexcept {
  if (z1.n_elem > soft_threshold::kVecSoftthreshInplaceSparse * z1.n_nonzero) {
    // This branch is better for large sparse vectors
    return soft_threshold::SoftThresholdSparse(z1, c, z2, gamma);
  } else {
    // This branch is better for all other vectors
    return soft_threshold::SoftThresholdDense(z1, c, z2, gamma);
  }
}

namespace soft_threshold {
//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1` that is *actually*
//! sparse.
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the thresholds.
inline arma::sp_vec SoftThresholdSparse(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                        const double gamma) noexcept {
  arma::uvec ind(z1.n_elem);
  arma::vec val(z1.n_elem);
  arma::uword nnz = 0;

  arma::sp_vec::const_iterator z1_read_iter = z1.begin();
  auto z2_read_iter = z2.begin();

  arma::uword i = 0;
  do {
    // The z1 iterator will only return the non-0 values. Therefore, we need to handle all
    // elements up the the current z1 iterator's index (`row`).
    const auto upper_bound = (z1_read_iter == z1.end() ? z1.n_elem : z1_read_iter.row());
    while (i < upper_bound) {
      // We know that z1 is 0 at index `i`
      const double tmp = c * (*z2_read_iter);
      if (tmp > gamma) {
        ind[nnz] = i;
        val[nnz++] = tmp - gamma;
      } else if (tmp < -gamma) {
        ind[nnz] = i;
        val[nnz++] = tmp + gamma;
      }
      // Advance the iterators for the non-sparse vectors
      ++z2_read_iter;
      ++i;
    }

    // If we are not yet at the end of the vectors, the z1 iterator points to a non-0 element in z1.
    if (i < z1.n_elem) {
      const double tmp = (*z1_read_iter) + c * (*z2_read_iter);
      if (tmp > gamma) {
        ind[nnz] = i;
        val[nnz++] = tmp - gamma;
      } else if (tmp < -gamma) {
        ind[nnz] = i;
        val[nnz++] = tmp + gamma;
      }
    }
    // Advance all the iterators. Both for the sparse vector z1 and the non-sparse vectors.
    ++z1_read_iter;
    ++z2_read_iter;
  } while (++i < z1.n_elem);

  // From the index and value vectors, create a sparse vector.
  if (nnz > 0) {
    return arma::sp_mat(ind.head_rows(nnz), arma::uvec {0, nnz}, val.head_rows(nnz), z1.n_elem, 1);
  } else {
    return arma::sp_vec(z1.n_elem);
  }
}

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1` that is
//! *actually dense*
//!
//! @param z a pointer to a vector of values to soft-threshold. Holds the thresholded values afterwards.
//! @param gamma the thresholds.
inline arma::sp_vec SoftThresholdDense(const arma::sp_vec& sz1, const double c, const arma::vec& z2,
                                       const double gamma) noexcept {
  arma::vec z1(sz1);
  const auto z1_end = z1.end();
  auto z2_read_iter = z2.cbegin();

  // Simple iteration is possible because we use dense vectors!
  for (auto z1_iter = z1.begin(); z1_iter != z1_end; ++z1_iter, ++z2_read_iter) {
    *z1_iter += c * (*z2_read_iter);
    if (*z1_iter > gamma) {
      *z1_iter -= gamma;
    } else if (*z1_iter < -gamma) {
      *z1_iter += gamma;
    } else {
      *z1_iter = 0;
    }
  }

  // Convert back to a sparse vector
  return arma::sp_vec(z1);
}

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma) for a sparse vector `z1` that is *actually*
//! sparse.
//!
//! @param gamma the thresholds.
inline arma::sp_vec SoftThresholdSparse(const arma::sp_vec& z1, const double c, const arma::vec& z2,
                                        const arma::vec& gamma) noexcept {
  arma::uvec ind(z1.n_elem);
  arma::vec val(z1.n_elem);
  arma::uword nnz = 0;

  arma::sp_vec::const_iterator z1_read_iter = z1.begin();
  auto z2_read_iter = z2.begin();
  auto gamma_read_iter = gamma.cbegin();

  arma::uword i = 0;
  do {
    // The z1 iterator will only return the non-0 values. Therefore, we need to handle all
    // elements up the the current z1 iterator's index (`row`).
    const auto upper_bound = (z1_read_iter == z1.end() ? z1.n_elem : z1_read_iter.row());
    while (i < upper_bound) {
      // We know that z1 is 0 at index `i`
      const double tmp = c * (*z2_read_iter);
      if (tmp > *gamma_read_iter) {
        ind[nnz] = i;
        val[nnz++] = tmp - *gamma_read_iter;
      } else if (tmp < -*gamma_read_iter) {
        ind[nnz] = i;
        val[nnz++] = tmp + *gamma_read_iter;
      }
      // Advance the iterators for the non-sparse vectors
      ++z2_read_iter;
      ++gamma_read_iter;
      ++i;
    }

    // If we are not yet at the end of the vectors, the z1 iterator points to a non-0 element in z1.
    if (i < z1.n_elem) {
      const double tmp = (*z1_read_iter) + c * (*z2_read_iter);
      if (tmp > *gamma_read_iter) {
        ind[nnz] = i;
        val[nnz++] = tmp - *gamma_read_iter;
      } else if (tmp < -*gamma_read_iter) {
        ind[nnz] = i;
        val[nnz++] = tmp + *gamma_read_iter;
      }
    }
    // Advance all the iterators. Both for the sparse vector z1 and the non-sparse vectors.
    ++z1_read_iter;
    ++gamma_read_iter;
    ++z2_read_iter;
  } while (++i < z1.n_elem);

  // From the index and value vectors, create a sparse vector.
  if (nnz > 0) {
    return arma::sp_mat(ind.head_rows(nnz), arma::uvec {0, nnz}, val.head_rows(nnz), z1.n_elem, 1);
  } else {
    return arma::sp_vec(z1.n_elem);
  }
}

//! Soft-thresholding sign(z1 + c * z2) * max(0, |z1 + c * z2| - gamma_j) for a sparse vector `z1` that is
//! *actually dense*
//!
//! @param gamma the thresholds.
arma::sp_vec SoftThresholdDense(const arma::sp_vec& sz1, const double c, const arma::vec& z2,
                                const arma::vec& gamma) noexcept {
  // Make a dense copy of the sparse vector `z1`. We know that `z1` is not that sparse
  // after all or small enough to make this operation faster than having to deal with a
  // sparse vector.
  arma::vec z1(sz1);
  const auto z1_end = z1.end();
  auto z2_read_iter = z2.cbegin();
  auto gamma_read_iter = gamma.cbegin();

  // Simple iteration is possible because we use dense vectors!
  for (auto z1_it = z1.begin(); z1_it != z1_end; ++z1_it, ++z2_read_iter, ++gamma_read_iter) {
    *z1_it += c * (*z2_read_iter);
    if (*z1_it > *gamma_read_iter) {
      *z1_it -= *gamma_read_iter;
    } else if (*z1_it < -*gamma_read_iter) {
      *z1_it += *gamma_read_iter;
    } else {
      *z1_it = 0;
    }
  }

  // Convert back to a sparse vector
  return arma::sp_vec(z1);
}
}  // namespace soft_threshold
}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_SOFT_THRESHOLD_HPP_
