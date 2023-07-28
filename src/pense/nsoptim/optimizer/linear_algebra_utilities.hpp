//
//  linear_algebra_utilities.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2019-11-14.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_LINEAR_ALGEBRA_UTILITIES_HPP_
#define NSOPTIM_OPTIMIZER_LINEAR_ALGEBRA_UTILITIES_HPP_

#include <memory>
#include <algorithm>
#include <limits>

#include "../armadillo.hpp"

#if !defined(ARMA_BLAS_CAPITALS)
#  define nsoptim_dtpsv dtpsv_
#else
#  define nsoptim_dtpsv DTPSV_
#endif

extern "C"
{
#if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
void nsoptim_dtpsv(const char* uplo, const char* trans, const char* diag, const arma::blas_int* n,
                   double* ap, double *x, const arma::blas_int* incx, arma::blas_len uplo_len,
                   arma::blas_len trans_len, arma::blas_len diag_len);
#else
void nsoptim_dtpsv(const char* uplo, const char* trans, const char* diag, const arma::blas_int* n,
                   double* ap, double *x, const arma::blas_int* incx);
#endif
}

namespace nsoptim {
namespace linalg {

//! Wrapper around the BLAS routine dtpsv.
inline void dtpsv(const char* uplo, const char* trans, const char* diag, const arma::blas_int* n,
                  double* ap, double *x, const arma::blas_int* incx) {
#if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
  nsoptim_dtpsv(uplo, trans, diag, n, ap, x, incx, 1, 1, 1);
#else
  nsoptim_dtpsv(uplo, trans, diag, n, ap, x, incx);
#endif
}

//! Solve a system of linear equations by using a Cholesky decomposition `chol`.
//!
//! Solves the system of equations `A x = b` for x using the Cholesky decomposition of `A = chol' . chol` in two steps:
//!  1. solve `chol' y = b` for `y`, and
//!  2. solve `chol x = y` for `x`.
//!
//! @param chol the Cholesky decomposition of matrix `A`.
//! @param b the target vector `b` on input and the solution `x` on output.
//! @return true if a solution is found, false otherwise.
inline bool SolveChol(const arma::mat& chol, arma::vec * const b) {
  // Instead of
  // return arma::solve(arma::trimatu(chol), arma::solve(arma::trimatl(arma::trans(chol)), x));
  // it is done manually to avoid the explicit transposition
  char uplo = 'U';
  char trans_n = 'N';
  char trans_y = 'T';
  char diag  = 'N';
  arma::blas_int n = arma::blas_int(chol.n_rows);
  arma::blas_int nrhs = 1;
  arma::blas_int info  = 0;

  arma::lapack::trtrs(&uplo, &trans_y, &diag, &n, &nrhs, chol.memptr(), &n, b->memptr(), &n, &info);
  if (info != 0) {
    return false;
  }

  arma::lapack::trtrs(&uplo, &trans_n, &diag, &n, &nrhs, chol.memptr(), &n, b->memptr(), &n, &info);
  if (info != 0) {
    return false;
  }

  return true;
}

//! Solve a system of linear equations `A x = b` by using conjugate gradient, starting from `x`.
//! @param A the symmetric positive definite design matrix
//! @param b the right-hand side of the equation
//! @param eps numerical tolerance for convergence
//! @param max_it maximum number of iterations for early stopping. If 0, no early stopping is performed
//! @param x the coefficient vector with a reasonably close initial guess on input and the solution on output.
//! @return the number of iterations required for convergence.
inline arma::uword SolveIterative(const arma::mat& A, const arma::vec& b, const double eps, arma::uword max_it,
                                  arma::vec * const x) {
  arma::uword it = 0;
  arma::vec resid = b - A * (*x);
  arma::vec step_dir = resid;
  double resid_norm = arma::norm(resid, 2);
  const double conv_threshold = eps * std::max(eps, resid_norm);

  if (max_it == 0) {
    max_it = b.n_elem;
  }

  while (true) {
    const arma::vec trans_step_dir = A * step_dir;
    const double step_size = resid_norm * resid_norm / arma::dot(step_dir, trans_step_dir);
    (*x) += step_size * step_dir;
    if (it % 4 == 0) {
      // at every 4th step, re-compute the residuals to avoid drift
      resid = b - A * (*x);
    } else {
      resid -= step_size * trans_step_dir;
    }
    if ((it++ >= max_it) || (arma::norm(resid, 2) < conv_threshold)) {
      break;
    }

    // update v
    const double prev_resid_norm = resid_norm;
    resid_norm = arma::norm(resid, 2);
    const double direction_update_sqrt = resid_norm / prev_resid_norm;
    step_dir = resid + (direction_update_sqrt * direction_update_sqrt) * step_dir;
  }

  return it;
}

//! Choleskey decomposition of a subset of the rows/columns of a symmetric positive-semidefinite matrix *A*.
class Cholesky {
 public:
  //! Initialize an *empty* Cholesky decomposition for matrix *matrix*.
  //!
  //! @param matrix the matrix to decompose.
  //! @param max_active the maximum number of "active" indices, i.e., the maximum size of the subset of rows/columns.
  Cholesky(const arma::mat& matrix, const arma::uword max_active) noexcept
    : gram_(matrix), max_active_(max_active), active_size_(0), active_cols_(max_active_),
      gram_decomp_packed_(new double[max_active * (max_active + 1) / 2]) {}

  //! Copy constructor, optionally resetting the decomposition to be empty.
  Cholesky(const Cholesky& other, const bool reset) noexcept
    : gram_(other.gram_), max_active_(other.max_active_),
      active_cols_(reset ? arma::uvec(max_active_) : other.active_cols_),
      gram_decomp_packed_(new double[max_active_ * (max_active_ + 1) / 2]) {
    if (!reset) {
      std::copy(other.gram_decomp_packed_.get(), other.gram_decomp_packed_.get() + max_active_ * (max_active_ + 1) / 2,
                gram_decomp_packed_.get());
    }
  }

  //! Copy constructor.
  Cholesky(const Cholesky& other) noexcept : Cholesky(other, false) {}

  //! Copy assignment.
  Cholesky& operator=(const Cholesky& other) {
    gram_ = other.gram_;
    max_active_ = other.max_active_;
    active_size_ = other.active_size_;
    active_cols_ = other.active_cols_;
    gram_decomp_packed_.reset(new double[max_active_ * (max_active_ + 1) / 2]);
    std::copy(other.gram_decomp_packed_.get(), other.gram_decomp_packed_.get() + max_active_ * (max_active_ + 1) / 2,
              gram_decomp_packed_.get());
    return *this;
  }

  //! Default move constructor.
  Cholesky(Cholesky&& other) = default;

  //! Default move assignment operator.
  Cholesky& operator=(Cholesky&& other) = default;

  //! Update the diagonal of the matrix and reset the decomposition.
  //!
  //! @param add value to add to the diagonal of the matrix.
  void UpdateMatrixDiagonal(const double add) noexcept {
    gram_.diag() += add;
    Reset();
  }

  //! Update the diagonal of the matrix and reset the decomposition.
  //!
  //! @param add values to add to the diagonal of the matrix.
  void UpdateMatrixDiagonal(const arma::vec& add) {
    gram_.diag() += add;
    Reset();
  }

  //! Reset the decomposition, but retaining the matrix.
  void Reset() noexcept {
    // Not much has to be done. Only the active_size_ needs to be reset to 0.
    active_size_ = 0;
  }

  //! Add a row/column of the matrix to the Cholesky decomposition.
  //!
  //! The row/column is added at the end.
  //! For example, if we start with an empty decomposition and add row/column #3 and then #1, the decomposition is
  //! for the permutated matrix *A_[(3, 1), (3, 1)]*!
  //!
  //! @param col the column index in the interval ``[0, matrix.n_cols) ``.
  //! @return ``true`` if the column was added, ``false`` if the column was not added because either it would make the
  //!         matrix singular and the Cholesky decomposition ill-defined or the gram matrix is at it's maximal size.
  bool Add(const arma::uword add) noexcept {
    const double sq_norm_new_x = gram_(add, add);
    const double norm_new_x = std::sqrt(sq_norm_new_x);

    if (active_size_ == 0) {
      // No active variables yet and the decomposition is empty.
      gram_decomp_packed_[0] = norm_new_x;
    } else if (active_size_ >= max_active_) {
      return false;
    } else {
      char upper = 'U';
      char trans_y = 'T';
      char diag_n = 'N';
      const arma::blas_int mat_size = arma::blas_int(active_size_);
      const arma::blas_int incx = 1;

      // Get a view to the next column (really only needed for easy computation of the norm below).
      double *next_column = &gram_decomp_col(active_size_);
      arma::vec l(next_column, active_size_, false, true);
      l = gram_.unsafe_col(add).elem(active_cols_.head(active_size_));

      // Solve the triangular system of linear equations
      dtpsv(&upper, &trans_y, &diag_n, &mat_size, gram_decomp_packed_.get(), l.memptr(), &incx);

      next_column += active_size_;  //< Now points to the diagonal element of this column.
      *next_column = sq_norm_new_x - arma::dot(l, l);

      // Check if the decomposition is singular (by machine precision)
      if (*next_column < std::numeric_limits<double>::epsilon()) {
        // Decomposition is singular. Don't add column!
        return false;
      } else {
        *next_column = std::sqrt(*next_column);
      }
    }

    active_cols_[active_size_++] = add;
    return true;
  }

  //! Drop one or more rows/columns from the decomposition.
  //!
  //! The indices to drop must be in terms of the added rows/columns.
  //! For example, if rows/columns were added in the following order *5, 3, 7, 2, 10* and one wants to remove
  //! column #7, one needs to drop index *2*!
  //!
  //! @param first iterator pointing to the first (i.e., **largest** index **in terms of added rows/colums**) to be
  //!              dropped.
  //! @param last iterator pointing one past the last (i.e., **smallest** index) to be dropped.
  template<typename InputIterator>
  void Drop(InputIterator first, InputIterator last) noexcept {
    static_assert(std::is_integral<typename InputIterator::value_type>::value,
                  "Iterator must point to an integral type");
    // It is assumed that the iterator range is in **descending** order!
    while (first != last) {
      const arma::uword drop = *first++;

      // If drop is the last column, we only need to reduce the active size to effectively
      // "remove" the last column of the decomposition.
      // Otherwise, more adjustments are necessary.
      if (drop < active_size_ - 1) {
        // Compute the Cholesky one-rank update of the lower part of the decomposition
        // A = L_{3, 3} L'_{3, 3} and update with L'_{2, 3} L_{2, 3}
        // L_{3,3} is the lower-right part of the upper-diagonal matrix L, i.e.,
        // L_{3,3} = ((L[drop + 1, drop + 1] ... L[drop + 1, active_size_ - 1]),
        //            (L[drop + 2, drop + 2] ... L[drop + 2, active_size_ - 1]),
        //            ...
        //            (L[active_size_ - 1, active_size_ - 1]))
        // L_{2,3} is the *drop*-th row in the decomposition with elements above index `drop`,
        // i.e., L_{2,3} = (L[drop, drop + 1], L[drop, drop + 2], ..., L[drop, active_size_ - 1]).
        double *out = &gram_decomp_col(drop);  // Start writing at the the dropped column.
        double *in = out + drop + 1;  //< now points to the beginning of the next column.
        for (auto k = drop + 1; k < active_size_; ++k) {
          double *drop_row = in + drop;  //< points row to drop in the k-th column
          // Copy first chunk, i.e., all rows of the column until the dropped row.
          out = std::copy(in, drop_row, out);
          // Copy second chunk, i.e., all rows of the column after the dropped row
          // and until the diagonal
          in += k;  //< now points to the diagonal element in the k-th column
          out = std::copy(drop_row + 1, in, out);

          // Compute new diagonal element.
          *out = std::sqrt((*in) * (*in) + (*drop_row) * (*drop_row));
          const double c = (*out) / (*in);
          const double s = (*drop_row) / (*in);

          // Compute off-diagonal elements and "update" update-vector.
          double *update = in;
          for (auto i = k + 1; i < active_size_; ++i) {
            update += i;
            drop_row += i;
            *update = (*update + s * (*drop_row)) / c;
            *drop_row = c * (*drop_row) - s * (*update);
          }
          ++out;
          ++in;  //< Now points to the beginning of the next column.
        }
        // Remove entry from `active_cols_`. Can not use `.shed_row()` because we don't want to change the vector's
        // size!
        std::copy(&active_cols_[drop + 1], &active_cols_[active_size_], &active_cols_[drop]);
      }
      --active_size_;
    }
  }

  //! Solve a linear system of equations of the form
  //! *B x = b* where B is a subset of matrix *A* with rows/columns added/dropped previously **in the order
  //! of adding/dropping**!
  //!
  //! @param b on input the right-hand-side of the equations, on output the solution to the system.
  void Solve(arma::vec* b) const {
    char upper = 'U';
    char trans_y = 'T';
    char trans_n = 'N';
    char diag_n = 'N';
    const arma::blas_int mat_size = arma::blas_int(active_size_);
    const arma::blas_int incx = 1;
    dtpsv(&upper, &trans_y, &diag_n, &mat_size, gram_decomp_packed_.get(), b->memptr(), &incx);
    dtpsv(&upper, &trans_n, &diag_n, &mat_size, gram_decomp_packed_.get(), b->memptr(), &incx);
  }

  const arma::mat& matrix() const noexcept {
    return gram_;
  }

  //! Get the rows/columns that are currently "active" in the order of the composition.
  const arma::subview_col<arma::uword> active() const noexcept {
    return active_cols_.head(active_size_);
  }

  //! Get the number of "active" rows/columns in the decomposition.
  arma::uword active_size() const noexcept {
    return active_size_;
  }

 private:
  //! Access the first element in column *column* of the Cholesky decomposition of the matrix.
  inline double& gram_decomp_col(arma::uword column) {
    return gram_decomp_packed_[column * (column + 1) / 2];
  }

  arma::mat gram_;
  arma::uword max_active_;
  arma::uword active_size_;
  arma::uvec active_cols_;
  std::unique_ptr<double[]> gram_decomp_packed_;
};

//! Define a proxy to compute elementwise products in-place for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template <typename T>
void InplaceElementwiseProduct(const arma::vec& rhs, T* lhs) {
  *lhs %= rhs;
}

//! Define a proxy to compute elementwise products in-place for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template <typename T>
void InplaceElementwiseProduct(const double rhs, T* lhs) {
  *lhs *= rhs;
}

//! Define a proxy to compute elementwise products for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template <typename T>
inline auto ElementwiseProduct(const T& lhs, const arma::vec& rhs) {
  return lhs % rhs;
}

//! Define a proxy to compute elementwise products for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template <typename T>
inline auto ElementwiseProduct(const T& lhs, const arma::sp_vec& rhs) {
  return lhs % rhs;
}

//! Define a proxy to compute elementwise products for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template <typename T>
inline auto ElementwiseProduct(const T& lhs, const double& rhs) {
  return lhs * rhs;
}

//! Define a proxy to compute elementwise products for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template<>
inline auto ElementwiseProduct<double>(const double& scalar, const arma::vec& rhs) {
  return scalar * rhs;
}

//! Define a proxy to compute elementwise products for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template<>
inline auto ElementwiseProduct<double>(const double& scalar, const arma::sp_vec& rhs) {
  return scalar * rhs;
}

//! Define a proxy to compute elementwise products for "any" type of left-hand-side and
//! a vector-type right-hand-side.
//! @param lhs the left-hand-side
//! @param rhs the right-hand-side vector
template<>
inline auto ElementwiseProduct<double>(const double& scalar, const double& rhs) {
  return scalar * rhs;
}
}  // namespace linalg
}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_LINEAR_ALGEBRA_UTILITIES_HPP_
