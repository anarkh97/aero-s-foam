#ifndef ROM_CHOLESKYUTILS_H
#define ROM_CHOLESKYUTILS_H

template <typename Scalar> class GenFullSquareMatrix;

namespace Rom {

// Replaces the upper-diagonal part of a symmetric matrix M by its upper Cholesky factor
// M -> U such that M = U^T * U
template <typename Scalar>
const GenFullSquareMatrix<Scalar> &
cholesky_factor_upper(GenFullSquareMatrix<Scalar> &);

// Replaces the lower-diagonal part of a symmetric matrix M by its lower Cholesky factor
// M -> L such that M = L * L^T
template <typename Scalar>
const GenFullSquareMatrix<Scalar> &
cholesky_factor_lower(GenFullSquareMatrix<Scalar> &);

// Solves the factorized linear system using the upper triangular part
// v -> U^{-1} U^{-T} v
template <typename Scalar>
const Scalar *
cholesky_solve_upper(const GenFullSquareMatrix<Scalar> &, Scalar *);

// Solves the factorized linear system using the lower triangular part
// v -> L^{-T} L^{-1} v
template <typename Scalar>
const Scalar *
cholesky_solve_lower(const GenFullSquareMatrix<Scalar> &, Scalar *);

// Replaces the upper triangular matrix by its inverse
// U -> U^{-1}
template <typename Scalar>
const GenFullSquareMatrix<Scalar> &
inverse_triangular_upper(GenFullSquareMatrix<Scalar> &);

// Replaces the lower triangular matrix by its inverse
// L -> L^{-1}
template <typename Scalar>
const GenFullSquareMatrix<Scalar> &
inverse_triangular_lower(GenFullSquareMatrix<Scalar> &);

} // end namespace Rom

#endif /* ROM_CHOLESKYUTILS_H */
