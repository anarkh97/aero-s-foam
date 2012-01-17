// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#ifndef ContactAlgorithm_h
#define ContactAlgorithm_h

#include <cmath>
#include <iostream>
#include "Contact_Defines.h"

namespace acme {

static const Real colinearity_tolerance = 0.99;

inline Real colinearityTolerance() 
{
  return colinearity_tolerance;
}

static const Real determinant_tolerance = 1.0e-30;

inline Real determinantTolerance() {
  return determinant_tolerance;
};

/*!
 * Compute the cross product of a and b.
 */
inline void Cross(const Real a[3], const Real b[3], Real c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/*!
 * Cross product where the first 3 real values represent the vector a
 * and the second set of three real values represents the vector b.
 */
inline void Cross(const Real dx1, const Real dy1, const Real dz1,
                  const Real dx2, const Real dy2, const Real dz2,
		  Real normal[3])
{
  normal[0] = dy1*dz2-dz1*dy2;
  normal[1] = dz1*dx2-dx1*dz2;
  normal[2] = dx1*dy2-dy1*dx2;
}

/*!
 * Compute the scalar triple product = (a X b).c
 */
inline Real ScalarTripleProduct(const Real a[3], const Real b[3], const Real c[3])
{
  return  (a[1]*b[2] - a[2]*b[1])*c[0]
        + (a[2]*b[0] - a[0]*b[2])*c[1]
        + (a[0]*b[1] - a[1]*b[0])*c[2];
}

namespace optimized {

  /*!
   * Compute the cross product of a and b.
   */
  inline void Cross(const Real a[3], const Real b[3], Real c[3])
  {
    const Real u1 = a[0];
    const Real u2 = a[1];
    const Real u3 = a[2];

    const Real v1 = b[0];
    const Real v2 = b[1];
    const Real v3 = b[2];

    const Real t1 = u1 - u2;
    const Real t2 = v2 + v3;
    const Real t3 = u1*v3;
    const Real t4 = t1*t2 - t3;

    c[0] = v2*(t1-u3) - t4;
    c[1] = u3*v1 - t3;
    c[2] = t4 - u2 *(v1 - t2);
  }

} // end namespace optimized

inline bool isZero(const Real a[3])
{
  return (   a[0] == 0.0
          && a[1] == 0.0
	  && a[2] == 0.0);
}

inline void Zero(Real a[3])
{
  a[0] = 0.0;
  a[1] = 0.0;
  a[2] = 0.0;
}

inline void Scale(Real a[3], const Real scale)
{
  a[0] *= scale;
  a[1] *= scale;
  a[2] *= scale;
}

inline void Copy(const Real source[3], Real destination[3])
{
  destination[0] = source[0];
  destination[1] = source[1];
  destination[2] = source[2];
}

/*!
 * Compute the sum of all elements of a 3D vector.
 */
inline Real Sum(const Real a[3])
{
  return (a[0]+a[1]+a[2]);
}

/*!
 * Compute the sum of the absolute value of all elements of a 3D vector.
 */
inline Real AbsSum(const Real a[3])
{
  return (std::fabs(a[0])+std::fabs(a[1])+std::fabs(a[2]));
}

/*!
 * Compute the dot product.
 */
inline Real Dot(const Real a[3])
{
  return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

/*!
 * Compute the dot product between two 3D vectors.
 */
inline Real Dot(const Real a[3], const Real b[3])
{
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

/*!
 * Compute the magnitude of a 3D vector.
 */
inline Real Magnitude(const Real a[3])
{
  return std::sqrt(Dot(a));
}

/*!
 * Compute the squared magnitude of a 3D vector.
 */
inline Real Magnitude2(const Real a[3])
{
  return Dot(a);
}

/*!
 * Compute the distance between vectors a and b.
 */
inline Real Distance(const Real a[3], const Real b[3])
{
  return std::sqrt((b[0]-a[0])*(b[0]-a[0]) +
              (b[1]-a[1])*(b[1]-a[1]) +
              (b[2]-a[2])*(b[2]-a[2]));
}

/*!
 * Compute the squared distance between vectors a and b.
 */
inline Real Distance2(const Real a[3], const Real b[3])
{
  return ((b[0]-a[0])*(b[0]-a[0]) +
          (b[1]-a[1])*(b[1]-a[1]) +
          (b[2]-a[2])*(b[2]-a[2]));
}

/*!
 * Normalize a 3D vector. v = v /||v||
 *
 * Return the magnitude of the input vector.
 */
inline Real Normalize(Real v[3])
{
  double mag = Magnitude(v);
  if ( mag != 0.0 ) {
    Real invMag = 1.0 / mag;
    v[0] *= invMag;
    v[1] *= invMag;
    v[2] *= invMag;
  }
  return mag;
}

/*!
 * Compute the area of the triangle formed from the 2 vectors a and b.
 */
inline Real Area(const Real a[3], const Real b[3])
{
  Real c[3];
  Cross(a, b, c);
  return 0.5*Magnitude(c);
}

inline Real det2x2(const Real A[2][2])
{
  // A = | a b |
  //   = | c d |
  return (A[0][0]*A[1][1] - A[1][0]*A[0][1]);
}

inline Real det2x2(Real a, Real b,
                   Real c, Real d)
{
  // A = | a b |
  //   = | c d |
  return (a*d - b*c);
}

inline void Solve_2x2_Matrix(const Real A[2][2], const Real x[2], Real y[2])
{
  Real a = A[0][0], b = A[0][1];
  Real c = A[1][0], d = A[1][1];

  Real det = det2x2(a,b,c,d);

  if ( det != 0.0 ) {
    Real inv_det = 1.0 / det;
    Real x0_det = x[0] * d - x[1] * b;
    Real x1_det = a * x[0] - c * x[1];
    y[0] = x0_det*inv_det;
    y[1] = x1_det*inv_det;
  } else {
    // error! singular matrix.
  }
}

inline void MatrixVectorProduct(const Real M[3][3], const Real x[3], Real y[3])
{ 
  Real M11 = M[0][0], M12 = M[0][1], M13 = M[0][2];
  Real M21 = M[1][0], M22 = M[1][1], M23 = M[1][2];
  Real M31 = M[2][0], M32 = M[2][1], M33 = M[2][2];

  Real x1 = x[0];
  Real x2 = x[1];
  Real x3 = x[2];

  y[0] = M11*x1 + M12*x2 + M13*x3;
  y[1] = M21*x1 + M22*x2 + M23*x3; 
  y[2] = M31*x1 + M32*x2 + M33*x3;
}

inline void SymMatrixVectorProduct(const Real M[3][3], const Real x[3], Real y[3])
{ 
  Real M11 = M[0][0], M12 = M[0][1], M13 = M[0][2];
  Real                M22 = M[1][1], M23 = M[1][2];
  Real                               M33 = M[2][2];

  Real x1 = x[0], x2 = x[1], x3 = x[2];

  y[0] = M11*x1 + M12*x2 + M13*x3;
  y[1] = M12*x1 + M22*x2 + M23*x3; 
  y[2] = M13*x1 + M23*x2 + M33*x3;
}

inline Real Invert_3x3_Sym_Matrix(const Real M [3][3], Real M_inv [3][3], const char *msg="")
{
  Real M11 = M[0][0], M12 = M[0][1], M13 = M[0][2];
  Real M21 = M[1][0], M22 = M[1][1], M23 = M[1][2];
  Real M31 = M[2][0], M32 = M[2][1], M33 = M[2][2];
  //
  // Cofactors
  //
  Real A11 = M22*M33 - M23*M32;
  Real A12 = M13*M32 - M12*M33;
  Real A13 = M12*M23 - M13*M22;
  Real A22 = M11*M33 - M13*M31;
  Real A23 = M13*M21 - M11*M23;
  Real A33 = M11*M22 - M12*M21;
  //
  // Determinant
  //
  Real det = M11*A11 + M12*A12 + M13*A13 ;
  if( std::fabs(det) < determinant_tolerance ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }
  Real det_inv = 1.0/det;
  M_inv[0][0] = A11*det_inv ;
  M_inv[1][0] = A12*det_inv ;
  M_inv[2][0] = A13*det_inv ;

  M_inv[0][1] = M_inv[1][0];
  M_inv[1][1] = A22*det_inv ;
  M_inv[2][1] = A23*det_inv ;

  M_inv[0][2] = M_inv[2][0];
  M_inv[1][2] = M_inv[2][1];
  M_inv[2][2] = A33*det_inv ;

  return det;
}

inline bool isEqual(const double x, const double y)
{
  const double epsilon = 1e-5; /* some small number such as 1e-5 */;
  return std::abs(x - y) <= epsilon * std::abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

inline bool isNormal(const Real a[3])
{
  //
  // KHP: Will floating point precision haunt me on this comparison?
  //
  // return ( std::fabs(Magnitude(a) - 1.0 ) < 1.0e-8 );
  //
  return isEqual( Magnitude(a), 1.0 );
}

inline void print(const Real a[3], const char *name="a")
{
  std::cout << name << "(" << a[0] << "," << a[1] << "," << a[2] << ")\n";
}

inline void printMag(const Real a[3], const char *name="a")
{
  std::cout << "(" << a[0] << "," << a[1] << "," << a[2] << "), ||" << name << "||= " << Magnitude(a) << std::endl;
}

inline Real Invert_2x2_Matrix(const Real M[2][2], Real M_inv[2][2], const char *msg="")
{
  Real a = M[0][0], b = M[0][1];
  Real c = M[1][0], d = M[1][1];

  Real det = a*d - b*c;
  if( std::fabs(det) < determinant_tolerance) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }

  Real inv_det = 1.0 / det;

  M_inv[0][0] =  d * inv_det; M_inv[0][1] = -b * inv_det;
  M_inv[1][0] = -c * inv_det; M_inv[1][1] =  a * inv_det;

  return det;
}

inline Real Invert_2x2_Sym_Matrix(const Real M[2][2], Real M_inv[2][2], const char *msg="")
{
  //
  // if ( M[0][1] != M[1][0] ) error, should have called regular
  // Invert_2x2_Matrix method.
  //
  Real a = M[0][0], b = M[0][1];
  Real              d = M[1][1];

  Real det = a*d - b*b;
  if( std::fabs(det) < determinant_tolerance ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }

  Real inv_det = 1.0 / det;

  M_inv[0][0] =  d * inv_det; M_inv[0][1] = -b * inv_det;
  M_inv[1][0] =  M_inv[0][1]; M_inv[1][1] =  a * inv_det;

  return det;
}

inline Real Invert_3x3Matrix(const Real M [3][3], Real M_inv [3][3], const char *msg="")
{
  //
  // Cofactors
  //
  Real cofM [3] [3];
  cofM[0][0] =   M[1][1]*M[2][2] - M[1][2]*M[2][1] ;
  cofM[1][0] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) ;
  cofM[2][0] =   M[0][1]*M[1][2] - M[0][2]*M[1][1] ;
  cofM[0][1] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) ;
  cofM[1][1] =   M[0][0]*M[2][2] - M[0][2]*M[2][0] ;
  cofM[2][1] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) ;
  cofM[0][2] =   M[1][0]*M[2][1] - M[1][1]*M[2][0] ;
  cofM[1][2] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) ;
  cofM[2][2] =   M[0][0]*M[1][1] - M[0][1]*M[1][0] ;
  //
  // Determinant
  //
  Real det = M[0][0]*cofM[0][0] + M[0][1]*cofM[0][1] + M[0][2]*cofM[0][2] ;
  if( std::fabs(det) < determinant_tolerance ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << msg << ", Determinant = 0" << std::endl;
#endif
    return 0;
  }
  Real det_inv = 1.0/det  ;
  M_inv[0][0] = cofM[0][0]*det_inv ;
  M_inv[1][0] = cofM[0][1]*det_inv ;
  M_inv[2][0] = cofM[0][2]*det_inv ;
  M_inv[0][1] = cofM[1][0]*det_inv ;
  M_inv[1][1] = cofM[1][1]*det_inv ;
  M_inv[2][1] = cofM[1][2]*det_inv ;
  M_inv[0][2] = cofM[2][0]*det_inv ;
  M_inv[1][2] = cofM[2][1]*det_inv ;
  M_inv[2][2] = cofM[2][2]*det_inv ;
  return det;
}

//
//  Compute ACME edge curvature from a pair of connected face normals
//  Curvature is an artifical construct where 0.0 is flat, and 2.0 is
//  fully opposed as:
//
//    curv = 1.0 - Theta = 1.0 - norm1.norm2
//
inline Real ComputeCurvatureFromNormals(const Real *norm1, const Real *norm2) {
  return (1.0 - Dot(norm1, norm2));
}




} // end namespace acme

#endif
