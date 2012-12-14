#ifndef __UTILITIES__HPP__
#define __UTILITIES__HPP__
#include <cmath>
#include <iostream>
#include <limits>
#include <complex>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>

const double epsilon2 = std::numeric_limits<double>::epsilon()*std::numeric_limits<double>::epsilon();

template<typename Scalar>
void
leftmult_rotvar(int num_nodes, int itrans, const Eigen::Matrix<Scalar,3,3> rotvar[3], 
                Eigen::Map<Eigen::Matrix<Scalar,18,18,Eigen::RowMajor> > &stiff)

/*****************************************************************
 *
 *  Purpose:
 *     Multiply | I  0 |*stiff or | I  0 |*stiff 
 *              | 0  V'|          | 0  V |
 *     where V is rotvar or variation of pseudovector rotation 
 *     with respect to instantanious rotation
 *
 *  Input:
 *     num_nodes : number of nodes
 *     itrans    : = 1 use transpose V
 *     rotvar    : rotvar[i][j][k] is component j k of rotation
 *                 gradient matrix for node i
 *     stiff     : stiffness
 *  Output:
 *     stiff     : stiff updated with the rotation gradient product
 *
 *  Coded by: Bjorn Haugen
 *******************************************************************/
{

   int i, j, inod, ndofs, ir;
   Scalar scr[3];

// Transform groups of 3 degrees of freedom

   ndofs = 6*num_nodes;

   if ( itrans == 1 ) {                       /* Use Transpose rotvar */
      for( inod=0; inod<num_nodes; inod++ ) {
         ir = 6*inod + 3; 
         for( j=0; j<ndofs; j++ ) {
            for( i=0; i<3; i++ ) {
                  scr[i] =  rotvar[inod](0,i)*stiff(ir +0,j)
                           +rotvar[inod](1,i)*stiff(ir +1,j)
                           +rotvar[inod](2,i)*stiff(ir +2,j);
            }
	    stiff(ir+0,j) = scr[0];
	    stiff(ir+1,j) = scr[1];
            stiff(ir+2,j) = scr[2];
         }
      }
   }

   else {                                     /*  Use rotvar */
      for( inod=0; inod<num_nodes; inod++ ) {
         ir = 6*inod + 3;
         for( j=0; j<ndofs; j++ ) {
            for( i=0; i<3; i++ ) {
                  scr[i] =  rotvar[inod](i,0)*stiff(ir +0,j)
                           +rotvar[inod](i,1)*stiff(ir +1,j)
                           +rotvar[inod](i,2)*stiff(ir +2,j);
            }
            stiff(ir+0,j) = scr[0];
            stiff(ir+1,j) = scr[1];
            stiff(ir+2,j) = scr[2];
         }
      }
   }

}

template<typename Scalar>
void rightmult_rotvar(int num_nodes, int itran, const Eigen::Matrix<Scalar,3,3> rotvar[3],
                      Eigen::Map<Eigen::Matrix<Scalar,18,18,Eigen::RowMajor> > &stiff)
/*****************************************************************
 *
 *  Purpose:
 *     Multiply stiff*| I  0 | or stiff*| I  0 |
 *                    | 0  V |          | 0  V'|
 *     where V is rotvar  or variation of pseudovector rotation 
 *     with respect to instantanious rotation
 *
 *  Input:
 *     num_nodes : number of nodes
 *     itran     : =1 user transpose rotvar
 *     rotvar    : rotvar[i][j][k] is compinent j k of rotation
 *                 gradient matrix for node i
 *     stiff     : stiffness
 *  Output:
 *     stiff     : stiff updated with the rotation gradient product
 *
 *  Coded by: Bjorn Haugen
 *******************************************************************/
{

   int i, j, ndofs, jnod, jr;
   Scalar scr[3];

// Transform groups of 3 degrees of freedom

   ndofs = 6*num_nodes;

   if( itran == 1 ) {                       /*  Use transpose rotvar */
      
      for( jnod=0; jnod<num_nodes; jnod++ ) {
         jr = 6*jnod +3;
         for( i=0; i<ndofs; i++ ) {
            for( j=0; j<3; j++ ) {
                  scr[j] =  stiff(i,jr +0)*rotvar[jnod](j,0)
                           +stiff(i,jr +1)*rotvar[jnod](j,1)
                           +stiff(i,jr +2)*rotvar[jnod](j,2);
            }
            stiff(i,jr+0) = scr[0];
            stiff(i,jr+1) = scr[1];
            stiff(i,jr+2) = scr[2];
         }
      }

   }

   else {                                     /*  Use rotvar */
      for( jnod=0; jnod<num_nodes; jnod++ ) {
         jr = 6*jnod +3;
         for( i=0; i<ndofs; i++ ) {
            for( j=0; j<3; j++ ) {
                  scr[j] =  stiff(i,jr +0)*rotvar[jnod](0,j)
                           +stiff(i,jr +1)*rotvar[jnod](1,j)
                           +stiff(i,jr +2)*rotvar[jnod](2,j);
            }
            stiff(i,jr+0) = scr[0];
            stiff(i,jr+1) = scr[1];
            stiff(i,jr+2) = scr[2];
         }
      }
   }

}

template<typename Scalar>
void pseudorot_var(const Eigen::Matrix<Scalar,3,1> &rvec, Eigen::Matrix<Scalar,3,3> &varmat)
/******************************************************************
 *
 *  Purpose: Compute the variation of a rotation pseudo vector r with
 *  respect to incremental rotation w where the rotations obey the rule
 *  
 *      R(r+w) = R(w)*R(r)
 *
 *  The variation is given as V = I -Spin/2 +eta*Spin*Spin
 *
 *  Input:
 *     rvec : initial rotation vector
 *
 *  Output:
 *     varmat: matrix of pseudovector variations with repspect to
 *             w = [eps, 0, 0], w = [0, eps, 0] and w = [0, 0, eps]
 *
 *  Coded by: Bjorn Haugen; adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int i, j;
   Scalar eta, th, sthh, cthh;
   Eigen::Matrix<Scalar,3,3> spin;

   Scalar th2 = rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];

// Check small value of theta^2 : th2 = 5e-6 gives small error between ex.
// and approximate solution

   if ( th2 < 5e-6 ) {
      eta = (1.0/12.0) + (1.0/720.0)*th2 + (1.0/30240.0)*(th2*th2);
   } else {
      using std::sqrt;
      th = sqrt( th2 );
      using std::sin;
      using std::cos;
      sthh = sin( 0.5*th );
      cthh = cos( 0.5*th );
      eta  = (sthh - 0.5*th*cthh)/( th2*sthh );
   }

// Compute the spin of the vector
   spin(0,0) = 0.0;
   spin(1,1) = 0.0;
   spin(2,2) = 0.0;

   spin(1,2) = -rvec[0];
   spin(0,2) =  rvec[1];
   spin(0,1) = -rvec[2];

   spin(2,1) =  rvec[0];
   spin(2,0) = -rvec[1];
   spin(1,0) =  rvec[2];

   // Compute first part of var = I - Spin/2

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) 
        varmat(i,j) = -0.5*spin(i,j);
      varmat(i,i) = 1.0;
   }

// Add the term eta*Spin*Spin
   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         varmat(i,j) += eta*(spin(i,0)*spin(0,j)
                            +spin(i,1)*spin(1,j)
                            +spin(i,2)*spin(2,j));
      }
   }

}

template<typename Scalar>
void pseudorot_2var(const Eigen::Matrix<Scalar,3,1> &r, const Eigen::Matrix<Scalar,3,1> &f, Eigen::Matrix<Scalar,3,3> &scndvar)
/******************************************************************
 *
 *  Purpose: Compute the variation of a contraction between
 *      Rmat'*f  where f is constant vector and Rmat is the
 *      rotation variation of a pseudovector r with respect
 *      to instantanious rotations
 *      
 *  Input:
 *     r : initial rotation vector ( vector to be varied );
 *     f : constant vector that is contracted with the varied
 *         rotation varitation matrix
 *
 *  Output:
 *     scndvar: matrix of the variation of the contraction with respect to
 *             w = [eps, 0, 0], w = [0, eps, 0] and w = [0, 0, eps]
 *
 *  Coded by: Bjorn Haugen; adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int i, j;
   Scalar eta, nu, th, sth, sthh, cthh, fac;
   Eigen::Matrix<Scalar,3,3> mat, sr, fr, spsq, fstvar;

   Scalar th2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

/* Check small value of theta^2 : th2 = 5e-6 gives small error between ex.
 * and approximate solution for eta
 */
 
   if ( th2 < 5e-6 ) {
      Scalar th4 = th2*th2;
      eta = (1.0/12.0)  + th2/720  + th4/30240;
      nu  = (1.0/360.0) + th2/7560 + th4/201600;
   } else {
      using std::sqrt;
      th = sqrt ( th2 );
      using std::sin;
      using std::cos;
      sth  = sin( th );
      sthh = sin( 0.5*th );
      cthh = cos( 0.5*th );
      eta  = (sthh - 0.5*th*cthh)/( th2*sthh );
      nu   = (th*(th + sth) - 8.0*sthh*sthh )/(4.0*(th2*th2)*sthh*sthh);
   }

/* Compute the spin of the rotation vectors */
   sr(0,0) = 0.0;
   sr(1,1) = 0.0;
   sr(2,2) = 0.0;

   sr(1,2) = -r[0];
   sr(0,2) =  r[1];
   sr(0,1) = -r[2];

   sr(2,1) =  r[0];
   sr(2,0) = -r[1];
   sr(1,0) =  r[2];


// Compute first part of mat =  -SpinF/2
   mat(0,0) = 0.0;
   mat(1,1) = 0.0;
   mat(2,2) = 0.0;

   mat(1,2) =  0.5*f[0];
   mat(0,2) = -0.5*f[1];
   mat(0,1) =  0.5*f[2];

   mat(2,1) = -0.5*f[0];
   mat(2,0) =  0.5*f[1];
   mat(1,0) = -0.5*f[2];

/* Add the terms  mat += eta*( (r'*f)I +r*f' -2f*r') */
   fac = r[0]*f[0] + r[1]*f[1] + r[2]*f[2];

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         spsq(i,j) = sr(i,0)*sr(0,j)
                    +sr(i,1)*sr(1,j)
                    +sr(i,2)*sr(2,j);
         fr(i,j) = f[i]*r[j];
         mat(i,j) += eta*( r[i]*f[j] - 2.0*f[i]*r[j] );
      }
      mat(i,i) += eta*fac;
   }

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         mat(i,j) += nu*( spsq(i,0)*fr(0,j)
                         +spsq(i,1)*fr(1,j)
                         +spsq(i,2)*fr(2,j) );
      }
   }

// Multiply with first variation
   pseudorot_var( r, fstvar );

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         scndvar(i,j) =  mat(i,0)*fstvar(0,j)
                        +mat(i,1)*fstvar(1,j)
                        +mat(i,2)*fstvar(2,j);
      }
   }

}

// PJSA: 5/2/2011 I observe much better results with the eigen version
// which uses algorithm from "Quaternion Calculus and Fast Animation",
// Ken Shoemake, 1987 SIGGRAPH course notes
#define USE_EIGEN_MAT_TO_QUAT_SCALAR

template<typename Scalar>
void mat_to_quat(const Eigen::Matrix<Scalar,3,3> &rten, Eigen::Matrix<Scalar,4,1> &q)
/*****************************************************************
 *  Compute the quaternion representation of a rotation tensor
 *  by means of equations from Richard A. Spurrier comments on
 *  "Singularity free Extraction of a Quaternion from a 
 *   Direction-Cosine Matrix" in J. of Rockets and Spacecraft 1978.
 *
 *  Input:
 *  rten:    rotation tensor
 *
 *  Output:
 *  q:    rotation quaternion ( Euler parameters )
 *
 *  Return:
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
#ifdef USE_EIGEN_MAT_TO_QUAT_SCALAR
      // using this version causes reparameterization of rotations at -2/3*pi,4/3*pi
      Eigen::Quaternion<Scalar> qq(rten);
      q[0] = qq.w();
      q[1] = qq.x();
      q[2] = qq.y();
      q[3] = qq.z();
#else
      // using this routine causes reparameterization of rotations at -1/2*pi,3/2*pi
      using std::sqrt;

      static int p[5] = {0,1,2,0,1};
      int    imax, i, j, k;
      Scalar trace;

      trace = rten(0,0) +rten(1,1) +rten(2,2);

      imax = 0;
      if ( rten(1,1) > rten(imax,imax) ) imax = 1;
      if ( rten(2,2) > rten(imax,imax) ) imax = 2;

      if ( trace > rten(imax,imax) ) {
         q[0] = sqrt( 1.0 +trace )/2.0;
         q[1] = ( rten(2,1) -rten(1,2) )/(4.0*q[0]);
         q[2] = ( rten(0,2) -rten(2,0) )/(4.0*q[0]);
         q[3] = ( rten(1,0) -rten(0,1) )/(4.0*q[0]);
      }
      else {
         i = p[imax];
         j = p[imax+1];
         k = p[imax+2];
         q[i+1] = sqrt( rten(i,i)/2.0 +(1.0 -trace)/4.0 );
         q[  0] = ( rten(k,j) -rten(j,k) )/(4.0*q[i+1]);
         q[j+1] = ( rten(j,i) +rten(i,j) )/(4.0*q[i+1]);
         q[k+1] = ( rten(k,i) +rten(i,k) )/(4.0*q[i+1]);
      }
/*
      // PJSA 2-12-09 this is a variant due to Markley (Journal of Guidance and Control Vol 31 No 2 pp 440-442, 2008)
      // using this version causes a reparameterization of the rotation at -pi,+pi
      if ( trace > rten(imax,imax) ) {
         q[0] = ( 1.0 +trace );
         q[1] = ( rten(2,1) -rten(1,2) );
         q[2] = ( rten(0,2) -rten(2,0) );
         q[3] = ( rten(1,0) -rten(0,1) );
      }
      else {
         i = p[imax];
         j = p[imax+1];
         k = p[imax+2];
         q[i+1] = 1.0 + rten(i,i) - rten(j,j) - rten(k,k);
         q[  0] = ( rten(k,j) -rten(j,k) );
         q[j+1] = ( rten(j,i) +rten(i,j) );
         q[k+1] = ( rten(k,i) +rten(i,k) );
      }
      q.normalize();
      //Scalar norm = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
      //for(int l=0; l<4; ++l) q[l] /= norm;
*/
/*
      // PJSA 3-12-09 method due to Bar-Itzhack (Journal of Guidance and Control Vol 23 No 6 pp 1085-1087, 2000)
      // which always gives the quaternion that corresponds to the orthogonal matrix closest to rten
      // this version causes a reparameterization of the rotation at -3/2*pi, 1/2*pi
      Eigen::Matrix<Scalar,4,4> A;
        A << (rten(0,0)-rten(1,1)-rten(2,2))/3.0, (rten(0,1)+rten(1,0))/3.0, (rten(0,2)+rten(2,0))/3.0, -(rten(2,1)-rten(1,2))/3.0,
              0.0, (rten(1,1)-rten(0,0)-rten(2,2))/3.0, (rten(1,2)+rten(2,1))/3.0, -(rten(0,2)-rten(2,0))/3.0,
              0.0, 0.0, (rten(2,2)-rten(0,0)-rten(1,1))/3.0, -(rten(1,0)-rten(0,1))/3.0,
              0.0, 0.0, 0.0, (rten(0,0)+rten(1,1)+rten(2,2))/3.0;

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar,4,4> > dec(A.template selfadjointView<Eigen::Upper>());
      Eigen::Matrix<Scalar,4,1> y = dec.eigenvectors().col(3);

      q[0] = y[3];
      q[1] = -y[0];
      q[2] = -y[1];
      q[3] = -y[2];
      if(q[0] < 0) for(int l=0; l<4; ++l) q[l] = -q[l];
*/
      //Scalar norm = std::sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
      //for(int l=0; l<4; ++l) q[l] /= norm;
#endif

  // always return the positive quaternion.
  // see: "INTERPOLATION OF ROTATIONAL VARIABLES IN NONLINEAR DYNAMICS OF 3D BEAMS" 
  // by Jelenic and Crisfield, Int. J. Numer. Meth. Engng. 43, 1193–1222 (1998)
  // section 3.5.2. for discussion of uniqueness of the extracted rotational pseudovector.
  // also appendix B1 of "A THREE-DIMENSIONAL FINITE-STRAIN ROD MODEL
  // PART II: COMPUTATIONAL ASPECTS" by Simo and Vu-Quoc, Computer Methods in Applied 
  // Mechanics and Engineering, Volume 58, Issue 1, October 1986

  if(q[0] < 0) {
    for(int l=0; l<4; ++l) q[l] = -q[l]; // this amounts to reparameterizing at -pi,pi
  }

}

//#define USE_MATRIX_LOGARITHM
template<typename Scalar>
void mat_to_vec(const Eigen::Matrix<Scalar,3,3> &rten, Eigen::Matrix<Scalar,3,1> &rvec)
/*****************************************************************
 *  Compute the rotation vector from a rotation tensor
 *
 *  Input:
 *  rten:    rotation tensor
 *
 *  Output:
 *  rvec:    rotation vector
 *
 *  Return:
 *****************************************************************/
{
#ifdef USE_MATRIX_LOGARITHM
  //another possibility is to use the matrix logarithm
  if(rten == Eigen::Matrix<Scalar,3,3>::Identity()) rvec = Eigen::Matrix<Scalar,3,1>::Zero();
  else {
    Eigen::Matrix<Scalar,3,3> skew;
    // method #1: use the matrix logarithm function in unsupported Eigen MatrixFunctions module
    skew = rten.log();
/*
    // method #2: use eigen decomposition (note: that rten is not self-adjoint)
    Eigen::EigenSolver<Eigen::Matrix<Scalar,3,3> > dec(rten);
    skew = (dec.eigenvectors() * dec.eigenvalues().array().log().matrix().asDiagonal() * dec.eigenvectors().adjoint()).real();
*/
    rvec[0] = skew(2,1);
    rvec[1] = skew(0,2);
    rvec[2] = skew(1,0);
  }
#else
      using std::sqrt;
      using std::asin;
      using std::acos;

      Eigen::Matrix<Scalar,4,1> q;
      Scalar th, coef;

      mat_to_quat<Scalar>( rten, q );

      Scalar sthh2 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];

      if( sthh2 <= epsilon2 ) {
        coef = Scalar(2.0);
      }
      else {
        Scalar sthh = sqrt ( sthh2 );
        if( sthh < 0.7 ) {
          sthh = sqrt( sthh2 );
          th = 2.0*asin( sthh );
        }
        else {
          Scalar cthh = q[0];
          th = 2.0*acos( cthh );
        }

        if( sthh > 1.0 ) {
          coef = th;
        }
        else {
          coef = th/sthh;
        }
      }

      rvec[0] = coef*q[1];
      rvec[1] = coef*q[2];
      rvec[2] = coef*q[3];
#endif
}

template<typename Scalar>
void vec_to_quat(const Eigen::Matrix<Scalar,3,1> &rvec, Eigen::Matrix<Scalar,4,1> &q)
/*****************************************************************
 *  Compute the quaterion representation from a rotation vector
 *
 *  Input:
 *  rvec:    rotation vector
 *
 *  Output:
 *  q :      quaterion representation of the rotation ( Euler parameters )
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
      using std::sqrt;
      using std::cos;
      using std::sin;
      Scalar coef;

      Scalar th2 = rvec.squaredNorm();

      if ( th2 <= epsilon2 ) {
        q[0] = Scalar(1);
        coef = Scalar(0.5);
      }
      else {
        if( th2 < 5e-6) {
          Scalar th4 = th2*th2;
          q[0] = 1 - th2/8 + th4/384;
          coef = 0.5 - th2/48 + th4/3840; 
        }
        else {
          using std::sqrt;
          Scalar th = sqrt(th2);
          using std::sin;
          using std::cos;
          q[0] = cos( 0.5*th );
          coef = sin(0.5*th)/th;
        }
      }

      q[1] = coef*rvec[0];
      q[2] = coef*rvec[1];
      q[3] = coef*rvec[2];
/* 
      // Compute norm of quaterion q
      Scalar norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
      if(norm == 0) {
        std::cerr << "WARNING: Skipping quaternion normalize( null )\n";
      }
      else {
        // Normalize quaterion q
        q[0] /= norm;
        q[1] /= norm;
        q[2] /= norm;
        q[3] /= norm;
      }
*/
      q = q.normalized(); // XXXX
      //q.normalize(); // TODO
}

#define USE_EIGEN_QUAT_TO_MAT_SCALAR

template<typename Scalar>
void quat_to_mat(const Eigen::Matrix<Scalar,4,1> &_q, Eigen::Matrix<Scalar,3,3> &rten)
/*****************************************************************
 *  Compute the rotation matrix from a quaternion representation
 *  of the rotation
 *
 *  Input:
 *  q:    rotation quaternion ( Euler parameters )
 *
 *  Output:
 *  rten:    rotation tensor
 *
 *  Return:
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
#ifdef USE_EIGEN_QUAT_TO_MAT_SCALAR
   Eigen::Quaternion<Scalar> q;
   q.w() = _q[0];
   q.x() = _q[1];
   q.y() = _q[2];
   q.z() = _q[3];
   q.normalize(); // TODO
   rten = q.toRotationMatrix();
#else
   // Compute norm of q
/*
   Scalar norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
   if(norm == 0) {
     std::cerr << "WARNING: Skipping quaternion normalize( null )\n";
   }
   else {
     // Normalize q
     q[0] /= norm;
     q[1] /= norm;
     q[2] /= norm;
     q[3] /= norm;
   }
*/
   Eigen::Matrix<Scalar,4,1> q = _q.normalized(); // XXXX

   rten(0,0) = 2.0*(q[1]*q[1] + q[0]*q[0]) - 1.0;
   rten(1,1) = 2.0*(q[2]*q[2] + q[0]*q[0]) - 1.0;
   rten(2,2) = 2.0*(q[3]*q[3] + q[0]*q[0]) - 1.0;

   rten(0,1) = 2.0*(q[1]*q[2] - q[3]*q[0]);
   rten(0,2) = 2.0*(q[1]*q[3] + q[2]*q[0]);
   rten(1,2) = 2.0*(q[2]*q[3] - q[1]*q[0]);

   rten(1,0) = 2.0*(q[2]*q[1] + q[3]*q[0]);
   rten(2,0) = 2.0*(q[3]*q[1] - q[2]*q[0]);
   rten(2,1) = 2.0*(q[3]*q[2] + q[1]*q[0]);

/*
   rten(0,0) = q[1]*q[1] - q[2]*q[2] - q[3]*q[3] + q[0]*q[0];
   rten(1,0) = 2.0*(q[1]*q[2] + q[3]*q[0]);
   rten(2,0) = 2.0*(q[1]*q[3] - q[2]*q[0]);
   rten(0,1) = 2.0*(q[1]*q[2] - q[3]*q[0]);
   rten(1,1) = -q[1]*q[1] + q[2]*q[2] - q[3]*q[3] + q[0]*q[0];
   rten(2,1) = 2.0*(q[2]*q[3] + q[1]*q[0]);
   rten(0,2) = 2.0*(q[1]*q[3] + q[2]*q[0]);
   rten(1,2) = 2.0*(q[2]*q[3] - q[1]*q[0]);
   rten(2,2) = -q[1]*q[1] - q[2]*q[2] + q[3]*q[3] + q[0]*q[0];
*/
#endif
}

//#define USE_MATRIX_EXPONENTIAL
template<typename Scalar>
void vec_to_mat(const Eigen::Matrix<Scalar,3,1> &rvec, Eigen::Matrix<Scalar,3,3> &rten)

/*****************************************************************
 *  Compute the rotation tensor from a rotation matrix
 *
 *  Input:
 *  rten:    rotation tensor
 *
 *  Output:
 *  rvec:    rotation vector
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
#ifdef USE_MATRIX_EXPONENTIAL
  // another possibility is to use the matrix exponential
  if(rvec == Eigen::Matrix<Scalar,3,1>::Zero()) rten = Eigen::Matrix<Scalar,3,3>::Identity();
  else {
    Eigen::Matrix<Scalar,3,3> skew;
    skew <<        0, -rvec[2],  rvec[1],
             rvec[2],        0, -rvec[0],
            -rvec[1],  rvec[0],        0;
    // method #1: use the matrix exponential function in unsupported Eigen MatrixFunctions module
    //rten = skew.exp();

    // method #2: use eigen decomposition (note: that skew is not self-adjoint)
    Eigen::EigenSolver<Eigen::Matrix<Scalar,3,3> > dec(skew);
    rten = (dec.eigenvectors() * dec.eigenvalues().array().exp().matrix().asDiagonal() * dec.eigenvectors().adjoint()).real();
  }
#else
   Eigen::Matrix<Scalar,4,1> q;

   vec_to_quat<Scalar>( rvec, q );

   quat_to_mat<Scalar>( q, rten );
#endif
}

template<typename Scalar>
void tangential_transf(const Eigen::Matrix<Scalar,3,1> &Psi, Eigen::Matrix<Scalar,3,3> &T)
{
  Scalar psi2 = Psi.squaredNorm();
  Eigen::Matrix<Scalar,3,3> Psiskew;
  Psiskew <<      0, -Psi[2],  Psi[1],
             Psi[2],       0, -Psi[0],
            -Psi[1],  Psi[0],       0;

  Scalar c1, c2, c3;
  if(psi2 < 5e-6) {
    Scalar psi4 = psi2*psi2;
    c1 = 1    - psi2/6   + psi4/120;  // + O(psi^6)
    c2 = 1/2. - psi2/24  + psi4/720;  // + O(psi^6)
    c3 = 1/6. - psi2/120 + psi4/5040; // + O(psi^6)
  }
  else {
    using std::sqrt;
    using std::sin;
    using std::cos;
    Scalar psi = sqrt(psi2);
    c1 = sin(psi)/psi;
    c2 = (1-cos(psi))/(psi2);
    c3 = (psi-sin(psi))/(psi*psi2);
  }
  T = c1*Eigen::Matrix<Scalar,3,3>::Identity() - c2*Psiskew + c3*(Psi*Psi.transpose());
}

template<typename Scalar>
void directional_deriv1(const Eigen::Matrix<Scalar,3,1> &Psi, const Eigen::Matrix<Scalar,3,1> &V,
                        Eigen::Matrix<Scalar,3,3> &C1)
{
  Scalar psi2 = Psi.squaredNorm();
  Eigen::Matrix<Scalar,3,3> Psiskew;
  Psiskew <<      0, -Psi[2],  Psi[1],
             Psi[2],       0, -Psi[0],
            -Psi[1],  Psi[0],       0;
  Eigen::Matrix<Scalar,3,3> Vskew;
  Vskew <<    0, -V[2],  V[1],
           V[2],     0, -V[0],
          -V[1],  V[0],     0;

  Scalar c1, c2, c3, c4, c5;
  if(psi2 < 5e-6) {
    Scalar psi4 = psi2*psi2;
    c1 = -1/3.  + psi2/30   - psi4/840;   // + O(psi^6)
    c2 = -1/12. + psi2/180  - psi4/6720;  // + O(psi^6)
    c3 = -1/60. + psi2/1260 - psi4/60480; // + O(psi^6)
    c4 = -1/2.  + psi2/24   - psi4/720;   // + O(psi^6)
    c5 = 1/6.   - psi2/120  + psi4/5040;  // + O(psi^6)
  }
  else {
    using std::sqrt;
    using std::cos;
    using std::sin;
    Scalar psi = sqrt(psi2);
    Scalar psi3 = psi*psi2;
    Scalar psi4 = psi*psi3;
    Scalar psi5 = psi*psi4;
    c1 = (psi*cos(psi)-sin(psi))/psi3;
    c2 = (psi*sin(psi)+2*cos(psi)-2)/psi4;
    c3 = (3*sin(psi)-2*psi-psi*cos(psi))/psi5;
    c4 = (cos(psi)-1)/psi2;
    c5 = (psi-sin(psi))/psi3;
  }
  C1 = c1*V*Psi.transpose() - c2*(Psiskew*V)*Psi.transpose()
       + c3*Psi.dot(V)*(Psi*Psi.transpose()) - c4*Vskew
       + c5*(Psi.dot(V)*Eigen::Matrix<Scalar,3,3>::Identity() + Psi*V.transpose());
}

template<typename Scalar>
void directional_deriv2(const Eigen::Matrix<Scalar,3,1> &Psi, const Eigen::Matrix<Scalar,3,1> &V,
                        Eigen::Matrix<Scalar,3,3> &C2)
{
  Scalar psi2 = Psi.squaredNorm();
  Eigen::Matrix<Scalar,3,3> Psiskew;
  Psiskew <<      0, -Psi[2],  Psi[1],
             Psi[2],       0, -Psi[0],
            -Psi[1],  Psi[0],       0;
  Eigen::Matrix<Scalar,3,3> Vskew;
  Vskew <<    0, -V[2],  V[1],
           V[2],     0, -V[0],
          -V[1],  V[0],     0;

  Scalar c1, c2, c3, c4, c5;
  if(psi2 < 5e-6) {
    Scalar psi4 = psi2*psi2;
    c1 = -1/3.  + psi2/30   - psi4/840;   // + O(psi^6)
    c2 = -1/12. + psi2/180  - psi4/6720;  // + O(psi^6)
    c3 = -1/60. + psi2/1260 - psi4/60480; // + O(psi^6)
    c4 = -1/2.  + psi2/24   - psi4/720;   // + O(psi^6)
    c5 = 1/6.   - psi2/120  + psi4/5040;  // + O(psi^6)
  }
  else {
    using std::sqrt;
    using std::cos;
    using std::sin;
    Scalar psi = sqrt(psi2);
    Scalar psi3 = psi*psi2;
    Scalar psi4 = psi*psi3;
    Scalar psi5 = psi*psi4;
    c1 = (psi*cos(psi)-sin(psi))/psi3;
    c2 = (psi*sin(psi)+2*cos(psi)-2)/psi4;
    c3 = (3*sin(psi)-2*psi-psi*cos(psi))/psi5;
    c4 = (cos(psi)-1)/psi2;
    c5 = (psi-sin(psi))/psi3;
  }
  C2 = c1*V*Psi.transpose() + c2*(Psiskew*V)*Psi.transpose()
       + c3*Psi.dot(V)*(Psi*Psi.transpose()) + c4*Vskew
       + c5*(Psi.dot(V)*Eigen::Matrix<Scalar,3,3>::Identity() + Psi*V.transpose());
}

template<typename Scalar>
void tangential_transf_dot(const Eigen::Matrix<Scalar,3,1> &Psi, const Eigen::Matrix<Scalar,3,1> &Psidot,
                           Eigen::Matrix<Scalar,3,3> &Tdot)
{
  Scalar psi2 = Psi.squaredNorm();
  Eigen::Matrix<Scalar,3,3> Psiskew;
  Psiskew <<      0, -Psi[2],  Psi[1],
             Psi[2],       0, -Psi[0],
            -Psi[1],  Psi[0],       0;

  Eigen::Matrix<Scalar,3,3> Psidotskew;
  Psidotskew <<         0, -Psidot[2],  Psidot[1],
                Psidot[2],          0, -Psidot[0],
               -Psidot[1],  Psidot[0],          0;

  Scalar c1, c2, c3, c4, c5;
  if(psi2 < 5e-6) {
    Scalar psi4 = psi2*psi2;
    c1 = -1/3.  + psi2/30   - psi4/840;   // + O(psi^6)
    c2 = -1/12. + psi2/180  - psi4/6720;  // + O(psi^6)
    c3 = -1/60. + psi2/1260 - psi4/60480; // + O(psi^6)
    c4 = -1/2.  + psi2/24   - psi4/720;   // + O(psi^6)
    c5 = 1/6.   - psi2/120  + psi4/5040;  // + O(psi^6)
  }
  else {
    using std::sqrt;
    using std::sin;
    using std::cos;
    Scalar psi  = sqrt(psi2);
    Scalar psi3 = psi*psi2;
    Scalar psi4 = psi*psi3;
    Scalar psi5 = psi*psi4;
    c1 = (psi*cos(psi)-sin(psi))/psi3;
    c2 = (psi*sin(psi)+2*cos(psi)-2)/psi4;
    c3 = (3*sin(psi)-2*psi-psi*cos(psi))/psi5;
    c4 = (cos(psi)-1)/psi2;
    c5 = (psi-sin(psi))/psi3;
  }
  Tdot = c1*Psi.dot(Psidot)*Eigen::Matrix<Scalar,3,3>::Identity() 
         - c2*Psi.dot(Psidot)*Psiskew + c3*Psi.dot(Psidot)*(Psi*Psi.transpose())
         + c4*Psidotskew + c5*(Psidot*Psi.transpose() + Psi*Psidot.transpose());
}

#endif
#endif
