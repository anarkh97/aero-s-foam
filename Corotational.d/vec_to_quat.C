#include <math.h>
#include <Corotational.d/utilities.h>
#include <iostream>

// const double epsilon = 1.0e-15;
const double epsilon = 1.0e-9;

typedef double Quat[4];
int EM_To_Q(double v[3], Quat q, int reparam);

void vec_to_quat( double rvec[3], double q[4] )
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
      double coef;

      // compute norm of rotation vector rvec
      double th = sqrt( rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2] );

      q[0] = cos( 0.5*th );

      if ( th < epsilon )
        coef = 0.5;
      else
        coef = sin(0.5*th)/th;

      q[1] = coef*rvec[0];
      q[2] = coef*rvec[1];
      q[3] = coef*rvec[2];
         
      // Compute norm of quaterion q
      double norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );

      // Normalize quaterion q
      q[0] /= norm;
      q[1] /= norm;
      q[2] /= norm;
      q[3] /= norm;
      //std::cerr << "1. q = " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << std::endl;

   if(q[0] < 0.0) for(int l=0; l<4; ++l) q[l] = -q[l]; // PJSA DEBUG
/*
  Quat quat;
  EM_To_Q(rvec, quat, 1);
  //std::cerr << "2. q = " << quat[0] << " " << quat[1] << " " << quat[2] << " " << quat[3] << std::endl;
  q[0] = quat[3];
  q[1] = quat[0];
  q[2] = quat[1];
  q[3] = quat[2];
*/  
}

