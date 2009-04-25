#include <math.h>
#include <Corotational.d/utilities.h>

// const double epsilon = 1.0e-15;
const double epsilon = 1.0e-9;

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

}

