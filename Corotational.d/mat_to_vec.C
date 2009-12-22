#include <math.h>
#include <Corotational.d/utilities.h>
#include <cmath>
#include <iostream>

const double epsilon = 1.0e-15;

void mat_to_vec(double rten[3][3] ,double rvec[3] )
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
      double th, q[4], coef;

      mat_to_quat( rten, q );

      double sthh = sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);

      double cthh = q[0];

      if ( sthh < 0.7 ) 
        th = 2.0*asin( sthh );
      else              
        th = 2.0*acos( cthh );

      if ( sthh < epsilon )
        coef = 2.0;
      else {
        if ( sthh > 1.0 ) sthh = 1.0;
        coef = th/sthh;
      }

      rvec[0] = coef*q[1];
      rvec[1] = coef*q[2];
      rvec[2] = coef*q[3];

}
