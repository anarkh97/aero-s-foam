#include <math.h>
#include <Corotational.d/utilities.h>

double form_rottensor( double rvec[3], double rten[3][3] )

/*****************************************************************
 *  Compute the rotation tensor based on a rotation vector
 *
 *  Input:
 *  rvec:   rotation vector
 *
 *  Output:
 *  rten:   rotation tensor
 *
 *  Return:
 *
 *  th = sqrt(rvec[0]*rvec[0] +rvec[1]*rvec[1] + rvec[2]*rvec[2])
 *     = rotation angle
 *****************************************************************/
{
      double q[4];

      vec_to_quat( rvec, q );

      quat_to_mat( q, rten );

      return sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]); 
}
