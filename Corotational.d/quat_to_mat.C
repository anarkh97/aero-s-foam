#include <stdio.h>
#include <math.h>
#include <Corotational.d/utilities.h>

void quat_to_mat( double q[4], double rten[3][3] )
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
/*
   // Compute norm of q
   double norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );

   // Normalize q
   q[0] /= norm;
   q[1] /= norm;
   q[2] /= norm;
   q[3] /= norm;

   rten[0][0] = 2.0*(q[1]*q[1] + q[0]*q[0]) - 1.0;
   rten[1][1] = 2.0*(q[2]*q[2] + q[0]*q[0]) - 1.0;
   rten[2][2] = 2.0*(q[3]*q[3] + q[0]*q[0]) - 1.0;

   rten[0][1] = 2.0*(q[1]*q[2] - q[3]*q[0]);
   rten[0][2] = 2.0*(q[1]*q[3] + q[2]*q[0]);
   rten[1][2] = 2.0*(q[2]*q[3] - q[1]*q[0]);

   rten[1][0] = 2.0*(q[2]*q[1] + q[3]*q[0]);
   rten[2][0] = 2.0*(q[3]*q[1] - q[2]*q[0]);
   rten[2][1] = 2.0*(q[3]*q[2] + q[1]*q[0]);
*/
   rten[0][0] = q[1]*q[1] - q[2]*q[2] - q[3]*q[3] + q[0]*q[0];
   rten[1][0] = 2.0*(q[1]*q[2] + q[3]*q[0]);
   rten[2][0] = 2.0*(q[1]*q[3] - q[2]*q[0]);
   rten[0][1] = 2.0*(q[1]*q[2] - q[3]*q[0]);
   rten[1][1] = -q[1]*q[1] + q[2]*q[2] - q[3]*q[3] + q[0]*q[0];
   rten[2][1] = 2.0*(q[2]*q[3] + q[1]*q[0]);
   rten[0][2] = 2.0*(q[1]*q[3] + q[2]*q[0]);
   rten[1][2] = 2.0*(q[2]*q[3] - q[1]*q[0]);
   rten[2][2] = -q[1]*q[1] - q[2]*q[2] + q[3]*q[3] + q[0]*q[0];
}
