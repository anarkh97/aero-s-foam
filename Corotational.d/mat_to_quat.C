#include <math.h>
#include <Corotational.d/utilities.h>
#include <Math.d/mathUtility.h>

void mat_to_quat( double rten[3][3], double q[4] )
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

      static int p[5] = {0,1,2,0,1};
      int    imax, i, j, k;
      double trace;

      trace = rten[0][0] +rten[1][1] +rten[2][2];

      imax = 0;
      if ( rten[1][1] > rten[imax][imax] ) imax = 1;
      if ( rten[2][2] > rten[imax][imax] ) imax = 2;

      if ( trace > rten[imax][imax] ) {
         q[0] = sqrt( 1.0 +trace )/2.0;
         q[1] = ( rten[2][1] -rten[1][2] )/(4.0*q[0]);
         q[2] = ( rten[0][2] -rten[2][0] )/(4.0*q[0]);
         q[3] = ( rten[1][0] -rten[0][1] )/(4.0*q[0]);
      }
      else {
         i = p[imax];
         j = p[imax+1];
         k = p[imax+2];
         q[i+1] = sqrt( rten[i][i]/2.0 +(1.0 -trace)/4.0 );
         q[  0] = ( rten[k][j] -rten[j][k] )/(4.0*q[i+1]);
         q[j+1] = ( rten[j][i] +rten[i][j] )/(4.0*q[i+1]);
         q[k+1] = ( rten[k][i] +rten[i][k] )/(4.0*q[i+1]);
      }
}
