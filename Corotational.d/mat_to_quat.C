#include <math.h>
#include <Corotational.d/utilities.h>
#include <Math.d/mathUtility.h>
#include <iostream>
#include <Utils.d/linkfc.h>

extern "C"      {
   void _FORTRAN(dsyev)(const char &, const char &, const int &,
                        double *, const int &, double *, double *,
                        const int &, int &);
}


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

/*
      // PJSA 2-12-09 this is a variant due to Markley (Journal of Guidance and Control Vol 31 No 2 pp 440-442, 2008)
      if ( trace > rten[imax][imax] ) {
         q[0] = ( 1.0 +trace );
         q[1] = ( rten[2][1] -rten[1][2] );
         q[2] = ( rten[0][2] -rten[2][0] );
         q[3] = ( rten[1][0] -rten[0][1] );
      }
      else {
         i = p[imax];
         j = p[imax+1];
         k = p[imax+2];
         q[i+1] = 1.0 + rten[i][i] - rten[j][j] - rten[k][k];
         q[  0] = ( rten[k][j] -rten[j][k] );
         q[j+1] = ( rten[j][i] +rten[i][j] );
         q[k+1] = ( rten[k][i] +rten[i][k] );
      }
      double norm = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
      for(int l=0; l<4; ++l) q[l] /= norm;
*/
/*
      // PJSA 3-12-09 method due to Bar-Itzhack (Journal of Guidance and Control Vol 23 No 6 pp 1085-1087, 2000)
      // which always gives the quaternion that corresponds to the orthogonal matrix closest to rten
      double A[16] = {
        (rten[0][0]-rten[1][1]-rten[2][2])/3.0, (rten[0][1]+rten[1][0])/3.0,            (rten[0][2]+rten[2][0])/3.0,            (rten[2][1]-rten[1][2])/3.0,
        0.0,                                    (rten[1][1]-rten[0][0]-rten[2][2])/3.0, (rten[1][2]+rten[2][1])/3.0,            (rten[0][2]-rten[2][0])/3.0,
        0.0,                                    0.0,                                    (rten[2][2]-rten[0][0]-rten[1][1])/3.0, (rten[1][0]-rten[0][1])/3.0,
        0.0,                                    0.0,                                    0.0,                                    (rten[0][0]+rten[1][1]+rten[2][2])/3.0 };

      int info;
      double w[4], work[11];
      _FORTRAN(dsyev)('V', 'L', 4, A, 4, w, work, 11, info);

      q[0] = A[15];
      q[1] = A[12];
      q[2] = A[13];
      q[3] = A[14];
*/

      // always return the positive quaternion
      //if(q[0] < 0.0) for(int l=0; l<4; ++l) q[l] = -q[l];

}
