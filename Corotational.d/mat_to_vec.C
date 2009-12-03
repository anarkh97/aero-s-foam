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

      //std::cerr << "here in mat_to_vec, sthh = " << sthh << ", cthh = " << cthh << ", th = " << th << std::endl;
      if ( sthh < epsilon )
        { coef = 2.0;  /*std::cerr << "here in mat_to_vec #1, sthh = " << sthh << std::endl;*/ }
      else {
        if ( sthh > 1.0 ) { std::cerr << "here in mat_to_vec #2, sthh = " << sthh << std::endl; sthh = 1.0; }
        coef = th/sthh;
      }

      rvec[0] = coef*q[1];
      rvec[1] = coef*q[2];
      rvec[2] = coef*q[3];

/*
// PJSA alternative method from www.euclideanspace.com
// heading = rotation about y axis, attitude = rotation about z axis, bank = rotation about x axis
// euler angle order: heading applied first, attitude applied second, bank applied last (231)
  double sqw = q[0]*q[0];
  double sqx = q[1]*q[1];
  double sqy = q[2]*q[2];
  double sqz = q[3]*q[3];
  double test = q[1]*q[2] + q[3]*q[0];
  double heading, attitude, bank;
  if (test > 0.499) { // singularity at north pole
    std::cerr << "singularity at north pole\n";
    heading = 2.0 * atan2(q[1],q[0]);
    attitude = M_PI/2;
    bank = 0;
  }
  else if (test < -0.499) { // singularity at south pole
    std::cerr << "singularity at south pole\n";
    heading = -2.0 * atan2(q[1],q[0]);
    attitude = -M_PI/2;
    bank = 0;
  }
  else {
    heading = atan2(2.0*q[2]*q[0]-2.0*q[1]*q[3], sqx-sqy-sqz+sqw);
    attitude = asin(2.0*test);
    bank = atan2(2.0*q[1]*q[0]-2.0*q[2]*q[3],-sqx+sqy-sqz+sqw);
  }
*/
  // another method from wikipedia using x-y-z convention
  //double phi = atan2(2.0*(q[0]*q[1]+q[2]*q[3]),1.0-2.0*(q[1]*q[1]+q[2]*q[2]));
  //double theta = asin(2.0*(q[0]*q[2]-q[3]*q[1]));
  //double psi = atan2(2.0*(q[0]*q[3]+q[1]*q[2]),1.0-2.0*(q[2]*q[2]+q[3]*q[3]));
  // another method from wikipedia using 3-1-3 convention
/*
  double qnorm = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  double q1 = q[0]/qnorm, q2 = q[1]/qnorm, q3 = q[2]/qnorm, q4 = q[3]/qnorm;
  double phi = atan2((q1*q3+q2*q4),(q2*q3-q1*q4));
  double theta = acos(-q1*q1-q2*q2+q3*q3+q4*q4);
  double psi = -atan2((q1*q3-q2*q4),(q2*q3+q1*q4));
*/
/*
        // it's not 123,231,312,132,213,321
        double q1 = q[1], q2 = q[2], q3 = q[3], q4 = q[0];
        double qnorm = sqrt(q1*q1+q2*q2+q3*q3+q4*q4);
        int EULER_order_out = 321, Euler_type;
        double psi, theta, phi;
        if (EULER_order_out==121) {
            psi=atan2((q1*q2+q3*q4),(q2*q4-q1*q3));
            theta=acos(q4*q4+q1*q1-q2*q2-q3*q3);
            phi=atan2((q1*q2-q3*q4),(q1*q3+q2*q4));
          Euler_type=2;
        } else if (EULER_order_out==232) {
            psi=atan2((q1*q4+q2*q3),(q3*q4-q1*q2));
            theta=acos(q4*q4-q1*q1+q2*q2-q3*q3);
            phi=atan2((q2*q3-q1*q4),(q1*q2+q3*q4));
          Euler_type=2;
        } else if (EULER_order_out==313) {
            psi=atan2((q1*q3+q2*q4),(q1*q4-q2*q3));
            theta=acos(q4*q4-q1*q1-q2*q2+q3*q3);
            phi=atan2((q1*q3-q2*q4),(q1*q4+q2*q3));
          Euler_type=2;
        } else if (EULER_order_out==131) {
            psi=atan2((q1*q3-q2*q4),(q1*q2+q3*q4));
            theta=acos(q4*q4+q1*q1-q2*q2-q3*q3);
            phi=atan2((q1*q3+q2*q4),(q3*q4-q1*q2));
          Euler_type=2;
        } else if (EULER_order_out==212) {
            psi=atan2((q1*q2-q3*q4),(q1*q4+q2*q3));
            theta=acos(q4*q4-q1*q1+q2*q2-q3*q3);
            phi=atan2((q1*q2+q3*q4),(q1*q4-q2*q3));
          Euler_type=2;
        } else if (EULER_order_out==323) {
            psi=atan2((q2*q3-q1*q4),(q1*q3+q2*q4));
            theta=acos(q4*q4-q1*q1-q2*q2+q3*q3);
            phi=atan2((q1*q4+q2*q3),(q2*q4-q1*q3));
          Euler_type=2;
        } else if (EULER_order_out==123) {
            psi=atan2(2.0*(q1*q4-q2*q3),(q4*q4-q1*q1-q2*q2+q3*q3));
            theta=asin(2.0*(q1*q3+q2*q4));
            phi=atan2(2.0*(q3*q4-q1*q2),(q4*q4+q1*q1-q2*q2-q3*q3));
            std::cerr << "psi     = " << psi     << ", theta   = " << theta   << ", phi      = " << phi      << std::endl;
          Euler_type=1;
        } else if (EULER_order_out==231) {
            psi=atan2(2.0*(q2*q4-q1*q3),(q4*q4+q1*q1-q2*q2-q3*q3));
            theta=asin(2.0*(q1*q2+q3*q4));
            phi=atan2(2.0*(q1*q4-q3*q2),(q4*q4-q1*q1+q2*q2-q3*q3));
            std::cerr << "phi     = " << phi     << ", psi     = " << psi     << ", theta    = " << theta    << std::endl;
        } else if (EULER_order_out==312) {
            psi=atan2(2.0*(q3*q4-q1*q2),(q4*q4-q1*q1+q2*q2-q3*q3));
            theta=asin(2.0*(q1*q4+q2*q3));
            phi=atan2(2.0*(q2*q4-q3*q1),(q4*q4-q1*q1-q2*q2+q3*q3));
            std::cerr << "theta   = " << theta   << ", phi     = " << phi     << ", psi      = " << psi      << std::endl;
          Euler_type=1;
        } else if (EULER_order_out==132) {
            psi=atan2(2.0*(q1*q4+q2*q3),(q4*q4-q1*q1+q2*q2-q3*q3));
            theta=asin(2.0*(q3*q4-q1*q2));
            phi=atan2(2.0*(q1*q3+q2*q4),(q4*q4+q1*q1-q2*q2-q3*q3));
            std::cerr << "psi     = " << psi     << ", phi     = " << phi     << ", theta    = " << theta    << std::endl;
          Euler_type=1;
        } else if (EULER_order_out==213) {
            psi=atan2(2.0*(q1*q3+q2*q4),(q4*q4-q1*q1-q2*q2+q3*q3));
            theta=asin(2.0*(q1*q4-q2*q3));
            phi=atan2(2.0*(q1*q2+q3*q4),(q4*q4-q1*q1+q2*q2-q3*q3));
            std::cerr << "theta   = " << theta   << ", psi     = " << psi     << ", phi      = " << phi      << std::endl;
          Euler_type=1;
        } else if (EULER_order_out==321) {
            psi=atan2(2.0*(q1*q2+q3*q4),(q4*q4+q1*q1-q2*q2-q3*q3));
            theta=asin(2.0*(q2*q4-q1*q3));
            phi=atan2(2.0*(q1*q4+q3*q2),(q4*q4-q1*q1-q2*q2+q3*q3));
            std::cerr << "phi     = " << phi     << ", theta   = " << theta   << ", psi      = " << psi      << std::endl;
          Euler_type=1;
        } else {
            std::cerr << "Error: Invalid output Euler angle order type (conversion string).\n";
        }
*/


  //std::cerr << "rvec[0] = " << rvec[0] << ", rvec[1] = " << rvec[1] << ", rvec[2]  = " << rvec[2]  << std::endl;
  //std::cerr << "bank    = " << bank    << ", heading = " << heading << ", attitude = " << attitude << std::endl;

}
