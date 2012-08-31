#include <cstdio>
#include <cmath>

#include "ControlInterface.h"

class MyControl : public ControlInterface {

  public:
  
    // initialization routine
    void init(double *displacement, double *velocity, double *acceleration);

    // control force
    void ctrl(double *displacement, double *velocity, double *acceleration, 
              double *force);
	      
    void usd_disp(double time, double *userDefineDisplacement, 
                  double *userDefineVelocity);

    void usd_forc(double time, double *userDefineForce);
    
};


// Define new control Object (i.e. control force function, user defined
// force function, user defined displacement function)

ControlInterface *controlObj = new MyControl();

// disp  = displacement
// vel   = velocity
// accel = acceleration

void
MyControl::ctrl(double *disp, double *vel, double *accel, double *force)
{

}


// init = initialization routine

void
MyControl::init(double *disp, double *vel, double *accel)
{

}

// usd_disp = user defined displacements
void
MyControl::usd_disp(double time, double *userDefineDisp, double *userDefineVel)
{
 fprintf(stderr,"Applying Time Dependent Prescribed Displacement, time = %f, 4*sin(time) = %f\n",time,4.0*sin(time));

 userDefineDisp[0] = 4.0*sin(time);
 userDefineDisp[1] = 4.0*sin(time);
 userDefineDisp[2] = 4.0*sin(time);
 userDefineVel[0] = 0.0;
 userDefineVel[1] = 0.0;
 userDefineVel[2] = 0.0;

/* 
 userDefineDisp[0] = 0.0;
 userDefineDisp[1] = 0.0;
 userDefineDisp[2] = 0.0;
 userDefineVel[0] = 4.0*sin(time);
 userDefineVel[1] = 4.0*sin(time);
 userDefineVel[2] = 4.0*sin(time);
*/
}

// usd_forc = user defined forces

void
MyControl::usd_forc(double time, double *usdForce)
{
}
