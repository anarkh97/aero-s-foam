#include <stdio.h>
#include <math.h>

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

const double F0 = 1.0;
const double omega = 1.0;
void
MyControl::ctrl(double *disp, double *vel, double *accel, double *force)
{
 fprintf(stderr,"Apply Control Force \n");
 force[0] = F0*cos(omega*disp[0]);
 force[1] = F0*cos(omega*disp[1]);
 force[2] = F0*cos(omega*disp[2]);
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
 // blank intentionally
}

// usd_forc = user defined forces

void
MyControl::usd_forc(double time, double *usdForce)
{


}
