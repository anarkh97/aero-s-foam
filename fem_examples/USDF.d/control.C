#include <cstdio>
#include <cmath>

#include "Control.d/ControlInterface.h"

class MyControl : public ControlInterface {

  public:
  
    // initialization routine
    void init(double *displacement, double *velocity, double *acceleration,
              SingleDomainDynamic *probDesc=0);

    // control force
    void ctrl(double *displacement, double *velocity, double *acceleration, 
              double *force, double time=0, SysState<Vector> *state=0, Vector *ext_f=0);
	      
    void usd_disp(double time, double *userDefineDisplacement, 
                  double *userDefineVelocity, double *userDefineAcc);

    void usd_forc(double time, double *userDefineForce);
    
};


// Define new control Object (i.e. control force function, user defined
// force function, user defined displacement function)

ControlInterface *controlObj = new MyControl();

// disp  = displacement
// vel   = velocity
// accel = acceleration

void
MyControl::ctrl(double *disp, double *vel, double *accel, double *force,
                double time, SysState<Vector> *state, Vector *ext_f)
{

}

// init = initialization routine

void
MyControl::init(double *disp, double *vel, double *accel, SingleDomainDynamic *probDesc)
{

}

// usd_disp = user defined displacements

void
MyControl::usd_disp(double time, double *userDefineDisp, double *userDefineVel,
                    double *userDefineAcc)
{
 // blank intentionally
}

// usd_forc = user defined forces

const double F0 = 1.0;
const double omega = 1.0;
void
MyControl::usd_forc(double time, double *usdForce)
{
 fprintf(stderr,"Applying Time Dependent Force\n");
 usdForce[0] = F0*cos(omega*time);
 usdForce[1] = F0*cos(omega*time);
 usdForce[2] = F0*cos(omega*time);
}
