#ifndef _CONTROL_INTERFACE_H_
#define _CONTROL_INTERFACE_H_

class ControlInterface {

  protected:
    double dt;	// time step size
    
  public:
  
    void setDt(double h) { dt = h; }

    // Control force initialization routine
    virtual void init(double *displacement, double *velocity, double *acc) = 0;

    // Control force routine
    virtual void ctrl(double *dis, double *vel, double *acc, double *f) = 0;

    // User defined force routine
    virtual void usd_forc(double time, double *userDefineForc) = 0;

    // User defined displacement routine
    virtual void usd_disp(double time, double *userDefineDisp,
                          double *userDefineVel) = 0;
};

#endif
