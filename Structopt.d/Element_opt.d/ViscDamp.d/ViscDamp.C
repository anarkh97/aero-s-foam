#ifdef STRUCTOPT

#include <cmath>
#include <cassert>

#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp.h>

const double ViscDamp::a      = 1.0;
const double ViscDamp::b      = 1.0;

// primary wave speed
double ViscDamp::V_p()
{
  const double E   = prop->E;
  const double nu  = prop->nu;
  const double rho = prop->rho;

  assert( rho > 0 );
  assert( nu >= 0 && nu < .5 );
  assert( E > 0 );

  return sqrt(E*(1.0-nu)/((1-2*nu)*(1+nu)*rho));
}

double ViscDamp::dV_p()
{
  const double E     = prop->E;
  const double nu    = prop->nu;
  const double rho   = prop->rho;
  const double dE    = gradprop->E;
  const double dnu   = gradprop->nu;
  const double drho  = gradprop->rho;

  return 0.5*(
	      dE*(1.0-nu)/((1-2*nu)*(1+nu)*rho)
	      - dnu*E/((1-2*nu)*(1+nu)*rho)
	      + dnu*2*(1+nu)*rho*E*(1.0-nu)/((1-2*nu)*(1+nu)*rho*(1-2*nu)*(1+nu)*rho)
	      - dnu*(1-2*nu)*rho*E*(1.0-nu)/((1-2*nu)*(1+nu)*rho*(1-2*nu)*(1+nu)*rho)
	      - drho*(1-2*nu)*(1+nu)*E*(1.0-nu)/((1-2*nu)*(1+nu)*rho*(1-2*nu)*(1+nu)*rho)
	      )/V_p();
}

// secondary wave speed
double ViscDamp::V_s()
{
  const double E   = prop->E;
  const double nu  = prop->nu;
  const double rho = prop->rho;

  assert( rho > 0 );
  assert( nu >= 0 && nu < .5 );
  assert( E > 0 );

  return sqrt(E/(2*(1+nu)*rho));
}

double ViscDamp::dV_s()
{
  const double E     = prop->E;
  const double nu    = prop->nu;
  const double rho   = prop->rho;
  const double dE    = gradprop->E;
  const double dnu   = gradprop->nu;
  const double drho  = gradprop->rho;

  return 0.5*(
	      dE/(2*(1+nu)*rho)
	      - dnu*2*rho*E/((2*(1+nu)*rho)*(2*(1+nu)*rho))
	      - drho*2*(1+nu)*E/((2*(1+nu)*rho)*(2*(1+nu)*rho))
	      )/V_s();  
}


double ViscDamp::c11(CoordSet& cs)
{
  const double rho   = prop->rho;
  return a*rho*V_p()*area(cs);
}

double ViscDamp::dc11(CoordSet& cs, CoordSet& dcs)
{
  const double rho   = prop->rho;
  const double drho  = gradprop->rho;
  const double ar    = area(cs);
  const double dar   = darea(cs, dcs);

  return
    a*drho*V_p()*ar + 
    a*rho*dV_p()*ar + 
    a*rho*V_p()*dar;
}

double ViscDamp::c00(CoordSet& cs)
{
  const double rho   = prop->rho;
  return b*rho*V_s()*area(cs);
}

double ViscDamp::dc00(CoordSet& cs, CoordSet& dcs)
{
  const double rho   = prop->rho;
  const double drho  = gradprop->rho;
  const double ar    = area(cs);
  const double dar   = darea(cs, dcs);

  return
    b*drho*V_s()*ar + 
    b*rho*dV_s()*ar + 
    b*rho*V_s()*dar;
}
#endif
