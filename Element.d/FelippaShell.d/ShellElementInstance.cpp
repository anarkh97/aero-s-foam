#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andescrd(int elm, double *x, double *y, double *z, double *rot,
           double *xlp, double *ylp, double *zlp, double &area);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesms(int elm, double *x, double *y, double *z, double *emass,
          double *gamma, double *grvfor, bool grvflg, double &totmas,
          bool masflg);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesstf(int elm, double *estiff, double *fint, double nu,
           double *x, double *y, double *z, double *v, int ctyp, int flag,
           double *ndtemps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesvms(int elm, int maxstr, double nu, double *globalx,
           double *globaly, double *globalz, double *v, double *stress,
           int ctyp, int strainflg, int surface, double *ndtemps);

template
double
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::equivstr(double sxx, double syy, double szz, double sxy);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::transform(double *lframe, double *gframe, double *str);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesups(int elm, double *state, double *X, double *Y, double *Z, double *_v);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesare(int elm, double *X, double *Y, double *Z, double &area);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesden(int elm, double *X, double *Y, double *Z, double &D);
#endif
