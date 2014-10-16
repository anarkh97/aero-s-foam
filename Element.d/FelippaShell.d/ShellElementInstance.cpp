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
::andesgf(int elm, double *x, double *y, double *z, double *_gravityForce,
          double *gamma, int gravflg, double rhoh);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesgfWRTcoord(int elm, double *x, double *y, double *z,
                  double *J, double *gamma, int gravflg, double rhoh,
                  int senMethod, double eps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesms(int elm, double *x, double *y, double *z, double *emass,
          double *gamma, double *grvfor, bool grvflg, double &totmas,
          bool masflg, double rhoh);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesmsWRTcoord(int elm, double *x, double *y, double *z, double *J,
                  double rhoh, int senMethod, double eps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesstf(int elm, double *estiff, double *fint, double nu,
           double *x, double *y, double *z, double *v, int ctyp, int flag,
           double *ndtemps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesstfWRTthick(int elm, double *destiffdthick, double nu,
                   double *x, double *y, double *z, int ctyp, int flag,
                   double *ndtemps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesstfWRTcoord(int elm, double *destiffdx[9], double E, double nu,
                   double rho, double eh, double Ta, double W, double *cFrame,
                   double *x, double *y, double *z, int ctyp, int flag,
                   int senMethod, double eps, double *ndtemps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesvms(int elm, int maxstr, double nu, double *x, double *y, double *z,
           double *v, double *stress, int ctyp, int strainflg, int surface,
           double *ndtemps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesvmsWRTdisp(int elm, double nu, double *x, double *y, double *z,
                  double *v, double *vmsWRTdisp, int ctyp, int surface,
                  double *ndtemps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesvmsWRTthic(int elm, double nu, double *x, double *y, double *z,
                  double *v, double *vmsWRTthic, int ctyp, int surface,
                  double *ndtemps);

template
void
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::andesvmsWRTcoord(int elm, double E, double nu, double rho, double eh,
                   double Ta, double W, double *cFrame, double *x,
                   double *y, double *z, double *v, double *vmsWRTcoord,
                   int ctyp, int surface, int senMethod, double eps,
                   double *ndtemps);

template
double
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::equivstr(double sxx, double syy, double szz, double sxy);

template
Eigen::Matrix<double,1,18>
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::equivstrSensitivityWRTdisp(double vms, double sxx, double syy, double szz, double sxy,
                             Eigen::Matrix<double,3,18> &dsigmadu);

template
double
ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>
::equivstrSensitivityWRTthic(double vms, double sxx, double syy, double szz, double sxy,
                             Eigen::Matrix<double,3,1> &dsigmadh);

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
