#ifdef USE_EIGEN3
#include <Element.d/Shell.d/ShellElementSemiTemplate.cpp>

template
void
ShellElementSemiTemplate<double>
::sands8(double *_xl,double *_yl, double *_zl,
               double e, double nu, double *_h,
               double *_v, double *_stress,
               bool strainFlg, int surface, double thrmStr);

template
void
ShellElementSemiTemplate<double>
::vms8WRTdisp(double *_xl,double *_yl, double *_zl,
              double e, double nu, double *_h,
              double *_v, double *_vmsWRTdisp,
              bool strainFlg, int surface, double thrmStr = 0);

#endif
