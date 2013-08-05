#ifndef _NONLINEAR_INFO_
#define _NONLINEAR_INFO_

#include <limits>

// 
// Nonlinear information used in solution algorithms: static, dynamic problems
// with Newton-Raphson iteration or Arclength method
//

struct NonlinearInfo {

   int updateK;         // number of Newton iterations between K updates in time/load step
   int stepUpdateK;     // number of time/load steps between K updates
   int kryflg;          // Krylov correction flag
   int initflg;         // initialization flag
   int reorthoflg;      // full reorthogonalization flag
   int maxiter;         // maximum # of Newton iterations
   int maxVec;          // maximum # of vectors to store
   int fitAlgShell;     // fit algorithm for shell elements
   int fitAlgBeam;      // fit algorithm for beam elements

   double tolRes;       // Newton iteration relative force residual tolerance
   double tolInc;       // Newton iteration relative displacement increment tolerance
   double absTolRes;    // Newton iteration absolute force residual tolerance
   double absTolInc;    // Newton iteration absolute displacement increment tolerance
   double dlambda;      // load step increment
   double maxLambda;    // maximum load step value to attain
   double lfactor;      // scaling factor to determine step size in the load
   int    extMin;       // lower iteration limit scaling in arclength 
   int    extMax;       // maximum iteration limit scaling in arclength

   bool unsymmetric;
   bool linesearch;
   bool failsafe;
   double failsafe_tol;

   bool linearelastic;

   NonlinearInfo() { setDefaults(); }

   void setDefaults() {
                     updateK     = 1; kryflg     =   0; initflg =   0;
                     reorthoflg  = 0; maxiter    = 100; maxVec  =   1;
                     fitAlgShell = 2; fitAlgBeam =   2; dlambda = 1.0;
                     tolRes = 1.0E-6; tolInc     = std::numeric_limits<double>::infinity();
                     absTolRes = 0.0; absTolInc  = std::numeric_limits<double>::infinity();
                     maxLambda = 1.0; lfactor    = 1.0; extMin  =   4;
                     extMax      = 6; unsymmetric = false; linesearch = false;
                     failsafe = false; failsafe_tol = std::numeric_limits<double>::epsilon();
                     stepUpdateK = 1; linearelastic = false; }

};

#endif
