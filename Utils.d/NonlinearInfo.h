#ifndef _NONLINEAR_INFO_
#define _NONLINEAR_INFO_

#include <limits>

// 
// Nonlinear information used in solution algorithms: static, dynamic problems
// with Newton-Raphson iteration or Arclength method
//

struct NonlinearInfo {

   int updateK;         // number of times to update K in NL analysis
   int kryflg;          // Krylov correction flag
   int initflg;         // initialization flag
   int reorthoflg;      // full reorthogonalization flag
   int maxiter;         // maximum # of Newton iterations
   int maxVec;          // maximum # of vectors to store
   int fitAlgShell;     // fit algorithm for shell elements
   int fitAlgBeam;      // fit algorithm for beam elements

   double tolRes;       // Newton iteration force residual tolerance
   double tolInc;       // Newton iteration displacement increment tolerance
   double dlambda;      // load step increment
   double maxLambda;    // maximum load step value to attain
   double lfactor;      // scaling factor to determine step size in the load
   int    extMin;       // lower iteration limit scaling in arclength 
   int    extMax;       // maximum iteration limit scaling in arclength

   bool unsymmetric;
   bool linesearch;
   bool failsafe;
   double failsafe_tol;

   NonlinearInfo() { updateK     = 1; kryflg     =   0; initflg =   0; 
                     reorthoflg  = 0; maxiter    = 100; maxVec  =   1; 
                     fitAlgShell = 2; fitAlgBeam =   2; dlambda = 1.0;
                     tolRes = 1.0E-6; tolInc     = std::numeric_limits<double>::infinity();
                     maxLambda = 1.0; lfactor    = 1.0; extMin  =   4;
                     extMax      = 6; unsymmetric = false; linesearch = false;
                     failsafe = false; failsafe_tol = std::numeric_limits<double>::epsilon(); }

};

#endif
