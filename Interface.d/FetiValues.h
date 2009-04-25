#ifndef CU_FETI_VALUES_H
#define CU_FETI_VALUES_H

class FetiValues
{
public:
  int number_iterations;      // number of iterations
  double primal_residual;     // primal residual
  double dual_residual;
  int number_reorthog_used;   // # reorthogonal vectors used
  double memory_use_this_sub; // total memory usage on this subdomain
  double time_in_feti;        // total time spent in FETI
  int feti_status;            // 0=success, 1=stagnation, 
                              // 2=iteration count exceeded, 3=other error
  int nrbms;

  FetiValues() { 
        number_iterations = 0;
        primal_residual = 0.0;
        dual_residual = 0.0;
        number_reorthog_used = 0;
        memory_use_this_sub = 0;
        time_in_feti = 0.0;
        feti_status = 0;
        nrbms = 0;
  };

}; // end of FetiValues definition

#endif
