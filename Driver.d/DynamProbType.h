#ifndef _DYNAMPROBTYPE_H_
#define _DYNAMPROBTYPE_H_

#include <Driver.d/EigenProbType.h>
#include <Problems.d/EigenDescr.h>
#include <Math.d/VectorSet.h>

class Structopt;
class DistrInfo;
class SingleInfo;
template <typename T> class SysState;

/*
template <class VecType> 
class SysState {
   VecType &d_n, &v_n, &a_n, &v_n_m;
 public:

   SysState(VecType &v): d_n(v), v_n(v), a_n(v), v_n_m(v) { }
   
   SysState(VecType &d, VecType &v, VecType &a, VecType &vm): 
     d_n(d), v_n(v), a_n(a), v_n_m(vm) { }
   
   SysState &operator=(const SysState &sys);  
   
   VecType &getDisp()      { return d_n; }
   VecType &getVeloc()     { return v_n; }
   VecType &getPrevVeloc() { return v_n_m; }
   VecType &getAccel()     { return a_n; }

   VecType &getDispConst()      const { return d_n; }
   VecType &getVelocConst()     const { return v_n; }
   VecType &getPrevVelocConst() const { return v_n_m; }
   VecType &getAccelConst()     const { return a_n; }

};
*/

template <class VecType,
          class ProblemDescriptor> 
class NewmarkWorkVec {

   // Type of Newmark (imlicit = 0 / explicit = 1 / quasistatic = -1)
   int typ;      

   VecType * d_n_p;
   VecType * o_n_p;
   VecType * v_n_p;
   VecType * a_n_p;
   VecType *   rhs;
   VecType * ext_f;
   VecType * d_n_h;
   VecType * v_n_h;
   VecType * Md_n_h;
   VecType * Cd_n_h;

   VecType * tmp1;

   union { 
      VecType * dnc;
      VecType * tmp2;
   };
   union {
      VecType * vnc;
      VecType * fint;
   };
   VecType * anc;

   public:
   
   NewmarkWorkVec(int _typ, ProblemDescriptor *probDesc);
   ~NewmarkWorkVec();

   NewmarkWorkVec &operator=(const NewmarkWorkVec &v);
      
   VecType & get_d_n_p()  { return *d_n_p; }
   VecType & get_o_n_p()  { return *o_n_p; }
   VecType & get_v_n_p()  { return *v_n_p; }
   VecType & get_a_n_p()  { return *a_n_p; }
   VecType & get_rhs()    { return *rhs;   }
   VecType & get_ext_f()  { return *ext_f; }
   VecType & get_d_n_h()  { return *d_n_h; }
   VecType & get_v_n_h()  { return *v_n_h; }
   VecType & get_Md_n_h() { return *Md_n_h;}
   VecType & get_Cd_n_h() { return *Cd_n_h;}
   VecType & get_tmp1()   { return *tmp1;  }
   VecType & get_dnc()    { return *dnc;   }
   VecType & get_vnc()    { return *vnc;   }
   VecType & get_anc()    { return *anc;   }
   VecType & get_tmp2()   { return *tmp2;  }
   VecType & get_fint()   { return *fint;  }

   VecType & get_d_n_pConst() const { return *d_n_p; }
   VecType & get_o_n_pConst() const { return *o_n_p; }
   VecType & get_v_n_pConst() const { return *v_n_p; }
   VecType & get_a_n_pConst() const { return *a_n_p; }
   VecType & get_rhsConst()   const { return *rhs;   }
   VecType & get_d_n_hConst() const { return *d_n_h; }
   VecType & get_v_n_hConst() const { return *v_n_h; }
   VecType & get_Md_n_hConst()const { return *Md_n_h;}
   VecType & get_Cd_n_hConst()const { return *Cd_n_h;}
   VecType & get_tmp1Const()  const { return *tmp1;  }
   VecType & get_dncConst()   const { return *dnc;   }
   VecType & get_vncConst()   const { return *vnc;   }
   VecType & get_ancConst()   const { return *anc;   }
   VecType & get_tmp2Const()  const { return *tmp2;  }
   VecType & get_fintConst()  const { return *fint;  }


};



template < 
     class DynOps, 
     class VecType, 
     class PostProcessor, 
     class ProblemDescriptor,
     class Scalar>
class DynamicSolver {

     ProblemDescriptor *probDesc;
     PostProcessor     *postProcessor;
     
     double beta, gamma, alphaf, alpham;
     double dt, tmax;
     double minVel, maxVel, delta;
     int    algType;
     
     // These are not used
     //double Wext, Wdmp, pWela, pWkin;

     int steadyFlag,steadyMin,steadyMax;
     double steadyTol;

     VecType * d_n;
     VecType * v_n;
     VecType * a_n;
     VecType * v_p;
     VecType * constForce;
     VecType * aeroForce;
     VecType * adjVec;
     
     SysState<VecType> * curState; 
     
     DynOps * dynOps; 

     NewmarkWorkVec<VecType,ProblemDescriptor> * workVec;

   public:

     DynamicSolver(ProblemDescriptor *PrbD) {
        probDesc = PrbD; 
        aeroAlg = -1;//copied from aeroFlag in Utils.d/SolverInfo.h - JFD
        dynOps = 0;
        postProcessor = 0;
     }
     //~DynamicSolver() { if(dynOps) delete dynOps; if(postProcessor) delete postProcessor; }
     ~DynamicSolver() { if(postProcessor) delete postProcessor; }

     void solve();

     void explicitNewmarkLoop(SysState<VecType>&,VecType&,
                              DynOps& dynOps, 
			      NewmarkWorkVec<VecType,ProblemDescriptor> &workVec,
                              double, double);
     void implicitNewmarkLoop(SysState<VecType>&,VecType&,
                              DynOps& dynOps, 
			      NewmarkWorkVec<VecType,ProblemDescriptor> &workVec,
                              double, double);
     void     quasistaticLoop(SysState<VecType>&, VecType&, DynOps& dynOps, 
			      NewmarkWorkVec<VecType,ProblemDescriptor> &workVec,
                              double, double, int =0);

     int checkSteadyState(double time, double step, double criteria=-1.0);

     VecType * getpDis()  { return d_n; }
     VecType * getpVel()  { return v_n; }
     VecType * getpAcc()  { return a_n; }

     DynOps  * getpOps()  { return dynOps; }

     int aeroAlg; 

     VecType * getaeroForce() { return aeroForce; }

};

#ifdef _TEMPLATE_FIX_
#include <Driver.d/DynamProbType.C>
#endif

#endif
