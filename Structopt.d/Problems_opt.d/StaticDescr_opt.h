#ifndef _STATIC_DESCR_OPT_H_
#define _STATIC_DESCR_OPT_H_

#include <Utils.d/MyComplex.h>
#include <Problems.d/StaticDescr.h>
#include <Driver.d/Domain.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/FullMatrix.h>

#include <Structopt.d/Driver_opt.d/Domain_opt.h>

#ifdef STRUCTOPT

template<class T, class VectorType, class SolverType>
class SingleDomainPostProcessor_opt :
  public SingleDomainPostProcessor< T, VectorType, SolverType >
{
public:
  SingleDomainPostProcessor_opt< T, VectorType, SolverType >(Domain *d, T *_bcx) : 
    SingleDomainPostProcessor< T, VectorType, SolverType >(d, _bcx) {}  

  SingleDomainPostProcessor_opt<T, VectorType, SolverType>(Domain *d, T *_bcx, 
							   StaticTimers *_times) :
    SingleDomainPostProcessor<T, VectorType, SolverType>(d, _bcx, _times) {}
  
  SingleDomainPostProcessor_opt<T, VectorType, SolverType>(Domain *d, T *_bcx, 
							   StaticTimers *_times,
							   SolverType *_solver) :
    SingleDomainPostProcessor<T, VectorType, SolverType>(d, _bcx, _times, _solver) {}
  
  void staticOutput(VectorType& solution, VectorType& rhs, double time);
};



template <class T, class VectorType, class SolverType>
class SingleDomainStatic_opt : public SingleDomainStatic< T, VectorType, SolverType >
{
private:
  std::auto_ptr<VectorType> tmpAdjV;

  void iniDensProj();
public:
  SingleDomainStatic_opt<T,VectorType,SolverType>(Domain_opt *d) :
    SingleDomainStatic<T,VectorType,SolverType>(d) {}
  
  void preoptProcess();
  void reBuild();
  //virtual void preProcess();
  //virtual void getRHS(VectorType &);
  void getPseudoLoad (VectorType &, VectorType &);
  double getAdjPseudoLoad(VectorType &, VectorType &);  
  FullSquareMatrix *getkelArray() { return this->kelArray;} 

  SingleDomainPostProcessor_opt<T,VectorType,SolverType> *getPostProcessor();
  int getNumLC();
  void setActiveLC(int);
  void densProjectVector(VectorType&);
};

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Problems_opt.d/StaticDescr_opt.C>
#endif

#endif
#endif
