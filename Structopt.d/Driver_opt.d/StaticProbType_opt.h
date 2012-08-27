#ifndef _STATICPROBTYPE_OPT_H_
#define _STATICPROBTYPE_OPT_H_

#include <cassert>

#include <Driver.d/StaticProbType.h>
#include <Utils.d/MathUtils.h>

#ifdef STRUCTOPT
class Structopt;

template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
class StaticSolver_opt :
  public StaticSolver < Scalar, OpSolver, VecType, PostProcessor, ProblemDescriptor, ComplexVecType >
{
private:
  VecType** solLC;
  VecType** rhsLC;
  VecType** gradLC;  
  double adjNorm;

public:
  StaticSolver_opt() :
    StaticSolver< Scalar, OpSolver, VecType, PostProcessor, ProblemDescriptor, ComplexVecType >(), 
    solLC(0), rhsLC(0), gradLC(0), adjNorm(0) {}
  explicit StaticSolver_opt(ProblemDescriptor* pd) : 
    StaticSolver< Scalar, OpSolver, VecType, PostProcessor, ProblemDescriptor, ComplexVecType >(pd), 
    solLC(0), rhsLC(0), gradLC(0), adjNorm(0) {}
  ~StaticSolver_opt() {}

  void optPreSolve();
  void optSolve();
  void optSolveDuds();     
  void optAdjSolve();     
  void optPostProcessor(double time);
  
  double optGetAdjointPseudo(int lc=0);
  
  VecType* getpsol   (int lc=0);
  VecType* getpgrad  (int lc=0);
  VecType* getpadj   (); 
};

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Driver_opt.d/StaticProbType_opt.C>
#endif

#endif // #ifdef STRUCTOPT
#endif
