#ifndef __DIST_DOM_OPT_HPP__
#define __DIST_DOM_OPT_HPP__

#include <Dist.d/DistDom.h>
#include <Structopt.d/Driver_opt.d/DecDomain_opt.h>

template<class Scalar>
class GenDistrDomain_opt : public GenDecDomain_opt<Scalar>, 
			   public GenDistrDomain<Scalar>
{
public:
  GenDistrDomain_opt(Domain_opt *domain) :
    GenDecDomain<Scalar>(domain),
    GenDecDomain_opt<Scalar>(domain),
    GenDistrDomain<Scalar>(domain) {}
  virtual ~GenDistrDomain_opt() {}

  void preProcess();

  void postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f,
		      GenDistrVector<Scalar> *aeroF = 0, int x = 0, MDDynamMat *dynOps = 0,
		      SysState<GenDistrVector<Scalar> > *distState = 0)
  { GenDistrDomain<Scalar>::postProcessing(u,f,aeroF,x,dynOps,distState); }

  void globalSum(int size, int* s) { this->communicator->globalSum(size, s); }
  void globalSum(int size, double* s) { this->communicator->globalSum(size, s); }
  void globalSum(int size, DComplex* s) { this->communicator->globalSum(size, s); }
  double globalMax(double d) { return this->communicator->globalMax(d); }
  double globalSum(double d) { return this->communicator->globalSum(d); }
};

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Dist_opt.d/DistDom_opt.C>
#endif


#endif
