#ifndef _MD_STATIC_DESCR_OPT_HPP_
#define _MD_STATIC_DESCR_OPT_HPP_

#ifdef STRUCTOPT

#include <Structopt.d/Driver_opt.d/DecDomain_opt.h>
#include <Structopt.d/Dist_opt.d/DistDom_opt.h>
#include <Paral.d/MDStatic.h>


template <class Scalar> 
class GenMultiDomainPostProcessor_opt : public GenMultiDomainPostProcessor<Scalar>
{
public:
  GenMultiDomainPostProcessor_opt(GenDecDomain<Scalar>* d, 
				  GenFetiSolver<Scalar> *s, 
				  StaticTimers* times=0) :
    GenMultiDomainPostProcessor<Scalar>(d, s, times) {}

    void staticOutput(GenDistrVector<Scalar>& solution, GenDistrVector<Scalar>& rhs, double time=0);
};

template<class Scalar>
class GenMultiDomainStatic_opt : public GenMultiDomainStatic<Scalar>
{
private:
  std::auto_ptr<GenDistrVector<Scalar> > tmpAdjV;
public:
  //#ifdef DISTRIBUTED
  //  typedef GenDistrDomain_opt<Scalar> DomainType_opt;
  //#else
  typedef GenDecDomain_opt<Scalar>   DomainType_opt;
  //#endif

  GenMultiDomainStatic_opt(DomainType_opt* domain);
  ~GenMultiDomainStatic_opt();
  virtual GenMultiDomainPostProcessor_opt<Scalar> *getPostProcessor()
    { return new GenMultiDomainPostProcessor_opt<Scalar>(this->decDomain, this->solver, this->times); }

  void preoptProcess();
  int getNumLC();
  void setActiveLC(int);
  void reBuild(); 
  void getPseudoLoad(GenDistrVector<Scalar>&, GenDistrVector<Scalar>&);
  double getAdjPseudoLoad(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj);

  void iniDensProj();
  void densProjectVector(GenDistrVector<Scalar>&);
};

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Paral_opt.d/MDStatic_opt.C>
#endif

#endif

#endif
