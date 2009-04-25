#ifndef __SUBDOMAIN_FACTORY_OPT_HPP__
#define __SUBDOMAIN_FACTORY_OPT_HPP__

#include <Driver.d/SubDomainFactory.h>
#include <Structopt.d/Driver_opt.d/SubDomain_opt.h>

template<class Scalar>
class GenSubDomainFactory_opt : public GenSubDomainFactory<Scalar>
{
 public:
  ~GenSubDomainFactory_opt() {}
  virtual
  GenSubDomain_opt<Scalar>* createSubDomain(Domain& d, int sn, int nNodes, int *nds,
					    int nElems, int *elems, int gn) const
  { return new GenSubDomain_opt<Scalar>(dynamic_cast<Domain_opt&>(d), sn, nNodes, nds, nElems, elems, gn); }

  virtual
  GenSubDomain_opt<Scalar>* createSubDomain(Domain & d, int sn, Connectivity &con, 
					    Connectivity &nds, int gn) const
  { return new GenSubDomain_opt<Scalar>(dynamic_cast<Domain_opt&>(d), sn, con, nds, gn); }

  virtual
  GenSubDomain_opt<Scalar>* createSubDomain(Domain& d, int sn, CoordSet* nodes, 
					    Elemset* elems, int *glNodeNums, 
					    int *glElemNums, int gn) const
  { return new GenSubDomain_opt<Scalar>(dynamic_cast<Domain_opt&>(d), sn, nodes, elems, 
					glNodeNums, glElemNums, gn); }
};

#endif
