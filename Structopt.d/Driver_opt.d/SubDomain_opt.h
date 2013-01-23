#ifndef __SUBDOMAIN_OPT_HPP__
#define __SUBDOMAIN_OPT_HPP__

#ifdef STRUCTOPT

#include <Driver.d/SubDomain.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#include <Structopt.d/Driver_opt.d/SubDomain_opt.h>

template<class Scalar>
class GenSubDomain_opt : public GenSubDomain<Scalar>, public Domain_opt
{
private:
  std::map<int, int> glToLocalElement;
  
public:
  GenSubDomain_opt(Domain_opt& d, int sn, Connectivity &con, Connectivity &nds, int gn) :
    Domain(d, con.num(gn), con[gn], nds.num(gn), nds[gn]), // virtual base first
    GenSubDomain<Scalar>(d, sn, con, nds, gn),
    Domain_opt()
  {}
    
  GenSubDomain_opt(Domain_opt& d, int sn, int nNodes, int *nds,
               int nElems, int *elems, int gn) :
    Domain(d, nElems, elems, nNodes, nds), // virtual base first
    GenSubDomain<Scalar>(d, sn, nNodes, nds, nElems, elems, gn),
    Domain_opt()
  {}

  GenSubDomain_opt(Domain_opt& d, int sn, CoordSet* nodes, Elemset* elems, 
		   int* glNodeNums, int* glElemNums, int gn) :
    Domain(d, elems, nodes), // virtual base first
    GenSubDomain<Scalar>(d, sn, nodes, elems, glNodeNums, glElemNums, gn),
    Domain_opt()
  {}

  void setUpData() { GenSubDomain<Scalar>::setUpData(); Domain_opt::setUpData(); }
  double getGradStrainEnergy(const Scalar* sol, const Scalar* grad, int* eleList=0, int listSize=0);
  double getStrainEnergy(const Scalar* sol, int* eleList=0, int listSize=0);
  void getGradduStrainEnergy(const Scalar* sol, Scalar* adj, int* eleList=0, int listSize=0);
  double getGradPartStrainEnergy(const Scalar* sol, int* eleList=0, int listSize=0);
  double getStructureMass(int* eleList, int listSize);
  double getGradStructureMass(int* eleList, int listSize);
  void buildDRHSforceDs(Scalar* adj);
  void buildDStiffnesDSmDU(Scalar* rhs, const Scalar* adj, int flag);
  Scalar DStiffnesDSmDUAdj(const Scalar* sol, const Scalar* adj);

  void makeGlobalToLocalElementMap();
  int globalToLocalElement(int glNum) const;

  void getMidPointMass(double& totmas, double* midPoint);

  Scalar getNodalDisp(const Scalar* sol, int node, int dispTyp, int weight, bool doublify);
  Scalar getGradNodalDisp(const Scalar* sol, const Scalar* grad, int glNode, int dispTyp, int weight, bool doublify);
  void getGradduNodalDisp(const Scalar* sol, Scalar* adj, int glNode, int dispTyp, int weight, bool doublify);

  Scalar getAverageDisp(const Scalar* sol, const int* nodes, const int* dtypes, int size, Connectivity& nodeToSub);  
  Scalar getGradAverageDisp(const Scalar* sol, const Scalar* grad, const int* nodes, const int* dtypes, 
			    int size, Connectivity& nodeToSub);
  void getGradduAverageDisp(const Scalar* sol, Scalar* adj, const int* nodes, const int* dtypes,
			    int size, Connectivity& nodeToSub);

  double getDisplevel(const Scalar* sol, int size, double powFac,
		      const int* nodeid, const int* distyp, 
		      const double* refVal, const Scalar* difVal, 
		      Connectivity& nodeToSub);
  double getGradDisplevel(const Scalar* sol, const Scalar* grad,
			  int size, double powFac,
			  const int* nodeid, const int* distyp, 
			  const double* refVal, 
			  const Scalar* difVal, const Scalar* gradDifVal,
			  Connectivity& nodeToSub);
};

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Driver_opt.d/SubDomain_opt.C>
#endif


#endif

#endif
