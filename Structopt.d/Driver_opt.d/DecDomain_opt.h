#ifndef __DECDOMAIN_OPT_HPP__
#define __DECDOMAIN_OPT_HPP__

#ifdef STRUCTOPT

#include <memory>
#include <boost/tuple/tuple.hpp>

#include <Driver.d/DecDomain.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#include <Structopt.d/Structopt_dec.h>
#include <Feti.d/DistrVector.h>

template<class Scalar>
class GenDecDomain_opt : virtual public GenDecDomain<Scalar>
{
protected:
  char * optinputfile;

  FILE * optunitout;        // opt - file
  FILE * optprotout;        // nlp - file

  int numAnalysis, numStaticAna;

public:
  GenDecDomain_opt(Domain_opt* d) : GenDecDomain<Scalar>(d) {}
  virtual ~GenDecDomain_opt() {}

  virtual void structoptSolve();
  //virtual void structoptInput();
  //virtual void openOutputFiles() { geoSource->openOutputFiles(); }

  virtual int getNumAnalysis()  { return dynamic_cast<Domain_opt*>(domain)->getNumAnalysis(); }
  virtual int getNumStaticAna() { return dynamic_cast<Domain_opt*>(domain)->getNumStaticAna(); }
  virtual void activateAnalysis (int,int=1);

  virtual int  probType() { return domain->probType(); }  
  virtual SolverInfo  &solInfo() { return domain->solInfo(); } 
  virtual MatrixTimers &getTimers() { return domain->getTimers(); }

  virtual int getNumLC() { return dynamic_cast<Domain_opt*>(domain)->getNumLC(); }
  virtual void setActiveLC(int lc) { dynamic_cast<Domain_opt*>(domain)->setActiveLC(lc); }

  virtual void buildOptGrad();
  virtual OptActInfo* buildOptInf();
  virtual void setoptInf(OptActInfo* info);
  //void getMidPointMass(double& totMass, double* midPoint) 
  //  { dynamic_cast<Domain_opt*>(domain)->getMidPointMass(totMass, midPoint); }
  virtual void zeroGrad();

  virtual double getStructureMass(int* eleList, int listSize);
  virtual void computeSubdStructureMass(int, int* eleList, int listSize, double* results);

  virtual double getGradStructureMass(int* eleList, int listSize);
  virtual void computeGradStructureMass(int, int* eleList, int listSize, double* results);

  virtual double getStrainEnergy(GenDistrVector<Scalar>& sol, int* eleList=0, int listSize=0);
  virtual void computeSubdStrainEnergy(int, GenDistrVector<Scalar>* sol, int* eleList, int listSize, double* results);

  virtual double getGradStrainEnergy(GenDistrVector<Scalar>&, GenDistrVector<Scalar>&,
			     int* eleList=0, int listSize=0);
  virtual void computeGradStrainEnergy(int, GenDistrVector<Scalar>*, GenDistrVector<Scalar>*,
			       int* eleList, int listSize, double* results);

  virtual void getGradduStrainEnergy(GenDistrVector<Scalar>& sol, 
			     GenDistrVector<Scalar>& adj, 
			     int* eleList, int listSize);
  virtual void computeGradduStrainEnergy(int iSub,
				 GenDistrVector<Scalar>* sol, 
				 GenDistrVector<Scalar>* adj, 
				 int* eleList, int listSize);

  virtual double getGradPartStrainEnergy(GenDistrVector<Scalar>& sol, int* eleList, int listSize);
  virtual void computeGradPartStrainEnergy(int, GenDistrVector<Scalar>* sol, int* eleList, int listSize, double* results);

  virtual void buildDRHSforceDs(GenDistrVector<Scalar>& adj);
  virtual void computeDRHSforceDs(int, GenDistrVector<Scalar>* adj);

  virtual void buildDStiffnesDSmDU(GenDistrVector<Scalar>&, GenDistrVector<Scalar>&, int flag=0);
  virtual void computeDStiffnesDSmDU(int iSub, GenDistrVector<Scalar>* rhs, GenDistrVector<Scalar>* sol, int flag);
  virtual Scalar DStiffnesDSmDUAdj(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj);
  virtual void computeDStiffnesDSmDUAdj(int, GenDistrVector<Scalar>* sol, GenDistrVector<Scalar>* adj, Scalar* results);
  virtual void preProcess();

  virtual void getMidPointMass(double& totmas, double* midPoint);
  virtual void computeMidPointMass(int, double* totmasses, double* midpoints);

  virtual Scalar getNodalDisp(GenDistrVector<Scalar>& sol, int node, int dispTyp, bool dblfy);
  virtual Scalar getGradNodalDisp(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad,
				  int node, int dispTyp, bool dblfy);
  virtual void getGradduNodalDisp(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj, 
				  int node, int dispType, bool dblfy);

  virtual Scalar getAverageDisp(GenDistrVector<Scalar>& sol, const int* nodes, const int* dtypes, int size);
  virtual void computeAverageDisp(int, GenDistrVector<Scalar>* sol, const int* nodes, const int* dtypes, int size,
				  Scalar* results);
  virtual Scalar getGradAverageDisp(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad,
				    const int* nodes, const int* dtypes, int size);
  virtual void computeGradAverageDisp(int iSub, GenDistrVector<Scalar>* sol, GenDistrVector<Scalar>* grad, 
				      const int* nodes, const int* dtypes, int size, Scalar* results);
  virtual void getGradduAverageDisp(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj,
				    const int* nodes, const int* dtypes, int size);
  virtual void computeGradduAverageDisp(int iSub, GenDistrVector<Scalar>* sol, GenDistrVector<Scalar>* adj, 
					const int* nodes, const int* dtypes, int size);

  virtual double getDisplevel(GenDistrVector<Scalar>& sol, 
			      int aFlag, int quoFlag, int size, const int* nodeid, const int* distyp, 
			      const double* refVal, const double* difVal, double powFac);
  // too many arguments!!! therefore, we use tuple here
  virtual void computeDisplevel(int iSub, GenDistrVector<Scalar>* sol,
				boost::tuple<int, double,
				const int*, const int*, const double*, const Scalar*> args,
				double* results);
  virtual double getGradDisplevel(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad, 
				  int aFlag, int quoFlag, int size, const int* nodeid, const int* distyp, 
				  const double* refVal, const double* difVal, double powFac);
  virtual void computeGradDisplevel(int iSub, 						    
				    GenDistrVector<Scalar>* sol,  GenDistrVector<Scalar>* grad, 
				    boost::tuple<int, double, 
				    const int*, const int*, const double*, 
				    const Scalar*, const Scalar*> args,
				    double* results);
  virtual void getGradduDisplevel(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj,
				  int aFlag, int quoFlag, int size, 
				  const int* nodeid, const int* distyp, 
				  const double* refVal, const double* difVal, double powFac);

  virtual void globalSum(int size, int* s) {/*do nothing*/}
  virtual void globalSum(int size, double* s) {/*do nothing*/}
  virtual void globalSum(int size, DComplex* s) {/*do nothing*/}
  virtual double globalMax(double d) { return d; }
  virtual double globalSum(double d) { return d; }

  void scalAdjVector(GenDistrVector<Scalar>& v, Scalar s);
  void execScalAdjVector(int iSub, GenDistrVector<Scalar>* v, Scalar s);


  // nodal density projection related
  void buildNodalDensity(const NodalDensityData& nodalDensData);
  void computeNodalDensity(int iSub, const NodalDensityData* nodalDensData);
  void buildDensProj(const DensityProjData& densProjData);
  void computeDensProj(int iSub, const DensityProjData* densProjData);
  void scaleDensProj();
  void computeScaleDensProj(int iSub, double pMax, double pMin);
  void densProjectVector(const DensityProjData& densProjData,  GenDistrVector<Scalar>& vec);
  void computeDensProjectVector(int iSub, const DensityProjData* densProjData,  GenDistrVector<Scalar>* vec);
  void densProjectVectorInv(const DensityProjData& densProjData, GenDistrVector<Scalar>& vec);
  void computeDensProjectVectorInv(int iSub, const DensityProjData* densProjData,  GenDistrVector<Scalar>* vec);
  void reduceNodalDens(const NodalDensityData& ndd);
  void reduceNodalDensOp(const NodalDensityData& ndd, int mySubs, double** subs);
  void reduceNodalDensPowSumOp(const NodalDensityData& ndd, int mySubs, double** subs);
  void reduceNodalDensMaxOp(const NodalDensityData& ndd, int mySubs, double** subs);

  //
  Domain* getDomain() { return this->domain; }
};


#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Driver_opt.d/DecDomain_opt.C>
#endif

#endif
#endif
