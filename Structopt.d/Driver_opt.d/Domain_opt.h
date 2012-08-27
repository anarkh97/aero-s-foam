#ifndef _DOMAIN_OPT_H_
#define _DOMAIN_OPT_H_

#include <map>
#include <memory>
#include <boost/scoped_array.hpp>
#include <cstdio>
#include <Driver.d/Domain.h>

#ifdef STRUCTOPT
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Driver_opt.d/GeoSource_opt.h>

class Optpro;
class Structopt;
class Structopt_sd;
class OptActInfo;

class Domain_opt : virtual public Domain
{
protected:
  // related to StressIntegral       
  double valvmint;    // stressintegral
  double valvol;      // volume

  // semi-analytical sensitivity analysis
  int semiSAflag;
  int structoptFlag;
  int reliabilityFlag;
  int optgradFlag;          // Flag set when arrays for analytical SA exist 

  int numAnalysis, numStaticAna;
  OptActInfo* actInf;

  std::auto_ptr<CoordSet> nodescopy;
  std::auto_ptr<CoordSet> gradnodes;     // all the nodes for sensitivity analysis

  boost::scoped_array<boost::scoped_array<BCond> > gradnbcLC;        // set of Neuman bc for each load case
  BCond*  gradnbc;                                     // set of gradients of Neuman bc for active load case
  boost::scoped_array<BCond*> nbcLC;            // set of Neuman bc for each load case

  boost::scoped_array<boost::scoped_array<ComplexBCond> > gradcnbcLC;
  ComplexBCond * gradcnbc;



  FILE * optunitout;        // opt - file
  FILE * optprotout;        // nlp - file

  FILE * relunitout;        // rel - file
  FILE * relprotout;        // rpt - file

  int numLC;                   // number of global load cases
  int actLC;                   // global id number of active load case
  boost::scoped_array<int>                 numAnaNmLC;     // number of load cases for each analysis
  boost::scoped_array<boost::scoped_array<int> >  mapAnaNmLC;     // maps analysisId related load case id to global loadcase id
  boost::scoped_array<int>                 numNmLC;        // number of Neuman bc for each load case

  ResizeArray<std::auto_ptr<SolverInfo> > analysisList;
  int activeAnalysis;


  boost::scoped_array<double> kelData;
  boost::scoped_array<double> gkelData;
  std::auto_ptr<Vector> elForce;
  std::auto_ptr<Vector> elKelD;
  std::auto_ptr<Vector> elAdj;
  std::auto_ptr<Vector> elGrad; 
  std::auto_ptr<Vector> gradelstress;
  std::auto_ptr<Vector> elRHS;
  std::auto_ptr<ComplexVector> c_elDisp;
  std::auto_ptr<ComplexVector> c_elKelD;
  std::auto_ptr<ComplexVector> c_elRHS;
  std::auto_ptr<ComplexVector> c_elGrad; 
  std::auto_ptr<ComplexVector> c_elAdj;

  int checkList(int,int*,int);

public:
  double densProj_minVal;
  double densProj_maxVal;
  boost::scoped_array<double> densProj;
  boost::scoped_array<double> nodalDens;

public:
  Domain_opt();
  virtual ~Domain_opt();
  
  template<class Scalar> void buildDRHSforceDs(GenVector<Scalar>& /*, GeomState* gs=0,int = 0*/);
  template<class Scalar> void buildDStiffnesDSmDU(GenVector<Scalar>&, GenVector<Scalar>&, const Scalar*, int flag=0);
  template<class Scalar> Scalar DStiffnesDSmDUAdj(GenVector<Scalar>&, GenVector<Scalar>&, const Scalar*);
  
  FILE * getOptUnit() { return optunitout; }
  FILE * getRelUnit() { return relunitout; }

  // Moved to the global scope; declared in Structopt_base.h
  //std::auto_ptr<Optpro>         optpro;
  //std::auto_ptr<Structopt_sd> structopt;  
  //std::auto_ptr<Optpro>         relpro;
  //std::auto_ptr<Structopt_sd> structrel;

  char * optinputfile;
  char * relinputfile;

  void copyNodes()         { if (nodescopy.get() == 0) { nodescopy.reset(nodes.copy()); } }
  CoordSet& getNodesCopy() { copyNodes(); return *nodescopy; }

  int getOptvarTyp();
  void setGradActiveLC(int lc);

  void setActiveLC(int,int*bc=0,double*bcx=0); 

  void setStructoptFlag (int flag) { structoptFlag = flag; }
  int  getStructoptFlag ()         { return structoptFlag; }
  
  void setReliabilityFlag (int flag) { reliabilityFlag = flag; }
  int  getReliabilityFlag ()         { return reliabilityFlag; }

  void postProcessing(Vector&, double*, Vector&, double);
  void postProcessing(ComplexVector&, DComplex*, ComplexVector&, double);

  int  getSemiSAflag()        { return semiSAflag; }
  void setSemiSAflag(int flag){ semiSAflag = flag; }
  int geteleNummap(int iele) { return geoSource->glToPackElem(iele); }

  void closeOutputFiles() { fclose(optunitout); fclose(optprotout); }

  int  getNumLC();
  int  nNeumann(int ana=0, int lc=0);
  BCond* getNBC(int ana=0, int lc=0); 

  int getNumAnalysis()  { return numAnalysis;  }
  int getNumStaticAna() { return numStaticAna; }
  void activateAnalysis (int,int=1);

  void reloptInitialize(Structopt_sd*);

  void buildOptGrad();
  OptActInfo* buildOptInf();  
  void setoptInf(OptActInfo* info) { actInf = info; }

  // function that returns gradients of coordinates and boundary conditions  
  CoordSet& getGradNodes()        { return *gradnodes; }
  BCond*    getGradNBC(int lc=0)  { return gradnbcLC[lc].get(); }

  double getStrainEnergy(Vector &sol, double *bcx, SparseMatrix *gStiff=0,  
			 FullSquareMatrix *kelArray=0, 
			 int* eleList=0, int listSize=0);
  double getGradStrainEnergy(Vector&,Vector&,double*,FullSquareMatrix* =0, 
			     int* eleList=0, int listSize=0);
  double getGradPartStrainEnergy(Vector&, double*, FullSquareMatrix* =0, 
				 int* eleList=0, int listSize=0 );
  void getGradduStrainEnergy(Vector&,Vector&,double*,FullSquareMatrix* = 0,
			     int* eleList=0, int listSize=0);

  double getStructureMass(int *, int);
  double getGradStructureMass(int *, int);

  double getMomOfInertia(int*, int, int);
  double getGradMomOfInertia(int*, int, int);

  double getDCmass(Vector &, double *,double &, int *, int);
  double getGradDCmass(Vector&, Vector&, double*, double&, int *, int);
  double getGradPartDCmass(Vector&, double*, double&, int *, int);
  void getGradduDCmass(Vector&, Vector&, double*, double&, int *, int);


  double getNodalStressStrain(Vector &sol ,double *bcx, int node, int stressIndex, int surface, int eleFlag=-1);  
  double getGradNodalStressStrain(Vector&, Vector&, double*, int, int, int, int eleFlag=-1);
  double getGradPartNodalStressStrain(Vector&, double*, int, int, int, int eleFlag=-1);
  void getGradduNodalStressStrain(Vector&,Vector&,double*,int,int,int, int eleFlag=-1);

  template<class Scalar> Scalar getNodalDisp(GenVector<Scalar>& sol, const Scalar* bcx,
					     int node, int dispTyp, bool dblfy);
  template<class Scalar> Scalar getGradNodalDisp(GenVector<Scalar>& sol, GenVector<Scalar>& grad,
						 const Scalar* bcx, int node, int dispTyp, bool dblfy);
  template<class Scalar> void getGradduNodalDisp(GenVector<Scalar>& sol, GenVector<Scalar>& adj, 
						 const Scalar* bcx, int node, int dispType, bool dblfy);
  template<class Scalar> Scalar getAverageDisp(GenVector<Scalar>& sol, 
					       const Scalar* bcx,
					       const int* nodes, const int* dtypes, int size);
  template<class Scalar> Scalar getGradAverageDisp(GenVector<Scalar>& sol, GenVector<Scalar>& grad,
						   const Scalar* bcx,
						   const int* nodes, const int* dtypes, int size);
  template<class Scalar> void getGradduAverageDisp(GenVector<Scalar>& sol, GenVector<Scalar>& adj,
						   const Scalar* bcx,
						   const int* nodes, const int* dtypes, int size);

  template<class Scalar> double getDisplevel(GenVector<Scalar>& sol, const Scalar* bcx,
					     int aFlag, int quoFlag, int size, 
					     const int* nodeid, const int* distyp, 
					     const double* refVal, const double* difVal, double powFac);
  template<class Scalar> double getGradDisplevel(GenVector<Scalar>& sol, GenVector<Scalar>& grad, 
						 const Scalar* bcx,
						 int aFlag, int quoFlag, int size, 
						 const int* nodeid, const int* distyp, 
						 const double* refVal, const double* difVal, double powFac);
  template<class Scalar> void getGradduDisplevel(GenVector<Scalar>& sol, GenVector<Scalar>& adj,
						 const Scalar* bcx,
						 int aFlag, int quoFlag, int size, 
						 const int* nodeid, const int* distyp, 
						 const double* refVal, const double* difVal, double powFac);
  
  //double getNodalDisp(Vector &sol ,double *bcx, int node, int dispTyp);
  //DComplex getNodalDisp(ComplexVector &sol ,DComplex *bcx, int node, int dispTyp);
  //double getGradNodalDisp(Vector&, Vector&, double*, int, int);
  //void getGradduNodalDisp(ComplexVector&, ComplexVector &, DComplex *, int , int);
  void getGradduNodalDisp(Vector&, Vector &, double*, int , int);
  void getGradduNodalDisp(Vector &sol, double *bcx, int node, int dispTyp, 
			  double* nodFac, int* nodLoc, int& numNod);
  void getGradduNodalDisp(ComplexVector &sol, DComplex *bcx, int node, int dispTyp, 
			  DComplex* nodFac, int* nodLoc, int& numNod);
 
  double getNodalInternalForce(Vector &sol,double *bcx,int node,int dispTyp);
  double getGradNodalInternalForce(Vector&, Vector&, double*, int, int);
  double getGradPartNodalInternalForce(Vector&, double*, int, int);
  void getGradduNodalInternalForce(Vector&, Vector &, double *, int , int);
  
  double getVariationNodalPos(int, int);
  double getGradVariationNodalPos(int, int);

  double getStressIntegral(Vector&,double*,double&,double&,int*,int,int,int);
  double getGradStressIntegral(Vector&,Vector&,double*,double&,double&,int*,int,int,int);
  void getGradduStressIntegral(Vector&,Vector&,double*,double&,double&,int*,int,int,int);
  double getGradPartStressIntegral(Vector&,Vector&,double*,double&,double&,int*,int,int,int);
  
  void zeroGrad();  
  void getMidPointMass(double&,double*);

  void reliabilitySolve();
  void reliabilityInput();
  
  void structoptSolve();
  void structoptInput();

  void buildDPressureForceDs(Vector&, GeomState* gs=0);

  void addAnalysis(int sact);
  void updateAnalysis();

  void setUpData();
  int setNeuman(int _numNmLC, BCond* _nbcLC);

  template<class Scalar> static void scalAdjVec(int len, Scalar* v, Scalar s);
  template<class Scalar> static Scalar unify(Scalar s, int disptype);

  // nodal density projection related
  void buildDensProj(const DensityProjData& densProjData);
  void scaleDensProj();
  void buildNodalDensity(const NodalDensityData& nodalDensData);
  virtual double densProjCoeff(int dof);
  virtual void densProjectStiffness(GenFullSquareMatrix<double>& kel, int num);
  virtual void densProjectStiffnessC(GenFullSquareMatrix<DComplex>& kel, int num);
  template<class Scalar> void densProjectVector(const DensityProjData& densProjData, Scalar* vec);
  template<class Scalar> void densProjectVectorInv(const DensityProjData& densProjData, Scalar* vec);


protected:
  template<class Scalar> GenVector<Scalar>* getElDisp();
  template<class Scalar> GenVector<Scalar>* getElKelD();
  template<class Scalar> GenVector<Scalar>* getElGrad();
  template<class Scalar> GenVector<Scalar>* getElAdj();
  template<class Scalar> GenVector<Scalar>* getElRHS();
};

template<> inline
double Domain_opt::unify(double, int) { return 1.0; }
template<> inline
DComplex Domain_opt::unify(DComplex s, int disptype)
{ 
  DComplex r;
  switch(disptype)
    {
    case 0:case 1:case 2:case 3:case 4:case 5:
      r = ScalarTypes::unify(s);
      break;
    case 6:
      r = ScalarTypes::unify(s)*ScalarTypes::id<DComplex>();
      break;
    default:
      assert(0);
      break;
    }
  return r;
}

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Driver_opt.d/Domain_opt.C>
#endif

#endif // #ifdef STRUCTOPT
#endif

