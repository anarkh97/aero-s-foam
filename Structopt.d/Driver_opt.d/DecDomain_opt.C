#include <Structopt.d/Driver_opt.d/DecDomain_opt.h>
#include <Structopt.d/Driver_opt.d/SubDomain_opt.h>

#include <Threads.d/PHelper.h>
#include <numeric>
#include <cfloat>

extern int yyoptparse(void);

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::structoptSolve() 
{
  optpro.reset(new Optpro(0));                 //Define Optimization Problem
  structopt.reset(new Structopt_dec<Scalar>(0, this, optpro.get()));              //Define Structural Opt. Interface
  
  dynamic_cast<Domain_opt*>(this->domain)->structoptInput();                          // Reading Optimiztion Input File 
  assert(dynamic_cast<Structopt_dec<Scalar>*>(structopt.get()) != 0);
  dynamic_cast<Structopt_dec<Scalar>*>(structopt.get())
    ->build(this, optpro.get());           //Build Structural Optimization    
                                             
  optpro->solve(structopt.get());                  //Solve Optimization problem

  dynamic_cast<Domain_opt*>(this->domain)->closeOutputFiles();
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::activateAnalysis (int sact, int nbFlag) 
{ 
  dynamic_cast<Domain_opt*>(domain)->activateAnalysis(sact, nbFlag);
  // go through all subdomains
  for(int i=0; i<this->getNumSub(); ++i)
    {
      this->getSubDomain(i)->solInfo() = domain->solInfo();
    }
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
double GenDecDomain_opt<Scalar>::getStrainEnergy(GenDistrVector<Scalar>& sol, int* eleList, int listSize)
{
  int numSub = this->getNumSub();
  double* results = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeSubdStrainEnergy, &sol, eleList, listSize, results);
  double result = std::accumulate(results, results + numSub, 0.0);
  globalSum(1, &result);
  return result;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeSubdStrainEnergy(int iSub, GenDistrVector<Scalar>* sol, int* eleList, int listSize, double* results)
{
  results[iSub] = dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getStrainEnergy(sol->subData(iSub), eleList, listSize);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDecDomain_opt<Scalar>::getStructureMass(int* eleList, int listSize)
{
  int numSub = this->getNumSub();
  double* results = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeSubdStructureMass, eleList, listSize, results);
  double result = std::accumulate(results, results + numSub, 0.0);
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeSubdStructureMass(int iSub, int* eleList, int listSize, double* results)
{
  results[iSub] = dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getStructureMass(eleList, listSize);
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::buildOptGrad()
{ 
  // temporary fix???
  dynamic_cast<Domain_opt*>(domain)->buildOptGrad(); 
  // go through all subdomains
  for(int i=0; i<this->getNumSub(); ++i)
    {
      dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(i))->buildOptGrad();
    }
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::zeroGrad()
{
  // temporary fix???
  dynamic_cast<Domain_opt*>(domain)->zeroGrad();
  // go through all subdomains
  for(int i=0; i<this->getNumSub(); ++i)
    {
      dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(i))->zeroGrad();
    }
}

//------------------------------------------------------------------------------
template<class Scalar>
OptActInfo* GenDecDomain_opt<Scalar>::buildOptInf()
{ 
  
  return dynamic_cast<Domain_opt*>(domain)->buildOptInf();
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::setoptInf(OptActInfo* info) 
{ 
  dynamic_cast<Domain_opt*>(domain)->setoptInf(info); 
  // go through all subdomains
  for(int i=0; i<this->getNumSub(); ++i)
    {
      dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(i))->setoptInf(info); 
    }
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::getGradduStrainEnergy(GenDistrVector<Scalar>& sol, 
						     GenDistrVector<Scalar>& adj, 
						     int* eleList, int listSize)
{
  execParal(this->getNumSub(), this, 
	    &GenDecDomain_opt<Scalar>::computeGradduStrainEnergy, &sol, &adj, eleList, listSize);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeGradduStrainEnergy(int iSub,
							 GenDistrVector<Scalar>* sol, 
							 GenDistrVector<Scalar>* adj, 
							 int* eleList, int listSize)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getGradduStrainEnergy(sol->subData(iSub), adj->subData(iSub), eleList, listSize);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDecDomain_opt<Scalar>::getGradPartStrainEnergy(GenDistrVector<Scalar>& sol, 
						       int* eleList, int listSize)
{
  int numSub = this->getNumSub();
  double* results = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeGradPartStrainEnergy, &sol, eleList, listSize, results);
  double result = std::accumulate(results, results + numSub, 0.0);
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeGradPartStrainEnergy(int iSub,
							   GenDistrVector<Scalar>* sol, 
							   int* eleList, int listSize, double* results)
{
  results[iSub] = dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getGradPartStrainEnergy(sol->subData(iSub), eleList, listSize);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::buildDRHSforceDs(GenDistrVector<Scalar>& adj)
{
  execParal(this->getNumSub(), this, &GenDecDomain_opt<Scalar>::computeDRHSforceDs, &adj);
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeDRHSforceDs(int iSub, GenDistrVector<Scalar>* adj)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    buildDRHSforceDs(adj->subData(iSub));
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDecDomain_opt<Scalar>::DStiffnesDSmDUAdj(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj)
{
  int numSub = this->getNumSub();
  Scalar* results = static_cast<Scalar*>(dbg_alloca(numSub*sizeof(Scalar)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeDStiffnesDSmDUAdj, &sol, &adj, results);
  Scalar result = std::accumulate(results, results + numSub, static_cast<Scalar>(0.0));
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeDStiffnesDSmDUAdj(int iSub, 
							GenDistrVector<Scalar>* sol, 
							GenDistrVector<Scalar>* adj,
							Scalar* results)
{
  results[iSub] = 
    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    DStiffnesDSmDUAdj(sol->subData(iSub), adj->subData(iSub));
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDecDomain_opt<Scalar>::getGradStructureMass(int* eleList, int listSize)
{
  int numSub = this->getNumSub();
  double* results = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeGradStructureMass, eleList, listSize, results);
  double result = std::accumulate(results, results + numSub, 0.0);
  globalSum(1, &result);
  return result;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeGradStructureMass(int iSub, int* eleList, int listSize, double* results)
{
  results[iSub] = 
    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getGradStructureMass(eleList, listSize);
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::preProcess()
{
  GenDecDomain<Scalar>::preProcess();
  assert(this->subToNode != 0);  
  if(this->nodeToSub == 0)
    { this->nodeToSub = this->subToNode->reverse(); }
  assert(this->nodeToSub != 0);  

  for(int i=0; i<this->getNumSub(); ++i)
    {
      //assert(dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(i)) != 0);      
      dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(i))->makeGlobalToLocalElementMap(); 
    }  
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::getMidPointMass(double& totmas, double* midPoint)
{
  int numSub = this->getNumSub();
  double* totmasses = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  double* midpoints = static_cast<double*>(dbg_alloca(numSub*3*sizeof(double)));  
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeMidPointMass, totmasses, midpoints);
  totmas = 0.0; midPoint[0] = midPoint[1] = midPoint[2] = 0.0;
  for(int i=0; i<numSub; ++i)
    {
      totmas += totmasses[i];
      midPoint[0] += midpoints[3*i];
      midPoint[1] += midpoints[3*i+1];
      midPoint[2] += midpoints[3*i+2];
    }
  globalSum(1, &totmas);
  globalSum(3, midPoint);
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeMidPointMass(int iSub, double* totmasses, double* midpoints)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getMidPointMass(totmasses[iSub], midpoints+3*iSub);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::buildDStiffnesDSmDU(GenDistrVector<Scalar>& rhs, GenDistrVector<Scalar>& sol, int flag)
{
  execParal(this->getNumSub(), this, &GenDecDomain_opt<Scalar>::computeDStiffnesDSmDU, &rhs, &sol, flag);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeDStiffnesDSmDU(int iSub, GenDistrVector<Scalar>* rhs, GenDistrVector<Scalar>* sol, int flag)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    buildDStiffnesDSmDU(rhs->subData(iSub), sol->subData(iSub), flag);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDecDomain_opt<Scalar>::getGradStrainEnergy(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad,
						     int* eleList, int listSize)
{
  int numSub = this->getNumSub();
  double* results = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeGradStrainEnergy, 
	    &sol, &grad, eleList, listSize, results);
  double result = std::accumulate(results, results + numSub, 0.0);
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeGradStrainEnergy(int iSub, GenDistrVector<Scalar>* sol, GenDistrVector<Scalar>* grad,
						       int* eleList, int listSize, double* results)
{
  results[iSub] = 
    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getGradStrainEnergy(sol->subData(iSub), grad->subData(iSub), eleList, listSize);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDecDomain_opt<Scalar>::getNodalDisp(GenDistrVector<Scalar>& sol, int node, int dispTyp,
					      bool dblfy)
{
  Scalar result = 0;
  int numSub = this->nodeToSub->num(node);
  for(int iSub=0; iSub<numSub; ++iSub)
    {
      int glbSubNo = (*(this->nodeToSub))[node][iSub];
      int locSubNo = (this->glSubToLocal)? this->glSubToLocal[glbSubNo] : glbSubNo;
      if(locSubNo<0) { continue; }
      result += dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(locSubNo))->
	getNodalDisp(sol.subData(locSubNo), node, dispTyp, numSub, dblfy);
    }
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::getGradduNodalDisp(GenDistrVector<Scalar>& sol, 
						  GenDistrVector<Scalar>& adj, 
						  int node, int dispTyp, bool dblfy)
{
  int numSub = this->nodeToSub->num(node);
  for(int iSub=0; iSub<numSub; ++iSub)
    {
      int glbSubNo = (*(this->nodeToSub))[node][iSub];
      int locSubNo = (this->glSubToLocal)? this->glSubToLocal[glbSubNo] : glbSubNo;
      if(locSubNo<0) { continue; }
      dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(locSubNo))->
	getGradduNodalDisp(sol.subData(locSubNo), adj.subData(locSubNo), node, dispTyp, numSub, dblfy);
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDecDomain_opt<Scalar>::getGradNodalDisp(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad,
						  int node, int dispTyp, bool dblfy)
{
  Scalar result = 0;
  int numSub = this->nodeToSub->num(node);
  for(int iSub=0; iSub<numSub; ++iSub)
    {
      int glbSubNo = (*(this->nodeToSub))[node][iSub];
      int locSubNo = (this->glSubToLocal)? this->glSubToLocal[glbSubNo] : glbSubNo;
      if(locSubNo<0) { continue; }
      result += dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(locSubNo))->
	getGradNodalDisp(sol.subData(locSubNo), grad.subData(locSubNo), 
			 node, dispTyp, numSub, dblfy);
    }
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDecDomain_opt<Scalar>::getAverageDisp(GenDistrVector<Scalar>& sol,
						const int* nodes, const int* dtypes,
						int size)
{
  int numSub = this->getNumSub();
  Scalar* results = static_cast<Scalar*>(dbg_alloca(numSub*sizeof(Scalar)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeAverageDisp, &sol, nodes, dtypes, size, results);
  Scalar result = std::accumulate(results, results + numSub, static_cast<Scalar>(0.0));
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeAverageDisp(int iSub, 
						  GenDistrVector<Scalar>* sol, 
						  const int* nodes, const int* dtypes,
						  int size,
						  Scalar* results)
{
  results[iSub] = 
    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getAverageDisp(sol->subData(iSub), nodes, dtypes, size, *(this->nodeToSub));
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDecDomain_opt<Scalar>::getGradAverageDisp(GenDistrVector<Scalar>& sol,
						    GenDistrVector<Scalar>& grad,
						    const int* nodes, const int* dtypes,
						    int size)
{
  int numSub = this->getNumSub();
  Scalar* results = static_cast<Scalar*>(dbg_alloca(numSub*sizeof(Scalar)));
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeGradAverageDisp, &sol, &grad,
	    nodes, dtypes, size, results);
  Scalar result = std::accumulate(results, results + numSub, static_cast<Scalar>(0.0));
  globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeGradAverageDisp(int iSub, 
						      GenDistrVector<Scalar>* sol, 
						      GenDistrVector<Scalar>* grad, 
						      const int* nodes, const int* dtypes,
						      int size,
						      Scalar* results)
{
  results[iSub] = 
    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getGradAverageDisp(sol->subData(iSub), grad->subData(iSub),
		       nodes, dtypes, size, *(this->nodeToSub));
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::getGradduAverageDisp(GenDistrVector<Scalar>& sol,
						    GenDistrVector<Scalar>& adj,
						    const int* nodes, const int* dtypes,
						    int size)
{
  int numSub = this->getNumSub();
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeGradduAverageDisp, &sol, &adj,
	    nodes, dtypes, size);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeGradduAverageDisp(int iSub, 
							GenDistrVector<Scalar>* sol, 
							GenDistrVector<Scalar>* adj, 
							const int* nodes, const int* dtypes,
							int size)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getGradduAverageDisp(sol->subData(iSub), adj->subData(iSub),
			 nodes, dtypes, size, *(this->nodeToSub));
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDecDomain_opt<Scalar>::getDisplevel(GenDistrVector<Scalar>& sol, 
					      int aFlag, int quoFlag, int size, 
					      const int* nodeid, const int* distyp, 
					      const double* refVal, const double* difVal, double powFac)
{
  int numSub = this->getNumSub();
  double* results = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  boost::scoped_array<Scalar> difVal1(new Scalar[size]);
  Scalar aveDisp;

  switch(aFlag)
    {
    case 0:
      std::copy(difVal, difVal+size, difVal1.get());
      break;
    case 1:
      aveDisp = getAverageDisp(sol, nodeid, distyp, size);
      std::fill(difVal1.get(), difVal1.get() + size, aveDisp);
      break;
    case 2:
      for(int i=0; i<size; ++i) 
	{ difVal1[i] = getNodalDisp(sol, static_cast<int>(difVal[i])-1, distyp[i], false); }
      break;
    default:
      assert(0);      
    }
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeDisplevel, 
	    &sol,
	    boost::make_tuple(size, powFac, nodeid, distyp, refVal, 
			      const_cast<const Scalar*>(difVal1.get())),
	    results);
  double quoVal = quoFlag ? size : 1.0;
  double result = std::accumulate(results, results + numSub, 0.0)/quoVal;
  globalSum(1, &result);
  return pow(result, 1.0/powFac);  
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeDisplevel(int iSub, GenDistrVector<Scalar>* sol,
				boost::tuple<int, double,
				const int*, const int*, const double*, const Scalar*> args,
				double* results)
{
  results[iSub] = 
    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getDisplevel(sol->subData(iSub), 
		 boost::get<0>(args), // size
		 boost::get<1>(args), // powFac
		 boost::get<2>(args), // nodeid
		 boost::get<3>(args), // distyp
		 boost::get<4>(args), // refVal
		 boost::get<5>(args), // difVal
		 (*this->nodeToSub));
    return;
}


//------------------------------------------------------------------------------
template<class Scalar>
double GenDecDomain_opt<Scalar>::getGradDisplevel(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad,
						  int aFlag, int quoFlag, int size, 
						  const int* nodeid, const int* distyp, 
						  const double* refVal, const double* difVal, double powFac)
{
  int numSub = this->getNumSub();
  double* results = static_cast<double*>(dbg_alloca(numSub*sizeof(double)));
  Scalar aveDisp, gradAveDisp;
  boost::scoped_array<Scalar> difVal1(new Scalar[size]);
  boost::scoped_array<Scalar> gradDifVal1(new Scalar[size]);

  switch(aFlag)
    {
    case 0:
      std::copy(difVal, difVal+size, difVal1.get());
      std::fill(gradDifVal1.get(), gradDifVal1.get() + size, 0);
      break;
    case 1:
      aveDisp = getAverageDisp(sol, nodeid, distyp, size);
      std::fill(difVal1.get(), difVal1.get() + size, aveDisp);
      gradAveDisp = getGradAverageDisp(sol, grad, nodeid, distyp, size);
      std::fill(gradDifVal1.get(), gradDifVal1.get() + size, gradAveDisp);
      break;
    case 2:
      for(int i=0; i<size; ++i) 
	{ 
	  difVal1[i]     = getNodalDisp(sol, static_cast<int>(difVal[i])-1, distyp[i], false);
	  gradDifVal1[i] = getGradNodalDisp(sol, grad, static_cast<int>(difVal[i])-1, distyp[i], false);
	}
      break;
    default:
      assert(0);      
    }
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeGradDisplevel, 
	    &sol, &grad,
	    boost::make_tuple(size, powFac, nodeid, distyp, refVal, 
			      const_cast<const Scalar*>(difVal1.get()),
			      const_cast<const Scalar*>(gradDifVal1.get())),
	    results);

  double quoVal = quoFlag ? size : 1.0;
  double sumVal = std::accumulate(results, results + numSub, 0.0)/quoVal;
  globalSum(1, &sumVal);
  double dl = getDisplevel(sol, aFlag, quoFlag, size, nodeid, distyp, refVal, difVal, powFac);
  return dl == 0?
    0 :
    pow(dl, 1.0 - powFac)/powFac*sumVal;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeGradDisplevel(int iSub, 						    
						    GenDistrVector<Scalar>* sol,  GenDistrVector<Scalar>* grad, 
						    boost::tuple<int, double, 
						    const int*, const int*, const double*, 
						    const Scalar*, const Scalar*> args,
						    double* results)
{
  results[iSub] = 
    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    getGradDisplevel(sol->subData(iSub), grad->subData(iSub),
		     boost::get<0>(args), // size
		     boost::get<1>(args), // powFac
		     boost::get<2>(args), // nodeid
		     boost::get<3>(args), // distyp
		     boost::get<4>(args), // refVal
		     boost::get<5>(args), // difVal
		     boost::get<6>(args), // gradDifVal
		     (*this->nodeToSub));
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::getGradduDisplevel(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj,
						  int aFlag, int quoFlag, int size, 
						  const int* nodeid, const int* distyp, 
						  const double* refVal, const double* difVal, double powFac)
{
  int numSub = this->getNumSub();
  Scalar aveDisp;
  GenDistrVector<Scalar> difAdj1(sol.info());
  GenDistrVector<Scalar> difAdj2(sol.info());

  switch(aFlag)
    {
    case 0:
      for(int i=0; i<size; ++i) 
	{	  
	  Scalar disVal  = (getNodalDisp(sol, nodeid[i], distyp[i], false) - difVal[i])/refVal[i];
	  difAdj1.zero();
	  getGradduNodalDisp(sol, difAdj1, nodeid[i], distyp[i], false);
   	  scalAdjVector(difAdj1,
			powFac*pow(ScalarTypes::doublify(disVal), powFac-1)/refVal[i]*
			Domain_opt::unify(disVal, distyp[i]));
	  adj += difAdj1;
	}
      break;
    case 1:
      assert(typeid(Scalar)!=typeid(DComplex)); // not implemented for complex
      aveDisp = getAverageDisp(sol, nodeid, distyp, size);
      difAdj2.zero();
      getGradduAverageDisp(sol, difAdj2, nodeid, distyp, size);
      for(int i=0; i<size; ++i) 
	{
	  Scalar disVal  = (getNodalDisp(sol, nodeid[i], distyp[i], false) - aveDisp)/refVal[i];
	  difAdj1.zero();
	  getGradduNodalDisp(sol, difAdj1, nodeid[i], distyp[i], false);
	  difAdj1 -= difAdj2;
	  scalAdjVector(difAdj1, powFac*pow(ScalarTypes::doublify(disVal), powFac-1)/refVal[i]*
			ScalarTypes::unify(disVal));
	  adj += difAdj1;
	}
      break;
    case 2:
      for(int i=0; i<size; ++i) 
	{ 
	  Scalar disVal  = (getNodalDisp(sol, nodeid[i], distyp[i], false) - 
			    getNodalDisp(sol, static_cast<int>(difVal[i])-1, distyp[i], false))/refVal[i];
	  difAdj1.zero();
	  getGradduNodalDisp(sol, difAdj1, nodeid[i], distyp[i], false);
	  difAdj2.zero();
	  int refNode = static_cast<int>(difVal[i])-1;
	  assert(refNode >= 0);
	  getGradduNodalDisp(sol, difAdj2, refNode, distyp[i], false);
	  difAdj1 -= difAdj2;
	  scalAdjVector(difAdj1, 
			powFac*pow(ScalarTypes::doublify(disVal), powFac-1)/refVal[i]*
			Domain_opt::unify(disVal, distyp[i]));
	  adj += difAdj1;
	}
      break;
    default:
      assert(0);      
    }
  double quoVal = quoFlag ? size : 1.0;
  double dl = getDisplevel(sol, aFlag, quoFlag, size, nodeid, distyp, refVal, difVal, powFac);
  if(dl == 0)
    { adj.zero(); }
  else
    { adj *= pow(dl, 1.0 - powFac)/powFac/quoVal; }
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::scalAdjVector(GenDistrVector<Scalar>& v, Scalar s)
{
  int numSub = this->getNumSub();
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::execScalAdjVector, &v, s);
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::execScalAdjVector(int iSub, GenDistrVector<Scalar>* v, Scalar s)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->
    scalAdjVec(v->subLen(iSub), v->subData(iSub), s);
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::buildNodalDensity(const NodalDensityData& nodalDensData)
{
  int numSub = this->getNumSub();
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeNodalDensity, &nodalDensData);
  reduceNodalDens(nodalDensData);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeNodalDensity(int iSub, const NodalDensityData* nodalDensData)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->buildNodalDensity(*nodalDensData);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::buildDensProj(const DensityProjData& densProjData)
{
  int numSub = this->getNumSub();
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeDensProj, &densProjData);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeDensProj(int iSub, const DensityProjData* densProjData)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->buildDensProj(*densProjData);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::scaleDensProj()
{
  int numSub = this->getNumSub();
  double pMax = -DBL_MAX;
  double pMin =  DBL_MAX;
  for(int iSub=0; iSub<numSub; ++iSub)
    {
      pMax = std::max(pMax, dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->densProj_maxVal);
      pMin = std::min(pMin, dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->densProj_minVal);
    }
  pMax =  this->globalMax(pMax);
  pMin = -this->globalMax(-pMin);

  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeScaleDensProj, pMax, pMin);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeScaleDensProj(int iSub, double pMax, double pMin)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->densProj_maxVal = pMax;
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->densProj_minVal = pMin;  
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))->scaleDensProj();
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::densProjectVector(const DensityProjData& densProjData,  GenDistrVector<Scalar>& vec)
{
  int numSub = this->getNumSub();
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeDensProjectVector, &densProjData, &vec);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeDensProjectVector(int iSub, const DensityProjData* densProjData,  GenDistrVector<Scalar>* vec)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))
    ->densProjectVector(*densProjData, vec->subData(iSub));
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::densProjectVectorInv(const DensityProjData& densProjData, GenDistrVector<Scalar>& vec)
{
  int numSub = this->getNumSub();
  execParal(numSub, this, &GenDecDomain_opt<Scalar>::computeDensProjectVectorInv, densProjData, &vec);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::computeDensProjectVectorInv(int iSub, const DensityProjData* densProjData,  GenDistrVector<Scalar>* vec)
{
  dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(iSub))
    ->densProjectVectorInv(*densProjData, vec->subData(iSub));
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::reduceNodalDens(const NodalDensityData& ndd)
{
  double** subs = static_cast<double**>(dbg_alloca(this->globalNumSub*sizeof(double**)));
  for(int iNode=0; iNode < this->nodeToSub->csize(); ++iNode)
    {
      const int numSub = this->nodeToSub->num(iNode);
      int mySubs = 0;
      if(numSub <= 1) { continue; }
      for(int iSub=0; iSub<numSub; ++iSub)
	{
	  const int glbSubNo = (*(this->nodeToSub))[iNode][iSub];
	  const int locSubNo = this->glSubToLocal? this->glSubToLocal[glbSubNo] : glbSubNo;
	  if(locSubNo<0) { continue; }
	  GenSubDomain_opt<Scalar>* subDomain = 
	    dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(locSubNo));
	  const int locNode = subDomain->globalToLocal(iNode);
	  if(locNode < 0) { continue; }
	  subs[mySubs] = subDomain->nodalDens.get() + locNode;
	  ++mySubs;
	}
      reduceNodalDensOp(ndd, mySubs, subs);
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::reduceNodalDensOp(const NodalDensityData& ndd,
						 int mySubs, double** subs)
{
  switch(ndd.nodDensFunc)
    {
    case 1:
      // nodalDens = (sum_i rho_i^p )^(1/p)
      reduceNodalDensPowSumOp(ndd, mySubs, subs);
      break;
    case 2:
      // nodalDen = 1/N*(sum_i rho_i^p )^(1/p)
      // not implemented
      assert(0);
      break;
    case 3:
      // nodalDens = max(connected nodes)
      reduceNodalDensMaxOp(ndd, mySubs, subs);
      break;
    default:
      assert(0);
      break;
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::reduceNodalDensMaxOp(const NodalDensityData& ndd,
						    int mySubs, double** subs)
{
  double buffer = -DBL_MAX;
  for(int i=0; i<mySubs; ++i)
    { buffer = std::max(*(subs[i]), buffer); }
  buffer = this->globalMax(buffer);
  for(int i=0; i<mySubs; ++i)
    { *(subs[i]) = buffer; }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain_opt<Scalar>::reduceNodalDensPowSumOp(const NodalDensityData& ndd,
						      int mySubs, double** subs)
{
  double buffer = 0;
  for(int i=0; i<mySubs; ++i)
    { buffer += (*(subs[i])>0)? pow(*(subs[i]), ndd.nodDensParam1) : 0; }
  buffer = this->globalSum(buffer);
  buffer = (buffer>0)? pow(buffer, 1.0/ndd.nodDensParam1) : 0;
  for(int i=0; i<mySubs; ++i)
    { *(subs[i]) = buffer; }
  return;
}
