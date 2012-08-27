#include <cassert>

#include <Structopt.d/Paral_opt.d/MDStatic_opt.h>

template<class Scalar>
GenMultiDomainStatic_opt<Scalar>::GenMultiDomainStatic_opt(DomainType_opt *d) :
  GenMultiDomainStatic<Scalar>()
{
  switch(d->solInfo().fetiInfo.version) {
  default:
  case FetiInfo::feti1:
    filePrint(stderr," ... FETI-1 is Selected             ...\n");
    break;
  case FetiInfo::feti2:
    filePrint(stderr," ... FETI-2 is Selected             ...\n");
    break;
  case FetiInfo::fetidp:
    if (!(d->solInfo().fetiInfo.dph_flag)) 
      filePrint(stderr," ... FETI-Dual/Primal is Selected   ...\n");
    else 
      filePrint(stderr," ... FETI-DPH is Selected           ...\n");
    break;
  }
  this->decDomain = dynamic_cast<GenDecDomain<Scalar>*>(d);
  this->times  = new StaticTimers;

  // debug: check static timers are initialized
  MatrixTimers &mt = d->getTimers();
}

template<class Scalar>
GenMultiDomainStatic_opt<Scalar>::~GenMultiDomainStatic_opt()
{  
  this->decDomain = 0;
  //Deleted in parent's destructor delete this->times;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenMultiDomainStatic_opt<Scalar>::preoptProcess()
{
  // decDomain->preprocess is called from structopt_dec!!! 
  // Makes renumbering, connectivities and dofsets
  // startTimerMemory(times->preProcess, times->memoryPreProcess);
  // decDomain->preProcess();
  // stopTimerMemory(times->preProcess, times->memoryPreProcess);

  iniDensProj();

  // Construct FETI solver
  this->times->getFetiSolverTime -= getTime();
  this->solver = this->decDomain->getFetiSolver();
  this->times->getFetiSolverTime += getTime();
}

//------------------------------------------------------------------------------
template<class Scalar>
int GenMultiDomainStatic_opt<Scalar>::getNumLC()
{
  return dynamic_cast<GenDecDomain_opt<Scalar>* >(this->decDomain)->getNumLC();
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenMultiDomainStatic_opt<Scalar>::setActiveLC(int lc)
{
  dynamic_cast<GenDecDomain_opt<Scalar>* >(this->decDomain)->setActiveLC(lc);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenMultiDomainStatic_opt<Scalar>::reBuild()
{
  iniDensProj();
  this->rebuildSolver();
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenMultiDomainStatic_opt<Scalar>::getAdjPseudoLoad(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj)
{
  Scalar part1 = 0;
  if(static_cast<GeoSource_opt*>(geoSource)->isRHSDsgnDep())
    {
      if(tmpAdjV.get()==0) { tmpAdjV.reset(new GenDistrVector<Scalar>(adj.info())); }
      tmpAdjV->zero();  
      // Derivative of the load vector for temperature, pressure and gravity
      dynamic_cast<GenDecDomain_opt<Scalar>* >(this->decDomain)->buildDRHSforceDs(*tmpAdjV);
      part1 = (*tmpAdjV)*adj;
    }
  //Derivative of the elastic stiffness
  Scalar part2 = dynamic_cast<GenDecDomain_opt<Scalar>* >(this->decDomain)->DStiffnesDSmDUAdj(sol,adj);
  return ScalarTypes::Real(part1 + part2);
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenMultiDomainPostProcessor_opt<Scalar>::staticOutput(GenDistrVector<Scalar>& sol, 
							   GenDistrVector<Scalar>& force,
							   double time)
{
  startTimerMemory(this->times->output, this->times->memoryOutput);
  dynamic_cast<GenDecDomain_opt<Scalar>*>(this->decDomain)->getDomain()->solInfo().setTimes(1.0e8, 1.0, 0.0);  
  dynamic_cast<typename GenMultiDomainStatic_opt<Scalar>::DomainType_opt*>(this->decDomain)
    ->postProcessing(sol,force,0, 0, static_cast<int>(time));
  stopTimerMemory(this->times->output, this->times->memoryOutput);
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenMultiDomainStatic_opt<Scalar>::getPseudoLoad(GenDistrVector<Scalar>& rhs,
                                                  GenDistrVector<Scalar>& sol)
{
  rhs.zero();
  // Derivative of the load vector for temperature, pressure and gravity  
  dynamic_cast<GenDecDomain_opt<Scalar>*>(this->decDomain)->buildDRHSforceDs(rhs);
  
  //Derivative of the elastic stiffness
  dynamic_cast<GenDecDomain_opt<Scalar>*>(this->decDomain)->buildDStiffnesDSmDU(rhs,sol);
}


//------------------------------------------------------------------------------
template<class Scalar>
void
GenMultiDomainStatic_opt<Scalar>::iniDensProj()
{
  //build density projector if needed
  const NodalDensityData& nodalDensData = static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData();
  const DensityProjData& densProjData   = static_cast<GeoSource_opt*>(geoSource)->getDensityProjData();
  if(densProjData.densProjFlag > 0)
    {
      assert(nodalDensData.nodDensFunc > 0);
      dynamic_cast<GenDecDomain_opt<Scalar>*>(this->decDomain)->buildNodalDensity(nodalDensData);
      dynamic_cast<GenDecDomain_opt<Scalar>*>(this->decDomain)->buildDensProj(densProjData);
      dynamic_cast<GenDecDomain_opt<Scalar>*>(this->decDomain)->scaleDensProj();
    }
}

//------------------------------------------------------------------------------
template<class Scalar>
void
GenMultiDomainStatic_opt<Scalar>::densProjectVector(GenDistrVector<Scalar>& v)
{
  const DensityProjData& densProjData = static_cast<GeoSource_opt*>(geoSource)->getDensityProjData();
  dynamic_cast<GenDecDomain_opt<Scalar>*>(this->decDomain)->densProjectVector(densProjData, v);
}
