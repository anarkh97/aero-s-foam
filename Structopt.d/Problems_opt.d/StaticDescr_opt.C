#include <Structopt.d/Problems_opt.d/StaticDescr_opt.h>

#ifdef STRUCTOPT

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic_opt<T, VectorType, SolverType>::preoptProcess()
{
  // Allocate space for the Static Timers
  this->times = new StaticTimers;
  
  startTimerMemory(this->times->preProcess, this->times->memoryPreProcess);
  
  // Makes renumbering, connectivities and dofsets
  domain->preProcessing();
  
  int numdof = domain->numdof();
  
  this->times->makeBCs -= getTime();
  boost::scoped_array<int> bc(new int[numdof]);
  this->bcx     = new T[numdof];
  
  // Make boundary conditions info
  domain->make_bc(bc.get(),this->bcx);
  this->times->makeBCs += getTime();
  
  // Now, call make_constrainedDSA(bc) to  built c_dsa 
  // that will incorporate all the boundary conditions info
  this->times->makeDOFs -= getTime();
  domain->make_constrainedDSA(bc.get());
  domain->makeAllDOFs();
  this->times->makeDOFs += getTime();
  
  // if we have initial displacements, we have to consider
  // the nonlinear tangent stiffness matrix instead of the
  // linear stiffness matrix. This is due to the prestress.
  this->kelArray  = 0;
  this->geomState = 0;
  this->allCorot  = 0;
  
  if(domain->numInitDisp6() > 0 && domain->solInfo().gepsFlg == 1) 
    {
      filePrint(stderr," ... Static Problem with Initial Displacement %d\n",
		domain->numInitDisp6());
      FullSquareMatrix *geomKelArray=0;
      domain->computeGeometricPreStress(this->allCorot, this->geomState, this->kelArray, this->times,
					geomKelArray);
    }
  
  stopTimerMemory(this->times->preProcess, this->times->memoryPreProcess);
  
  iniDensProj();

  // domain->template getSolverAndKuc<T>(solver, kuc, kelArray);
  dynamic_cast<Domain_opt*>(domain)->template getSolverAndKuc<T>(this->allOps, this->kelArray, 0);  //build Operator without factorization
  // PJSA 10-5-04 for freq sweep compatibility
  this->solver = this->allOps.sysSolver;                                 // also need M, Muc, C and Cuc (in allOps)
  this->kuc = this->allOps.Kuc;

  Rbm *rigidBodyModes = 0;
  int useProjector=domain->solInfo().filterFlags;
  if(useProjector) {
    std::cout << " ... RBMfilter Requested            ..." << std::endl;
    rigidBodyModes = domain->constructRbm();
    this->projector_prep(rigidBodyModes);
  }
  
  int useHzemFilter = domain->solInfo().hzemFilterFlag;
  if(useHzemFilter) 
    {
      std::cout << " ... HZEMfilter Requested           ..." << std::endl;
      rigidBodyModes = domain->constructHzem();
      this->projector_prep(rigidBodyModes);
    }

  
  return;
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic_opt<T, VectorType, SolverType>::reBuild()
{
  if (this->kuc) { this->kuc->zeroAll(); }
  
  this->allOps.Kuc       = this->kuc;
  this->allOps.sysSolver = this->solver;
  
  Rbm *rbm = 0; 
  
  iniDensProj();

  if (this->kelArray) { domain->getElasticStiff(this->kelArray); }
  
  domain->rebuildOps(this->allOps, 1.0, 0.0, 0.0, rbm, this->kelArray);
}


template<class T, class VectorType, class SolverType>
void
SingleDomainStatic_opt<T, VectorType, SolverType>::getPseudoLoad(VectorType &rhs,
								 VectorType &sol)
{
  rhs.zero();  
  // Derivative of the load vector for temperature, pressure and gravity  
  dynamic_cast<Domain_opt*>(domain)->buildDRHSforceDs(rhs);
  
  //Derivative of the elastic stiffness
  dynamic_cast<Domain_opt*>(domain)->buildDStiffnesDSmDU(rhs,sol,this->bcx);
}

template<class T, class VectorType, class SolverType>
double
SingleDomainStatic_opt<T, VectorType, SolverType>::getAdjPseudoLoad(VectorType &sol,
								    VectorType &adj)
{
  T part1 = 0;
  if(static_cast<GeoSource_opt*>(geoSource)->isRHSDsgnDep())
    {
      if(tmpAdjV.get()==0) { tmpAdjV.reset(new VectorType(adj.size())); }
      tmpAdjV->zero();  
      // Derivative of the load vector for temperature, pressure and gravity
      dynamic_cast<Domain_opt*>(domain)->buildDRHSforceDs(*tmpAdjV);
      part1 = (*tmpAdjV)*adj;
    }
  //Derivative of the elastic stiffness
  T part2 = dynamic_cast<Domain_opt*>(domain)->DStiffnesDSmDUAdj(sol,adj,this->bcx);
  return ScalarTypes::Real(part1 + part2);
}


template<class T, class VectorType, class SolverType>
void
SingleDomainPostProcessor_opt< T, VectorType, SolverType >::staticOutput(VectorType &sol, 
									 VectorType &force,
									 double time)
{
  this->times->output -= getTime();
  domain->solInfo().setTimes(1.0e8, 1.0, 0.0);
  dynamic_cast<Domain_opt*>(domain)->postProcessing(sol,this->bcx,force,time);
  this->times->output += getTime();
}

template<class T, class VectorType, class SolverType>
int
SingleDomainStatic_opt<T, VectorType, SolverType>::getNumLC()
{
  return dynamic_cast<Domain_opt*>(domain)->getNumLC();
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic_opt<T, VectorType, SolverType>::setActiveLC(int lc)
{
  dynamic_cast<Domain_opt*>(domain)->setActiveLC(lc);
}

template<class T, class VectorType, class SolverType>
SingleDomainPostProcessor_opt<T,VectorType,SolverType> *
SingleDomainStatic_opt<T, VectorType, SolverType>::getPostProcessor()
{
  return new SingleDomainPostProcessor_opt<T,VectorType,SolverType>(this->domain,this->bcx,this->times,this->solver);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic_opt<T, VectorType, SolverType>::iniDensProj()
{
  //build density projector if needed
  const NodalDensityData& nodalDensData = static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData();
  const DensityProjData& densProjData   = static_cast<GeoSource_opt*>(geoSource)->getDensityProjData();
  if(densProjData.densProjFlag > 0)
    {
      assert(nodalDensData.nodDensFunc > 0);
      dynamic_cast<Domain_opt*>(this->domain)->buildNodalDensity(nodalDensData);
      dynamic_cast<Domain_opt*>(this->domain)->buildDensProj(densProjData);
      dynamic_cast<Domain_opt*>(this->domain)->scaleDensProj();
    }
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic_opt<T, VectorType, SolverType>::densProjectVector(VectorType& v)
{
  const DensityProjData& densProjData = static_cast<GeoSource_opt*>(geoSource)->getDensityProjData();
  dynamic_cast<Domain_opt*>(domain)->densProjectVector(densProjData, v.data());  
}

#endif
