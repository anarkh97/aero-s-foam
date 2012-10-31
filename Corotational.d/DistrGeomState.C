#include <Corotational.d/DistrGeomState.h>
#include <Feti.d/DistrVector.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/DecDomain.h>
#include <Threads.d/PHelper.h>
#include <Corotational.d/TemperatureState.h>
#ifdef USE_MPI
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

DistrGeomState::DistrGeomState(DecDomain* domain)
{
 // Number of subdomains
 numSub = domain->getNumSub();

 // array of pointers to the subdomain's geometry states
 gs = new GeomState*[numSub];

 // array of data structures to store lagrange multipliers 
 mu = new std::map<std::pair<int,int>,double>[numSub];
 lambda = new std::vector<double>[numSub];

 // parallel execution of subdomain geometry state construction
 execParal(numSub,this,&DistrGeomState::makeSubGeomStates,domain); 
}

// Subdomain geom state construction
void
DistrGeomState::makeSubGeomStates(int isub, DecDomain *domain)
{
 SubDomain *sd = domain->getSubDomain(isub);
 if(sd->solInfo().soltyp == 2)
   gs[isub] = new TemperatureState( *sd->getDSA(), *sd->getCDSA(), sd->getNodes() );
 else
   gs[isub] = new GeomState( *sd->getDSA(), *sd->getCDSA(), sd->getNodes(), &sd->getElementSet() );
}

DistrGeomState::DistrGeomState(const DistrGeomState &g2)
{
  numSub = g2.getNumSub();
  gs = new GeomState*[numSub];
  mu = new std::map<std::pair<int,int>,double>[numSub];
  lambda = new std::vector<double>[numSub];
  execParal(numSub,this,&DistrGeomState::subCopyConstructor,g2);
}

void
DistrGeomState::subCopyConstructor(int isub, const DistrGeomState &g2)
{
  TemperatureState* ts;
  if(ts = dynamic_cast<TemperatureState*>(g2[isub])) gs[isub] = new TemperatureState(*ts);
  else gs[isub] = new GeomState(*(g2[isub]));

  mu[isub] = g2.mu[isub];
  lambda[isub] = g2.lambda[isub];
}

DistrGeomState::~DistrGeomState()
{
  for(int i=0; i<numSub; ++ i) delete gs[i];
  delete [] gs;
  delete [] mu;
  delete [] lambda;
}

// Subdomain update
void
DistrGeomState::subUpdate(int isub, DistrVector &v)
{
  StackVector vec(v.subData(isub), v.subLen(isub));
  gs[isub]->update(vec);
}

void
DistrGeomState::subSetVelocity(int isub, DistrVector &v, DistrVector &a)
{
  StackVector vsub(v.subData(isub), v.subLen(isub));
  StackVector asub(a.subData(isub), a.subLen(isub));
  gs[isub]->setVelocity(vsub, asub);
}

void
DistrGeomState::subStep_update(int isub, DistrVector &v_n, DistrVector &a_n, 
                               double &delta, DistrGeomState &ss,
                               double beta, double gamma, double alphaf, double alpham, bool zeroRot)
{
 StackVector vel(v_n.subData(isub),v_n.subLen(isub));
 StackVector acc(a_n.subData(isub),a_n.subLen(isub));
 GeomState &step = *ss[isub];
 gs[isub]->midpoint_step_update(vel,acc,delta,step,beta,gamma,alphaf,alpham,zeroRot);
}

void
DistrGeomState::midpoint_step_update(DistrVector &veloc_n, DistrVector &accel_n,
                                     double &delta, DistrGeomState &ss,
                                     double beta, double gamma, double alphaf, double alpham, bool zeroRot)
{
 execParal9R(numSub,this,&DistrGeomState::subStep_update,veloc_n,accel_n,delta,ss,beta,gamma,alphaf,alpham,zeroRot);
}

void
DistrGeomState::subInc_get(int isub,DistrVector &inc_vec, DistrGeomState &ss, bool zeroRot)
{
 StackVector v(inc_vec.subData(isub),inc_vec.subLen(isub));
 gs[isub]->get_inc_displacement(v,*ss[isub],zeroRot);
}

void
DistrGeomState::get_inc_displacement(DistrVector &inc_vec, DistrGeomState &ss, bool zeroRot)
{
 execParal3R(numSub,this,&DistrGeomState::subInc_get, inc_vec, ss, zeroRot);
}

void
DistrGeomState::subTot_get(int isub, DistrVector &tot_vec)
{
 StackVector v(tot_vec.subData(isub), tot_vec.subLen(isub));
 gs[isub]->get_tot_displacement(v);
}

void
DistrGeomState::get_tot_displacement(DistrVector &tot_vec)
{
 execParal1R(numSub, this, &DistrGeomState::subTot_get, tot_vec);
}

void
DistrGeomState::subInterp(int isub, double &alpha, DistrGeomState &u,
                          DistrGeomState &un)
{
 GeomState &uR  = *u[isub];
 GeomState &unR = *un[isub];
 gs[isub]->interp(alpha, uR, unR );
}

void
DistrGeomState::interp(double alpha, DistrGeomState &u, DistrGeomState &un)
{
  execParal3R(numSub, this, &DistrGeomState::subInterp, alpha, u, un);
}

void
DistrGeomState::subDiff(int isub, DistrGeomState &unp, DistrVector &un)
{
 StackVector u(un.subData(isub), un.subLen(isub));
 GeomState &unpR = *unp[isub];
 gs[isub]->diff(unpR, u);
}

void
DistrGeomState::diff(DistrGeomState &unp, DistrVector &un)
{
  execParal2R(numSub, this, &DistrGeomState::subDiff, unp, un);
}

void
DistrGeomState::update(DistrVector &v)
{
  execParal1R(numSub, this, &DistrGeomState::subUpdate, v);
}

void
DistrGeomState::setVelocity(DistrVector &v, DistrVector &a)
{
  execParal2R(numSub, this, &DistrGeomState::subSetVelocity, v, a);
}

DistrGeomState &
DistrGeomState::operator=(DistrGeomState &unp)
{
  execParal1R(numSub, this, &DistrGeomState::subCopy, unp);
  return *this;
}

void
DistrGeomState::subCopy(int isub, DistrGeomState &unp)
{
  GeomState &unpR = *unp[isub];
  *(gs[isub]) = unpR; 

  mu[isub] = unp.mu[isub];
  lambda[isub] = unp.lambda[isub];
}

int 
DistrGeomState::getTotalNumElemStates()
{
  int ret = 0;
  for(int i=0; i<numSub; ++i) ret += gs[i]->getTotalNumElemStates();
#ifdef USE_MPI
  ret = structCom->globalSum(ret);
#endif
  return ret;
}
