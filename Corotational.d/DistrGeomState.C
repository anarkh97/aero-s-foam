#include <Corotational.d/DistrGeomState.h>
#include <Feti.d/DistrVector.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/DecDomain.h>
#include <Threads.d/PHelper.h>


DistrGeomState::DistrGeomState(DecDomain* domain)
{
 // Number of subdomains
 numSub = domain->getNumSub();

 // array of pointers to the subdomain's geometry states
 gs = new GeomState*[numSub];

 // parallel execution of subdomain geometry state construction
 execParal(numSub,this,&DistrGeomState::makeSubGeomStates,domain); 
}

// Subdomain geom state construction
void
DistrGeomState::makeSubGeomStates(int isub, DecDomain *domain)
{
 SubDomain *sd = domain->getSubDomain(isub);
 gs[isub] = new GeomState( *sd->getDSA(), *sd->getCDSA(), sd->getNodes() );
}

DistrGeomState::DistrGeomState(const DistrGeomState &g2)
{
  numSub = g2.getNumSub();
  gs = new GeomState*[numSub];
  execParal(numSub,this,&DistrGeomState::subCopyConstructor,g2);
}

void
DistrGeomState::subCopyConstructor(int isub, const DistrGeomState &g2)
{
  gs[isub] = new GeomState(*(g2[isub]));
}

// Subdomain update
void
DistrGeomState::subUpdate(int isub, DistrVector &v)
{
 StackVector vec(v.subData(isub), v.subLen(isub) );
 gs[isub]->update(vec);
}

void
DistrGeomState::subStep_update(int isub, DistrVector &v_n, double &delta, 
                               DistrGeomState &ss)
{
 StackVector vel(v_n.subData(isub),v_n.subLen(isub));
 GeomState &step = *ss[isub];
 gs[isub]->midpoint_step_update(vel,delta,step);
}

void
DistrGeomState::midpoint_step_update(DistrVector &veloc_n, double &delta, 
                                     DistrGeomState &ss)
{
 execParal3R(numSub,this,&DistrGeomState::subStep_update,veloc_n,delta,ss);
}

void
DistrGeomState::subInc_update(int isub,DistrVector &inc_vec, DistrGeomState &ss, bool zeroRot)
{
 StackVector v(inc_vec.subData(isub),inc_vec.subLen(isub));
 gs[isub]->get_inc_displacement(v,*ss[isub],zeroRot);
}

void
DistrGeomState::get_inc_displacement(DistrVector &inc_vec, DistrGeomState &ss, bool zeroRot)
{
 execParal3R(numSub,this,&DistrGeomState::subInc_update, inc_vec, ss, zeroRot);
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
 execParal3R(numSub,this,&DistrGeomState::subInterp,alpha,u,un);
}

void
DistrGeomState::subDiff(int isub, DistrGeomState &unp, DistrVector &un)
{
 StackVector u(un.subData(isub), un.subLen(isub) );
 GeomState &unpR = *unp[isub];
 gs[isub]->diff(unpR, u);
}

void
DistrGeomState::diff(DistrGeomState &unp, DistrVector &un)
{
 execParal2R(numSub,this,&DistrGeomState::subDiff,unp,un);
}


void
DistrGeomState::update(DistrVector &v)
{
   execParal1R(numSub, this,&DistrGeomState::subUpdate, v);
}

DistrGeomState &
DistrGeomState::operator=(DistrGeomState &unp)
{
  execParal1R(numSub, this,&DistrGeomState::subCopy, unp);
  return *this;
}

void
DistrGeomState::subCopy(int isub, DistrGeomState &unp)
{
  GeomState &unpR = *unp[isub];
  *(gs[isub]) = unpR; 
}

