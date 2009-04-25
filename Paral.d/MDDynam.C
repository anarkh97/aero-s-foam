#include <iostream>
#include <Driver.d/Domain.h>
#include <Driver.d/DynamProbType.h>
#include <Paral.d/MDDynam.h>
#include <Threads.d/Paral.h>
#include <Driver.d/Dynam.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Paral.d/MDOp.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/Vector.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/SGISparseMatrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/Skyline.d/SGISky.h>
#include <Timers.d/GetTime.h>
#include <Control.d/ControlInterface.h>
#include <Threads.d/PHelper.h>
#include <Paral.d/GenMS.h>
#include <Paral.d/SubDOp.h>
#ifdef DISTRIBUTED
#include <Utils.d/DistHelper.h>
#include <Dist.d/DistDom.h>
#endif
#include <Driver.d/DecDomain.h>

#include <Hetero.d/DistFlExchange.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/SolverInfo.h>

class IntFullM;

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
               DistrVector *_v1, DistrVector*_v2, double c, CuCSparse **_Kuc)
{
 v1  = _v1;
 v2  = _v2;
 f   = _f;
 c1  = c;
 sd  = _sd;
 Kuc = _Kuc;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                  DistrVector *_v1, DistrVector*_v2, 
                  double c, double *_userDefDisps, DistrGeomState *_geomState)
{
 v1 = _v1;
 v2 = _v2;
 f  = _f;
 c1 = c;
 sd = _sd;

 userDefDisps = _userDefDisps;
 geomState = _geomState;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                  DistrVector *_v1)
{
 v1 = _v1;
 f  = _f;
 sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd)
{
  f = _f;
 sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                  DistrVector *_v1, DistrVector*_v2, DistrVector *_v3)
{
 v1 = _v1;
 v2 = _v2;
 v3 = _v3;
 f  = _f;
 sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                  DistrVector *_v1, DistrVector*_v2, DistrVector *_v3,  DistrVector *_v4)
{
 v1 = _v1;
 v2 = _v2;
 v3 = _v3;
 v4 = _v4;
 f  = _f;
 sd = _sd;
}


void
MultiDomainOp::runFor(int isub)
{
  (this->*f)(isub);
}

void
MultiDomainOp::computeExtForce(int isub)
{
  // Get the pointer to the part of the vector f corresponding to subdomain isub
  StackVector localF(v1->subData(isub),v1->subLen(isub));

  // Get the pointer to the part of the vector cnst_f corresponding to subdomain isub
  StackVector cnst_f(v2->subData(isub),v2->subLen(isub));

  // IF YOUR CODE CRASHES WITHIN THIS CALL, IT MEANS MultiDomainOp MUST BE UPDATED FOR PrevFrc
  int *userDataMap = sd[isub]->getUserDispDataMap();
  PrevFrc dummy(0);
  sd[isub]->computeExtForce4(dummy, localF, cnst_f, 0, c1, sd[isub]->getKuc(), userDefDisps, userDataMap);
 
  // build pressure forces (by convention pressure forces are follower for nonlinear dynamics)
  if(sd[isub]->pressureFlag() && domain->solInfo().isNonLin()) sd[isub]->buildPressureForce<double>(localF, (*geomState)[isub]);
}

void
MultiDomainOp::getConstForce(int isub)
{
 // Get the pointer to the part of the vector f correspoding to subdomain sNum
 StackVector f(v1->subData(isub),v1->subLen(isub));
 f.zero();

 // build gravity forces
 if(sd[isub]->gravityFlag()) sd[isub]->buildGravityForce<double>(f);

 // build pressure forces (by convention pressure forces are non-follower for linear dynamics)
 if(sd[isub]->pressureFlag() && !domain->solInfo().isNonLin()) sd[isub]->buildPressureForce<double>(f);
}

void
MultiDomainOp::makeAllDOFs(int isub)
{
  sd[isub]->makeAllDOFs();
}

void
MultiDomainOp::getInitState(int isub)
{
 StackVector disp( v1->subData(isub),v1->subLen(isub));
 StackVector veloc(v2->subData(isub),v1->subLen(isub));
 StackVector accel(v3->subData(isub),v1->subLen(isub));
 StackVector   v_p(v4->subData(isub),v1->subLen(isub));

 sd[isub]->initDispVeloc(disp, veloc, accel, v_p);
}

void
MultiDomDynPostProcessor::setPostProcessor(DistFlExchanger *exchanger) 
{
  distFlExchanger = exchanger;
}

void
MultiDomDynPostProcessor::setUserDefs(double **disps, double **vels) 
{
  usrDefDisps = disps;
  usrDefVels = vels;
}

void
MultiDomDynPostProcessor::dynamOutput(int tIndex, MDDynamMat &dynOps, DistrVector &distForce, 
                                      DistrVector *distAeroF, SysState<DistrVector>& distState)
{
  if(!times) times = new StaticTimers;
  startTimerMemory(times->output, times->memoryOutput);

  // PJSA 4-15-08: update bcx for time dependent prescribed displacements and velocities (previously done in computeExternalForce2)
  ControlLawInfo *claw = geoSource->getControlLaw();
  ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
  if(claw && claw->numUserDisp) {
    double t = double(tIndex)*domain->solInfo().dt;
    double *userDefineDisp = new double[claw->numUserDisp];
    double *userDefineVel  = new double[claw->numUserDisp];
    //cerr << "getting usdd at time " << t << " for dynamOutput\n";
    userSupFunc->usd_disp(t,userDefineDisp,userDefineVel);
    paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC, userDefineDisp, userDefineVel);
    delete [] userDefineDisp; delete [] userDefineVel;
  }

  decDomain->postProcessing(distState.getDisp(), distForce, 0.0, distAeroF, tIndex, &dynOps, &distState); 
  stopTimerMemory(times->output, times->memoryOutput);

  SolverInfo& sinfo = domain->solInfo();
  if(sinfo.aeroFlag >= 0)
    if(tIndex != sinfo.initialTimeIndex)
      distFlExchanger->sendDisplacements(distState, usrDefDisps, usrDefVels);
}

