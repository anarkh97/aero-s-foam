#include <Timers.d/GetTime.h>
#include <Math.d/mathUtility.h>
#include <unistd.h>
#include <Paral.d/SubDOp.h>
#include <Driver.d/SysState.h>

//-------------------------------------------------------------------------------------------
template<class VecType>
SysState<VecType> & SysState<VecType>::operator=(const SysState<VecType> &v2)
{
   d_n = v2.getDispConst();
   v_n = v2.getVelocConst();
   a_n = v2.getAccelConst();
   v_n_p = v2.getPrevVelocConst();

   return *this;

}

//-------------------------------------------------------------------------------------------

template <class VecType,
          class ProblemDescriptor> 
NewmarkWorkVec<VecType,ProblemDescriptor>::NewmarkWorkVec(int _typ,ProblemDescriptor *probDesc)
{
   typ=_typ;

   switch (typ) {

     case -1:
       d_n_p  = new VecType( probDesc->solVecInfo() );
       o_n_p  = new VecType( probDesc->solVecInfo() );
       ext_f  = new VecType( probDesc->solVecInfo() );
       rhs    = new VecType( probDesc->solVecInfo() );
       break;
     case 0:
       tmp1   = new VecType( probDesc->solVecInfo() );
       tmp2   = new VecType( probDesc->solVecInfo() );
       fint   = new VecType( probDesc->solVecInfo() );
       ext_f  = new VecType( probDesc->solVecInfo() );
       v_n_h  = new VecType( probDesc->solVecInfo() );
       break;
     case 1:
       d_n_p  = new VecType( probDesc->solVecInfo() );
       v_n_p  = new VecType( probDesc->solVecInfo() );
       a_n_p  = new VecType( probDesc->solVecInfo() );
       rhs    = new VecType( probDesc->solVecInfo() );
       ext_f  = new VecType( probDesc->solVecInfo() );
       d_n_h  = new VecType( probDesc->solVecInfo() );
       v_n_h  = new VecType( probDesc->solVecInfo() );
       Md_n_h = new VecType( probDesc->solVecInfo() );
       Cd_n_h = new VecType( probDesc->solVecInfo() );
       tmp1   = new VecType( probDesc->solVecInfo() );
       dnc    = new VecType( probDesc->bcInfo() );
       vnc    = new VecType( probDesc->bcInfo() );
       anc    = new VecType( probDesc->bcInfo() );
       break;
   }
}

//-----------------------------------------------------------------------------
template <class VecType,
          class ProblemDescriptor>
NewmarkWorkVec<VecType,ProblemDescriptor> & NewmarkWorkVec<VecType,ProblemDescriptor>:: operator=(const NewmarkWorkVec<VecType,ProblemDescriptor> &v)
{

     if (typ!=v.typ) {
           
         switch (typ) {
                                                                                                  
	     case -1:
       		delete  d_n_p;
       		delete  o_n_p;
       		delete  ext_f;
       		delete  rhs;
       		break;
     	     case 0:
       		delete  tmp1;
       		delete  tmp2;
       		delete  fint;
       		break;
     	     case 1:
       		delete  d_n_p;
       		delete  v_n_p;
       		delete  a_n_p;
       		delete  rhs;
       		delete  ext_f;
       		delete  d_n_h;
       		delete  v_n_h;
       		delete  Md_n_h;
       		delete  Cd_n_h;
       		delete  tmp1;
       		delete  dnc;
       		delete  vnc;
       		delete  anc;
       		break;
      	 }
         typ=v.typ;
     }

     switch (typ) {
                                                                                                  
         case -1:
            d_n_p  = new VecType();
            *d_n_p = v.get_d_n_pConst();
            o_n_p  = new VecType();
            *o_n_p = v.get_o_n_pConst();
            ext_f  = new VecType();
            *ext_f = v.get_ext_fConst();
            rhs    = new VecType();
            *rhs   = v.get_rhsConst();
            break;
         case 0:
            tmp1  = new VecType();
            *tmp1 = v.get_tmp1Const();
            tmp2  = new VecType();
            *tmp2 = v.get_tmp2Const();
            fint  = new VecType();
            *fint = v.get_fintConst();
            break;
         case 1:
            d_n_p  = new VecType();
            *d_n_p = v.get_d_n_pConst();
            v_n_p  = new VecType();
            *v_n_p = v.get_v_n_pConst();
            a_n_p  = new VecType();
            *a_n_p = v.get_a_n_pConst();
            rhs    = new VecType();
            *rhs   = v.get_rhsConst();
            ext_f  = new VecType();
            *ext_f = v.get_ext_fConst();
            d_n_h  = new VecType();
            *d_n_h = v.get_d_n_hConst();
            v_n_h  = new VecType();
            *v_n_h = v.get_v_n_hConst();
            Md_n_h = new VecType();
            *Md_n_h = v.get_Md_n_hConst();
            Cd_n_h = new VecType();
            *Cd_n_h = v.get_Cd_n_hConst();
            tmp1   = new VecType();
            *tmp1  = v.get_tmp1Const();
            dnc    = new VecType();
            *dnc   = v.get_dncConst();
            vnc    = new VecType();
            *vnc   = v.get_vncConst();
            anc    = new VecType();
            *anc   = v.get_ancConst();
            break;
   }

}
//------------------------------------------------------------------------------

template <class VecType, 
          class ProblemDescriptor> 
NewmarkWorkVec<VecType,ProblemDescriptor>::~NewmarkWorkVec()
{
   switch (typ) {
   
     case -1:
       delete  d_n_p;
       delete  o_n_p;
       delete  ext_f;
       delete  rhs;
       break;
     case 0:
       delete  tmp1;
       delete  tmp2;
       delete  fint;
       break;
     case 1:
       delete  d_n_p;  
       delete  v_n_p;  
       delete  a_n_p;  
       delete  rhs;    
       delete  ext_f; 
       delete  d_n_h;  
       delete  v_n_h;
       delete  Md_n_h;
       delete  Cd_n_h;
       delete  tmp1;
       delete  dnc;
       delete  vnc;
       delete  anc;
       break;
   }
}

//------------------------------------------------------------------------------

template<
     class DynOps,             // Data Structure for K, C, M and dynMat
     class VecType,            // Vector type used
     class PostProcessor,      
     class ProblemDescriptor,
     class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::solve()
{
   // Construct renumbering, connectivities, dofsets
   // and boundary conditions. Suppose to be common for
   // all the different dynamic schemes.

   probDesc->preProcess();

   postProcessor = probDesc->getPostProcessor();

   // Allocate vectors for displacement, velocity, 
   //                      acceleration and last velocity
   d_n = new VecType( probDesc->solVecInfo() );
   v_n = new VecType( probDesc->solVecInfo() );
   a_n = new VecType( probDesc->solVecInfo() );
   v_p = new VecType( probDesc->solVecInfo() );

   // Get time loop information 
   aeroAlg = probDesc->getAeroAlg();
   if(aeroAlg >= 0 || probDesc->getThermoeFlag() >= 0) {
     probDesc->computeTimeInfo(); // computes a new tmax if necessary
   }
   probDesc->getTimes( dt, tmax );

   // Set up initial conditions
   curState = new SysState<VecType>( *d_n, *v_n, *a_n, *v_p);
   probDesc->getInitState( *curState );

   // The aeroPreProcess is done later now, so that the correct
   // time-step is passed to the fluid in the explicit case
   // However, if we are doing a ping-pong, first we run aeroPreProcess and return 
   if(aeroAlg == 1 || aeroAlg == 8) {
     probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p );
     return;
   }

   aeroForce = (aeroAlg >= 0 || probDesc->getThermoeFlag() >= 0) ? new VecType(probDesc->solVecInfo()) : 0;

   // Build time independent forces i.e. gravity force, pressure force
   constForce = new VecType( probDesc->solVecInfo() );
   //probDesc->getConstForce( *constForce ); PJSA moved to after buildOps

   // Get SteadyState Flag and Parameters
   probDesc->getSteadyStateParam(steadyFlag, steadyMin, steadyMax, steadyTol);

   // Get Time Integration Scheme
   algType = probDesc->getTimeIntegration();
   if(aeroAlg == 10 && algType != 1) {
      fprintf(stderr, "WARNING: B0 AERO type only valid for quasi-static.  Running QUASISTATICS for 1 iteration\n");
      algType = 1;
   }
   
   double timeLoop = -getTime();
   switch ( algType )
   {
     // Newmark
     case 0:

       // Get Newmark Parameter  
       probDesc->getNewMarkParameters( beta, gamma, alphaf, alpham );
       
       // ... Newmark Beta == 0 -> CENTRAL DIFFERENCES ALGORITHM
       // ... Newmark Beta != 0 -> NEWMARK ALGORITHM
       //     This is now the Generalized Alpha Method, of which
       //     the NEWMARK algorithm is a subset.

       if(beta == 0.0) {

         // Defining Working Arrays   
         workVec = new NewmarkWorkVec<VecType,ProblemDescriptor>(0,probDesc);
 
         if(domain->solInfo().order != 1 && gamma != 0.5) {
           filePrint(stderr," ### WATCH gamma for explicit Newmark is set to 0.5 ###\n");
           gamma=0.5;
         }

         // Build Necessary Operators (K, M)
         bool fourthOrder = probDesc->getDomain()->solInfo().modifiedWaveEquation;
         if (fourthOrder) 
           dynOps = probDesc->buildOps(1.0, dt*gamma, dt*dt*beta);
         else
           dynOps = probDesc->buildOps(1.0, 0.0, 0.0);

         probDesc->getConstForce( *constForce );

         // Check stability time step
         if(domain->solInfo().stable) probDesc->computeStabilityTimeStep(dt, *dynOps);

         if(aeroAlg == 20) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_n ); //e Se Eq. 51 of C.Farhat et al. IJNME(2010) Robust and provably ... 
         else
           if(aeroAlg >= 0) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p );
         if(probDesc->getThermoeFlag() >= 0) probDesc->thermoePreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getAeroheatFlag() >= 0) probDesc->aeroHeatPreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getThermohFlag() >= 0) probDesc->thermohPreProcess(*d_n, *v_n, *v_p);
         

         // Time Integration Loop 
         explicitNewmarkLoop( *curState, *constForce, *dynOps, *workVec, dt, tmax);
       } 
       else {

         if(aeroAlg >= 0) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p );
         if(probDesc->getThermoeFlag() >= 0) probDesc->thermoePreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getAeroheatFlag() >= 0) probDesc->aeroHeatPreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getThermohFlag() >= 0) probDesc->thermohPreProcess(*d_n, *v_n, *v_p);

         // Defining Working Arrays   
         workVec = new NewmarkWorkVec<VecType,ProblemDescriptor>(1,probDesc);

         if(domain->solInfo().order == 1) { // heat
           // build K, M and dynK = (M + (dt*gamma)*K
           dynOps = probDesc->buildOps(1, 0, dt*gamma);
         }
         else { // mech and acou
           // Build K, M, C and dynK = ((1-alpham)/(1-alphaf))*M + (dt*gamma)*C + (dt*dt*beta)*K
           dynOps = probDesc->buildOps(((1.0-alpham)/(1.0-alphaf)), dt*gamma, dt*dt*beta);
         }

         probDesc->getConstForce( *constForce );

         // Time Integration Loop 
         implicitNewmarkLoop( *curState, *constForce, *dynOps, *workVec, dt, tmax);
       }
       break;
       
     // Quasi-Static
     case 1:
       cerr << "here in Driver.d/DynamProbType.C #1\n";
       if(aeroAlg >= 0) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p );
       if(probDesc->getThermoeFlag() >= 0) probDesc->thermoePreProcess(*d_n, *v_n, *v_p);
       if(probDesc->getAeroheatFlag() >= 0) probDesc->aeroHeatPreProcess(*d_n, *v_n, *v_p);
       if(probDesc->getThermohFlag() >= 0) probDesc->thermohPreProcess(*d_n, *v_n, *v_p);

       probDesc->getQuasiStaticParameters(maxVel, delta);
       
       // Defining Working Arrays   
       workVec = new NewmarkWorkVec<VecType,ProblemDescriptor>(-1,probDesc);

       // Build Necessary Operators (only K!)
       dynOps = probDesc->buildOps(0,0,1); 
 
       probDesc->getConstForce( *constForce );
       
       // Quasi-Static Loop
       quasistaticLoop( *curState, *constForce, *dynOps, *workVec, dt, tmax, aeroAlg);     
   }
   timeLoop += getTime();
   probDesc->printTimers(dynOps, timeLoop);
   

   // Delete arrays
   delete curState;
   delete workVec;

   delete d_n;
   delete v_n;
   delete a_n;
   delete constForce;
   delete v_p;

}

// -----------------------------------------------------------------------------//
//
// ... Quasi-Static Time Loop
//
// -----------------------------------------------------------------------------

template< class DynOps,        class VecType, 
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::quasistaticLoop(SysState<VecType>& curState, VecType& constForce,
                  DynOps& dynOps, 
		  NewmarkWorkVec<VecType,ProblemDescriptor>& workVec,
                  double dt, double tmax, int aeroFlg)
{ 
   filePrint(stderr, " ... Quasistatic loop               ... \n");
   if (aeroFlg == 10)  steadyMax = 1;
   // get initial displacements

   VecType &d_n = curState.getDisp();
   
   // Initialize some vectors 

   VecType  &d_n_p = workVec.get_d_n_p();
   VecType    &rhs = workVec.get_rhs();
   VecType  &ext_f = workVec.get_ext_f();

   // Initialize some parameters
  // Get initial Time
   int tIndex;
   int initIndex;
   double initialTime = 0.0;

   probDesc->getInitialTime(initIndex, initialTime);
   double initExtForceNorm = probDesc->getInitialForceNorm();
   tIndex = initIndex;

   int iSteady  = 0;

   double forceRef;

   double relaxFac = maxVel;
   
   // Output state of model

   postProcessor->dynamOutput( tIndex, dynOps, ext_f, aeroForce, curState );

   //-----------------------------------------------------------------------
   // ... BEGIN MAIN TIME-LOOP
   //-----------------------------------------------------------------------

   double totalTime = -getTime();

  for (tIndex = tIndex+1; tIndex <= steadyMax; tIndex++) {

    // ... call projector for RBMs in case of rbmfilter level 2
    if (probDesc->getFilterFlag() == 2) probDesc->project( d_n );

    // ... compute external force
    probDesc->computeExtForce2( curState, ext_f, constForce, tIndex, (double)tIndex*delta, aeroForce);

    // ... build force reference norm 
    if (tIndex==initIndex+1) {
      if (initExtForceNorm == 0.0)
        forceRef=ext_f.norm();
      else 
        forceRef = initExtForceNorm;
      if(verboseFlag) filePrint(stderr, " ... Initial Force: %8.2e ...\n", forceRef);
    }

    // ... build internal force 
    probDesc->getInternalForce(d_n, rhs, (double)tIndex*delta);

    // ... check for convergence
    double relres = 0.0;
    if (forceRef != 0.0)  relres = norm(rhs-ext_f)/forceRef;
    else {
      relres = norm(rhs-ext_f);
      filePrint(stdout, " ... WARNING: Reference External Force is zero, Relative residual is absolute error norm ...\n");
    }

    if(relres <= steadyTol && delta == 0) iSteady = 1;

    if(aeroAlg >= 0 || probDesc->getAeroheatFlag() >= 0) {
      filePrint(stderr," ... Pseudo-Step = %d  Rel. Res. = %10.4e ...\n",tIndex, relres);
  
      // command communication with fluid
      if(tIndex == steadyMax && !iSteady) { 
        probDesc->processLastOutput();
        if(aeroFlg != 10) {
          postProcessor->dynamOutput( tIndex, dynOps, ext_f, aeroForce, curState );
          probDesc->cmdCom(1); 
          break;
        }
      }
      else
        iSteady = probDesc->cmdCom(iSteady);
    }

    // ... stop quasi-transient simulation if converged
    if(iSteady) {
      filePrint(stderr," --------------------------------------\n");
      filePrint(stderr," ... Quasistatic Analysis Converged After %d Steps ...\n",tIndex);
      filePrint(stderr," --------------------------------------\n");
      probDesc->processLastOutput();
      postProcessor->dynamOutput( tIndex, dynOps, ext_f, aeroForce, curState );
      break; 
    }
    else if (tIndex == steadyMax)
      probDesc->processLastOutput();

    // ... save load vector
    rhs=ext_f;

    // ... solve System for current load
    dynOps.dynMat->reSolve( rhs );

    // ... compute displacement increment;
    d_n_p.linC(rhs, -1.0, d_n);

    // ... apply relaxation factor
    d_n_p *= relaxFac;

    // ... update solution
    d_n   += d_n_p;
//   filePrint(stderr, " ... ||u|| = %e\n", d_n.norm());

    // ... output current solution and send displacements to fluid
    postProcessor->dynamOutput( tIndex, dynOps, ext_f, aeroForce, curState );
  }

  if (!iSteady && aeroAlg != 10) {
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," ... Quasistatic Analysis Did Not Converge After %d Steps ...\n",tIndex);
    filePrint(stderr," --------------------------------------\n");
  }

  // ... output CPU time spent in quasi-static loop
  totalTime += getTime();
#ifdef PRINT_TIMERS
  filePrint(stderr," ... Total Loop Time = %.2e s   ...\n",totalTime/1000.0);
#endif

}

// -----------------------------------------------------------------------------//
// ... GENERAL IMPLICIT NEWMARK TIME LOOP
//
// DESCRIPTION: this routine implements the implicit Generalized Alpha method 
//              for time integration solution of the second-order differential equation:
//              M*(d^2u/dt^2) + C*(du/dt) + K*u = fext(u)
//
// WARNINGS   : 1) Second-order accuracy is obtained if and only if gamma = 1/2
//              2) Unconditionally stable when 2*beta >= gamma >= 1/2
//              3) This function does not work for beta == 0.0 ... we use the
//                 explicitNewmarkLoop function in this case
//
//  See J. Chung and G.M. Hulbert, "A time integration algorithm for structural
//  dynamics with improved numerical dissipation: the Generalized Alpha Method",
//  Journal of Applied Mechanics, no 60, 1993, pp 371-375.
//
//  This is a more general method that uses two parameters, alphaf and alpham to 
//  control numerical disspation.  In particular, it allows the user to limit
//  high frequency disspation while minimizing unwanted low frequency disspation.
//  The Newmark algorithm becomes a particular case of this method, achieved with
//  alphaf = alpham = 0.
// -----------------------------------------------------------------------------
template< class DynOps,        class VecType, 
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::implicitNewmarkLoop(SysState<VecType>& curState, VecType& constForce,
                      DynOps& dynOps, 
		      NewmarkWorkVec<VecType,ProblemDescriptor>& workVec,
                      double dt, double tmax)
{
   filePrint(stderr, " ... Implicit Newmark Time Integration Scheme: beta = %4.2f, gamma = %4.2f, alphaf = %4.2f, alpham = %4.2f ...\n",beta,gamma,alphaf,alpham);

   int parity = 0;
   SysState<VecType> *bkState = 0;
   // Allocate backup state for A5 algorithm
   if(aeroAlg == 5) {
     VecType d_bk(probDesc->solVecInfo());
     VecType v_bk(probDesc->solVecInfo());
     VecType a_bk(probDesc->solVecInfo());
     VecType v_p_bk(probDesc->solVecInfo());
     bkState = new SysState<VecType>(d_bk, v_bk, a_bk, v_p_bk);
   }

   VecType &d_n = curState.getDisp();
   VecType &v_n = curState.getVeloc();
   VecType &a_n = curState.getAccel();
   VecType &v_p = curState.getPrevVeloc();

   // project initial displacements in case of rbmfilter
   if(probDesc->getFilterFlag() > 0) {
      probDesc->project( d_n );
      probDesc->project( v_n );
   }

   // Some vectors for implicit Newmark time loop
   VecType  &d_n_p = workVec.get_d_n_p();
   VecType  &v_n_p = workVec.get_v_n_p();
   VecType  &a_n_p = workVec.get_a_n_p();
   VecType   & rhs = workVec.get_rhs();
   VecType  &ext_f = workVec.get_ext_f();
   VecType  &d_n_h = workVec.get_d_n_h();
   VecType  &v_n_h = workVec.get_v_n_h();
   VecType &Md_n_h = workVec.get_Md_n_h();
   VecType &Cd_n_h = workVec.get_Cd_n_h();
   VecType   &tmp1 = workVec.get_tmp1();

   // New vector added for prescribed boundary conditions
   //VecType &tmp1 = workVec.get_tmp1();
   //tmp1.zero();

   // These vectors have length equal to the number of constraints
   // ( number of dirichlet boundary conditions )
   //VecType &dnc = workVec.get_dnc();
   //VecType &vnc = workVec.get_vnc();
   //VecType &anc = workVec.get_anc();

   // zero these vectors
   //dnc.zero();
   //vnc.zero();
   //anc.zero();

   // Get initial time and time index
   int n = 0;
   double t = 0.0;
   probDesc->getInitialTime(n, t);

   // Compute the external force at the initial time: fext^0
   probDesc->computeExtForce2(curState, ext_f, constForce, -1, t, aeroForce, gamma, alphaf);

   // Compute the initial acceleration: a^0 = M^{-1}(fext^0 - Ku^0 - Cv^0)
   if(domain->solInfo().iacc_switch && dynOps.Msolver) {
     if(domain->solInfo().order == 1) {
       if(verboseFlag) filePrint(stderr," ... Computing initial first time derivative of temperature ...\n");
       dynOps.K->mult(d_n, tmp1);
       v_n.linC(ext_f, -1.0, tmp1);
       dynOps.Msolver->reSolve(v_n);
     }
     else {
       if(verboseFlag) filePrint(stderr," ... Computing initial acceleration ...\n");
       dynOps.K->mult(d_n, tmp1);
       a_n.linC(ext_f, -1.0, tmp1);
       if(dynOps.C) {
         dynOps.C->mult(v_n, tmp1);
         a_n -= tmp1;
       }
       dynOps.Msolver->reSolve(a_n);
       if(probDesc->getFilterFlag() == 2) probDesc->project(a_n);
     }
   }

   // Output initial state of model
   postProcessor->dynamOutput(n, dynOps, ext_f, aeroForce, curState);

   // ... BEGIN MAIN TIME-LOOP
   double totalTime = -getTime();
   char ch[4] = { '|', '/', '-', '\\' };

   for( ; t < tmax-0.01*dt; t += dt) {
     if(aeroAlg < 0) {
       filePrint(stderr,"\r  %c  Time Integration Loop: t = %9.3e, %3d%% complete ",
                 ch[int((totalTime + getTime())/250.)%4], t+dt, int((t+dt)/(tmax-0.01*dt)*100));
     }

     // ... For Aeroelastic A5 Algorithm, Do restore and backup here
     if(aeroAlg == 5) probDesc->a5StatusRevise(parity, curState, *bkState);

     // Mode decomposition of displacement
     if(probDesc->getModeDecompFlag()) probDesc->modeDecomp(t, n, d_n);

     // ... Construct force vector, includes time-independent constant force

     // ... Compute external force at time t+dt*(1-alphaf)
     // (if alphaf=0.5, external force is at time t+0.5*dt)
     probDesc->computeExtForce2(curState, ext_f, constForce, n, 
                                t+(dt*(1-alphaf)), aeroForce, gamma, alphaf);

     if(domain->solInfo().order == 1) { // heat (XXXX CURRENTLY ONLY IMPLEMENTED FOR alpham = alphaf = 1/2)
       // Solve for temperature: d^{n+1/2} = (M + gamma*dt*K)^{-1}(gamma*dt*f^{n+1/2} + M*(d^n+dt/2*(1-2*gamma)*v^n))
       d_n_h.linC(d_n, dt/2.0*(1.0-2.0*gamma), v_n);
       dynOps.M->mult( d_n_h, Md_n_h );
       rhs.linC( Md_n_h, gamma*dt, ext_f );
       dynOps.dynMat->reSolve( rhs );

       //if(probDesc->getHzemFlag() == 2) probDesc->tempProject(rhs);

       // Extrapolate temperature solution to t^{n+1} : d^{n+1} = 2*d^{n+1/2} - d^n
       d_n.linC(2.0, rhs, -1.0, d_n);

       // Compute the first time derivative of temperature at t^{n+1}: v^{n+1} = 2/(gamma*dt)*(d^{n+1/2 - d^n) - (1-gamma)/(gamma)*v^n
       v_n_p.linC(2.0/(gamma*dt), d_n, -2.0/(gamma*dt), rhs);
       if(gamma != 1.0) v_n_p.linAdd(-(1.0-gamma)/gamma, v_n);
     }
     else { // mech, acou
       // ... Construct R.H.S. vector
       // ... d_n_h = ((1-alpham)/(1-alphaf))*d_n 
       //           + dt*(1-alpham)*v_n 
       //           + dt*dt*((1-alpham)/2 - beta)*a_n
       // (if alphaf=alpham=0.5,beta=0.25: d_n_h = d_n + dt*0.5*v_n + zero*a_n)

       // First: d_n_h = ((1-alpham)/(1-alphaf))*d_n + dt*(1-alpham)*v_n 
       d_n_h.linC( ((1.0-alpham)/(1.0-alphaf)), d_n, dt*(1.0-alpham), v_n );

       // Second: d_n_h = d_n_h + dt*dt*((1-alpham)/2 - beta)*a_n
       d_n_h.linC( d_n_h, dt*dt*(0.5*(1.0-alpham)-beta), a_n );

       // Third: Multiply by Mass Matrix M
       dynOps.M->mult( d_n_h, Md_n_h );

       // Accumulate in rhs vector: rhs = Md_n_h + beta*dt*dt*ext_f
       rhs.linC( Md_n_h, beta*dt*dt, ext_f );

       if(dynOps.C) {
         // ... d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma(1-alphaf))*v_n
         // ...       - dt*dt*dt*0.5(1-alphaf)*(2*beta - gamma)*a_n
         // (if alphaf=alpham=gamma=0.5,beta=0.25: d_n_h = dt*0.5*d_n - zero*v_n - zero*a_n)

         // ... d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma(1-alphaf))*v_n
         d_n_h.linC( dt*gamma, d_n, dt*dt*(gamma*(1.0-alphaf)-beta), v_n );

         // ... d_n_h = d_n_h - dt*dt*dt*0.5(1-alphaf)*(2*beta - gamma)*a_n
         d_n_h.linC( d_n_h, dt*dt*dt*0.5*(alphaf-1.0)*(2.0*beta - gamma), a_n );

         // Multiply by Damping Matrix C
         dynOps.C->mult( d_n_h, Cd_n_h );
  
         // Accumulate in rhs vector 
         rhs += Cd_n_h;

       }
       // add prescribed boundary condition contributions to rhs vector
       // and update accelerations at prescribed dofs
       // probDesc->addPrescContrib( dynOps.M12, dynOps.C12, dnc, vnc, anc, tmp1, t);
       //rhs += tmp1;
       //tmp1.zero();

       dynOps.dynMat->reSolve( rhs ); // Now rhs contains d_(n+1-alphaf)

       // call projector for RBMs in case of rbmfilter level 2
       if (probDesc->getFilterFlag() == 2) probDesc->project( rhs );

       // one time step forward
       // d_n_p = 1/(1-alphaf)*[d_(n+1-alphaf)-alphaf*d_n] = d_(n+1)
       d_n_p.linC(rhs,(-1.0*alphaf),d_n);
       d_n_p *= (1.0/(1-alphaf));
   
       // a_n_p = 1/(dt^2*beta)*[d_(n+1)-d_n] - 1/(dt*beta)*v_n + (1-1/(2*beta))*a_n = a_(n+1)
       a_n_p = d_n_p;
       a_n_p -= d_n;
       a_n_p *= (1.0/(dt*dt*beta));
       a_n_p.linAdd( (-1.0/(dt*beta)), v_n, (1.0-1.0/(2.0*beta)), a_n );

       // v_n_h = gamma/(beta*dt)*[d_(n+1-alphaf) - d_n] + (1.0-(1.0-alphaf)*gamma/beta)*v_n + dt*(1.0-alphaf)*(2.0*beta-gamma)/(2*beta)*a_n
       v_n_h =  rhs;
       v_n_h -= d_n;
       v_n_h *= (gamma/(beta*dt));
       v_n_h.linAdd( (1.0-(gamma*(1.0-alphaf)/beta)), v_n,
                     (dt*0.5*(1.0-alphaf)*(2.0*beta-gamma)/beta), a_n);

       // v_n_p = 1/(1-alphaf)*(v_n_h - alphaf*v_n) = v_(n+1)
       v_n_p.linC(v_n_h,(-1.0*alphaf),v_n);
       v_n_p *= (1.0/(1-alphaf));
     
       // Now swap v_n_p -> v_n and d_n_p -> d_n
       v_p = v_n;
       d_n.swap( d_n_p );
       a_n.swap( a_n_p );
     }
     v_n.swap( v_n_p );

     // Increment time index
     n++;

     // ... current state is replaced by predicted value 
     // NOTE: current state is modified here only for output
     // it will be restored to the backup state at the
     // beginning of the next iteration
     if(!parity && aeroAlg == 5) {
       curState.getDisp().linC(0.5, curState.getDisp(), 0.5, bkState->getDisp());
       curState.getVeloc().linC(0.5, curState.getVeloc(), 0.5, bkState->getVeloc());
     }
     else {
       // FORCE PRINTING AT LAST ITERATION
       if(t+1.01*dt > tmax)  probDesc->processLastOutput();
     }

     postProcessor->dynamOutput(n, dynOps, ext_f, aeroForce, curState);

     // ... For A5 Algorithm, do one time step back if necessary
     // add n so that time index is wound back as well as t
     if(aeroAlg == 5) probDesc->a5TimeLoopCheck( parity, t, dt );

   }
   if(aeroAlg < 0)
     filePrint(stderr,"\r ... Time Integration Loop: t = %9.3e, 100%% complete ...\n", t);

   totalTime += getTime();
#ifdef PRINT_TIMERS
   if(verboseFlag) filePrint(stderr," ... Total Loop Time = %.2e s   ...\n",totalTime/1000.0);
#endif
}

// -----------------------------------------------------------------------------
//
// ... CENTRAL DIFFERENCE TIME LOOP
//
// DESCRIPTION: this routine implements the explicit central difference time-integrator
//              for the solution of the second-order differential equation:
//              (a) LINEAR mech or acou: M*(d^2u/dt^2) + C*(du/dt) + K*u = fext(u)
//              (b) NONLINEAR mech: M*(d^2u/dt^2) + C*(du/dt) + fint(u) = fext(u)
//
// WARNINGS:    1. The mass matrix is automatically LUMPED (in main.C) unless MRATIO is set to 1
//              2. Viscous damping is supported, but to keep the scheme explicit the equilibrium 
//                 condition is expressed as M*a^{n+1} + C*v^{n+1/2} + K*u^{n+1} = fext^{n+1}
//                 where v^{n+1/2} = v^n + dt/2*a^n
//              3. Contact/tied surfaces are supported (using ACME) but LMPCs are not supported
//              4. Velocity and/or acceleration controls (ACTUATORS) are not strictly correct since we
//                 use v^n and a^n to compute fext^{n+1}
//              5. USDD needs to be checked
//
// -----------------------------------------------------------------------------
 
template< class DynOps,        class VecType, 
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::explicitNewmarkLoop(SysState<VecType>& curState, VecType& constForce, 
                      DynOps& dynOps, 
                      NewmarkWorkVec<VecType,ProblemDescriptor>& workVec,
                      double dt, double tmax)
{
  // XXXX this hasn't been generalized for first order yet
  filePrint(stderr, " ... Explicit Newmark Time Integration Scheme: beta = %4.2f, gamma = %4.2f, alphaf = %4.2f, alpham = %4.2f ...\n",0.0,0.5,0.0,0.0);

  int parity = 0;
  SysState<VecType> *bkState = 0;
  // Allocate backup state for A5 algorithm
  if(aeroAlg == 5) {
    VecType d_bk(probDesc->solVecInfo());
    VecType v_bk(probDesc->solVecInfo());
    VecType a_bk(probDesc->solVecInfo());
    VecType v_p_bk(probDesc->solVecInfo());
    bkState = new SysState<VecType>(d_bk, v_bk, a_bk, v_p_bk);
  }

  // Vectors we will need to use
  VecType &v_n = curState.getVeloc();
  VecType &v_p = curState.getPrevVeloc();
  VecType &d_n = curState.getDisp();
  VecType &a_n = curState.getAccel();
  VecType &v_n_h = workVec.get_v_n_h();
  VecType &fext = workVec.get_ext_f();
  VecType &fint = workVec.get_fint();
  VecType &tmp1 = workVec.get_tmp1();
  VecType &tmp2 = workVec.get_tmp2();
  VecType v_h_p(probDesc->solVecInfo());
  v_h_p = 0.0;

  // project initial displacements in case of rbmfilter
  if(probDesc->getFilterFlag() > 0) {
    probDesc->project(d_n);
    probDesc->project(v_n);
  }

  // Initialize time index (n) and time (t^n)
  int n = 0;
  double t = 0.0;
  probDesc->getInitialTime(n, t);

  // Get initial external force vector fext^0
  probDesc->computeExtForce2(curState, fext, constForce, n, t, aeroForce, 0.5, 0.0);

  // Compute the initial internal forces fint^0
  if(domain->solInfo().isNonLin()) probDesc->getInternalForce(d_n, fint, t); else // PJSA 3-31-08
  dynOps.K->mult(d_n, fint);

  // Compute the initial acceleration a^0 = M^{-1}(fext^0 - fint^0 - C*v^0)
  // XXXX this should be the initial velocity for first order problems
  if(verboseFlag) filePrint(stderr," ... Computing initial acceleration ...\n");
  if(dynOps.C) {
    dynOps.C->mult(v_n,tmp2);
    fint.linC(fint,1.0,tmp2);
  }
  a_n.linC(1.0, fext, -1.0, fint);
  dynOps.dynMat->reSolve(a_n);
  if(domain->tdenforceFlag() || domain->solInfo().penalty) { // Contact corrector step: a^0 += M^{-1}*Fctc
    tmp1.linC(dt, v_n, 0.5*dt*dt, a_n); tmp1 += d_n; // predicted displacement d^1 = d^0 + dt*v^0 + dt*dt/2*a^0
    probDesc->getContactForce(tmp1, tmp2);
    dynOps.dynMat->reSolve(tmp2);
    a_n += tmp2;
  }
  if(probDesc->getFilterFlag() == 2) probDesc->project(a_n);

  // Compute velocity at first half-time station: v^{1/2} = v^0 + dt/2*a^0
  v_n_h.linC(1.0, v_n, 0.5*dt, a_n);

  // Output the state at t^0: d^0, v^0, a^0, fext^0
  postProcessor->dynamOutput(n, dynOps, fext, aeroForce, curState);

  // for using the modified wave equation in the acoustic problemType
  double coeff = 0.0;
  double coeff2 = -dt*dt/(probDesc->getDomain()->solInfo().modifiedWaveEquationCoef);
  bool fourthOrder = probDesc->getDomain()->solInfo().modifiedWaveEquation;

  // ... BEGIN MAIN TIME-LOOP
  double totalTime = -getTime();
  char ch[4] = { '|', '/', '-', '\\' };

  //double t1 = 0, t2 = 0, t3 = 0, t4 = 0;

  for( ; t < tmax-0.01*dt; t += dt) {

    if(aeroAlg < 0) {
      filePrint(stderr,"\r  %c  Time Integration Loop: t = %9.3e, %3d%% complete ",
                ch[int((totalTime + getTime())/250.)%4], t, int(t/(tmax-0.01*dt)*100));
    }
    //t1 -= getTime(); 

    if (fourthOrder) {
      // this is as in the previous release of the FEM code (before august 28th 2008)
      //  the variables have been change to match the one declared above.
      // the alpha rayleigh damping is probably neglectedó
      //
      // External force
      probDesc->computeExtForce2(curState, fext, constForce, n+1, t+dt, aeroForce, 0.5, 0.0);
      d_n.linAdd(dt,v_n_h);
      // Internal force
      dynOps.K->mult(d_n, fint);
      // 4th order loop see Cohen et al, Finite Elements in Analysis and Design 16(1994) pp 329-336
      //"Higher-order finite elems with mass lumping for teh 1D wave equation"
      dynOps.dynMat->reSolve(fint);//ne peut pas aller a l'order 6 ... dynMat contient dt*C/2 !
      tmp1.linC(1.0,d_n,coeff2,fint);
      dynOps.K->mult(tmp1,fint);
      if (dynOps.C) {
        tmp1.linC(1.0,fext,-1.0,fint);
        dynOps.C->mult(v_n_h,tmp2);
        tmp1.linC(2.0,tmp1,-1.0,tmp2);
        dynOps.M->mult(v_n_h,tmp2);
        coeff=0.5*dt;
        //v_p=v_n_h;//temporary
        v_n_h.linC(1.0,tmp2,coeff,tmp1);
        dynOps.dynMat->reSolve(v_n_h);
        //coeff=1/dt;
        //a_n.linC(coeff,v_n_h,-coeff,v_p);
      } else {
        coeff=dt;
        a_n.linC(1.0,fext,-1.0,fint);
        dynOps.dynMat->reSolve(a_n);
        v_n_h.linC(1.0,v_n_h,coeff,a_n);
      }
      //coeff=0.5*dt;
      //v_p=v_n;
      //v_n.linC(1.0,v_n_h,coeff,a_n);
      
      // Update the time index
      n += 1;
      // Output the state at t^n: d^n, v^n, a^n, fext^n
      postProcessor->dynamOutput(n, dynOps, fext, aeroForce, curState);

    }
    else {
      if(aeroAlg == 5) probDesc->a5StatusRevise(parity, curState, *bkState);

      // Mode decomposition of displacement
      if(probDesc->getModeDecompFlag()) probDesc->modeDecomp(t, n, d_n);

      // Update the displacement at t^(n+1): d^{n+1} = d^n + dt*v^{n+1/2}
      d_n.linAdd(dt, v_n_h);

      // C0: Send predicted displacement at t^{n+1.5} to fluid
      //if(aeroAlg == 20) probDesc->aeroSend(t+dt, d_n, v_n, a_n, v_n_h); // note: v_n, a_n haven't been updated yet!!!
      if(aeroAlg == 20) probDesc->aeroSend(t+dt, d_n, v_n_h, a_n, v_h_p);

      // Compute the external force at t^{n+1}
      //t2 -= getTime();
      probDesc->computeExtForce2(curState, fext, constForce, n+1, t+dt, aeroForce, 0.5, 0.0);
      //t2 += getTime();

      // Compute the internal force at t^{n+1}
      //t3 -= getTime();
      if(domain->solInfo().isNonLin()) probDesc->getInternalForce(d_n,fint,t+dt);
      else {
        dynOps.K->mult(d_n, fint);
      }
      //t3 += getTime();

      // Compute the acceleration at t^{n+1}: a^{n+1} = M^{-1}(fext^{n+1}-fint^{n+1}-C*v^{n+1/2})
      if(dynOps.C) {
         dynOps.C->mult(v_n_h, tmp1);
         fint.linAdd(1.0, tmp1);
      }
      a_n.linC(1.0, fext, -1.0, fint);
      //t4 -= getTime();
      dynOps.dynMat->reSolve(a_n);
      //t4 += getTime();
      if(domain->tdenforceFlag() || domain->solInfo().penalty) { // Contact corrector step
        tmp1.linC(dt, v_n_h, dt*dt, a_n); tmp1 += d_n; // predicted displacement d^{n+2} = d^{n+1} + dt*(v^{n+1/2} + dt*a^{n+1})
        probDesc->getContactForce(tmp1, tmp2);
        dynOps.dynMat->reSolve(tmp2);
        a_n += tmp2;
      }
      if(probDesc->getFilterFlag() == 2) probDesc->project(a_n);

      // Update the velocity at t^{n+1}: v^{n+1} = v^{n+1/2}+dt/2*a^n
      v_p = v_n;
      v_n.linC(1.0, v_n_h, 0.5*dt, a_n);

      // Update the time index
      n += 1;

      // ... current state is replaced by predicted value
      // NOTE: current state is modified here only for output
      // it will be restored to the backup state at the
      // beginning of the next iteration
      if(!parity && probDesc->getAeroAlg() == 5) {
        curState.getDisp().linC(0.5, curState.getDisp(), 0.5, bkState->getDisp());
        curState.getVeloc().linC(0.5, curState.getVeloc(), 0.5, bkState->getVeloc());
      }
      else {
        if(t+1.01*dt > tmax) probDesc->processLastOutput(); // force printing at last iteration
      }

      // Output the state at t^n: d^n, v^n, a^n, fext^n
      postProcessor->dynamOutput(n, dynOps, fext, aeroForce, curState);

      // Compute midpoint velocity: v^{n+1/2} = v^{n-1/2} + dt*a^n
      v_h_p = v_n_h;
      v_n_h.linAdd(dt, a_n);
  
      // ... For A5 Algorithm, do one time step back if necessary
      // add n so that time index is wound back as well as t
      if(aeroAlg == 5) probDesc->a5TimeLoopCheck(parity, t, dt);
    } 
    //t1 += getTime();
    //filePrint(stderr,"\r  %c  Time Integration Loop: t = %9.3e, %3d%% complete %e %e %e %e",
    //          ch[int((totalTime + getTime())/250.)%4], t, int(t/(tmax-0.01*dt)*100), t1/1000/n, t2/1000/n, t3/1000/n, t4/1000/n);
  }
  if(aeroAlg < 0)
    filePrint(stderr,"\r ... Time Integration Loop: t = %9.3e, 100%% complete ...\n", t);

  totalTime += getTime();
#ifdef PRINT_TIMERS
  if(verboseFlag) filePrint(stderr," ... Total Loop Time = %.2e s   ...\n",totalTime/1000.0);
#endif
}

//------------------------------------------------------------------------------

template<
     class DynOps,             // Data Structure for K, C, M and dynMat
     class VecType,            // Vector type used
     class PostProcessor,      
     class ProblemDescriptor,
     class Scalar> 

int
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar >
::checkSteadyState(double time, double step, double criteria)
{
  int nstep = (int)(time / step);
  
  filePrint(stderr," ... Pseudo-Step = %d  Rel. Res. = %10.4e\n",nstep,criteria);

  //if ( nstep < steadyMin ) return 0;
  if(nstep > steadyMax) return 2;

  return (criteria < steadyTol) ? 1 : 0;
}
