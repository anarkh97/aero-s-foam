#ifndef _TIMESLICE_C_
#define _TIMESLICE_C_

#include <Pita.d/TimeSlice.h>
#include <Pita.d/DistrTimeDecompSolver.h>
#include <Driver.d/DynamProbType.h>
#include <Pita.d/PitaTimers.h>

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::LinTimeSlice()
{
   ptdynOps = 0;
   probDesc = 0;
   postProcessor      = 0;
   ptTimeDecompSolver = 0;

   numITA = 0;

   convergence = false;
   
   ptaeroForce  = 0;
   ptconstForce = 0;
   ptworkVec    = 0;
   ptcurState   = 0;

   d = 0; v = 0; a = 0; vp = 0;
}

//------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::print()
{
 cout << "Time Slice no: " << SliceRank << endl;
 cout << "Ti_initial: " <<Ti_initial<< ", Ti_final: " << Ti_final << endl;
}

//------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::
    setProblemParam(DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> *_pt,
                        ProblemDescriptor *_probDesc, PostProcessor *_postProcessor)
{
   probDesc      = _probDesc;
   postProcessor = _postProcessor;

   ptTimeDecompSolver = _pt;
}

//-------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::
     setTimeParam(double _Ti_initial, double _Ti_final, double _dt, int _tIndex, int _numTimeSteps)
{
 
  Ti_initial = _Ti_initial;
  Ti_final   = _Ti_final;
  dt = _dt;

  TimeIndex = _tIndex;
  numSteps  = _numTimeSteps;
}

//-------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::
     setParam(int _SliceRank, int _locRank, int _numTSinmyCPU)
{
   SliceRank    = _SliceRank;
   locRank      = _locRank;
   numTSinmyCPU = _numTSinmyCPU;
}

//----------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::~LinTimeSlice()
{
    if (d) { delete d; d = 0 ; } if (v)  { delete v; v = 0; }
    if (a) { delete a; a = 0 ; } if (vp) { delete vp; vp = 0; }
                                                                                                                                                                       
    if (ptcurState)    { delete ptcurState;   ptcurState   = 0; }
    if (ptaeroForce)   { delete ptaeroForce;  ptaeroForce  = 0; }
    if (ptconstForce)  { delete ptconstForce; ptconstForce = 0; }
    if (ptworkVec)     { delete ptworkVec;    ptworkVec    = 0; }
}

//-----------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::changeParam(bool ck, int numsave,double dtsave)
{

 if (ck) {
   this->dt *= this->numSteps; 
   this->numSteps = 1;
 } else {
   this->dt = dtsave;
   this->numSteps = numsave;
 }

}

//-----------------------------------------------------------------

template<
    class DynOps, 
    class VecType, 
    class PostProcessor, 
    class ProblemDescriptor,
    class InfoSize> 
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::initialization(VecType *_ptconstForce, 
     VecType *_aeroForce,double _alphaf, double _alpham, double _beta, double _gamma, DynOps *_ptdynOps)
{  
    d = new VecType( this->probDesc->solVecInfo() );
    v = new VecType( this->probDesc->solVecInfo() );
    a = new VecType( this->probDesc->solVecInfo() );
    vp = new VecType( this->probDesc->solVecInfo() );
    ptcurState = new SysState<VecType>(*d,*v,*a,*vp);
 
    ptworkVec  = new NewmarkWorkVec<VecType,ProblemDescriptor> (1, this->probDesc);
    
    ptconstForce  = new VecType(this->probDesc->solVecInfo());
    *ptconstForce = *_ptconstForce;

    if(_aeroForce) {
       ptaeroForce  = new VecType(this->probDesc->solVecInfo()) ; 
       *ptaeroForce = *_aeroForce;
    }

    this->alphaf = _alphaf;
    this->alpham = _alpham;
    this->beta   = _beta;
    this->gamma  = _gamma;

    this->ptdynOps = _ptdynOps;
}

//------------------------------------------------------------------

template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::getInitCond(bool step0, 
           bool ck, bool bi, int myCPU, int cycleTSid_0, StateSet<VecType,DynOps,InfoSize> *tmp)
{

     // in a PITA mode
     VecType &tmp1  = ptworkVec->get_tmp1();
     VecType &dnc   = ptworkVec->get_dnc();
     VecType &vnc   = ptworkVec->get_vnc();

     VecType &d_n = ptcurState->getDisp();
     VecType &v_n = ptcurState->getVeloc();
     VecType &v_p = ptcurState->getPrevVeloc();      

     v_p.zero();                           // not for aero

     tmp1.zero();
     dnc.zero(); vnc.zero();              
                                         
     //d_n & v_n
     if (bi) {
       d_n.zero(); v_n.zero();
     }
     if (step0) {
       *ptcurState = *(this->ptTimeDecompSolver->getptCurState());
     }
     if (ck){ 
       d_n=tmp->getDisp(0);
       v_n=tmp->getVel(0); 
     } 
     if (!bi && !step0 && !ck) {
       this->ptTimeDecompSolver->get_Seed(d_n,v_n,this->SliceRank-cycleTSid_0);
     }

}

//-------------------------------------------------------------------------
template<
    class DynOps, 
    class VecType, 
    class PostProcessor, 
    class ProblemDescriptor,
    class InfoSize> 
void LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::
     solve_ITA_linearDynam(bool bi, bool step0, bool ck, int myCPU, 
     StateSet<VecType,DynOps,InfoSize> *tmp, int cycleTSid_0)
{

    /*!!!!!!!! we are in the case where alpham = 0.5, beta = 0.25, gamma = 0.5 !!!!!!!!
      ( (1-alpham)/2 - beta ) == 0 and ( 2beta - gamma ) == 0
      an has no influence on the computation of dn and vn. 
      So we don't need to calculate an or anp
    */
    
    /*This method is used to solve ITA in each case: 
      -resolution on the coarse grid to compute Ck and Step0 (ck or step0) 
      -resolution on the fine grid to compute Bi             (bi)
      -resolution on the fine grid on each TS and at each    (!ck && !step0 && !bi)
       iteration 
    */
   
    //get Yki
    PitaTimers &Ptimes = this->ptTimeDecompSolver->getTimers();

    //this timer was usefull when we had to compute 
    //an = M-1(F-Dv-Ku) 
    Ptimes.getInitCondtime -= getTime();
    getInitCond(step0, ck, bi, myCPU, cycleTSid_0, tmp); 
    Ptimes.getInitCondtime += getTime();

    //need for computation
    double t = this->Ti_initial;
    int tI = this->TimeIndex;
    
    VecType  &d_n = ptcurState->getDisp();
    VecType  &v_n = ptcurState->getVeloc();
    VecType  &v_p = ptcurState->getPrevVeloc();

    VecType  &d_n_p = ptworkVec->get_d_n_p();
    VecType  &v_n_p = ptworkVec->get_v_n_p();
    VecType    &rhs = ptworkVec->get_rhs();
    VecType  &ext_f = ptworkVec->get_ext_f();
    VecType  &d_n_h = ptworkVec->get_d_n_h();
    VecType  &v_n_h = ptworkVec->get_v_n_h();
    VecType &Md_n_h = ptworkVec->get_Md_n_h();
    VecType &Cd_n_h = ptworkVec->get_Cd_n_h();

    VecType &tmp1 = ptworkVec->get_tmp1();
    //VecType &dnc  = ptworkVec->get_dnc();
    //VecType &vnc  = ptworkVec->get_vnc();
    //VecType &anc  = ptworkVec->get_anc();
    
    double *pt_dt = &(this->dt);
    double *pt_t = &(t);

    //TStime == time of resolution without computation added
    //by Pita. To check if TStime is equivalent to the time
    //of resolution of the sequential method on a same time 
    //interval

    if (this->SliceRank==0 && !bi && !ck && !step0 && this->numITA == 0) {
       Ptimes.TStime =- getTime();
    }

    bool performOutput = !ck && !step0 && !bi;

    if (performOutput) {
      this->postProcessor->openOutputFilesForPita(this->SliceRank);
      this->postProcessor->pitaDynamOutput(tI,*(this->ptdynOps),ext_f,ptaeroForce,*ptcurState,this->SliceRank,*pt_t);
    }
    
    for (int i=0; i < this->numSteps; i++) {

      // Mode decomposition of displacement
      if (this->probDesc->getModeDecompFlag() ) {
        this->probDesc->modeDecomp(t,tI,d_n); 
      }

       // ... For Aeroelastic A5 Algorithm, Do restore and backup here
       //if(aeroAlg == 5) this->probDesc->a5StatusRevise(parity, *ptcurState, *bkState);

        // ... Put initial seed  
        // ... Construct force vector, includes time-independent constant force
        // ... Compute external force at time t+dt*(1-alphaf)
        // ... time modification if ( dt != sinfo.dt ) 

       if (this->SliceRank==0 && !bi && !ck && !step0 && this->numITA==0) Ptimes.TStime += getTime();
       Ptimes.StoreFineGridTime -= getTime();
       if (step0) tmp->put_end_StateSet(d_n,v_n);
       Ptimes.StoreFineGridTime += getTime();
       if (this->SliceRank==0 && !bi && !ck && !step0 && this->numITA==0) Ptimes.TStime -= getTime();
            
       //if (!ck) probDesc->computeExtForce(*ptcurState, ext_f,*ptconstForce,tI,t+(dt*(1-alphaf)),ptaeroForce,gamma,alphaf,pt_dt);
       if (!ck) this->probDesc->computeExtForce2(*ptcurState, ext_f,*ptconstForce,tI,t+(this->dt*(1-this->alphaf)),ptaeroForce,this->gamma,this->alphaf,pt_dt); //HB: due to changes from templating DynamPbDescr ...
       else ext_f.zero();

       // ... Construct R.H.S. vector
       // ... d_n_h = ((1-alpham)/(1-alphaf))*d_n + dt*(1-alpham)*v_n + dt*dt*((1-alpham)/2 - beta)*a_n
       // ... d_n_h = ((1-alpham)/(1-alphaf))*d_n + dt*(1-alpham)*v_n  
       
       // First: d_n_h = ((1-alpham)/(1-alphaf))*d_n + dt*(1-alpham)*v_n 
       d_n_h.linC( ((1.0-this->alpham)/(1.0-this->alphaf)), d_n, this->dt*(1.0-this->alpham), v_n );

       // Second: d_n_h = d_n_h + dt*dt*((1-alpham)/2 - beta)*a_n
       //d_n_h.linC( d_n_h, dt*dt*(0.5*(1.0-alpham)-beta), a_n );

       // Third: Multiply by Mass Matrix M
       (*(this->ptdynOps)).M->mult( d_n_h, Md_n_h );

       // Accumulate in rhs vector: rhs = Md_n_h + beta*dt*dt*ext_f
       rhs.linC( Md_n_h, this->beta*this->dt*this->dt, ext_f );
                                                                                     
       if( (*(this->ptdynOps)).C ) {
          // ... d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma(1-alphaf))*v_n
          // ...       - dt*dt*dt*0.5(1-alphaf)*(2*beta - gamma)*a_n

          // ... d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma(1-alphaf))*v_n
          d_n_h.linC( this->dt*this->gamma, d_n, this->dt*this->dt*(this->gamma*(1.0-this->alphaf)-this->beta), v_n );

          // ... d_n_h = d_n_h - dt*dt*dt*0.5(1-alphaf)*(2*beta - gamma)*a_n
          //d_n_h.linC( d_n_h, dt*dt*dt*0.5*(alphaf-1.0)*(2.0*beta - gamma), a_n );

          // Multiply by Damping Matrix C
          (*(this->ptdynOps)).C->mult( d_n_h, Cd_n_h );
          
          // Accumulate in rhs vector 
          rhs += Cd_n_h;
       }
       // add prescribed boundary condition contributions to rhs vector
       // and update accelerations at prescribed dofs
      /* if( !step0 && !ck){                                //time modification if ( dt != sinfo.dt )
         probDesc->addPrescContrib( ptdynOps.M12, ptdynOps.C12, dnc, vnc, anc, tmp1, t);
       }else{
         probDesc->addPrescContrib( ptdynOps.M12, ptdynOps.C12, dnc, vnc, anc, tmp1, t, pt_dt);
       }
       */
       rhs += tmp1;
       tmp1.zero();

       if (ck || step0) (*(this->ptdynOps)).dynMat_Dt->reSolve( rhs ); // Now rhs contains (d_n+1-alphaf)
       else (*(this->ptdynOps)).dynMat->reSolve( rhs );

       //   call projector for RBMs in case of rbmfilter level 2
       if (this->probDesc->getFilterFlag() == 2) this->probDesc->project( rhs );

       // one time step forward
       d_n_p.linC(rhs,(-1.0*this->alphaf),d_n);
       d_n_p *= (1.0/(1-this->alphaf));

       v_n_h =  rhs;
       v_n_h -= d_n;                                  // v_n_h currently is d_(n+1-alphaf) - d_n
       v_n_h *= (this->gamma/(this->beta*this->dt));
       v_n_h.linAdd( (1.0-(this->gamma*(1.0-this->alphaf)/this->beta)), v_n );

       v_n_p.linC(v_n_h,(-1.0*this->alphaf),v_n);
       v_n_p *= (1.0/(1-this->alphaf));
   
       // Now swap v_n_p -> v_n and d_n_p -> d_n
       v_p = v_n;
       v_n.swap( v_n_p );
       d_n.swap( d_n_p );
       
       if (this->SliceRank==0 && !bi && !ck && !step0 && this->numITA==0) Ptimes.TStime += getTime();
       //store data if it's the last loop
       //this timer is used to check that the put_end_StateSet or
       //replace_end_StateSet method don't take a long time
       Ptimes.StoreFineGridTime -= getTime();
       if ( i==(this->numSteps-1) && !step0 ){
          if (ck) tmp->replace_end_StateSet(d_n,v_n);
          else    tmp->put_end_StateSet(d_n,v_n); 
       }
       Ptimes.StoreFineGridTime += getTime();
       if (this->SliceRank==0 && !bi && !ck && !step0 && this->numITA==0) Ptimes.TStime -= getTime();

       // ... For A5 Algorithm, do one time step back if necessary
       // add tIndex so that tIndex is wound back as well as t
       //if(aeroAlg == 5) this->probDesc->a5TimeLoopCheck(parity, t, this->dt);

       // Increment time index and time
       tI ++;
       t+=this->dt;

       // Output current solution 
       if (performOutput) {  
          this->postProcessor->pitaDynamOutput(tI,*(this->ptdynOps),ext_f,ptaeroForce,*ptcurState,this->SliceRank,*pt_t);
       }
    
    }

    if (performOutput) {
      this->postProcessor->closeOutputFiles();
    }

    if (this->SliceRank==0 && !bi && !ck && !step0 && this->numITA==0) {
        Ptimes.TStime += getTime();
        cout<<"TS time "<<endl;
    }
}
#endif
