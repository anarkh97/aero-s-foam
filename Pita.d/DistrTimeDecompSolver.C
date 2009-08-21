#ifndef _DISTRTIMEDECOMPSOLVER_C_
#define _DISTRTIMEDECOMPSOLVER_C_

#include <Pita.d/DistrTimeDecompSolver.h>

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::getnumCPU()
{
   numCPU = structCom->numCPUs();
   myCPU  = structCom->myID();
   
   numTotalActiveTS    = numCPU*numTSperCycleperCPU;
   x = numTSperCycleperCPU+1;
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::exchange_InfoCPU()
{
  timeCom->template allGatherv<int>(infoCPU+position_Info[myCPU],numdata_Info[myCPU],infoCPU,numdata_Info,position_Info);
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::exchange_PropValue()
{
   //Exchange the vectors yi(Ti+1). These vectors have to be stored in ptPropState by id's ascending order.
   //During the allGatherv method, datas are stored by MPI Process's origin. After allGatherv, we have to 
   //store datas in a correct order. 

   //To avoid problem due to exchanging a big number of datas, we first send the displacements vectors then
   //the velocity vectors.

   //ptPropState is a StateSet object. Vectors are stored at the end of the stack, at the position 
   //(numVectors-1). If we want to store a vector before the end of the stack, we have to change the cursor
   //ie numVectors. To do that, we use adjustnumVectors.

   int pos1 = cycleTSid[0];
   int pos2 = cycleTSid[numTotalActiveTS-1];
   int pos3, cpu, cursor;
   bool b=true;

   if (active) {
     //stores all the displacement vectors in an arrray
     for (int i=0; i<infoCPU[myCPU*x]; i++){ 
       (ptProptmp->getdispState(i)).putIn(buf,i*probsize,probsize);
     }
   }else{ 
     //if isn't active anymore, the MPI Process sends 0.0
     buf[0]=0.0;
   }
   
   //exchanges displacement 
   timeCom->template allGatherv<double>(buf,numdata[myCPU],ptBuffeur,numdata,position);

   if (active){       
     //if a MPI Process is active, it stores all the TS's value in ptPropState by id's 
     //ascending order starting position pos1 = first active TS. The propoagation of 
     //the inactive TS are kept, that's why the staring position isn't 0.
     ptPropState->adjustnumVectors(pos1);
     for (int i=pos1; i<pos2+1; i++){
       //find the CPU of this TS
       cpu=(*TStoCPU)[i][0];
       //find the position of this TS in its CPU's last cycle's active TS
       cursor=0;
       pos3=0;
       while ( b && cursor<infoCPU[cpu*x] ) {
         cursor++;
         if ( infoCPU[cpu*x+cursor]==i )  b=false; 
       }     
       //put in PropState this TS's displacement. During the exchange, this displacement 
       //was stored in ptBuffeur at the adress:
       //ptBuffeur + beginning of the CPU + (cursor-1)*size
       ptPropState->put_end_StateSet(ptBuffeur+position[cpu]+(cursor-1)*probsize,probsize,*ptZero);
       //as we have only the displacements vectors, we store them with a velocity=0
       b=true;
     }
   }

   if (active) {
      //stores all the velocity vectors in an arrray
      for (int i=0; i<infoCPU[myCPU*x]; i++){
        (ptProptmp->getvelState(i)).putIn(buf,i*probsize,probsize);
      }
   }else{
      buf[0]=0.0;
   }

   //exchanges velocity 
   timeCom->template allGatherv<double>(buf,numdata[myCPU],ptBuffeur,numdata,position);

   if (active){                                           
     //displacements vectors have already been stored
     for (int i=pos1; i<pos2+1; i++){
       ptPropState->adjustnumVectors(i+1);
       cpu=(*TStoCPU)[i][0];
       cursor=0;
       pos3=0;
       while ( b && cursor<infoCPU[cpu*x] ) {
         cursor++;
         if ( infoCPU[cpu*x+cursor]==i )  b=false;
       }
       ptPropState->add_end_StateSet(*ptZero,ptBuffeur+position[cpu]+(cursor-1)*probsize,probsize);
       b=true;
     }
   }
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::exchange_BiValue()
{
   //this method is used only one time, that's why arrays are created and deleted here
   //same kind of exchanging technicals than exchange_PropValue()
   int cursor=0;

   int *numdata_Bi  = new int [numCPU];
   int *position_Bi = new int [numCPU];

   for (int i=0; i<numCPU; i++){
     numdata_Bi[i]  = CPUtoTS->num(i)*probsize;
     position_Bi[i] = cursor;
     cursor += numdata_Bi[i];
   }
  
   double *buf_Bi       = new double[numdata_Bi[myCPU]];
   double *buf_Bi2      = new double[numdata_Bi[myCPU]];
   double *ptBuffeur_Bi = new double[nTS*probsize];
                                                                               
   //exchanges displacement and velocity
   for (int i=0; i<CPUtoTS->num(myCPU); i++){
     (ptBfine->getdispState(i)).putIn(buf_Bi,i*probsize,probsize);
     (ptBfine->getvelState(i)).putIn(buf_Bi2,i*probsize,probsize);
   }
   
   timeCom->template allGatherv<double>(buf_Bi,numdata_Bi[myCPU],ptBuffeur_Bi,numdata_Bi,position_Bi);
 
   ptBfine->adjustnumVectors(0);
   for (int i=0; i<numCPU; i++){
     for (int j=0; j<CPUtoTS->num(i); j++){
       ptBfine->put_end_StateSet(ptBuffeur_Bi+position_Bi[i]+j*probsize,probsize,*ptZero);
     }
   }

   timeCom->template allGatherv<double>(buf_Bi2,numdata_Bi[myCPU],ptBuffeur_Bi,numdata_Bi,position_Bi);
  
   cursor=1;
   for (int i=0; i<numCPU; i++){
     for (int j=0; j<CPUtoTS->num(i); j++){
       ptBfine->adjustnumVectors(cursor);
       ptBfine->add_end_StateSet(*ptZero,ptBuffeur_Bi+position_Bi[i]+j*probsize,probsize);
       cursor++;
     }
   }

   delete[] numdata_Bi;
   delete[] position_Bi;
   delete[] buf_Bi;
   delete[] buf_Bi2;
   delete[] ptBuffeur_Bi;
    
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::buildTSandConnect()
{

   //Tfinal modification so that Tfinal = k*Dt and nTS computation
   double tmp = Tinitial;
   int count = 0;
   if (myCPU==0) cout<<" Tfinal "<<Tfinal<<endl;
   while(tmp<Tfinal){
    tmp+=Dt;
    count++;
   }
 
   Tfinal = tmp;
   nTS    = count;
   
   if (myCPU==0) cout<<"nTS "<<nTS<<" Tfinal "<<Tfinal<<" Dt "<<Dt<<endl;                                                           
    
   //If the seed Yi0 are computed on a coarser mesh, checks that the data are
   //correctly read from the input file
   if (nvStep0) {
    int numDofPitaDisp6 = geoSource->getNumPitaIDis6();
    int numDofPitaVel6  = geoSource->getNumPitaIVel6();
    int numTSPitaDisp6  = geoSource->getNumTSPitaIDis6();
    int numTSPitaVel6   = geoSource->getNumTSPitaIVel6();

    if ( numDofPitaDisp6 != numDofPitaVel6 || numTSPitaDisp6 != nTS || numTSPitaVel6 !=nTS){
      fprintf(stderr," ... Error in PitaDisp6 or PitaVel6 ... \n");
      cout<<"numDofPitaDisp6 "<<numDofPitaDisp6<<", numDofPitaVel6 "<<numDofPitaVel6<<endl;
      cout<<"numTSPitaDisp6 "<<numTSPitaDisp6<<" , numTSPitaVel6 "<<numTSPitaVel6<<", nTS "<<nTS<<endl;
      exit(-1);
    }
   }
   
   //exit if we are using too many CPUs 
   if (numCPU*numTSperCycleperCPU>nTS) {
     fprintf(stderr," *** ERROR: number of CPU > number of TS ***\n");
     exit(-1);
   }
   active = true;
 
   //Connectivity and duplicate outputfiles to get one per TS 
   getConnect();
   numTS = CPUtoTS->num(myCPU);
   geoSource->duplicateFilesForPita(CPUtoTS->num(myCPU), CPUtoTS->operator[](myCPU));

   //initialization of array
   arrayInitialization();

   //build TS
   double Ti_initial,Ti_final;
   int tIndex, id, nTimeSteps, locRank, numTSinmyCPU;

   arrayTimeSlices = new LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> [numTS];
  
   nTimeSteps = Jratio;
   for (int i=0; i<numTS; i++){

     id = (*CPUtoTS)[myCPU][i];
     tIndex = id * Jratio;
     Ti_initial = Tinitial + id * Dt;
     Ti_final = Ti_initial + Dt;
      
     numTSinmyCPU = CPUtoTS->num(myCPU); 
     locRank = getlocRank(id);

     arrayTimeSlices[i].setProblemParam(this,probDesc,postProcessor);    
     arrayTimeSlices[i].setTimeParam(Ti_initial,Ti_final,dt,tIndex,nTimeSteps);
     arrayTimeSlices[i].setParam(id,locRank,numTSinmyCPU);
   }

   nTimeSteps = nTS; 
   ptCoarseGrid = new LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>;
   ptCoarseGrid->setProblemParam(this,probDesc,postProcessor);
   ptCoarseGrid->setTimeParam(Tinitial,Tfinal,Dt,InitTimeIndex,nTimeSteps);
   ptCoarseGrid->setParam(0,0,1);

} 

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::getConnect()
{
 
   //TStoCPU and CPUtoTS
   int *pt1 = new int[nTS+1];
   int *target1 = new int[nTS];    
   
   //TS->CPU 
   //TS are distributed on the differents MPI Process by group of size numTSperCycleperCPU 
   //the numTSperCycleperCPU first TS     -> CPU0
   //the numTSperCycleperCPU nest ones TS -> CPU1

   int cursor1=0;
   int cursor2=0;
   for(int i=0; i<nTS; i++){
     pt1[i]=i;
     target1[i]=cursor1%numCPU;
     cursor2++;
     if (cursor2==numTSperCycleperCPU) {
       cursor1++;
       cursor2=0;
     }
   }
   pt1[nTS]=nTS;

   TStoCPU = new Connectivity (nTS,pt1,target1);
   CPUtoTS=TStoCPU->reverse();
   CPUtoTS->sortTargets();

   if (myCPU==0) printConnect(CPUtoTS);

}

//-------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::vector_Initialize(bool first)
{
   // d_n,v_n,a_n,v_p put to 0 is !step0
   // project d_n,v_n
   // tmp1,anc,dnc put to 0

   if (!first) { //d_n,v_n,a_n,vp put to 0 + extract user supplied initial disp and vel
     if (!nvStep0) probDesc->getInitState(*ptcurState);   
     else probDesc->getPitaState(*ptcurState,0);
   }

   // project initial displacements in case of rbmfilter
   if (probDesc->getFilterFlag() > 0) {
      probDesc->project(ptcurState->getDisp());
      probDesc->project(ptcurState->getVeloc());
   }

   // New vector added for prescribed boundary conditions
   ptworkVec->get_tmp1().zero();

   // These vectors have length equal to the number of constraints
   // ( number of dirichlet boundary conditions )
   ptworkVec->get_dnc().zero();
   ptworkVec->get_vnc().zero();
   ptworkVec->get_anc().zero();
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::normRefcomputation()
{
   //normDisp and normVel ==  max of the jump
   double normDisp, normVeloc;
   int stop = 0;

   while (stop<ptJumpState->getnumVectors()){
      normDisp  = ptJumpState->getDisp(stop).norm();
      normVeloc = ptJumpState->getVel(stop).norm();
      if (normDisp>normDisp_ref && normVeloc>normVeloc_ref) {
        normDisp_ref  = normDisp;
        normVeloc_ref = normVeloc;
      }
      stop++;
   }

   //if (myCPU==0) cout << "!!! normRef !!! " << normDisp_ref << " " << normVeloc_ref << endl;
}

//-----------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::newInitCond(int cycle)
{
   //Yi = yi-1(Ti) + Cki or Yi = Yi0

   int firstTS = cycleTSid[0];
   int lastTS  = cycleTSid[numTotalActiveTS-1];

   ptSeedState->adjustnumVectors();

   for (int i=firstTS; i<=lastTS; i++){
     if (cycle>0 && i>0 && TSflag[i-1]){ 
       ptSeedState->put_end_StateSet(*ptPropState,i-1);
       ptSeedState->add_end_StateSet(*ptCk,i);
     }else{
       ptSeedState->put_end_StateSet(*ptStep0,i);
     }
   } 
  
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::jumpEvaluation()
{
  int imin = cycleTSid[0];
  int pos;

  ptJumpState->adjustnumVectors(imin);

  //TSflage[i]=true to remember this TS has already been activated.
  //This, to be able to use Jump(i) and Ck(i) to initialize TSi+1.

  if (imin==0) { 
     ptJumpState->put_end_StateSet(*ptZero,*ptZero); 
     TSflag[0]=true;
     imin++; 
  }

  //Jump(i) = yi-1(Ti) - Yi
  for (int i=imin; i<=cycleTSid[numTotalActiveTS-1]; i++){
    if (!TSflag[i]) TSflag[i]=true;
    pos = i-1;
    ptJumpState->put_end_StateSet(*ptPropState,pos);
    pos = i-cycleTSid[0];
    ptJumpState->sub_end_StateSet(*ptSeedState,pos);
  }
  
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::stopIfJumpSmall()
{
  bool check = true;  
  int cursor = 0; double normVeloc = 0;
  int stop   = 0; double normDisp  = 0; 

  //firstTS = lowest TS of the cycle
  //TSid    = lowest TS of the MPI Process
  int firstTS = cycleTSid[0]; 
  int TSid  = arrayTimeSlices[locidTS].getSliceRank();
   
  //check if the MPI Process can stop ie if its TS have converged. 
  //check as long as cursor=number of TS checked < number of 
  //active TS = infoCPU[myCPU*x] and stop as soon as one of its 
  //TS hasn't converged.
  while (check==true && cursor<infoCPU[myCPU*x]) {
    stop = 0;
    for (int i=firstTS; i<=TSid; i++){
      normDisp  = (ptJumpState->getDisp(i)).norm();
      normVeloc = (ptJumpState->getVel(i)).norm();
      if ( (normDisp<min(atol,rtol*normDisp_ref)) && (normVeloc<min(atol,rtol*normVeloc_ref)) ) stop++;
    }
    if ( stop==TSid-firstTS+1 ){
      //a TS converges if its jumps and all the previous ones obtained during the same cycle are small. 
      //we consider that a TS that isn't active anymore has already converged. Indeed, as these TS aren't 
      //active anymore (because numITA==kiter), we can't improve their jumps which may be not small. The 
      //current active TSs can't convergence towards a better solution, even if we apply one ITA more to 
      //them. So even if the inactive TS haven't really converged, if an actual TS has a small jump, we 
      //considerate it as a convergent TS. 
      arrayTimeSlices[locidTS+cursor].converge();
      filePrint(stderr,"... convergence of TS %d, iter %d ... \n",arrayTimeSlices[locidTS+cursor].getSliceRank(),
                            arrayTimeSlices[locidTS+cursor].getnumITA()-1);
      if (locidTS+cursor==nTS-1){
        //filePrint(stderr," --------------------------------------\n");
        //filePrint(stderr,"      ... Convergence of PITA ...\n") ;
        //filePrint(stderr," --------------------------------------\n");
        active = false;
        numTotalActiveTS=0;
      }
      //this TS converges so can check the next TS
      firstTS=TSid+1;
      TSid=firstTS;
      cursor++;
    }else{ 
      check=false; 
    }
 }

}

//------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::arrayInitialization()
{

   ptBuffeur = new double[numTotalActiveTS*probsize];
   buf       = new double[numTSperCycleperCPU*probsize];

   cycleTSid = new int[numTotalActiveTS];
   infoCPU   = new int[x*numCPU];

   numdata   = new int[numCPU];
   position  = new int[numCPU];
   
   numdata_Info  = new int[numCPU];
   position_Info = new int[numCPU];
   
   TSflag = new bool[nTS];

   //contain the number and id of each CPU's active TS
   //infoCPU[x*i]   = number of active TS
   //infoCPU[x*i+j] = id of the jth active TS of the ith MPI Process

   for (int i=0; i<numCPU; i++){
     infoCPU[x*i]=numTSperCycleperCPU;
     for (int j=1; j<x; j++){
       infoCPU[x*i+j]=(*CPUtoTS)[i][j-1];
     }
   }
    
   //for infoCPU exchange 
   for (int i=0; i<numCPU; i++){
     numdata_Info[i]=x;
     position_Info[i]=x*i;
   }

   //if a TS has already been activated
   for (int i=0; i<nTS; i++) TSflag[i]=false;
}

//------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::getActiveTSidAndBuffeur()
{     
     //use infoCPU to fill cycleTSid, numdata and position 

     //cycleTSid contains id of the active TS

     //numdata[i]  = number of data the ith MPI Process will send to the other ones
     //            = number of TS*size of a vector

     //position[i] = position of the first data send by this MPI Process in the total 
     //              buffeur

     //if a MPI Process isn't active anymore, it will send only one data.
     
     int cursor  = 0;   //number of active TS
     int cursor2 = 0;   //total number of data that will be exchanged

     for (int i=0; i<numCPU; i++){
       if ( infoCPU[x*i]!=0 ) {
          for (int j=1; j<infoCPU[x*i]+1; j++){
            cycleTSid[cursor]=infoCPU[x*i+j];
            cursor++;
          }
          numdata[i]  = infoCPU[x*i]*probsize;
          position[i] = cursor2;
          cursor2 += numdata[i];
       }else{
          numdata[i]  = 1;
          position[i] = cursor2;
          cursor2++;
       }
     }
                                                                                                                                 
     if (numTotalActiveTS>0) stable_sort(cycleTSid,cycleTSid+numTotalActiveTS);
 
}

//-----------------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::myNextActiveTS()
{
  //fills infoCPU ie determines which TS will be active during the next cycle
  for (int i=0; i<numCPU; i++){
    if (i!=myCPU) infoCPU[x*i]=0;
    for (int j=1; j<numTSperCycleperCPU+1; j++)   
      infoCPU[x*i+j]=-1;
  }
                                                                                                                                    
  if (active) {
    bool firstid = true; int cursor=0; int cursor2=0;
    bool cond1, cond2, cond3;
    int newlocid;
    
    /*a TS isn't activated anymore if 
        -the number of ITA applyed == kiter
        -it has already converged
        -the number of ITA applyed == TSid+1

      cursor2 = next cycle's number of TS
      locidTS = local rank of the last cycle's lowest TS = position 
                in arrayTimeSlices of the lowest active TS
      locidTS+cursor  = id of the TS checked < numTS
    */ 
    
    while (cursor2<infoCPU[x*myCPU] && locidTS+cursor<numTS ) {
      cond1 = ( arrayTimeSlices[locidTS+cursor].getnumITA()==kiter );
      cond2 = ( arrayTimeSlices[locidTS+cursor].getnumITA()==arrayTimeSlices[locidTS+cursor].getSliceRank()+1 );
      cond3 = ( arrayTimeSlices[locidTS+cursor].getconvergence()==true );
      //if one of this cond=true, don't activate this TS
      //and check the next one
      if ( cond1 || cond2 || cond3 ){
        cursor++;
      }else{
        //else activate this TS and store its id in infoCPU
        if (firstid) { 
           //if it's the first one stored, we have to update locidTS
           newlocid = locidTS+cursor;  
           firstid  = false; 
        }
        cursor2++;
        infoCPU[x*myCPU+cursor2]=(*CPUtoTS)[myCPU][locidTS+cursor];
        cursor++;
      }
    }
    if (!firstid) locidTS = newlocid;
    infoCPU[x*myCPU]=cursor2;
                                                                                                                                    
    if (cursor2==0){
      //no TS will be activated during the next cycle, MPI Propcess will be inactive
      //and will send only one data during Propagation exchanges.
      active = false;
      delete[] buf;
      buf = new double[1];
    }
  }

}

//-----------------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::Ck_computation(int cycle)
{
   //build the Base
   Ptimes->basetime -= getTime();
   ptSeedState->buildBases(*ptPropState, *ptSk, *ptAmplSk, *ptdynOps, *ptKd_Mv, *ptBfine, numTotalActiveTS, cycleTSid);
   Ptimes->basetime += getTime();
                   
   //ck computation       
   VecType tmpd_FinePropCK(info);   // (GFita)J Ck == propagation on the fine grid
   VecType tmpv_FinePropCK(info);   
     
   VecType tmpd_Ckort(info);        // Ckort == (I-Pk)Ck   
   VecType tmpv_Ckort(info);

   VecType zero(info); zero.zero();

   double ps=0.0;   

   //we compute Ck only for this cycle's Active TS. That's why  
   // i >= cycleTSid[0]. As the Ck computation depends on Jump(i-1) 
   //and Ck(i-1), we can compute Ck for the first TS that hasn't 
   //been activated yet, this to use for this TS better initial   
   //conditions than Step0.

   int imin = cycleTSid[0];
   int imax; 
   if (cycleTSid[numTotalActiveTS-1]<nTS-1) imax=cycleTSid[numTotalActiveTS-1]+1; 
   else imax=cycleTSid[numTotalActiveTS-1]; 

  
   ptCk->adjustnumVectors(imin);
   //C0=0
   if (imin==0)  { ptCk->put_end_StateSet(zero,zero); imin++; }

   for (int i=imin; i<imax+1; i++){
  
     ptCk->put_end_StateSet(*ptCk,i-1);              //Cki = Cki-1 + jump Ti-1
     ptCk->add_end_StateSet(*ptJumpState,i-1);                  

     tmpd_FinePropCK.zero(); tmpv_FinePropCK.zero();
     tmpd_Ckort=ptCk->getDisp(i);
     tmpv_Ckort=ptCk->getVel(i);
                                                                                                
     //(GFita)J PkCk = (GFita)J [Sumj <Ck|Skj>Skj]= Sumj [<Ck|Skj> (GFita)J Skj] = Sumj [<Ck|Skj> AmplSkj]
     Ptimes->CkfineProjtime -= getTime();
     for (int j=0; j<ptSk->getnumVectors(); j++){

       ps = ptCk->getDisp(i)*ptKd_Mv->getDisp(j)+ ptCk->getVel(i)*ptKd_Mv->getVel(j);

       tmpd_FinePropCK.linAdd(ps,ptAmplSk->getDisp(j));
       tmpv_FinePropCK.linAdd(ps,ptAmplSk->getVel(j));

       tmpd_Ckort.linAdd(-1*ps,ptSk->getDisp(j));
       tmpv_Ckort.linAdd(-1*ps,ptSk->getVel(j));
     }
     ptCk->replace_end_StateSet(tmpd_FinePropCK,tmpv_FinePropCK);
     Ptimes->CkfineProjtime += getTime();

     //ITA on the coarse grid applyed to Ckortho ie (I-Pk)Ck
     if (CkCoarse) {
       
       bool ck=true,step0=true,bi=true;
       Ptimes->CkcoarseProjtime -= getTime();

       StateSet<VecType,DynOps,InfoSize> Ctmp(1,info);
       Ctmp.put_end_StateSet(tmpd_Ckort,tmpv_Ckort);
       int numsave=0; double dtsave=0.0; 

       arrayTimeSlices[0].getParam(numsave,dtsave);
       arrayTimeSlices[0].changeParam(ck);
       arrayTimeSlices[0].solve_ITA_linearDynam(!bi,!step0,ck,myCPU,&Ctmp,0);
       arrayTimeSlices[0].changeParam(!ck,numsave,dtsave);

       Ptimes->CkcoarseProjtime += getTime();

       ptCk->add_end_StateSet(Ctmp,0);
     }

   }
}

//-----------------------------------------------------------------------------------------------------

template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::solve_PITA_linearDynam()
{
   // Get myCPU and numCPU == number of CPU groups.
   // Each CPU group solves the problem on a given set of TS. 
   getnumCPU();

   //PitaTotalTime == solver totalTime without PitaprintFile
   Ptimes->PitaTotaltime -= getTime();

   probDesc->preProcess();
   postProcessor=probDesc->getPostProcessor();
 
   info = probDesc->solVecInfo();

   probsize = info.totLen();
   if (myCPU == 0) cout << "probsize " << probsize << endl;

   ptZero = new VecType(info);
   ptZero->zero();

   ptconstForce = new VecType(info); 
   
   //Allocate vectors for displacement, velocity,
   //acceleration and last velocity
   VecType d_n( probDesc->solVecInfo() );
   VecType v_n( probDesc->solVecInfo() );
   VecType a_n( probDesc->solVecInfo() );
   VecType v_p( probDesc->solVecInfo() );
                                                                                                                            
   // Set up initial conditions
   ptcurState = new SysState<VecType>( d_n, v_n, a_n, v_p);
   probDesc->getInitState( *ptcurState );

   if (!nvStep0) probDesc->getInitState(*ptcurState);
   else probDesc->getPitaState(*ptcurState,0);

   aeroAlg = probDesc->getAeroAlg();
   if(aeroAlg >= 0) probDesc->aeroPreProcess(d_n,v_n,a_n,v_p);

   // if we are doing a ping-pong, return
   if(aeroAlg == 1 || aeroAlg == 8) return;

   if(probDesc->getThermoeFlag() >= 0) probDesc->thermoePreProcess(d_n, v_n, v_p);

   ptaeroForce = (aeroAlg >= 0 || probDesc->getThermoeFlag() >= 0) ? new VecType(info) : 0;   //otherwise =0

   // Build time independent forces i.e. gravity force, pressure force
   probDesc->getConstForce(*ptconstForce);

   // Get SteadyState Flag and Parameters
   int steadyFlag,steadyMin,steadyMax;                                       //useless for implicit
   double steadyTol;
   probDesc->getSteadyStateParam(steadyFlag, steadyMin, steadyMax, steadyTol);

   // Get Time Integration Scheme
   int algType = probDesc->getTimeIntegration();

   switch (algType)
   {
     // Newmark
     case 0:

       // Get Newmark Parameter

       double beta, gamma, alphaf, alpham;
       probDesc->getNewMarkParameters(beta,gamma,alphaf,alpham );

       // ... Newmark Beta == 0 -> CENTRAL DIFFERENCES ALGORITHM
       // ... Newmark Beta != 0 -> NEWMARK ALGORITHM
       //     This is now the Generalized Alpha Method, of which
       //     the NEWMARK algorithm is a subset.

       if(beta == 0) {

         ptworkVec=new NewmarkWorkVec<VecType,ProblemDescriptor> (0,probDesc);

         // Build Necessary Operators (K, M)
         //ptdynOps = probDesc->buildOps(1.0, 0.0, 1.0);

         //Time Integration Loop
         filePrint(stderr," ... Running explicit Newmark ...\n");
         //explicitNewmarkLoop( *ptcurState, *ptconstForce, ptdynOps_dt, ptworkVec, dt, Tfinal,);
       }
       else {
//----------------------------------------------------------------------
//----------------------beginning of PITA-------------------------------

         implicitNewmark_PITA(alphaf,alpham,beta,gamma);

//---------------------End of PITA--------------------------------------
//----------------------------------------------------------------------
       }
       break;

      // Quasi-Static
     case 1:

       double maxVel, delta;
       probDesc->getQuasiStaticParameters(maxVel, delta);

       ptworkVec=new NewmarkWorkVec<VecType,ProblemDescriptor> (-1,probDesc);

       // Build Necessary Operators (only K!)
       //ptdynOps = probDesc->buildOps(0,0,1);

       //if(ptdynOps->dynMat->numRBM() > 0)
       //    fprintf(stderr, " ... Number of RBM(s)     =   %d     ...\n", ptdynOps->dynMat->numRBM());
       // Quasi-Static Loop
       //quasistaticLoop( *ptcurState, *ptconstForce, *ptdynOps, *ptworkVec, dt, Tfinal, aeroAlg);
   }

}
                                                                                                                                                     
//---------------------------------------------------------------------------------------------

template< class DynOps,
          class VecType,
          class PostProcessor,
          class ProblemDescriptor,
          class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::implicitNewmark_PITA(double alphaf,double alpham,double beta,double gamma,int optFlag=0,int* stcopt=0)
{
   // PITA is valid only if alpham=1/2, beta=1/4 and gamma=1/2. In such a case,
   // in the implicit Newmark computation (Pita.d/TimeSlice.C), the terms that 
   // multiply a_n are null. 
   // So, we don't need to compute an=M-1(F-Kd-Dv) when we are updating each 
   // Time Slice Initial Condition.
     
   if (alpham!=0.5 || beta!=0.25 || gamma!=0.5){
      fprintf(stderr,"!!!  PITA isn't valid for these Newmark parameters  !!! ");
      exit(-1);
   }

   // don't need to initialize ext_f. Indeed, ext_f is put to 0 before each computeExt 
   // computation.

   //int parity = 0;
   SysState<VecType> *bkState = 0;
   // Allocate backup state for A5 algorithm
   if(aeroAlg == 5) {
     VecType d_bk(probDesc->solVecInfo());
     VecType v_bk(probDesc->solVecInfo());
     VecType a_bk(probDesc->solVecInfo());
     VecType v_p_bk(probDesc->solVecInfo());
     bkState = new SysState<VecType>(d_bk, v_bk, a_bk, v_p_bk);
   }

   //build workVec
   ptworkVec = new NewmarkWorkVec<VecType,ProblemDescriptor> (1,probDesc);

   //Build K, M, C and dynK = ((1-alpham)/(1-alphaf))*M + (dt*gamma)*C + (dt*dt*beta)*K
   // if cr = true, ptdynOps = buildPitaOps else ptdynOps = buildOps
   Ptimes->dynOpstime -= getTime();

   bool cr  = CkCoarse || !nvStep0;
   ptdynOps = probDesc->buildPitaOps( ((1.0-alpham)/(1.0-alphaf)), dt*gamma, dt*dt*beta,
                                 ((1.0-alpham)/(1.0-alphaf)), Dt*gamma, Dt*Dt*beta, cr);

   //if(ptdynOps->dynMat->numRBM() > 0)
   //    fprintf(stderr, " *** WARNING: %d ZEM(s) Found!  \n", ptdynOps->dynMat->numRBM());

   Ptimes->dynOpstime += getTime();
   if (myCPU==0) cout<<"DynOpstime "<<Ptimes->dynOpstime/1000.0<<endl;
   
   vector_Initialize(true);

   //build TimeSlices and initialization
   Ptimes->buildTStime -= getTime();
   buildTSandConnect();

   for(int i = 0; i<numTS; i++) 
     arrayTimeSlices[i].initialization(ptconstForce,ptaeroForce,alphaf,alpham,beta,gamma,ptdynOps);

   ptCoarseGrid->initialization(ptconstForce,ptaeroForce,alphaf,alpham,beta,gamma,ptdynOps);
   Ptimes->buildTStime += getTime();

   // maxnumITA is used only to make sure it will stop, 
   // has no influence on the code 
   int maxnumITA = nTS * kiter - ((kiter - 1) * kiter) / 2;
   //cout << " !! maxnumITA !! " << maxnumITA << endl;

   int baseSize = std::max(maxnumITA, 2 * probsize);

   /*
    Vectors disp and vel are stored in a StateSet object.  

    ptSeedState = store the initial conditions of the active TS. These 
                  vectors don't need to be kept at the end of a cycle. 
                  So SeedState's maximal size = numTotalActiveTS

    ptBfine, ptStep0, ptPropState, ptJumpState, ptCk have a size = nTS

    Indeed, we need to keep Bi and Step0 values during all the program. 
    we keep Bi values to compute the base and Step0 values to use them as 
    initial conditions for a TS whose previous TS hasn't been activated 
    yet.

    We keep the propagation, jump and Ck values of the inactive TS too 
    in order to compute new Ck.

    Values stored in ptBfine, ptStep0, ptPropState, ptJumpState, ptCk
    are stored by id. TS0 -> pos=0, TS1 -> pos=1...

    Values stored in ptSeedState are stored at the position 
    pos = TSid - cycleTSid[0] 

   */

   //Initialization af all StateSet

   ptSeedState  = new StateSet<VecType,DynOps,InfoSize>(numTotalActiveTS,info);

   ptBfine     = new StateSet<VecType,DynOps,InfoSize>(nTS,info);
   ptStep0     = new StateSet<VecType,DynOps,InfoSize>(nTS,info);
   ptPropState = new StateSet<VecType,DynOps,InfoSize>(nTS,info);
   ptJumpState = new StateSet<VecType,DynOps,InfoSize>(nTS,info);
   ptCk        = new StateSet<VecType,DynOps,InfoSize>(nTS,info);

   ptSk      = new StateSet<VecType,DynOps,InfoSize>(baseSize,info);
   ptAmplSk  = new StateSet<VecType,DynOps,InfoSize>(baseSize,info);
   ptKd_Mv   = new StateSet<VecType,DynOps,InfoSize>(baseSize,info);

   ptProptmp  = new StateSet<VecType,DynOps,InfoSize>(numTSperCycleperCPU,info);

   //beginning of the algorithm
   double totalTime = -getTime();

   bool bi=true; 
   bool ck=true;
   bool step0=true; 

   //Bifine and exchange 
   Ptimes->Bitime -= getTime();
   if (NoForce) {
     //If we don't have any force, don't compute Bi.
     for (int i=0; i<nTS; i++) 
       ptBfine->put_end_StateSet(*ptZero,*ptZero);
   } else if (ConstForce) {
     //When forces are not time-dependent, the Bi's are all equal 
     //Hence we need to compute only the first term. 
     arrayTimeSlices[0].solve_ITA_linearDynam(bi,!step0,!ck,myCPU,ptBfine,0);
     for (int i=1; i<nTS; i++)
       ptBfine->put_end_StateSet(ptBfine->getDisp(0),ptBfine->getVel(0));  
   }else{
     for (int i=0; i<CPUtoTS->num(myCPU); i++)
       arrayTimeSlices[i].solve_ITA_linearDynam(bi,!step0,!ck,myCPU,ptBfine,0);
     exchange_BiValue();
   }
   Ptimes->Bitime += getTime();
   
   //step0
   Ptimes->step0time -= getTime();
   if (!nvStep0) {
     //cout<<" !!!!!!!!!!!!!!!! old step0 !!!!!!!!!!!!!!!!!!!!"<<endl;
     vector_Initialize();    
     ptCoarseGrid->solve_ITA_linearDynam(!bi,step0,!ck,myCPU,ptStep0,0);
   }else{
     //cout<<" !!!!!!!!!!!!!!!! new step0 !!!!!!!!!!!!!!!!!!!!"<<endl;
     ptStep0->adjustnumVectors();
     for (int l=0; l<nTS; l++){
       probDesc->getPitaState(*ptcurState,l);
       ptStep0->put_end_StateSet(ptcurState->getDisp(),ptcurState->getVeloc());
     }
   }
   Ptimes->step0time += getTime();

   //ITA on the fine grid ( don't need to re-initialize ext_f because are set to 0 in computeExtForce)
   vector_Initialize();                         

   VecType dispJump(info); 
   VecType velJump(info);
   
   int cycle = 0; locidTS = 0; 
   
   //beginning of the loop
   while (numTotalActiveTS != 0 && cycle < maxnumITA){

     //print for each MPI Process the number of active TS and their id 
     if (myCPU==0) printInfoCPU(); 
     
     //fill this cycle's cycleTSid and buffeurs 
     getActiveTSidAndBuffeur();

     //get init cond                            
     Ptimes->newCondtime -= getTime();
     if (active) newInitCond(cycle);
     Ptimes->newCondtime += getTime();
      
     //ITA and store Prop value in buffeur
     if (active) {
       //solveITA == time used to solve the problem
       //on the fine grid (without Bi computation)
       Ptimes->solveITA -= getTime();
       ptProptmp->adjustnumVectors();
       for (int i=0; i<infoCPU[myCPU*x]; i++){
         arrayTimeSlices[locidTS+i].solve_ITA_linearDynam(!bi,!step0,!ck,myCPU,ptProptmp,cycleTSid[0]);
         arrayTimeSlices[locidTS+i].addnumITA();
       }
       Ptimes->solveITA += getTime();
     }

     //Exchange buffeur 
     if (numCPU>0) {
       Ptimes->exchangeProptime -= getTime();
       exchange_PropValue();  
       Ptimes->exchangeProptime += getTime();
     }

     //jump
     Ptimes->jump -= getTime();
     if (active) jumpEvaluation();
     Ptimes->jump += getTime();

     //norm of reference
     if (cycle==0) normRefcomputation();
 
     //stop if all jumps are small
     Ptimes->compNorm -= getTime();
     if (active) stopIfJumpSmall();
     Ptimes->compNorm += getTime();

     //Ck
     Ptimes->CkTotaltime -= getTime();
     if (active) Ck_computation(cycle); 
     Ptimes->CkTotaltime += getTime();

     //number and TSid of the next cycle for this MPIProcess 
     Ptimes->newTSid -= getTime();
     myNextActiveTS();
     Ptimes->newTSid += getTime();

     //exchange this info with other MPIProcess
     if (numCPU>0) {
       Ptimes->exchange_info -= getTime();
       exchange_InfoCPU();
       Ptimes->exchange_info += getTime();
     }

     // new numTotalActiveTS 
     numTotalActiveTS = 0;
     for (int i=0; i<numCPU; i++){
       if (infoCPU[x*i] != 0) numTotalActiveTS += infoCPU[x*i];
     }

     cycle++;
   }

   Ptimes->PitaTotaltime += getTime();

   if ( ! optFlag ) {
     totalTime += getTime();
     filePrint(stderr," --------------------------------------\n");
     filePrint(stderr," ... Total Loop Time      = %.2e s\n",totalTime/1000.0);
     fflush(stderr);

   if (myCPU==0) {
     double timeLoopTime = 0.0; // TODO meaningful value
     probDesc->printTimers(ptdynOps, timeLoopTime);
   }
   }
  
   Ptimes->printPitaTimesFile();
}

#endif
