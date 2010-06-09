#ifndef _STATESET_C_
#define _STATESET_C_

#include <Pita.d/StateSet.h>
#include <Math.d/SparseMatrix.h>

class Communicator;
extern Communicator* structCom;

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
StateSet<VecType,DynOps,InfoSize>::StateSet()
{
   len = 0;
   maxnumVectors = 0;
   numVectors  = currentVector = 0;

   dispSet = NULL;
   velSet  = NULL;
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
StateSet<VecType,DynOps,InfoSize>::StateSet(int maxnum, const InfoSize &inf)
{

  maxnumVectors = 0;
  numVectors = currentVector = 0;
  
  info = inf;  
  len = info.totLen();
  
  setMaxNumVectors(maxnum);
 
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
StateSet<VecType,DynOps,InfoSize>::~StateSet()
{
  emptySet();
}

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::copySet(const StateSet &ss)
{

  emptySet();
  maxnumVectors = ss.maxnumVectors;
  numVectors    = ss.numVectors;  
  len           = ss.len;
  info          = ss.info;
  currentVector = ss.currentVector;    
 
  if (maxnumVectors > 0)
  {
  
    VecType **newDispSet = new VecType *[maxnumVectors];
    VecType **newVelSet = new VecType *[maxnumVectors];
 
    for (int i = 0; i < numVectors; ++i)
    {
      dispSet[i] = new VecType(ss.dispSet[i]);
      velSet[i] = new VecType(ss.velSet[i]); 
    }
    for (int i = numVectors; i < maxnumVectors; ++i)
    {
      dispSet[i] = NULL;
      velSet[i] = NULL;
    }

  }

}

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::mergeSet(const StateSet &ss)
{

  if (ss.len != len)
  {
    fprintf(stderr, "Error in StateSet::mergeSet : Lengths do not match");
    return;
  }

  if (ss.numVectors > 0)
  {
    int newNumVectors = numVectors + ss.numVectors;
    if (newNumVectors > maxnumVectors)
      setMaxNumVectors(newNumVectors);
    for (int i = 0; i < ss.numVectors; ++i)
      addState(ss.getDisp(i), ss.getVel(i)); 
    numVectors = newNumVectors;
  }

  currentVector = numVectors;

}

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::emptySet()
{

   if (maxnumVectors > 0)
   {

     for (int i = 0; i < numVectors; ++i)
       delete dispSet[i];
     delete[] dispSet;
     dispSet = NULL;
     
     for (int i = 0; i < numVectors; ++i)
       delete velSet[i];
     delete[] velSet;
     velSet = NULL;

     maxnumVectors = 0;
     numVectors = currentVector = 0;
  
   }

}

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::newSet(int maxnum, const InfoSize &inf)
{

  emptySet();
  info = inf;
  len = info.totLen();
  setMaxNumVectors(maxnum);

}


//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::setMaxNumVectors(int n)
{
  if (maxnumVectors == n)
    return;

  if (maxnumVectors > 0)
  {
   
    VecType **newDispSet = new VecType *[n];
    VecType **newVelSet = new VecType *[n];
    
    int newNumVectors = std::min(n, numVectors);
    int newCurrentVector = std::min(n, currentVector); 
    for (int i = 0; i < newNumVectors; ++i)
    {
      newDispSet[i] = dispSet[i];
      newVelSet[i] = velSet[i];
    }
    for (int i = newNumVectors; i < numVectors; ++i)
    {
      if(dispSet[i]) delete dispSet[i];
      if(velSet[i]) delete velSet[i];
    } 
   
    delete[] dispSet;
    delete[] velSet;

    maxnumVectors = n;
    numVectors = newNumVectors;
    currentVector = newCurrentVector;
    dispSet = newDispSet;
    velSet  = newVelSet; 

  }
  else
  {
    maxnumVectors = n;
    if (maxnumVectors > 0)
    {
      dispSet = new VecType*[maxnumVectors];
      for (int i = 0; i < maxnumVectors; ++i)
        dispSet[i] = NULL;
      
      velSet  = new VecType*[maxnumVectors];
      for (int i = 0; i < maxnumVectors; ++i)
        velSet[i] = NULL;
    }
  }
}

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::setState(const VecType &disp, const VecType &vel, int i)
{
  if (i >= numVectors)
    fprintf(stderr," .... Error in StateSet::setState, out of range .... \n");
  else
  {
    if (disp.size() != len || vel.size() != len)
      fprintf(stderr," .... Error in StateSet::setState, lengths do not match .... \n");
    else
    {
      *(dispSet[i]) = disp;
      *(velSet[i]) = vel;
    }
  }

}


//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::adjustnumVectors(int i) 
{
 
  if (i > numVectors)
   {
     fprintf(stderr," . ... Error in StateSet::adjustnumVectors : out of range (%d/%d)\n", i, numVectors);
     exit(1);
   }
  else
    currentVector = i;
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::put_end_StateSet(const VecType &disp, const VecType &vel)
{

  if ((disp.size()!=len)||(vel.size()!=len))
  {
    fprintf(stderr," .... Error in Put_end_StateSet, size pb .... \n");
    return;
  }
  if (currentVector >= maxnumVectors)
  {
    fprintf(stderr," .... Error in Put_end_StateSet, numvect problem .... \n");
    return;
  }

  if (currentVector >= numVectors)
  {
    currentVector = numVectors;
    dispSet[numVectors] = new VecType(info);
    velSet[numVectors]  = new VecType(info);
    ++numVectors;
  }

  *(dispSet[currentVector]) = disp;
  *(velSet[currentVector]) = vel;

  ++currentVector;
  
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::put_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i)
{
  if (S.getSize()!=len)        fprintf(stderr," .... Error in put_end_StateSet, length problem .... \n");
  if (i>=S.getnumVectors()) fprintf(stderr," .... Error in put_end_StateSet, index problem .... \n");
  if (currentVector>=maxnumVectors) fprintf(stderr," .... Error in put_end_StateSet, numvect pb .... \n");

  if (currentVector >= numVectors)
  {
    currentVector = numVectors;
    dispSet[numVectors] = new VecType(info);
    velSet[numVectors]  = new VecType(info);
    ++numVectors;
  }

  *(dispSet[currentVector]) = S.getDisp(i);
  *(velSet[currentVector]) = S.getVel(i);

  ++currentVector;
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::put_end_StateSet(double *array, int num, const VecType &vel)
{
  if (len!=num)                  fprintf(stderr," .... Error in put_end_StateSet, wrong array size .... \n");
  if (currentVector == maxnumVectors) fprintf(stderr," .... Error in Put_end_StateSet, numvect problem .... \n");

  if (currentVector >= numVectors)
  {
    currentVector = numVectors;
    dispSet[numVectors] = new VecType(info);
    velSet[numVectors]  = new VecType(info);
    ++numVectors;
  }

  dispSet[currentVector]->getDataFrom(array,num);
  *(velSet[currentVector])=vel;

  ++currentVector;
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::add_end_StateSet(const VecType &disp, const VecType &vel)
{
   if (disp.size() != len || vel.size() != len )
   {
     fprintf(stderr, " .... Error in add_end_StateSet, lengths do not match ... \n");
     return;
   }

   *(dispSet[currentVector-1]) += disp;
   *(velSet[currentVector-1]) += vel;
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::add_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i)
{
   if ( S.getSize()!=len  || i>=S.getnumVectors() )
   {
      fprintf(stderr, " .... Error in add_end_StateSet, lengths do not match .... \n");
      return;
   }

   *(dispSet[currentVector-1]) += S.getDisp(i);
   *(velSet[currentVector-1]) += S.getVel(i);
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::add_end_StateSet(const VecType &disp, double *array, int num)
{
  if (len!=num || len!=disp.size())
  {
    fprintf(stderr," .... Error in add_end_StateSet, wrong array size .... \n");
    return;
  }

  (velSet[currentVector-1])->addDataFrom(array,num);
  *(dispSet[currentVector-1])+=disp;
}

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::sub_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i)
{
   if (S.getSize()!=len)        fprintf(stderr," .... Error in add_end_StateSet, length problem .... \n");
   if (i>=S.getnumVectors()) fprintf(stderr," .... Error in add_end_StateSet, numVectors problem .... \n");
                                                                                                                                                
   *(dispSet[currentVector-1]) -= S.getDisp(i);
   *(velSet[currentVector-1]) -= S.getVel(i);
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::replace_end_StateSet(const VecType &disp, const VecType &vel)
{
  setState(disp, vel, currentVector - 1);
}
                                                                                                                                      
//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::replace_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i)
{
  setState(S.getDisp(i), S.velDisp(i), currentVector - 1);
}

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::replace_StateSet(int k, const VecType &disp, const VecType &vel)
{
  setState(disp, vel, k);
}

//-------------------------------------------------------------------------------
template < 
  class VecType, 
  class DynOps, 
  class InfoSize >
void StateSet<VecType,DynOps,InfoSize>::buildBases(StateSet<VecType,DynOps,InfoSize> &Prop, 
           StateSet<VecType,DynOps,InfoSize> &Sk, StateSet<VecType,DynOps,InfoSize> &AmplSk, DynOps& dynOps, 
           StateSet<VecType,DynOps,InfoSize> &Kd_Mv, StateSet<VecType,DynOps,InfoSize> &Bi, int numActiveTS, 
           int *cycleTSid)
{
  /* 
    this method is applyed to SeedState Yi=yi(Ti). 
    Prop = yi(Ti+1) = (Gfine)J Yi + Bi with  Gfine = amplification matrix on the fine time grid.
   
    This method builds a base and its amplification using an Iterative Modified Grahm 
    Schmidt (IMGS) algorithm:
         .Sk is the base obtained with the Seed vectors ie Yi=yi(Ti)
         .AmplSk is the Amplification of the base Sk. ie AmplSki = (Gfine)J Ski
          It corresponds to the propagation vectors ie yi(Ti+1) and the data Bi
    
    tmpd == temporary vector for disp 
    tmpv == temporary vector for vel
    
    _S == for Seed base;  _AmplS == for the Amplification base
    
    The metric used here is Q = (M 0 ; 0 K).  
    
    <a|b> = dot product = bT Q a

    (Yi)i         ->  (Si)i     = ( Yi - <Yi|Sj>Sj ) / norm    
    (Pi)i, (Bi)i  ->  (AmplSi)i = ( Ampl(Yi) - <Yi|Sj> Ampl(Sj)) / norm
                                = ( (yi(T+1) - Bi) - <Yi|Sj> AmplSj )   

    For each given vector a, <a|Si> = aT Q Si = a_vT (M Si_v) + a_dT (K Si_d).
    The vectors (M Si_v) and (K Si_d) are temporarily stored to improve performance.

  */


  double atol = 1e-10;
  double rtol = 1e-4;

  double ps = 0.0;

  double norm0, norm1, norm;
  norm = norm0 = norm1 = 0.0;

  int cursor = cycleTSid[0];

  VecType tmpd_S(info);
  VecType tmpv_S(info);
  VecType tmpd_AmplS(info);
  VecType tmpv_AmplS(info);

  VecType d(info); d.zero(); 
  VecType v(info); v.zero();

  int maxit = 1;

  int count = 0;
 
  for (int i = 0; i < numActiveTS; i++){
     if (Sk.getnumVectors() < 2 * len) {

       tmpd_S = getDisp(i);          //tmpd_S=Yi
       tmpv_S = getVel(i);

       tmpd_AmplS.linC(1, Prop.getDisp(cursor+i), -1, Bi.getDisp(cursor+i));
       tmpv_AmplS.linC(1, Prop.getVel(cursor+i), -1, Bi.getVel(cursor+i));
  
       dynOps.K->mult(tmpd_S,d); dynOps.M->mult(tmpv_S,v);

       norm0 = sqrt(tmpd_S*d + tmpv_S*v);
       norm1 = norm0;
 
       //IMGS
       for (int iter = 0; iter < maxit; iter++){              
         for (int j = 0; j < Sk.getnumVectors(); j++){
 
            ps = (tmpd_S * Kd_Mv.getDisp(j)) + (tmpv_S * Kd_Mv.getVel(j));                       
  
            //Yi - <Yi|Sj>Sj   
            tmpd_S.linAdd(-1*ps, Sk.getDisp(j));        
            tmpv_S.linAdd(-1*ps, Sk.getVel(j));             

            //(yi(T+1)-Bi) - <Yi|Sj> AmplSj
            tmpd_AmplS.linAdd(-1*ps, AmplSk.getDisp(j));    
            tmpv_AmplS.linAdd(-1*ps, AmplSk.getVel(j));
         }

         dynOps.K->mult(tmpd_S,d); dynOps.M->mult(tmpv_S,v);
         norm=sqrt(tmpd_S*d + tmpv_S*v);

         if(norm>0.5*norm1) break;                  //IMGS
         norm1=norm;
       }
       //cout<<"norm0, norm , getDisp(i)[86] "<<norm0<<" "<<norm<<" "<<getDisp(i)[86]<<endl;

       //norm<1.5*norm0 is used to avoid case where norm>>norm0. This appears when
       //the base isn't built correctly and "expolde" 
       if (norm>max(atol,rtol*norm0) && norm<1.5*norm0){

         Kd_Mv.put_end_StateSet((1.0/norm)*d,(1.0/norm)*v);

         Sk.put_end_StateSet((1.0/norm)*tmpd_S,(1.0/norm)*tmpv_S);
         AmplSk.put_end_StateSet((1.0/norm)*tmpd_AmplS,(1.0/norm)*tmpv_AmplS);

       } else { ++count; }

     }
  }
     //cout << "# vectors : Final / Initial : " << Sk.getnumVectors() << " / " << count << endl;
 
}

template <
  class VecType,
  class DynOps,
  class InfoSize >
void StateSet < VecType, DynOps, InfoSize >::buildOrtho(StateSet<VecType,DynOps,InfoSize> &final, const DynOps &dynOps)
{
  StateSet temp;
  return;
}


#endif
