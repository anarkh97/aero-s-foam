#include <Feti.d/Feti.h>

template<>
void
GenFetiSolver<DComplex>::reSendInterfaceRBM(int iSub)
{
  fprintf(stderr, "WARNING: GenFetiSolver<DComplex>::reSendInterfaceRBM(...) not implemented \n");
}

template<>
void
GenFetiSolver<double>::reSendInterfaceRBM(int iSub)
{
 // sd[iSub]->sendRBMs(fetiOps[iSub]->numRBM, fetiOps[iSub]->locRBMs, 
 //                    fetiOps[iSub]->locInterfRBMs);
 sd[iSub]->extractInterfRBMs(fetiOps[iSub]->numRBM, fetiOps[iSub]->locRBMs,
                             fetiOps[iSub]->locInterfRBMs);

 GenCoarseSet<double> &thisSet = fetiOps[iSub]->control->cset[sd[iSub]->subNum()];

 thisSet.numGs = fetiOps[iSub]->numRBM;
 thisSet.locGs = fetiOps[iSub]->locInterfRBMs;

 if(QGisLocal == 0) { // If there is a Q
   if(fetiOps[iSub]->rbm) {  // In the dynamic case, Q is Fi
     sd[iSub]->multMFi(fetiOps[iSub]->solver, thisSet.locGs, 
                       thisSet.locQGs, fetiOps[iSub]->numRBM);
   } 
   else {
     for(int i = 0; i < fetiOps[iSub]->numRBM; ++i)
       sd[iSub]->multKbb(thisSet.locGs + i*sd[iSub]->interfLen(),
                         thisSet.locQGs +i*sd[iSub]->interfLen());
   }
 } 
 else {
   fetiOps[iSub]->control->cset[sd[iSub]->subNum()].locQGs 
         = fetiOps[iSub]->locInterfRBMs; // no preconditioning
 }

 if(isFeti2 && (isDynamic == 0)) {
   if(QGisLocal == 0)
     thisSet.locFGs = thisSet.locQGs;
   else
     sd[iSub]->multMFi(fetiOps[iSub]->solver, thisSet.locGs, 
                       thisSet.locFGs, fetiOps[iSub]->numRBM);
 }
}
