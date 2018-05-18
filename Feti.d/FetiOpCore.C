#include <Feti.d/Feti.h>
#include <Feti.d/FetiOp.h>

#if defined(WINDOWS) || defined(MACOSX)
 #include <cfloat>
#else
 #include <climits>
#endif
#include <float.h>

template<>
void
GenFetiOp<DComplex>::sendInterfRBM(FSCommPattern<DComplex> *rbmPat)
{
  fprintf(stderr, "WARNING: GenFetiOp<DComplex>::sendInterfRBM() not implemented  \n");
}

template<> 
void
GenFetiOp<double>::sendInterfRBM(FSCommPattern<double> *rbmPat)
{
 int isUsual = 1;

 locInterfRBMs = new double[numRBM*sd->interfLen()];

 // sd->sendRBMs(numRBM, locRBMs, locInterfRBMs);
 sd->extractInterfRBMs(numRBM, locRBMs, locInterfRBMs);
 sd->sendInterfRBMs(numRBM, locInterfRBMs, rbmPat);

 GenCoarseSet<double> &thisSet = control->cset[sd->localSubNum()];
 thisSet.numGs = numRBM;
 thisSet.locGs = locInterfRBMs;

 if(QGisLocal == 0) { // If there is a Q
   int i;
   thisSet.locQGs = new double[numRBM*sd->interfLen()];
   if(rbm) {  // In the dynamic case, Q is Fi
     // multiple rhs version of multFi
     sd->multMFi(solver, thisSet.locGs, thisSet.locQGs, numRBM);
   } 
   else {
     if(control->nQ == 4) {
       for(i = 0; i< numRBM; ++i)
         sd->multDiagKbb(thisSet.locGs  + i*sd->interfLen(),
                         thisSet.locQGs + i*sd->interfLen());
     } 
     else {
       for(i = 0; i < numRBM; ++i)
          sd->multKbb(thisSet.locGs  + i*sd->interfLen(),
                      thisSet.locQGs + i*sd->interfLen());
     }
   }
 }
 else if(isUsual)
   control->cset[sd->localSubNum()].locQGs = locInterfRBMs;// no preconditioning
 else {
   DofSetArray *dsa = sd->c_dsa;
   int numDofs = dsa->size();
   double *mask = (double *) dbg_alloca(sizeof(double)*numDofs);
   int numNodes = dsa->numNodes();
   int idof,inode;
   for(idof = 0; idof < numDofs; ++idof)
     mask[idof] = 1;
   for(inode = 0; inode < numNodes; ++inode) {
     idof = dsa->locate(inode,DofSet::Xrot);
     if(idof >= 0) mask[idof] = 0;
     idof = dsa->locate(inode,DofSet::Yrot);
     if(idof >= 0) mask[idof] = 0;
     idof = dsa->locate(inode,DofSet::Zrot);
     if(idof >= 0) mask[idof] = 0;
   }
   
   auto weight = sd->weight;
   auto allBoundDofs = (*sd->scomm->sharedDOFs)[0];
   thisSet.locQGs = new double[numRBM*sd->interfLen()];
   int interfaceLen = sd->interfLen();
   int i,j;
   double *scale = sd->scaling;
   double minScale = DBL_MAX, maxScale = DBL_MIN;
   for(i = 0; i < numRBM; ++i) {
     for(j=0; j< interfaceLen; ++j) {
       thisSet.locQGs[i*interfaceLen+j] = mask[allBoundDofs[j]]*
       locInterfRBMs[i*interfaceLen+j]/(weight[allBoundDofs[j]] -1);
       if(scale[j] > maxScale) maxScale = scale[j];
       if(scale[j] < minScale) minScale = scale[j];
     }
   }
   fprintf(stderr,"Min Max scale %e %e\n",minScale, maxScale);
 }

 if (isFeti2 && isDynamic == 0) {
   if (QGisLocal == 0) {
     thisSet.locFGs = thisSet.locQGs;
   } else {
     thisSet.locFGs = new double[numRBM*sd->interfLen()];
     sd->multMFi(solver, thisSet.locGs, thisSet.locFGs, numRBM);
   }
 }
}
