#include <Utils.d/dbg_alloca.h>
#include <stdlib.h>

#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/Skyline.d/SGISky.h>
#include <Solvers.d/UFront.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Driver.d/Domain.h>
template<>
void
Domain::makeFrontalOps(AllOps<DComplex> &ops, double Kcoef, double Mcoef, double Ccoef,
                       Rbm *rbm, FullSquareMatrix *kelArray)
{
  fprintf(stderr, "WARNING: Domain::makeFrontalOps(...) is not implmented for DComplex type \n");
}


template<>
void
Domain::makeFrontalOps(AllOps<double> &ops, double Kcoef, double Mcoef, double Ccoef,
                       Rbm *rbm, FullSquareMatrix *kelArray)
{
 if(matrixTimers) matrixTimers->assemble -= getTime();

 int makeMass = Mcoef != 0 || ops.M != 0 || ops.C != 0;

 // Rayleigh Damping coefficients: C = alpha*M + beta*K
 double alpha = sinfo.alphaDamp;
 double  beta = sinfo.betaDamp;

 int size    = sizeof(double)*maxNumDOFs*maxNumDOFs;
 double *marray = (double *) dbg_alloca(size);
 double *karray = (double *) dbg_alloca(size);

 FullSquareMatrix kel(1,karray);
 FullSquareMatrix mel(1,marray);

 int *nodeOrder = renumb.order;

 int i, j, inode, iele, ele;

 // Flags elements that have already been introduced to the front
 int *hasAppeared =(int*) dbg_alloca(sizeof(int)*numele);

 int numdofs   = numUncon();
 int maxfrsize = 0;
 int curfrsize = 0;

 int *appears = (int*) dbg_alloca(sizeof(int)*numnodes);

 for(i=0; i<numnodes; ++i)
   appears[i] = 0;

 for(i=0; i<numele; ++i)
   hasAppeared[i] = 0;

 // Loop over the node by order of elimination to find the elements
 // that need to be introduced in the front before the node can be eliminated
 // Also determine the maximum front size
 for(inode = 0; inode < numnodes && nodeOrder[inode] >= 0; ++inode) {
   int node = nodeOrder[inode];
   for(iele = 0; iele < nodeToElem->num(node); ++iele) {
     int ele = (*nodeToElem)[node][iele];
     if(hasAppeared[ele] == 0) { // insert the element now
        for(j = 0; j < elemToNode->num(ele); ++j) {
           if(appears[(*elemToNode)[ele][j]] == 0) {
              appears[(*elemToNode)[ele][j]] = 1;
              curfrsize += c_dsa->weight((*elemToNode)[ele][j]);
           }
        }
        hasAppeared[ele] =1;
     }
   }
   if(curfrsize > maxfrsize)
      maxfrsize = curfrsize;
   curfrsize -= c_dsa->weight(node);
 }

 // ... Call the Unrolled Frontal's constructor
 UFront *fr = constructFrontal(maxfrsize, rbm);

 for(i=0; i<numele; ++i)
   hasAppeared[i] = 0;

// Add discrete mass contribution to the Mass Matrix
// Three matrices need to be changed. Mass matrix itself,
// Damping Matrix and K tilda. The latter has to be on the fly, so we create an
// array of discrete masses for each DOF

 double *dMass = (double *) dbg_alloca(sizeof(double)*numdofs);
 for(i = 0; i < numdofs; ++i)
   dMass[i] = 0.0;

 if(makeMass) {
   DMassData *current = firstDiMass;
   while(current != 0) {
     int cdof = c_dsa->locate(current->node, (1 << current->dof));
     if(cdof >= 0) {
        if(ops.M) ops.M->addDiscreteMass(cdof, current->diMass);
        if(ops.C) {
           ops.C->addDiscreteMass(cdof, alpha*current->diMass);
           dMass[cdof] += Ccoef*alpha*current->diMass;
        }
        dMass[cdof] += Mcoef*current->diMass;
     }
     current = current->next;
   }
 }


 for(inode = 0; inode < numnodes && nodeOrder[inode] >= 0; ++inode) {
   int node = nodeOrder[inode];

   int cmflg = 1;
   if ((sinfo.probType == SolverInfo::Dynamic)
      &&(solInfo().newmarkBeta == 0.0))
      cmflg = 0;
   for(iele = 0; iele < nodeToElem->num(node); ++iele) {
      ele = (*nodeToElem)[node][iele];
      if(hasAppeared[ele] == 0) {
         Element * c_ele = packedEset[ele] ;
         if(c_ele == 0)
           continue ;
         if(matrixTimers) matrixTimers->formTime -= getTime();

         if(kelArray)
           kel = kelArray[iele];
         else
           kel = c_ele->stiffness(nodes,karray) ;

         if(matrixTimers) matrixTimers->formTime += getTime();
	 if(makeMass) {
           mel = c_ele->massMatrix(nodes,marray,cmflg);
         }

         if(ops.K)   ops.K->add(kel,(*allDOFs)[ele]);
         if(ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[ele]);
         if(ops.M)   ops.M->add(mel,(*allDOFs)[ele]);
         if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[ele]);
         if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[ele]);

         int dim = kel.dim();
         int i,j;
         if(makeMass)
           for(i = 0; i < dim; ++i)
             for(j = 0; j < dim; ++j) {
               double m  = mel[i][j];
               double k  = kel[i][j];
	       mel[i][j] = alpha*m + beta*k;
               kel[i][j] = Kcoef*k + Ccoef*mel[i][j] + Mcoef*m;
             }
	 else
	   for(i = 0; i < dim; ++i)
             for(j = 0; j < dim; ++j)
	       kel[i][j] *= Kcoef;

         fr->addkel(kel,(*allDOFs)[ele]);

         // Damping matrix
         if(ops.C)   ops.C->add(mel,(*allDOFs)[ele]);
         if(ops.Cuc) ops.Cuc->add(mel,(*allDOFs)[ele]);

         hasAppeared[ele] = 1;
      }
   }
   // first dof of this node
   int firstdof = c_dsa->firstdof(node);
   int idof;
   for(idof = 0; idof < c_dsa->weight(node); ++idof) {
       int dof = firstdof + idof;
       fr->addToDiag(dof,dMass[dof]);
       fr->elim(dof) ;
   }
 }

 fr->finishUpdate();

 ops.sysSolver = fr;

 if(matrixTimers) matrixTimers->assemble += getTime();
}

SGISky *
Domain::constructSGISkyMatrix(Rbm *rbm)
{
  return new SGISky(nodeToNode,dsa,c_dsa,sinfo.trbm,rbm);
}

UFront *
Domain::constructFrontal(int maxFrontSize, Rbm *rbm)
{
  return new UFront(c_dsa, maxFrontSize, sinfo.trbm, rbm);
}

