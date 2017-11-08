#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>
#include <Utils.d/SolverInfo.h>
#include <Utils.d/OutputInfo.h>
#include <Utils.d/resize_array.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/Vector.h>
#include <Solvers.d/CRSolver.h>
#include <Solvers.d/BCGSolver.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <HelmAxi.d/FourierDescrip.h>
#include <HelmAxi.d/AxiHElem.h>
#include <Driver.d/GeoSource.h>

FourierStatic::FourierStatic(Domain *d, FourierHelmBCs *fouBC) {
  domain = d;
  bcs = fouBC;
  hParams = new GenHelmParam;
  hParams->kappaHelm = geoSource->kappa();
}


void FourierStatic::make_bc(int *bc) {
 int i;
 int maxDof = domain->numdof();
 DofSetArray *tab = domain->getDSA();

 for (i=0; i<maxDof; ++i)
   bc[i]  = BCFREE;

// Set the Complex Neuman boundary conditions
 for (i=0; i<bcs->numNeu; ++i) {
   int dof=tab->locate(bcs->neuBC[i].nnum1, (1 << 7));
   bc[dof] = BCLOAD;
   dof=tab->locate(bcs->neuBC[i].nnum2, (1 << 7));
   bc[dof] = BCLOAD;
   if (bcs->neuBC[i].nnum3>-1) {
     dof=tab->locate(bcs->neuBC[i].nnum3, (1 << 7));
     bc[dof] = BCLOAD;
   }
 }

// Set the Complex Dirichlet boundary condtions
 for (i=0; i<bcs->numDir; ++i) {
   int dof=tab->locate(bcs->dirBC[i].nnum, (1 << 7));
   if(bc[dof] == BCFIXED) {
     fprintf(stderr,"WARNING: check input, found repeated AXIHDIR "
                    "(node %d, dof %d)\n",bcs->dirBC[i].nnum,1);
   }
   bc[dof] = BCFIXED;
 }
}


ComplexVector FourierStatic::make_bcxC(int mode) {
 int i;
 ComplexVector V(domain->numdof(),DComplex(0,0));
 DofSetArray *tab = domain->getDSA();
 CoordSet &noeuds = domain->getNodes();

// Set the Complex Dirichlet boundary conditions
 for (i=0; i<bcs->numDir; ++i) {
   int dof = tab->locate(bcs->dirBC[i].nnum, (1 << 7));
   if (dof < 0) continue;
   Node nd = noeuds.getNode(bcs->dirBC[i].nnum);
   V[dof] = bcs->Dir2Four(mode,i,hParams->kappaHelm,nd.x,nd.y);
 }

 return (V);
}


void FourierStatic::preProcess() {

 BCond *buffer;
 int i;

 buffer = new BCond[bcs->numDir];
 for (i=0; i<bcs->numDir; ++i) {
    buffer[i].nnum = bcs->dirBC[i].nnum;
    buffer[i].dofnum = 7;
    buffer[i].val = 0.0;
 }

 domain->preProcessing();

 int numdof = domain->numdof();
 int *bc = (int *) dbg_alloca(sizeof(int)*numdof);

 make_bc(bc);

 domain->makeCDSA(bcs->numDir,buffer);

 ConstrainedDSA* c_tab = domain->getCDSA();
 DofSetArray *tab = domain->getDSA();

 delete[] buffer;
 domain->makeAllDOFs();
 //SolverInfo &sol = domain->solInfo();

 if (bcs->numDir>0)
   kuc = new CuCComplexSparse(domain->getNodeToNode(),tab,c_tab);
 else
   kuc = 0;

}


int
FourierStatic::solVecInfo() {
  return solver->dim();
}


ComplexSolver*
FourierStatic::getSolver() {
  return solver;
}


int
FourierStatic::getModes() {
  return bcs->numModes;
}


int
FourierStatic::getSlices() {
  return bcs->numSlices;
}


void FourierStatic::constructK() {

 SolverInfo &sinfo = domain->solInfo();

 switch (sinfo.solvercntl->type) {
   default:
   case 0:
     switch (sinfo.solvercntl->subtype) {
       default :
          fprintf(stderr," ... WARNING: invalid direct method ...\n");
       case 1 :
          fprintf(stderr," ... Skyline Solver Selected        ...\n");
          solver = constructSkyMatrixC();
          break;
       case 3 :
          fprintf(stderr," ... Complex Sparse Solver Selected ...\n");
          solver = constructBLKSparseMatrixC();
          break;
     }
     break;
//   case 1:
//      SparseMatrix *SpM;
//      SpM = constructBLKSparseMatrixC();
//      makeComplexK(mode, SpM);
//      switch(sinfo.iterType) {
//          default:
//          case 4:
//            fprintf(stderr," ... Using Bi-Conjugate Gradient    ...\n");
//            solver = new ComplexBCGSolver(sinfo.maxit, sinfo.tol, SpM);
//            break;
//          case 5:
//            fprintf(stderr," ... Using Conjugate Residual       ...\n");
//            solver = new ComplexCRSolver(sinfo.maxit, sinfo.tol, SpM);
//            break;
//      }
//     break;
 }

}


void FourierStatic::BuildKs(int mode) {

 SolverInfo &sinfo = domain->solInfo();
 ComplexSparseMatrix *spm;

 switch (sinfo.solvercntl->type) {
   default:
   case 0:
     switch (sinfo.solvercntl->subtype) {
       default :
       case 1 :
         spm = (SkyMatrixC*) solver;
         makeComplexK(mode,spm);
//         solver = (SkyMatrixC*) spm;
         break;
       case 3 :
         spm = (BLKSparseMatrixC*) solver;
         makeComplexK(mode,spm);
//         solver = (BLKSparseMatrixC*) spm;
         break;
     }
     break;
//   case 1:
//      SparseMatrix *SpM = solver->getA();
//      makeComplexK(mode, SpM);
//      solver->setA(SpM);
//      break;
 }

}


SkyMatrixC *
FourierStatic::constructSkyMatrixC() {

 SolverInfo &sinfo = domain->solInfo();
 ConstrainedDSA* c_tab = domain->getCDSA();
 DofSetArray *tab = domain->getDSA();

 SkyMatrixC *K = new SkyMatrixC(domain->getNodeToNode(),tab,c_tab,sinfo.solvercntl->trbm);

 return K;

}


BLKSparseMatrixC *
FourierStatic::constructBLKSparseMatrixC() {

 SolverInfo &sinfo = domain->solInfo();
 ConstrainedDSA* c_tab = domain->getCDSA();
 DofSetArray *tab = domain->getDSA();

 BLKSparseMatrixC *K = new BLKSparseMatrixC(domain->getNodeToNode(),
                           tab, c_tab, sinfo.solvercntl->trbm, *sinfo.solvercntl);

 return K;

}


void
FourierStatic::makeComplexK(int mode, ComplexSparseMatrix *K) {

 DofSetArray *tab = domain->getDSA();

 double pi=4*atan(1.0);
 int maxDof = domain->maxNumDOF();
 double *pointer = (double *) dbg_alloca(sizeof(double)*maxDof*maxDof);
 double a1=2*pi;
 double a2=0;
 double kappa=hParams->kappaHelm;

 K->zeroAll();
 if (kuc) kuc->zeroAll();

 if (mode>0) {
    a1 = pi;
    a2 = pi*mode*mode;
 }

 int iele;
 int maxele = domain->numElements();
 CoordSet &noeuds = domain->getNodes();
 Connectivity *alltab = domain->getAllDOFs();
 Elemset &table = domain->getElementSet();
 for (iele=0; iele<maxele; ++iele) {
    AxiHElement *elem = dynamic_cast<AxiHElement *>(table[iele]);
    if(elem == 0)  {
       int one=1;
       fprintf(stderr,"Element chosen non axisymmetric. Aborting\n");
       exit(one);
    }

    FullSquareMatrix kel = elem->stiffness(noeuds,pointer);

    kel *= a1;

    K->add(kel,(*alltab)[iele]);
    if (kuc) kuc->add(kel,(*alltab)[iele]);

    kel.zero();
    kel = elem->stiffteta(noeuds,pointer);

    kel *= a2;

    K->add(kel,(*alltab)[iele]);
    if (kuc) kuc->add(kel,(*alltab)[iele]);

 }

 int i;
 //int MaxNodes = domain->numNode();

 if (bcs->numSomm>0) {

   int somDofs = bcs->somBC[0]->numDofs();
   int *somP = (int *) dbg_alloca(sizeof(int)*somDofs);
   DComplex *buffSom = (DComplex *) dbg_alloca(sizeof(DComplex)*somDofs*somDofs);

   for (i=0; i<bcs->numSomm; ++i) {
     FullSquareMatrixC ks = bcs->somBC[i]->turkelMatrix(noeuds,kappa,mode,
                            buffSom);
     bcs->somBC[i]->dofs(*tab, somP);
     ks *= -a1;
     K->add(ks,somP);
     if (kuc) kuc->add(ks,somP);
   }

 }

}


void FourierStatic::buildRHS(ComplexVector &force, int mode) {

  int i;
  ConstrainedDSA *c_tab = domain->getCDSA();
  DofSetArray *tab = domain->getDSA();
  CoordSet &noeuds = domain->getNodes();

  bcxC = make_bcxC(mode);

  // ... COMPUTE EXTERNAL FORCE FROM COMPLEX NEUMAN BC
  for (i=0; i < bcs->numNeu; ++i) {

   double *x = (double *) dbg_alloca(3*sizeof(double));
   double *y = (double *) dbg_alloca(3*sizeof(double));

   int dof1 = c_tab->locate(bcs->neuBC[i].nnum1, (1 << 7));
   if (dof1 < 0)
      continue;
   Node nd = noeuds.getNode(bcs->neuBC[i].nnum1);
   x[0] = nd.x;
   y[0] = nd.y;

   int dof2 = c_tab->locate(bcs->neuBC[i].nnum2, (1 << 7));
   if (dof2 < 0)
      continue;
   nd = noeuds.getNode(bcs->neuBC[i].nnum2);
   x[1] = nd.x;
   y[1] = nd.y;

   int dof3 = (bcs->neuBC[i].nnum3==-1) ? -1 :
               c_tab->locate(bcs->neuBC[i].nnum3, (1 << 7));
   if (dof3>=0) {
     nd = noeuds.getNode(bcs->neuBC[i].nnum3);
     x[2] = nd.x;
     y[2] = nd.y;
   }

   double nx = 0.0;
   double ny = 0.0;

   domain->getNormal2D(bcs->neuBC[i].nnum1,bcs->neuBC[i].nnum2,nx,ny);

   DComplex *f = (DComplex*) dbg_alloca(sizeof(DComplex)*3);
   f[0] = DComplex(0.0, 0.0);
   f[1] = DComplex(0.0, 0.0);
   f[2] = DComplex(0.0, 0.0);

   if (dof3==-1) {
     bcs->IntNeu2(mode,i,hParams->kappaHelm,x,y,nx,ny,f);
     if(dof1 >= 0) force[dof1] += f[0];
     if(dof2 >= 0) force[dof2] += f[1];
   }
   else {
     bcs->IntNeu3(mode,i,hParams->kappaHelm,x,y,nx,ny,f);
     if(dof1 >= 0) force[dof1] += f[0];
     if(dof2 >= 0) force[dof2] += f[1];
     if(dof3 >= 0) force[dof3] += f[2];
   }

  }

  ComplexVector Vc(bcs->numDir, 0.0);

  // CONSTRUCT THE NON-HOMOGENEOUS COMPLEX DIRICHLET BC VECTOR
  for (i=0; i< bcs->numDir; ++i) {
    int dof = tab->locate(bcs->dirBC[i].nnum, (1 << 7));
    if (dof < 0) continue;
    int dof1 = c_tab->invRCN(dof);
    if (dof1 >= 0)
       Vc[dof1] = bcxC[dof];
  }

  // ... PERFORM MULTIPLICATION TO GET NON-HOMOGENEOUS FORCE:
  //                    Fnh -= [Kuc]*[Vc]
  if (kuc) kuc->multSubtract(Vc,force);

}


ComplexVector FourierStatic::MergeSolBCs(ComplexVector sol) {

  int iNode;
  int maxNode = domain->numNode();
  DofSetArray *tab = domain->getDSA();
  ConstrainedDSA *c_tab = domain->getCDSA();
  ComplexVector u(maxNode,0.0);

  for (iNode=0; iNode<maxNode; ++iNode) {
     int loc = c_tab->locate(iNode, (1 << 7));
     int loc1= tab->locate(iNode, (1 << 7));
     if (loc >= 0)
          u[iNode] = sol[loc];
     else if (loc1 >= 0)
             u[iNode] = bcxC[loc1];
     else
             u[iNode] = 0.0;
  }

  return u;
}


void FourierStatic::postProcessing() {

 // Print the cpu times

}


void FourierStatic::Reconstruction3D(ComplexVector *globalsol) {

   double pi = 4*atan(1.0);
   double *c = new double[bcs->numSlices];
   double *s = new double[bcs->numSlices];
   double *coeff = new double[bcs->numModes+1];

   DComplex result;
   int i, j, k;
   int iInfo, number;
   int maxNode = domain->numNode();

   for (i=0; i<bcs->numSlices; ++i) {
      double angle = 2*pi*i/bcs->numSlices;
      c[i] = cos(angle);
      s[i] = sin(angle);
   }

   fprintf(stderr," ... Output Complex Solution        ...\n");

   //domain->openOutputFiles();
   geoSource->openOutputFiles();

// ... PRINT THE REAL PART OF HELMHOLTZ SOLUTION

   //number = domain->numOutputFiles();
   number = geoSource->getNumOutInfo();
   //ResizeArray <OutputInfo> &sortie = domain->getOINFO();
   OutputInfo *sortie = geoSource->getOutputInfo();

   for (iInfo=0; iInfo<number; ++iInfo) {
       //int w = sortie[iInfo].width;
       //int p = sortie[iInfo].precision;

       fprintf(sortie[iInfo].filptr,"%d\n",bcs->numSlices*maxNode);
       fprintf(sortie[iInfo].filptr,"0.0000000000000\n");

       for (i=0; i<bcs->numSlices; ++i) {
           for (j=0; j<maxNode; ++j) {
               result = globalsol[0][j];
               for (k=1; k<=bcs->numModes; ++k) {
                  int angle = (k*i)%bcs->numSlices;
                  result = result + c[angle]*globalsol[2*k-1][j];
                  result = result + s[angle]*globalsol[2*k][j];
               }

//               fprintf(sortie[iInfo].filptr,"% *.*e\n",w,p,real(result) );
               fprintf(sortie[iInfo].filptr,"% 1.15e\n",real(result) );
           }
        }

        for (k=0; k<=bcs->numModes; ++k) {
           coeff[k] = 0;
           for (j=0; j<maxNode; ++j) {
             coeff[k] += real(globalsol[2*k][j])*real(globalsol[2*k][j]);
             coeff[k] += imag(globalsol[2*k][j])*imag(globalsol[2*k][j]);
             if (k>0) {
              coeff[k] += real(globalsol[2*k-1][j])*real(globalsol[2*k-1][j]);
              coeff[k] += imag(globalsol[2*k-1][j])*imag(globalsol[2*k-1][j]);
             }
           }
        }

        geoSource->outputHeader(iInfo);

        fprintf(sortie[iInfo].filptr,"%d\n",bcs->numSlices*maxNode);
        fprintf(sortie[iInfo].filptr,"0.0000000000000\n");
        for (i=0; i<bcs->numSlices; ++i) {
           for (j=0; j<maxNode; ++j) {
               result = globalsol[0][j];
               for (k=1; k<=bcs->numModes; ++k) {
                  int angle = (k*i)%bcs->numSlices;
                  result = result + c[angle]*globalsol[2*k-1][j];
                  result = result + s[angle]*globalsol[2*k][j];
               }
//               fprintf(sortie[iInfo].filptr,"% *.*e\n",w,p,imag(result) );
               fprintf(sortie[iInfo].filptr,"% 1.15e\n",imag(result) );
           }
        }

        fflush(sortie[iInfo].filptr);
   }

   for (iInfo=0; iInfo<number; ++iInfo) {
      if ((sortie[iInfo].interval == 1) &&
          (sortie[iInfo].type == (OutputInfo::Helmholtz))) {
        fclose(sortie[iInfo].filptr);
      }
   }
}

