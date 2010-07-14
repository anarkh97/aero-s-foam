/*
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/dofset.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Math.d/matrix.h>
#include <Utils.d/Connectivity.h>
#include <Threads.d/Paral.h>
#include <Solvers.d/Rbm.h>
#include <Threads.d/Paral.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
// #include <Utils.d/MathUtils.h>

//#define TESTRECTCG
//#define TEST_R2


extern "C" {
 extern void _FORTRAN(zgeqrf)(const int &m, const int &n, ComplexD* a,
                        const int& lda, ComplexD *tau,
                        ComplexD *work, const int &lwork, int &info);
}

extern "C" {
extern void _FORTRAN(zgemv)(const char &, const int &, const int &, const ComplexD &,
                 ComplexD *, const int &, ComplexD *, const int &,
                 const ComplexD &, ComplexD *, const int &);
}

extern "C" {
extern void _FORTRAN(dsyevx)(const char &,const char &,const char &, const int &,
                 double *, int &, double &, double &, int &, int &, double &,
                 int &, double *, double *, int &, double *, int &, int *,
                 int *, int &);
}

extern "C" {
  void _FORTRAN(zgemm)(const char &, const char &, const int &,const int &,
                  const int &, const ComplexD &, ComplexD *, const int &,
                  ComplexD *, const int &, const ComplexD &, ComplexD *,
                  const int &);
}

template<class Scalar>
void
GenSubDomain<Scalar>::updateLocalMatrices(GenSparseMatrix<Scalar> *K, int *dofs,
                                          FullSquareMatrix *reEl, FullSquareMatrix *imEl)
{
 if(reEl) {
   if(K) K->add(*reEl, dofs);
   if(Kcc) Kcc->add(*reEl, dofs);
   if(Krc) Krc->add(*reEl, dofs);
   if(Kuc) Kuc->add(*reEl, dofs);
   if(Kbb) Kbb->add(*reEl, dofs);
   if(Kib) Kib->add(*reEl, dofs);
   if(KiiSparse) KiiSparse->add(*reEl, dofs);
 }
 if(imEl) {
   if(K) K->addImaginary(*imEl, dofs);
   if(Kcc) Kcc->addImaginary(*imEl, dofs);
   if(Krc) Krc->addImaginary(*imEl, dofs);
   if(Kuc) Kuc->addImaginary(*imEl, dofs);
   if(Kbb) Kbb->addImaginary(*imEl, dofs);
   if(Kib) Kib->addImaginary(*imEl, dofs);
   if(KiiSparse) KiiSparse->addImaginary(*imEl, dofs);
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::updateLocalDampingMatrices(int *dofs,
  FullSquareMatrix *reEl, FullSquareMatrix *imEl, double ss, int n)
{
 if(reEl) {
// RT: 11/17/2008
   *reEl /= pow(ss,n); 
   if(C_deriv && C_deriv[n-1]) C_deriv[n-1]->add(*reEl, dofs);
   if(Cuc_deriv && Cuc_deriv[n-1]) Cuc_deriv[n-1]->add(*reEl, dofs);
 }
 if(imEl) {
 // RT: 11/17/2008
   *imEl /= pow(ss,n);
   if(C_deriv && C_deriv[n-1]) C_deriv[n-1]->addImaginary(*imEl, dofs);
   if(Cuc_deriv && Cuc_deriv[n-1]) Cuc_deriv[n-1]->addImaginary(*imEl, dofs);
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleLocalSommer(GenSparseMatrix<Scalar> *K)
{
 // this is the same function as Domain::assembleGlobalSommer but calls updateLocalMatrices 
 // instead of updateGlobalMatrices
 if(numSommer > 0) {
   //filePrint(stderr," ... Assemble Sommer els, sub %2d, numSommer = %d\n",subNumber,numSommer); //HB
   if(sommer[0]->dim() == 3) {
     if(sommerfeldType==1) sommerfeldType = 3;
     if(sommerfeldType==2) sommerfeldType = 4;
   }

   if(subNumber == 0) {
     if(sommerfeldType == 1)
       fprintf(stderr, " with 1st order Bayliss-Turkel");
     else if(sommerfeldType == 2)
       fprintf(stderr, " with 2nd order Bayliss-Turkel");
     else if(sommerfeldType == 3)
       fprintf(stderr, " with 1st order 3D Bayliss-Turkel");
     else if(sommerfeldType == 4)
       fprintf(stderr, " with 2nd order 3D Bayliss-Turkel");
   }

   if(sommerfeldType == 1 || sommerfeldType == 2)
     getCurvatures(this);
   if(curvatureFlag != 1) {
     if(sommerfeldType == 3 || sommerfeldType == 4)
       getCurvatures3D(this);
   }

   double *v = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
   double *vbt = 0;
   // This loops adds the contribution of the terms emanating
   // from the non-reflecting boundary conditions
   int i;

   for(i=0; i<numSommer; i++) {
     ComplexD *bt2Matrix = 0;
     ComplexD **bt2nMatrix = 0;
     int *dofs = sommer[i]->dofs(*dsa);
     FullSquareMatrix ms = sommer[i]->sommerMatrix(nodes,v);
     FullSquareMatrix mm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
     double kappa = sommer[i]->el->getProperty()->kappaHelm; // PJSA 1-15-2008

     ms.multiply(mm,kappa);  // mm = kappa*ms

     // "Zero-order" term
     updateLocalMatrices(K,dofs,0,&mm);

     double psi; // curvature of the boundary
     double HH, KK;
     // 1st order Bayliss-Turkel boundary condition
     if(sommerfeldType == 1 ) {
       psi = curvatures[i];
       HH = psi/2.0;
       ms.multiply(mm,-HH);
       updateLocalMatrices(K,dofs,&mm,0);
     }

     // 2nd order Bayliss-Turkel boundary condition
     else if(sommerfeldType == 2) {
       vbt = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
       FullSquareMatrix ks = sommer[i]->turkelMatrix(nodes,vbt);
       psi = curvatures[i];
       // fprintf(stderr,"curvature %f\n",psi);

       ComplexD cm = ComplexD(-0.5*psi,0.0) + psi*psi/8.0/ComplexD(psi,-kappa);
       ComplexD cs = -0.5/ComplexD(psi,-kappa);
       double c1,c2,c3,c4;
       c1 = imag(cm);
       c2 = real(cm);
       c3 = imag(cs);
       c4 = real(cs);

       FullSquareMatrix mm1(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
       // mm = ms*c2 + ks*c4;
       ms.multiply(mm,c2);
       ks.multiply(mm1,c4);
       mm += mm1;
       updateLocalMatrices(K,dofs,&mm,0);
       // mm = ms*c1 + ks*c3;
       ms.multiply(mm,c1);
       ks.multiply(mm1,c3);
       mm += mm1;
       updateLocalMatrices(K,dofs,0,&mm);
     }

     // 1st order 3D Bayliss-Turkel boundary condition
     else if(sommerfeldType == 3) {
       if(curvatureFlag != 1) {
         HH = 0.0;
         for(int iNode=0;iNode<somElemToNode->num(i);iNode++) {
           int iNodeNumber = (*somElemToNode)[i][iNode];
           int iSommerNode = nodeToSommerNodeMap[iNodeNumber];
           HH += curvaturesH[iSommerNode];
         }
         HH /= somElemToNode->num(i);
       }
       else {
         HH = 1.0/curvatureConst1;
       }
       ms.multiply(mm,-HH);
       updateLocalMatrices(K,dofs,&mm,0);
     }

     // 2nd order 3D Bayliss-Turkel boundary condition
     else if(sommerfeldType == 4) {
       HH = 0.0;
       KK = 0.0;
       int iNode;
       if(curvatureFlag != 1) {
         for(iNode=0;iNode<somElemToNode->num(i);iNode++) {
           int iNodeNumber = (*somElemToNode)[i][iNode];
           int iSommerNode = nodeToSommerNodeMap[iNodeNumber];
           HH += curvaturesH[iSommerNode];
           KK += curvaturesK[iSommerNode];
         }
         HH /= somElemToNode->num(i);
         KK /= somElemToNode->num(i);
       }
       else {
         HH = 0.5*(1.0/curvatureConst1+1.0/curvatureConst1);
         KK = 1.0/curvatureConst1*1.0/curvatureConst1;
       }
       ComplexD ii=ComplexD(0.0, 1.0);
       ComplexD cm = -ii/2.0/kappa/(1.0+ii*2.0*HH/kappa)*(KK-HH*HH) - ComplexD(HH,0.0);

       double c1,c2;
       c1 = imag(cm);
       c2 = real(cm);

       if(curvatureFlag != 2) {
         ms.multiply(mm,c2);
         updateLocalMatrices(K,dofs,&mm,0);
         ms.multiply(mm,c1);
         updateLocalMatrices(K,dofs,0,&mm);
       }

       bt2Matrix = new DComplex [ms.dim()*ms.dim()*sizeof(DComplex)];
       for(int j=0; j<ms.dim()*ms.dim(); ++j) bt2Matrix[j] = 0.0;  // PJSA
       if(solInfo().doFreqSweep) {
         int N = solInfo().nFreqSweepRHS - 1; // number of derivatives
         bt2nMatrix = new DComplex * [N];
         for(int n=1; n<=N; ++n) {
           bt2nMatrix[n-1] = new DComplex [ms.dim()*ms.dim()*sizeof(DComplex)];
           for(int j=0; j<ms.dim()*ms.dim(); ++j) bt2nMatrix[n-1][j] = 0.0;  // PJSA
         }
       }
       FullSquareMatrix ksRe(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
       FullSquareMatrix ksIm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));

       if(curvatureFlag == 2) {
         double HHH[3];
         double KKK[3];
         for(iNode=0;iNode<somElemToNode->num(i);iNode++) {
           int iNodeNumber = (*somElemToNode)[i][iNode];
           int iSommerNode = nodeToSommerNodeMap[iNodeNumber];
           HHH[iNode] = curvaturesH[iSommerNode];
           KKK[iNode] = curvaturesK[iSommerNode];
         }

         sommer[i]->sommerMatrixEllipsoid(nodes, kappa, HHH, KKK, bt2Matrix);

         int kDof=0;
         int iDof, jDof;
         for(iDof=0; iDof<ms.dim(); iDof++) {
           for(jDof=0; jDof<ms.dim(); jDof++) {
             ksRe[iDof][jDof] = real(bt2Matrix[kDof]);
             ksIm[iDof][jDof] = imag(bt2Matrix[kDof]);
             kDof++;
           }
         }

         updateLocalMatrices(K,dofs,&ksRe,&ksIm);
       }

       double curv_e[3];
       double curv_f[3];
       double curv_g[3];
       double tau1[3][3];
       double tau2[3][3];
       if(curvatureFlag != 1) {
         int n1 = (*somElemToNode)[i][0];
         int n2 = (*somElemToNode)[i][1];
         int n3 = (*somElemToNode)[i][2];
         int nn[3];
         nn[0] = nodeToSommerNodeMap[n1];
         nn[1] = nodeToSommerNodeMap[n2];
         nn[2] = nodeToSommerNodeMap[n3];
         double *curv_normal[3];
         int iiNode;
         for(iiNode=0; iiNode<3; iiNode++) {
           int curvNode = nn[iiNode];
           curv_e[iiNode]      = curvatures_e[curvNode];
           curv_f[iiNode]      = curvatures_f[curvNode];
           curv_g[iiNode]      = curvatures_g[curvNode];
           curv_normal[iiNode] = curvatures_normal[nn[iiNode]];
           getTau(curv_normal[iiNode], tau1[iiNode], tau2[iiNode]);

           // This is for experiments with exact curvatures for ellipsoid
           if(curvatureFlag == 2) {
             tau1[iiNode][0] = curvatures_tau1[curvNode][0];
             tau1[iiNode][1] = curvatures_tau1[curvNode][1];
             tau1[iiNode][2] = curvatures_tau1[curvNode][2];
             tau2[iiNode][0] = curvatures_tau2[curvNode][0];
             tau2[iiNode][1] = curvatures_tau2[curvNode][1];
             tau2[iiNode][2] = curvatures_tau2[curvNode][2];
           }
         }
       }
       else {
         curv_e[0] = curv_e[1] = curv_e[2] = 1.0/curvatureConst1;
         curv_g[0] = curv_g[1] = curv_g[2] = 1.0/curvatureConst1;
         curv_f[0] = curv_f[1] = curv_f[2] = 0.0;
         int nds[3],iiNode;
         sommer[i]->nodes(nds);
         for(iiNode=0; iiNode<3; iiNode++) {
           double nrmal[3];
           Node nd = nodes.getNode(nds[iiNode]);
           double l = sqrt(nd.x*nd.x+nd.y*nd.y+nd.z*nd.z);
           nrmal[0] = nd.x/l;
           nrmal[1] = nd.y/l;
           nrmal[2] = nd.z/l;
           getTau(nrmal, tau1[iiNode], tau2[iiNode]);
         }
       }
       if(curvatureFlag != 1) {
         sommer[i]->BT2(nodes, curv_e, curv_f, curv_g, tau1, tau2, kappa, bt2Matrix);
         if(solInfo().doFreqSweep) {
           int N = solInfo().nFreqSweepRHS - 1; // number of derivatives
           for(int n=1; n<=N; ++n) 
             sommer[i]->BT2n(nodes, curv_e, curv_f, curv_g, tau1, tau2, kappa, bt2nMatrix[n-1], n);
         }
       }
       else {
         sommer[i]->sphereBT2(nodes, curvatureConst1 , kappa, bt2Matrix);
       }
       if(curvatureFlag == 2) {
         sommer[i]->ellipsoidBT2(nodes, curvatureConst2, curvatureConst1, kappa, bt2Matrix);
       }

       int kDof=0;
       int iDof, jDof;
       for(iDof=0; iDof<ms.dim(); iDof++) {
         for(jDof=0; jDof<ms.dim(); jDof++) {
           ksRe[iDof][jDof] = real(-1.0/(2.0*ii*kappa)*bt2Matrix[kDof]);
           ksIm[iDof][jDof] = imag(-1.0/(2.0*ii*kappa)*bt2Matrix[kDof]);
           kDof++;
         }
       }

       updateLocalMatrices(K,dofs,&ksRe,&ksIm);
     }

     if(solInfo().doFreqSweep) {
       switch(sommerfeldType) {
         case 0:
         case 1:
         case 3:
           updateLocalDampingMatrices(dofs,0,&ms,
             real(sommer[i]->el->getProperty()->soundSpeed),1);
                // zero- and first-order sommerfeld
           break;
         case 4:
           computeSommerDerivatives(HH, KK, curvatureFlag, dofs, ms, bt2nMatrix, kappa); // 3D second-order sommerfeld
           break;
         case 2:
         default:
           cerr << " *** ERROR: Sommerfeld type " << sommerfeldType << " is not supported for frequency sweep \n";
           break;
       }
     } 
     if(bt2Matrix) delete [] bt2Matrix;
     if(bt2nMatrix) {
       int N = solInfo().nFreqSweepRHS - 1; // number of derivatives
       for(int n=1; n<=N; ++n) if(bt2nMatrix[n-1]) delete [] bt2nMatrix[n-1];
       delete [] bt2nMatrix;
     }
     delete [] dofs;
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::dumpMiscData(int handle) 
{
  // Dump local node to global node mapping
  write(handle,&numnodes,sizeof(int));        // number of nodes in this subdomain
  write(handle,glNums,numnodes*sizeof(int));  // global number of each node in this subdomain

  int *dofs_array = new int[numnodes*4];

  int i;
  for(i=0;i<numnodes;i++) {
    int hLoc  = c_dsa->locate(i, DofSet::Helm);
    int xLoc  = c_dsa->locate(i, DofSet::Xdisp);
    int yLoc  = c_dsa->locate(i, DofSet::Ydisp);
    int zLoc  = c_dsa->locate(i, DofSet::Zdisp);
    dofs_array[4*i] = hLoc;
    dofs_array[4*i+1] = xLoc;
    dofs_array[4*i+2] = yLoc;
    dofs_array[4*i+3] = zLoc;
  }

  // degree of freedom number for each node in this subdomain
  // for each node four integers are written, for pressure, x-displacement,
  // y-displacement, and z-displacement, -1 indicates that the degree of freedom
  // is not present
  write(handle,dofs_array,numnodes*4*sizeof(int));

  delete[] dofs_array;
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeSommerDerivatives(double HH, double KK, int curvatureFlag, int *dofs, FullSquareMatrix &ms,
                                               DComplex **bt2nMatrix, double kappa)
{
  // PJSA 5-26-05
  // this function is for 3D second-order sommerfeld
  FullSquareMatrix mm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
  ComplexD ii=ComplexD(0.0, 1.0);
  int N = solInfo().nFreqSweepRHS - 1;  // number of derivatives
  for(int n=1; n<=N; ++n) {
    DComplex cm = (n==1) ? ii : 0;
    cm -= ((KK-HH*HH)*ii/2.0 * pow(-1.0,n)*double(DFactorial(n))/pow((kappa+ii*2.0*HH),n+1));
    if(curvatureFlag != 2) {
      ms.multiply(mm,real(cm)); // mm = real part of cm*ms
      updateLocalDampingMatrices(dofs,&mm,0,
         sqrt(geoSource->shiftVal())/kappa,n);
      ms.multiply(mm,imag(cm)); // mm = imaginary part of cm*ms
      updateLocalDampingMatrices(dofs,0,&mm,
          sqrt(geoSource->shiftVal())/kappa,n);
    }
    else cerr << " *** WARNING: 3D 2nd order Sommerfeld with curvatureFlag 2 is not supported for frequency sweep \n";
  }

// PJSA DEBUG    
//  GenFullSquareMatrix<DComplex> bt2(ms.dim(),(DComplex*)dbg_alloca(ms.dim()*ms.dim()*sizeof(DComplex)));
//  int kDof=0;
//  for(int iDof=0; iDof<ms.dim(); iDof++) {
//    for(int jDof=0; jDof<ms.dim(); jDof++) {
//      bt2[iDof][jDof] = bt2Matrix[kDof];
//      kDof++;
//    }
//  }
//  FullSquareMatrix ksRe(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
//  FullSquareMatrix ksIm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
//
//  GenFullSquareMatrix<DComplex> bt2n(bt2,1.0);
//  for(int n=1; n<=N; ++n) {
//    GenFullSquareMatrix<DComplex> copy(bt2n, 1.0);
//    copy.multiply(bt2,bt2n);  // bt2n = bt2n*bt2
//    DComplex cm = -pow(1.0/(2.0*ii*kappa),n+1)*pow(-2.0*ii,n)*double(DFactorial(n));
//    for(int iDof=0; iDof<ms.dim(); iDof++) {
//      for(int jDof=0; jDof<ms.dim(); jDof++) {
//        DComplex cmbt = cm*bt2n[iDof][jDof];
//        ksRe[iDof][jDof] = real(cmbt);
//        ksIm[iDof][jDof] = imag(cmbt);
//      }
//    }
//    updateLocalDampingMatrices(dofs,&ksRe,&ksIm,
//        sqrt(geoSource->shiftVal())/kappa,n);
//  }
//
  FullSquareMatrix ksRe(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
  FullSquareMatrix ksIm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
  for(int n=1; n<=N; ++n) {
    DComplex cm = -pow(1.0/(2.0*ii*kappa),n+1)*pow(-2.0*ii,n)*double(DFactorial(n));
    int kDof = 0;
    for(int iDof=0; iDof<ms.dim(); iDof++) {
      for(int jDof=0; jDof<ms.dim(); jDof++) {
        DComplex cmbt = cm*bt2nMatrix[n-1][kDof];
        ksRe[iDof][jDof] = real(cmbt);
        ksIm[iDof][jDof] = imag(cmbt);
        kDof++;
      }
    }
    updateLocalDampingMatrices(dofs,&ksRe,&ksIm,
         sqrt(geoSource->shiftVal())/kappa,n);
  }
}
*/

template<class Scalar>
void
GenSubDomain<Scalar>::addSommer(SommerElement *ele)
{
 ele->dom = this;
 sommer[numSommer++] = ele;
 //if(sinfo.ATDARBFlag != -2.0) packedEset.elemadd(numele++,ele);  // XDEBUG
 if(sinfo.ATDARBFlag != -2.0)
  { ele->renum(glToLocalNode); packedEset.elemadd(numele++,ele); }
}

