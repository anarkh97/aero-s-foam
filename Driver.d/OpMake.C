#include <typeinfo>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#ifdef USE_EIGEN3
#include <Math.d/EiSparseMatrix.h>
#endif
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/SGISparseMatrix.h>
#include <Math.d/Skyline.d/BlockSky.h>
#include <Math.d/Skyline.d/SGISky.h>
#include <Solvers.d/PCGSolver.h>
#include <Solvers.d/UFront.h>
#include <Solvers.d/CRSolver.h>
#include <Solvers.d/BCGSolver.h>
#include <Solvers.d/GmresSolver.h>
#include <Solvers.d/Spooles.h>
#include <Solvers.d/Mumps.h>
#include <Solvers.d/GoldfarbIdnani.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Driver.d/GeoSource.h>
#include <Solvers.d/MappedAssembledSolver.h>
#include <Driver.d/Dynam.h>
#include <Sfem.d/Sfem.h>
#include <Sfem.d/SfemBlockMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/DiagMatrix.h>
#include <Feti.d/DistrVector.h>
#include <Hetero.d/FlExchange.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/MpcElement.d/MpcElement.h>
#include <Utils.d/MFTT.h>
#include <Utils.d/MathUtils.h>
#include <Element.d/Sommerfeld.d/TrianglePressureBC.h>
#include <Element.d/Sommerfeld.d/QuadPressureBC.h>

#include <Rom.d/GaussNewtonSolver.h>
#include <Rom.d/GappyProjectionSolver.h>
#include <Rom.d/GalerkinProjectionSolver.h>

extern Sfem* sfem;
extern int verboseFlag;
#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

/** abstract method of operator assembly */
template<class Scalar, class OpList>
void
Domain::assembleSparseOps(OpList &ops) {
  int size = sizeof(Scalar)*maxNumDOFs*maxNumDOFs;
  Scalar *karray = (Scalar *) dbg_alloca(size);
  Scalar *marray = (Scalar *) dbg_alloca(size);
  for(int iele=0; iele < numele; ++iele) {
    // Skip all irrelevant elements
    if(packedEset[iele]->isSommerElement()) continue;

    if(ops.needStiff()) {
      GenFullSquareMatrix<Scalar> kel = packedEset[iele]->stiffness(nodes, karray);
      ops.addStiffness((*allDOFs)[iele], kel);
    }
    if(ops.needMass()) {
      GenFullSquareMatrix<Scalar> kel = packedEset[iele]->massMatrix(nodes, karray);
      ops.addMass((*allDOFs)[iele], kel);
    }
  }

  if(ops.needMass()) {
    // Add discrete mass contribution to the Mass Matrix
    // Three matrices need to be changed. Mass matrix itself,
    // Damping Matrix and K tilda.
    Scalar m;
    DMassData *current = firstDiMass;
    while(current != 0) {
      int dof = dsa->locate(current->node, (1 << current->dof));
      int jdof = (current->jdof > -1) ?// PJSA 10-9-06 for off-diagonal mass terms eg. products of inertia I21, I31, I32
         dsa->locate(current->node, (1 << current->jdof)) : dof;

      ops.addMass(dof, jdof, current->diMass);
      current = current->next;
    }
  }
}

template<class Scalar>
void
Domain::makeSparseOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
		      double Ccoef, GenSparseMatrix<Scalar> *mat,
                      FullSquareMatrix *kelArray, FullSquareMatrix *melArray)
{
 if(matrixTimers) matrixTimers->memoryForm -= memoryUsed();
 int makeMass = Mcoef != 0 || ops.M != 0 || ops.C != 0;

 // Rayleigh damping coefficients: C = alpha*M + beta*K
 double alphaDamp = sinfo.alphaDamp, alpha;
 double  betaDamp = sinfo.betaDamp, beta;

 int size = sizeof(double)*maxNumDOFs*maxNumDOFs;

 double *karray = new double[maxNumDOFs*maxNumDOFs];
 double *marray = new double[maxNumDOFs*maxNumDOFs];

 FullSquareMatrix kel, mel, cel;

 double mratio = geoSource->getMRatio(); // 0 for lumped, 1 for consistent mass matrix

 complex<double> *kcarray = new complex<double>[4*maxNumDOFs*maxNumDOFs];
 complex<double> *mcarray = new complex<double>[4*maxNumDOFs*maxNumDOFs];

 FullSquareMatrixC kcel(1, kcarray);
 FullSquareMatrixC mcel(1, mcarray);

 bool isShifted = geoSource->isShifted();
 bool isDamped = (alphaDamp != 0.0) || (betaDamp != 0.0) || packedEset.hasDamping();

 double *izarray = (isShifted || isDamped) ? static_cast<double *>(dbg_alloca(size)) : 0;
 FullSquareMatrix izel(1, izarray); // izel stores the imaginary part of impedence matrix zel
 double omega, omega2;

 int sizeC = sizeof(DComplex)*maxNumDOFs*maxNumDOFs;
 DComplex *karrayC = new DComplex[maxNumDOFs*maxNumDOFs];
 DComplex *marrayC = new DComplex[maxNumDOFs*maxNumDOFs];
 FullSquareMatrixC kelC(1, karrayC);
 FullSquareMatrixC melC(1, marrayC);

 if(sinfo.isCoupled) computeCoupledScaleFactors();

 GenSubDomain<Scalar> *subCast = dynamic_cast<GenSubDomain<Scalar>*>(this);
 if(!subCast) 
   if(sinfo.ATDARBFlag >= 0.0 || sinfo.ATDROBalpha != 0.0) checkSommerTypeBC(this);

 bool zeroRot = (sinfo.zeroRot && sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0);
 int *dofType = (zeroRot) ? dsa->makeDofTypeArray() : 0;

 if(sinfo.farfield) { addSBoundNodes(); makeKss(this); } // for Farfield output (TODO check with Radek)

 int iele;
 int *m_dofs;

 // LOOP over elements except for sommer
 for(iele = 0; iele < numele; ++iele) {
   StructProp *prop = packedEset[iele]->getProperty();
   if(packedEset[iele]->isSommerElement()) continue;
   bool isComplexF = (prop && prop->fp.PMLtype != 0);
   alpha = (packedEset[iele]->isDamped()) ? prop->alphaDamp : alphaDamp;
   beta = (packedEset[iele]->isDamped()) ? prop->betaDamp : betaDamp;
   m_dofs = (sinfo.type == 0 && subCast) ? (*domain->getAllDOFs())[subCast->getGlElems()[iele]] : (*allDOFs)[iele];
   complex<double> kappa2 = packedEset[iele]->helmCoefC();
   omega2 = geoSource->shiftVal();

   // Form element real and complex stiffness matrices in kel and kcel
   if(matrixTimers) matrixTimers->formTime -= getTime();
   if(isComplexF) {
     kcel = packedEset[iele]->stiffness(nodes, kcarray);
   }
   else {
     if(kelArray) kel.copy(kelArray[iele]);
     else kel = packedEset[iele]->stiffness(nodes, karray);
     this->densProjectStiffness(kel, iele);
   }
   if(sinfo.isCoupled) {
      if(isStructureElement(iele)) kel *= cscale_factor2;
      else if(packedEset[iele]->isFsiElement()) kel *= cscale_factor;
   }
   if(matrixTimers) matrixTimers->formTime += getTime();

   // Assemble element real and complex stiffness matrices in ops.K
   if(matrixTimers) matrixTimers->assemble -= getTime();
   if(isComplexF || (imag(kappa2) != 0)) {
     if(ops.K) ops.K->add(kcel,(*allDOFs)[iele]);
   }
   else {
     if(ops.K) ops.K->add(kel,(*allDOFs)[iele]);
     if(!isShifted && ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]);
     if(!isShifted && ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);
     if(packedEset[iele]->isConstraintElement()) { // XXXX
       if(sinfo.isNonLin() && Mcoef == 1 && Kcoef == 0 && Ccoef == 0 && sinfo.newmarkBeta != 0) {
         kel = packedEset[iele]->stiffness(nodes, karray);
         if(mat) mat->add(kel,(*allDOFs)[iele]);
       }
       else {
         if(ops.Msolver) ops.Msolver->add(kel,(*allDOFs)[iele]);
       }
     }
   }
   if(matrixTimers) matrixTimers->assemble += getTime();

   int i, j, dim = kel.dim();
   // Form element real and complex mass matrices in mel and mcel
   if(matrixTimers) matrixTimers->formTime -= getTime();
   if(makeMass || isShifted) {
     if(isComplexF) {
       mcel = packedEset[iele]->massMatrix(nodes, mcarray);
     }
     else {
       if(melArray) mel.copy(melArray[iele]);
       else mel = packedEset[iele]->massMatrix(nodes, marray, mratio);
       this->densProjectStiffness(mel, iele);
     }
     if(sinfo.isCoupled) { 
       if(isStructureElement(iele)) mel *= cscale_factor2;
       else if(packedEset[iele]->isFsiElement()) mel *= coupledScaling;
     }
   }

   if(isShifted) {
     if(isDamped && isStructureElement(iele)) {
       omega = sqrt(omega2);
       if(packedEset[iele]->hasDamping()) {
	 izel = packedEset[iele]->dampingMatrix(nodes, izarray);
	 this->densProjectStiffness(izel, iele);
	 for(i = 0; i < dim; ++i)
	   for(j = 0; j < dim; ++j)
             izel[i][j] = -omega*(izel[i][j]+beta*kel[i][j] + alpha*mel[i][j]);
       }
       else if(isDamped) {
	 izel.setSize(mel.dim());
	 for(i = 0; i < dim; ++i)
	   for(j = 0; j < dim; ++j)
             izel[i][j] = -omega*(beta*kel[i][j] + alpha*mel[i][j]);
       }
     }
     if(isComplexF) {
       int dim = kcel.dim();
       for(i = 0; i < dim; ++i)
         for(j = 0; j < dim; ++j)
           kcel[i][j] -= kappa2*mcel[i][j];
     } 
     else if(imag(kappa2) != 0.0) {
       kcel.setSize(kel.dim());
       for(i = 0; i < dim; ++i)
         for(j = 0; j < dim; ++j)
           kcel[i][j] = complex<double>(kel[i][j])-kappa2*mel[i][j];
     }
     else {
       double o2 = (isFluidElement(iele)) ? real(kappa2) : omega2;
       for(i = 0; i < dim; ++i)
         for(j = 0; j < dim; ++j) {
           if(isDamped && isStructureElement(iele))
             izel[i][j] = -omega*(beta*kel[i][j] + alpha*mel[i][j]);
           kel[i][j] -= o2*mel[i][j];
         }
     }
     if(sinfo.doFreqSweep && isFluidElement(iele)) {
       if(isComplexF) {
         mcel *= kappa2/omega2;
       }
       else if (imag(kappa2) != 0) {
         mcel.setSize(mel.dim());
         for(i = 0; i < dim; ++i)
           for(j = 0; j < dim; ++j)
             mcel[i][j] = kappa2/omega2*mel[i][j];
       }
       else {
         mel *= real(kappa2)/omega2;
       }
     }
   }
   else {
     // do something for coupled eigen/dynamics, set fluidCelerity = 1500
     if(makeMass && sinfo.isCoupled && isFluidElement(iele)) mel /= (1500.0 * 1500);
     if(sinfo.acoustic == 1 && prop) {
       double c = prop->ss; //speed of sound
       mel /= (c*c);
     }
   }
   if(matrixTimers) matrixTimers->formTime += getTime();

   // Assemble element real and complex mass matrices in ops.M, ops.Muc, ops.Mcc and ops.Msolver
   if(matrixTimers) matrixTimers->assemble -= getTime();
   if(isComplexF || (imag(kappa2) != 0)) {
     if(ops.M)   ops.M->add(mcel,(*allDOFs)[iele]);
     if(ops.Muc) ops.Muc->add(mcel,(*allDOFs)[iele]);
     if(ops.Mcc) ops.Mcc->add(mcel,(*allDOFs)[iele]);
     if(ops.Msolver) ops.Msolver->add(mcel,(*allDOFs)[iele]);
   }
   else {
     if(ops.M)   ops.M->add(mel,(*allDOFs)[iele]);
     if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[iele]);
     if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[iele]);
     if(ops.Msolver) ops.Msolver->add(mel,(*allDOFs)[iele]);
   }
   if(matrixTimers) matrixTimers->assemble += getTime();

   // Form the element impedance matrix in kel and the element damping matrix in mel
   if(matrixTimers) matrixTimers->formTime -= getTime();
   if(!isShifted) {
     if(makeMass) {
       for(i = 0; i < dim; ++i)
         for(j = 0; j < dim; ++j) {
           double m  = mel[i][j];
           double k  = kel[i][j];
           mel[i][j] = (zeroRot && (dofType[ (*allDOFs)[iele][i] ] == 1 
                        || dofType[ (*allDOFs)[iele][j] ] == 1)) ? 0 : alpha*m + beta*k; // mel is now the damping element matrix
           kel[i][j] = Kcoef*k + Ccoef*mel[i][j] + Mcoef*m;
         }
     }
     else {
       for(i = 0; i < dim; ++i)
         for(j = 0; j < dim; ++j)
           kel[i][j] *= Kcoef;
     }
   }
   if(matrixTimers) matrixTimers->formTime += getTime();

   // Assemble the element impedance matrix kel in mat, ops.Kuc and ops.spp
   if(matrixTimers) matrixTimers->assemble -= getTime();
   if(isComplexF || (imag(kappa2) != 0)) {
     if(mat) mat->add(kcel,m_dofs);
     if(isShifted && ops.Kuc) ops.Kuc->add(kcel,(*allDOFs)[iele]);
     if(isShifted && ops.Kcc) ops.Kcc->add(kcel,(*allDOFs)[iele]);
     if(ops.spp) ops.spp->add(kcel,(*allDOFs)[iele]);
   }
   else {
     if(mat) mat->add(kel,m_dofs);
     if(isShifted && ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]); // note: Kuc is [K-omega2*M]_{uc} for IMPE (TODO check eigen)
     if(isShifted && ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);
     if(ops.spp) ops.spp->add(kel,(*allDOFs)[iele]);
     if(Kss) Kss->add(kel,(*allDOFs)[iele]); // for farfield output (TODO: check with Radek)
     if(isShifted && isDamped && isStructureElement(iele)) {
       if(mat) mat->addImaginary(izel,m_dofs);
       if(ops.Kuc) ops.Kuc->addImaginary(izel,(*allDOFs)[iele]);
       if(ops.Kcc) ops.Kcc->addImaginary(izel,(*allDOFs)[iele]);
       if(ops.spp) ops.spp->addImaginary(izel,(*allDOFs)[iele]);
     }
   }

   // Assemble the element damping matrix mel in ops.C and ops.C_deriv
   if(isShifted) {
     if(isDamped && isStructureElement(iele)) {
       izel /= omega;
       if(ops.C_deriv && ops.C_deriv[0]) ops.C_deriv[0]->addImaginary(izel,(*allDOFs)[iele]);
       if(ops.Cuc_deriv && ops.Cuc_deriv[0]) ops.Cuc_deriv[0]->addImaginary(izel,(*allDOFs)[iele]);
     }
   }
   else {
     if(ops.C) ops.C->add(mel,(*allDOFs)[iele]);
     if(ops.Cuc) ops.Cuc->add(mel,(*allDOFs)[iele]);
   }
   if(matrixTimers) matrixTimers->assemble += getTime();
 }

 // now deal with complex elements
 if(typeid(Scalar)==typeid(DComplex)) {
   // we assume they have no damping
   for(iele = 0; iele < numele; ++iele) {
     StructProp *prop = packedEset[iele]->getProperty();
     if(packedEset[iele]->isSommerElement()) continue;
     if(!packedEset[iele]->isComplex()) continue;
     m_dofs = (sinfo.type == 0 && subCast) ? (*domain->getAllDOFs())[subCast->getGlElems()[iele]] : (*allDOFs)[iele];

     if(matrixTimers) matrixTimers->formTime -= getTime();
     kelC = packedEset[iele]->complexStiffness(nodes, karrayC);
     this->densProjectStiffnessC(kelC, iele);
     if(sinfo.isCoupled) {  // coupled scaling
       if(isStructureElement(iele)) kelC *= cscale_factor2;
       else if(packedEset[iele]->isFsiElement()) kelC *= cscale_factor;
     }

     int i, j, dim = kel.dim();
     if(makeMass || isShifted) {
       melC = packedEset[iele]->complexMassMatrix(nodes, marrayC, mratio);
       this->densProjectStiffnessC(melC, iele);
       if(sinfo.isCoupled) {
         if(isStructureElement(iele)) melC *= cscale_factor2;
         else if(packedEset[iele]->isFsiElement()) melC *= coupledScaling;
       }
     }
     if(isShifted) {
       omega2 = (isFluidElement(iele)) ? packedEset[iele]->helmCoef() : geoSource->shiftVal();
       if(isStructureElement(iele)) {
         omega = sqrt(omega2);
       }
       for(i = 0; i < dim; ++i)
         for(j = 0; j < dim; ++j)
           kelC[i][j] -= omega2*melC[i][j];
       if(sinfo.doFreqSweep && isFluidElement(iele))
         melC *= packedEset[iele]->helmCoef()/geoSource->shiftVal();
     }
     else {
       // do something for coupled eigen/dynamics, set fluidCelerity = 1500
       if(makeMass && sinfo.isCoupled && isFluidElement(iele)) melC /= (1500.0 * 1500);
       if(sinfo.acoustic == 1 && prop) {
         double c = prop->ss;//speed of sound
         melC /= (c*c);
       }
     }
     if(matrixTimers) matrixTimers->formTime += getTime();

     if(matrixTimers) matrixTimers->assemble -= getTime();
     if(ops.K)   ops.K->add(kelC,(*allDOFs)[iele]);
     if(ops.Kuc) ops.Kuc->add(kelC,(*allDOFs)[iele]);
     if(ops.Kcc) ops.Kcc->add(kelC,(*allDOFs)[iele]);
     if(ops.M)   ops.M->add(melC,(*allDOFs)[iele]);
     if(ops.Muc) ops.Muc->add(melC,(*allDOFs)[iele]);
     if(ops.Mcc) ops.Mcc->add(melC,(*allDOFs)[iele]);
     if(ops.Msolver) ops.Msolver->add(melC,(*allDOFs)[iele]);
     if(mat) mat->add(kelC,m_dofs);
     if(ops.spp) ops.spp->add(kelC,(*allDOFs)[iele]);
     if(matrixTimers) matrixTimers->assemble += getTime();
   }
 }

 if(sinfo.ATDARBFlag >= 1.0) {
   getCurvatures3Daccurate(this);
 }
 int passage = 0;
 bool Complx = false;//Test for stabilising absorbing condition of high order - JF

 // LOOP over the sommer elements
 if((sinfo.ATDARBFlag == 1.5) && Complx) { // kel goes for DComplex, mat and ops.K(s) too
   for(iele = 0; iele < numele; ++iele) {
     if(!packedEset[iele]->isSommerElement() || sinfo.ATDARBFlag == -2.0) continue;
     m_dofs = (sinfo.type == 0 && subCast) ? (*domain->getAllDOFs())[subCast->getGlElems()[iele]] : (*allDOFs)[iele];
     //add massMatrix
     mel.zero();
     mel = packedEset[iele]->massMatrix(nodes, marray,mratio);
     if (passage==0) {
       int si = mel.dim();
       cel.setSize(si);
       passage=1;
     }
     if(ops.M)   ops.M->add(mel,(*allDOFs)[iele]);
     if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[iele]);
     if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[iele]);
     if(ops.Msolver) ops.Msolver->add(mel,(*allDOFs)[iele]);
     mel *= Mcoef;
     if(mat) mat->add(mel,m_dofs);
     if(ops.spp) ops.spp->add(mel,(*allDOFs)[iele]);

     //add dampingMatrix
     cel.zero();
     cel = packedEset[iele]->dampingMatrix(nodes, marray);
     if(ops.C) ops.C->add(cel,(*allDOFs)[iele]);
     if(ops.Cuc) ops.Cuc->add(cel,(*allDOFs)[iele]);
     cel *= Ccoef;
     if(mat) mat->add(cel,m_dofs);
     if(ops.spp) ops.spp->add(cel,(*allDOFs)[iele]);

     //add stiffness
     kel.zero();
     kel = packedEset[iele]->stiffness(nodes, karray);
     if(ops.K) ops.K->add(kel,(*allDOFs)[iele]);
     if(ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]);
     if(ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);
     kel *= Kcoef;
     if(mat) mat->add(kel,m_dofs);
     if(ops.spp) ops.spp->add(kel,(*allDOFs)[iele]);
     //add Imaginary part of stiffness
     kel.zero();
     kel = packedEset[iele]->imStiffness(nodes, karray);
     if(ops.K) ops.K->addImaginary(kel,(*allDOFs)[iele]);
     if(ops.Kuc) ops.Kuc->addImaginary(kel,(*allDOFs)[iele]);
     if(ops.Kcc) ops.Kcc->addImaginary(kel,(*allDOFs)[iele]);
     kel *= Kcoef;
     if(mat) mat->addImaginary(kel,m_dofs);
     if(ops.spp) ops.spp->addImaginary(kel,(*allDOFs)[iele]);
   }
 }
 else {
   for(iele = 0; iele < numele; ++iele) {
     if(!packedEset[iele]->isSommerElement() || sinfo.ATDARBFlag == -2.0) continue;
     m_dofs = (sinfo.type == 0 && subCast) ? (*domain->getAllDOFs())[subCast->getGlElems()[iele]] : (*allDOFs)[iele];

     // add massMatrix
     mel.zero();
     mel = packedEset[iele]->massMatrix(nodes, marray, mratio);
     if (passage==0) {
       int si = mel.dim();
       cel.setSize(si);
       passage=1;
     }
     if(ops.M)   ops.M->add(mel,(*allDOFs)[iele]);
     if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[iele]);
     if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[iele]);
     if(ops.Msolver) ops.Msolver->add(mel,(*allDOFs)[iele]);
     mel *= Mcoef;
     if(mat) mat->add(mel,m_dofs);
     if(ops.spp) ops.spp->add(mel,(*allDOFs)[iele]);

     // add dampingMatrix
     cel.zero();
     cel = packedEset[iele]->dampingMatrix(nodes, marray);
     if(ops.C) ops.C->add(cel,(*allDOFs)[iele]);
     if(ops.Cuc) ops.Cuc->add(cel,(*allDOFs)[iele]);
     cel *= Ccoef;
     if(mat) mat->add(cel,m_dofs);
     if(ops.spp) ops.spp->add(cel,(*allDOFs)[iele]);

     // add stiffness
     kel.zero();
     kel = packedEset[iele]->stiffness(nodes, karray);
     if(ops.K) ops.K->add(kel,(*allDOFs)[iele]);
     if(ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]);
     if(ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);
     kel *= Kcoef;
     if(mat) mat->add(kel,m_dofs);
     if(ops.spp) ops.spp->add(kel,(*allDOFs)[iele]);
   }
 }

 if(makeMass || isShifted) {
   // Add discrete mass contribution to the Mass Matrix
   // Three matrices need to be changed. Mass matrix itself,
   // Damping Matrix and K tilda.
   Scalar m;
   DMassData *current = firstDiMass;
   double omega2 = geoSource->shiftVal();
   double omega = sqrt(geoSource->shiftVal()); 
   while(current != 0) {
     // PJSA: modified all addDiscreteMass functions to accept dsa dof rather than cdsa dof (much safer)
     // TODO what about Kuc, Cuc?
     int dof = dsa->locate(current->node, (1 << current->dof));
     int m_dof = (sinfo.type == 0 && subCast) ? domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->dof)) : dof; 
     if(dof == -1) 
       { current = current->next; continue; }
     if(current->jdof > -1) { // PJSA 10-9-06 for off-diagonal mass terms eg. products of inertia I21, I31, I32
       int jdof = dsa->locate(current->node, (1 << current->jdof));
       int m_jdof = (sinfo.type == 0 && subCast) ? domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->jdof)) : jdof;
       if(jdof == -1) 
         { current = current->next; continue; }
       if(isShifted) {
         // note: for coupled problems we are assuming lumped mass is structural not fluid
         double mass = sinfo.isCoupled ? current->diMass*cscale_factor2 : current->diMass;
         double m_real = -omega2*mass;
         double m_imag = isDamped ? -omega*alpha*mass : 0.0; 
         ScalarTypes::initScalar(m, m_real, m_imag);
         if(mat) mat->add(m_dof, m_jdof, m);
         if(ops.spp) ops.spp->add(dof, jdof, m);
         if(isDamped) {
           ScalarTypes::initScalar(m, 0.0, alpha*mass);
           if(ops.C_deriv && ops.C_deriv[0]) ops.C_deriv[0]->add(dof, jdof, m);
         }
         if(ops.M) ops.M->add(dof, jdof, mass);
         if(ops.Msolver) ops.Msolver->add(dof, jdof, mass);
       }
       else {
         if(ops.M) ops.M->add(dof,jdof,current->diMass);
         if(ops.Msolver) ops.Msolver->add(dof,jdof,current->diMass);
         if(ops.C) ops.C->add(dof,jdof,alpha*current->diMass);
         double mass = Mcoef*current->diMass;
         if (ops.C) mass += Ccoef*alpha*mass;
         if(mat) mat->add(m_dof,m_jdof,mass);
         if(ops.spp) ops.spp->add(dof,jdof,mass);
       }
     }
     else {
       if(isShifted) {  // add discrete mass contributions to global matrices for frequency response analysis
         double mass = sinfo.isCoupled ? current->diMass*cscale_factor2 : current->diMass;
         double m_real = -omega2*mass;
         double m_imag = isDamped ? -omega*alpha*mass : 0.0;
         ScalarTypes::initScalar(m, m_real, m_imag);
         if(mat) mat->addDiscreteMass(m_dof, m);
         if(ops.spp) ops.spp->addDiscreteMass(dof, m);
         if(isDamped) {
           ScalarTypes::initScalar(m, 0.0, alpha*mass);
           if(ops.C_deriv && ops.C_deriv[0]) ops.C_deriv[0]->addDiscreteMass(dof, m);
         }
         if(ops.M) ops.M->addDiscreteMass(dof, mass);
         if(ops.Msolver) ops.Msolver->addDiscreteMass(dof, mass);
       }
       else {
         if(ops.M) ops.M->addDiscreteMass(dof, current->diMass);
         if(ops.Msolver) ops.Msolver->addDiscreteMass(dof, current->diMass);
         if(ops.C) ops.C->addDiscreteMass(dof, alpha*current->diMass);
         double mass = Mcoef*current->diMass;
         if (ops.C) mass += Ccoef*alpha*mass;
         if(mat) mat->addDiscreteMass(m_dof, mass);
         if(ops.spp) ops.spp->addDiscreteMass(dof, mass);
       }
     }
     current = current->next;
   }
 }

 //ADDED FOR HEV PROBLEM, EC, 20070820
 if((sinfo.HEV) && (probType() == SolverInfo::Modal)) {
   Mff = new GenBLKSparseMatrix<double>(nodeToNodeFluid, dsaFluid, c_dsaFluid, sinfo.trbm, sinfo.sparse_renum, 0);
   Mff->zeroAll();

   for(iele = 0; iele < geoSource->numElemFluid(); ++iele) {
     //element matrix already premultiplied by fluid density
     kel = (*(geoSource->getPackedEsetFluid()))[iele]->stiffness(nodes,karray);
     Mff->add(kel,(*allDOFsFluid)[iele]);
   }
   fprintf(stderr," ... Factoring Fluid Mass Matrix ...\n");
   //double trbmtemp = sinfo.trbm;
   //sinfo.setTrbm(1e-16);
   Mff->factor();
   //sinfo.setTrbm(trbmtemp);
   fprintf(stderr," ... Fluid Mass Matrix Factored...\n");

   if(sinfo.addedMass == 1) {

     //get fluid density from first element since HOMOGENEOUS fluid is assumed
     double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;

     double** C_NZrows = getCMatrix();
     double** FCtinv_NZcol = new double* [domain->nuNonZero];

     int nnp = domain->numUnconFluid();

     for (int i=0; i < domain->nuNonZero; ++i)  {
       FCtinv_NZcol[i] = new double[nnp];
       for (int jp=0; jp < nnp; ++jp)  {
          C_NZrows[i][jp] = rhoFluid * C_NZrows[i][jp];
          FCtinv_NZcol[i][jp] = C_NZrows[i][jp];
       }
     }

     //obtain non zero parts of F^-1C^T
     Mff->reSolve(domain->nuNonZero,FCtinv_NZcol);

     //initialize the added mass operator (Ma_NZ only contains the non-zero parts)
     FullSquareMatrix Ma_NZ(domain->nuNonZero,0);
     Ma_NZ.zero();

     //obtain the added mass operator
     for (int i=0; i<domain->nuNonZero; ++i)  {
       for (int j=i; j<domain->nuNonZero; ++j)  {
         for (int k=0; k<domain->npNonZero[i]; ++k)  {
           int K = domain->pmap[i][k];
             Ma_NZ[i][j] += C_NZrows[i][K]*FCtinv_NZcol[j][K];
         }
       }
     }
     for (int i=0; i<domain->nuNonZero; ++i)  {
       for (int j=i; j<domain->nuNonZero; ++j)  {
         Ma_NZ[j][i] = Ma_NZ[i][j];
       }
     }

     fprintf(stderr," ... HEV Problem: Added mass obtained ...\n");

     //add added mass to structural mass for EVP!
     ops.M->add(Ma_NZ,domain->umap_add);
     fprintf(stderr," ... HEV Problem: Added mass contribution included in structural mass ...\n");

     if(isShifted) {
       omega2 = geoSource->shiftVal();
       Ma_NZ *= (-omega2*Kcoef);
       if (ops.K) ops.K->add(Ma_NZ,domain->umap_add);
       if (mat)   mat->add(Ma_NZ,domain->umap_add);
     if (ops.spp)   ops.spp->add(Ma_NZ,domain->umap_add);
       fprintf(stderr," ... HEV Problem: Added mass contribution included in shifted opertor ...\n");
     }
   }

 }

 delete [] karray;
 delete [] marray;
 delete [] kcarray;
 delete [] mcarray;
 delete [] karrayC;
 delete [] marrayC;

 if(sinfo.isAcousticHelm()) assembleSommer<Scalar>(mat, &ops);

 if(sinfo.ATDROBalpha != 0.0) assembleATDROB<Scalar>(mat, &ops,Kcoef);

 if(matrixTimers) matrixTimers->memoryForm += memoryUsed();
}

template<class Scalar>
GenDBSparseMatrix<Scalar> *
Domain::constructDBSparseMatrix(DofSetArray *dof_set_array, Connectivity *cn)
{
 if(dof_set_array == 0) dof_set_array = c_dsa;
 if(cn == 0)
   return new GenDBSparseMatrix<Scalar>(nodeToNode, dsa, c_dsa);
 else
   return new GenDBSparseMatrix<Scalar>(cn, dsa, c_dsa);
}

template<class Scalar>
GenEiSparseMatrix<Scalar> *
Domain::constructEiSparseMatrix(DofSetArray *c_dsa, Connectivity *nodeToNode, bool flag)
{
#ifdef USE_EIGEN3
  if(c_dsa == 0) c_dsa = Domain::c_dsa;
  if(nodeToNode == 0) nodeToNode = Domain::nodeToNode;
  if(sinfo.mpcDual) {
    Connectivity *nodeToNodeG = nodeToNode;
    DofSetArray *g_dsa = new ConstrainedDSA(*dsa, *Domain::c_dsa);
    typename WrapEiSparseMat<Scalar>::CtorData baseArg(nodeToNodeG, dsa, g_dsa);
    return new GoldfarbIdnaniQpSolver<WrapEiSparseMat<Scalar>, Scalar>(baseArg, Domain::c_dsa, sinfo.goldfarb_tol);
  }
  else {
   return new GenEiSparseMatrix<Scalar>(nodeToNode, dsa, c_dsa, flag);
  }
#else
 cerr << "USE_EIGEN3 is not defined\n";
#endif
}

template<class Scalar>
GenNBSparseMatrix<Scalar> *
Domain::constructNBSparseMatrix()
{
 return new GenNBSparseMatrix<Scalar>(nodeToNode, c_dsa);
}

template<class Scalar>
GenCuCSparse<Scalar> *
Domain::constructCuCSparse(DofSetArray *dof_set_array)
{
 if(dof_set_array == 0) dof_set_array = c_dsa;
 GenCuCSparse<Scalar> *sp = 0;
 if((dsa->size() - dof_set_array->size()) == 0)
   sp = 0;
 else
   sp = new GenCuCSparse<Scalar>(nodeToNode, dsa, dof_set_array);
 return sp;
}

template<class Scalar>
GenCuCSparse<Scalar> *
Domain::constructCCSparse(DofSetArray *dof_set_array)
{
  if(dof_set_array == 0)
    dof_set_array = c_dsa;
  return new GenCuCSparse<Scalar>(nodeToNode, dsa, dof_set_array->getConstrndNum(),
		                  dof_set_array->getConstrndNum());
}

template<class Scalar>
GenSkyMatrix<Scalar> *
Domain::constructSkyMatrix(DofSetArray *DSA, Rbm *rbm)
{
  if(DSA==0) DSA=c_dsa;
  if(!sinfo.getDirectMPC())
    return new GenSkyMatrix<Scalar>(nodeToNode, DSA, sinfo.trbm, rbm);
  else {
    if(nodeToNodeDirect) delete nodeToNodeDirect;
    nodeToNodeDirect = prepDirectMPC();
    DOFMap *baseMap = new DOFMap[dsa->size()];
    DOFMap *eqMap = new DOFMap[DSA->size()];
    // TODO Examine when DSA can be different from c_dsa
    ConstrainedDSA *MpcDSA = makeMaps(dsa, c_dsa, baseMap, eqMap);
    typename WrapSkyMat<Scalar>::CtorData baseArg(nodeToNodeDirect, MpcDSA, sinfo.trbm, rbm);
    int nMappedEq = DSA->size();
    return
      new MappedAssembledSolver<WrapSkyMat<Scalar>, Scalar>(baseArg, dsa->size(), baseMap,
          nMappedEq, eqMap, c_dsa);
  }
}

template<class Scalar>
GenBlockSky<Scalar> *
Domain::constructBlockSky(DofSetArray *DSA)
{
  if(DSA==0) DSA=c_dsa;
  return new GenBlockSky<Scalar>(nodeToNode, DSA, sinfo.trbm);
}

template<class Scalar>
GenBLKSparseMatrix<Scalar> *
Domain::constructBLKSparseMatrix(DofSetArray *DSA, Rbm *rbm)
{
  if(DSA == 0) DSA = c_dsa;
  if(sinfo.newmarkBeta == 0.0 && sinfo.acoustic) {
    if(sinfo.ATDARBFlag >- 1.0) {
      return new GenBLKSparseMatrix<Scalar>(nodeToNode_sommer->modify(), dsa, DSA, sinfo.trbm, sinfo.sparse_renum, rbm);
    }
    else {
      int numN=numNodes();
      // import a diagonal connectivity for the mass matrix
      Connectivity *connForMas = 0;
      connForMas = new Connectivity(numN);
      //connForMass->print();
      connForMas = connForMas->modify();
      GenBLKSparseMatrix<Scalar> * retur =  new GenBLKSparseMatrix<Scalar>(connForMas, dsa, DSA, sinfo.trbm, sinfo.sparse_renum, rbm);
      delete connForMas;
      return retur;
    }
  }
  else {
    if(!sinfo.getDirectMPC())
      return new GenBLKSparseMatrix<Scalar>(nodeToNode, dsa, DSA, sinfo.trbm, sinfo.sparse_renum, rbm);
    else {
      if(nodeToNodeDirect) delete nodeToNodeDirect;
      nodeToNodeDirect = prepDirectMPC();
      DOFMap *baseMap = new DOFMap[dsa->size()];
      DOFMap *eqMap = new DOFMap[DSA->size()];
      // TODO Examine when DSA can be different from c_dsa
      ConstrainedDSA *MpcDSA = makeMaps(dsa, c_dsa, baseMap, eqMap);
      typename WrapSparseMat<Scalar>::CtorData
        baseArg(nodeToNodeDirect, dsa, MpcDSA, sinfo.trbm, sinfo.sparse_renum, /*rbm*/ (Rbm*)NULL); // TODO consider rbm issue
      int nMappedEq = DSA->size();
      return new MappedAssembledSolver<WrapSparseMat<Scalar>, Scalar>(baseArg, dsa->size(), baseMap, nMappedEq, eqMap, c_dsa);
    }
  }
}

template<class Scalar>
GenSGISparseMatrix<Scalar> *
Domain::constructSGISparseMatrix(Rbm *rbm)
{
  return new GenSGISparseMatrix<Scalar>(nodeToNode, dsa, c_dsa, rbm);
}

template<class Scalar>
GenSGISparseMatrix<Scalar> *
Domain::constructSGISparseMatrix(int subNumber, Rbm *rbm)
{
  return new GenSGISparseMatrix<Scalar>(nodeToNode, dsa, c_dsa, rbm, subNumber);
}

template<class Scalar>
GenPCGSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar> > *
Domain::constructPCGSolver(GenSparseMatrix<Scalar> *K, Rbm *rbm)
{
  return new GenPCGSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar> >
                                (K, sinfo.precond, sinfo.maxit, sinfo.tol, sinfo.maxvecsize, rbm);
}

template<class Scalar>
GenSpoolesSolver<Scalar> *
Domain::constructSpooles(ConstrainedDSA *DSA, Rbm *rbm)
{
  if(DSA == 0) DSA = c_dsa;
  if(!sinfo.getDirectMPC())
    return new GenSpoolesSolver<Scalar>(nodeToNode, dsa, DSA, rbm);
  else {
    if(nodeToNodeDirect) delete nodeToNodeDirect;
    nodeToNodeDirect = prepDirectMPC();
    DOFMap *baseMap = new DOFMap[dsa->size()];
    DOFMap *eqMap = new DOFMap[DSA->size()];
    // TODO Examine when DSA can be different from c_dsa
    ConstrainedDSA *MpcDSA = makeMaps(dsa, c_dsa, baseMap, eqMap);
    typename WrapSpooles<Scalar>::CtorData baseArg(nodeToNodeDirect, dsa, MpcDSA, rbm);
    int nMappedEq = DSA->size();
    return
      new MappedAssembledSolver<WrapSpooles<Scalar>, Scalar>(baseArg, dsa->size(), baseMap,
          nMappedEq, eqMap, c_dsa);
  }
}

template<class Scalar>
GenMumpsSolver<Scalar> *
Domain::constructMumps(ConstrainedDSA *DSA, Rbm *rbm, FSCommunicator *com)
{
  if(DSA == 0) DSA = c_dsa;
  if(!sinfo.getDirectMPC())
    return new GenMumpsSolver<Scalar>(nodeToNode, dsa, DSA, com);
  else {
    if(nodeToNodeDirect) delete nodeToNodeDirect;
    nodeToNodeDirect = prepDirectMPC();
    DOFMap *baseMap = new DOFMap[dsa->size()];
    DOFMap *eqMap = new DOFMap[DSA->size()];
    // TODO Examine when DSA can be different from c_dsa
    ConstrainedDSA *MpcDSA = makeMaps(dsa, c_dsa, baseMap, eqMap);
    typename WrapMumps<Scalar>::CtorData baseArg(nodeToNodeDirect, dsa, MpcDSA, com);
    int nMappedEq = DSA->size();
    return
      new MappedAssembledSolver<WrapMumps<Scalar>, Scalar>(baseArg, dsa->size(), baseMap,
          nMappedEq, eqMap, c_dsa);
  }
}

template<class Scalar>
Rom::GenGaussNewtonSolver<Scalar> *
Domain::constructGaussNewtonSolver()
{
  return new Rom::GenGaussNewtonSolver<Scalar>(nodeToNode, dsa, c_dsa);
}

template<class Scalar>
Rom::GenGalerkinProjectionSolver<Scalar> *
Domain::constructGalerkinProjectionSolver()
{
  return new Rom::GenGalerkinProjectionSolver<Scalar>(nodeToNode, dsa, c_dsa);
}

template<class Scalar>
Rom::GenGappyProjectionSolver<Scalar> *
Domain::constructGappyProjectionSolver()
{
  return new Rom::GenGappyProjectionSolver<Scalar>(nodeToNode, dsa, c_dsa);
}

template<class Scalar>
void
Domain::buildOps(AllOps<Scalar> &allOps, double Kcoef, double Mcoef, double Ccoef,
                 Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray, bool factorize)
{
 if(matrixTimers) matrixTimers->memorySolve -= memoryUsed();

 if(allOps.sysSolver) delete allOps.sysSolver;
 GenSolver<Scalar> *systemSolver = 0;
// RT: 032010 based on Phil's input
// if(!sinfo.inpc) {
//   if (allOps.sysSolver) delete allOps.sysSolver;
//   allOps.sysSolver = 0;
// }

 //GenSparseMatrix<Scalar> *spm = 0;
 allOps.spm = 0;
 SfemBlockMatrix<Scalar> *sfbm = 0;
 int L = 1;
 int ndim = 0;
 if(sinfo.inpc||sinfo.noninpc) {
   L =  sfem->getL();
   int n = domain->numUncon();
   int P = sfem->getP();
   ndim = sfem->getndim();
   int output_order = sfem->getoutput_order();
   if(sinfo.inpc)  sfbm = new SfemBlockMatrix<Scalar>(L,n,P,ndim,output_order);
 }

// for(int i=0; i<L; ++i) {
 for(int i=0; i<ndim+1; ++i) {
   if(sinfo.inpc)  domain->setNewProperties(i);
   if(sinfo.noninpc && i>0) break;
   switch(sinfo.type) {
    default:
      fprintf(stderr," *** WARNING: Solver not Specified  ***\n");
    case 0:
      makeStaticOpsAndSolver<Scalar>(allOps, Kcoef, Mcoef, Ccoef,
                                     systemSolver, allOps.spm, rbm, kelArray, melArray); // also used for eigen
      break;
    case 1:
      makeDynamicOpsAndSolver<Scalar>(allOps, Kcoef, Mcoef, Ccoef,
                                      systemSolver, allOps.spm, rbm, kelArray, melArray);
      break;
   }
   if(sinfo.inpc) {
     sfbm->setKi(allOps.spm,i);
     if(i==0 && sinfo.precond==3) {
       GenBLKSparseMatrix<Scalar> *prec_solver = constructBLKSparseMatrix<Scalar>(c_dsa, rbm);
       prec_solver->zeroAll();
       AllOps<Scalar> allOps_tmp;
       makeSparseOps<Scalar>(allOps_tmp,Kcoef,Mcoef,Ccoef,prec_solver,kelArray,melArray);
       prec_solver->factor();
       sfbm->setMeanSolver(prec_solver);
     }
     if(i==0) {
       allOps.rhs_inpc = new GenVector<Scalar>(domain->numUncon());
       domain->template buildRHSForce<Scalar>(*allOps.rhs_inpc, allOps.Kuc);
     }
     if(allOps.Kuc) allOps.Kuc->zeroAll();
   }
 }

 // Set allOps pointer to system solver
 if(!sinfo.inpc) allOps.sysSolver = systemSolver;
 else {
   switch(sinfo.iterType) {
     default:
     case 0 : {
       filePrint(stderr," ... CG Solver is Selected           ...\n");
       allOps.sysSolver = new GenPCGSolver<Scalar, GenVector<Scalar>, SfemBlockMatrix<Scalar> >(sfbm, sinfo.precond, sinfo.maxit,
                                                                                                sinfo.tol, sinfo.maxvecsize);
       break;
     }
     case 4: {
       filePrint(stderr," ... Bi-CG Solver is Selected       ...\n");
       allOps.sysSolver = new GenBCGSolver<Scalar, GenVector<Scalar>, SfemBlockMatrix<Scalar>, GenSolver<Scalar> >(sinfo.maxit, sinfo.tol, sfbm);
       break;
     }
     case 5: {
       filePrint(stderr," ... CR Solver is Selected          ...\n");
       allOps.sysSolver = new GenCRSolver<Scalar, GenVector<Scalar>, SfemBlockMatrix<Scalar>, GenSolver<Scalar> >(sinfo.maxit, sinfo.tol, sfbm);
       break;
     }
   }
 }

 if(factorize)
   {
     // Time system matrix factorization
     if(matrixTimers) matrixTimers->factor -= getTime();

     if(systemSolver) systemSolver->factor();

     if(matrixTimers) matrixTimers->factor += getTime();
     if(matrixTimers) matrixTimers->memorySolve += memoryUsed();
   }
}

template<class Scalar>
void
Domain::rebuildOps(AllOps<Scalar> &allOps, double Kcoef, double Mcoef, double Ccoef,
                   Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray,
                   bool factorize)
{
 GenSolver<Scalar> *systemSolver;
 GenSparseMatrix<Scalar> *spm;

 switch(sinfo.type) {

  case 0:

     switch( sinfo.subtype ) {
       default:
         fprintf(stderr,"Driver.d/OpMake.C - rebuildOps, sinfo.subtype = %d not implemented, going for sky\n",sinfo.subtype);
       case 0: { //case 1:
         spm = (GenSkyMatrix<Scalar>*)allOps.sysSolver;
         spm->zeroAll();
         makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
         systemSolver  = (GenSkyMatrix<Scalar>*) spm;
       }
       break;
       case 5: { //case 2:
	 makeFrontalOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,rbm,kelArray,melArray);
         systemSolver = allOps.sysSolver;
       }
       break;
       case 1: { //case 3:
         spm = (GenBLKSparseMatrix<Scalar>*)allOps.sysSolver;
         spm->zeroAll();
         makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
         systemSolver   = (GenBLKSparseMatrix<Scalar>*) spm;
       }
       break;
       case 2: { //case 4:
         spm = (GenSGISparseMatrix<Scalar>*)allOps.sysSolver;
         spm->zeroAll();
         makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
         systemSolver   = (GenSGISparseMatrix<Scalar>*) spm;
       }
       break;
       case 3: { //case 5:
#ifdef NO_COMPLEX
	 spm = dynamic_cast<SGISky*>(allOps.sysSolver);
	 spm->zeroAll();
	 makeSparseOps(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
	 systemSolver   = dynamic_cast<SGISky*>(spm);
#else
	 fprintf(stderr,"ERROR: templated SGISkyMatrix class is not implemeted \n");
	 assert(0);
#endif
       }
       break;
       case 8: { //case 8:
         spm =(GenSpoolesSolver<Scalar>*)allOps.sysSolver;
         spm->zeroAll();
         makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
         systemSolver   = (GenSpoolesSolver<Scalar>*) spm;
       }
       break;
     }
     break;

  case 1:
    fprintf(stderr," Error in rebuildOps: rebuild for this solver\n");
    fprintf(stderr,"                      not implemented ...exit\n");
    exit(-1);
    break;
 }

 if(factorize)
   {
     systemSolver->factor();
   }
}

template<class Scalar>
void
Domain::getSolverAndKuc(GenSolver<Scalar> *&solver, GenSparseMatrix<Scalar> *&kuc,
                        FullSquareMatrix *kelArray, bool factorize)
{
 // ... Declare an AllOps structure
 AllOps<Scalar> allOps;

 // ... Call necessary Operator's constructors
 allOps.Kuc = constructCuCSparse<Scalar>();

 Rbm *rbm = 0;
 // ... Construct geometric rigid body modes if necessary
 if(sinfo.rbmflg) rbm = constructRbm();

 // ... Construct "thermal rigid body mode" if necessary
 if(sinfo.hzemFlag) {
   rbm = constructHzem();
   if(rbm->numRBM() == 0) rbm = 0;
 }

 // ... Construct "Sloshing rigid body mode" if necessary
 // ADDED FOR SLOSHING PROBLEM, EC, 20070723
 if(sinfo.slzemFlag) {
   fprintf(stderr," ... constructSlzem called in OpMake.C -- THIS CALL HAS NOT BEEN DEBUGGED ...\n");
   rbm = constructSlzem();
   if(rbm->numRBM() == 0) rbm = 0;
 }

 // ... Build stiffness matrix K and Kuc
 buildOps<Scalar>(allOps, 1.0, 0.0, 0.0, rbm, kelArray, factorize);

 // ... Return with solver and Kuc
 solver = allOps.sysSolver;
 kuc    = allOps.Kuc;
}

template<class Scalar>
void
Domain::getSolverAndKuc(AllOps<Scalar> &allOps, FullSquareMatrix *kelArray, bool factorize)
{
 // PJSA 10-5-04: new version, this one can pass all sys matricies back via allOps
 // not just K and Kuc. This is required for frequency sweep analysis

 // ... Call necessary Operator's constructors
 allOps.Kuc = constructCuCSparse<Scalar>();
 allOps.Kcc = constructCCSparse<Scalar>();

 Rbm *rbm = 0;
 // ... Construct geometric rigid body modes if necessary
 if(sinfo.rbmflg) rbm = constructRbm();

 // ... Construct "thermal rigid body mode" if necessary
 if(sinfo.hzemFlag) {
   rbm = constructHzem();
   if(rbm->numRBM() == 0) rbm = 0;
 }

 // ... Construct "Sloshing rigid body mode" if necessary
 // ADDED FOR SLOSHING PROBLEM, EC, 20070723
 if(sinfo.slzemFlag) {
   fprintf(stderr," ... constructSlzem called in OpMake.C -- THIS CALL HAS NOT BEEN DEBUGGED ...\n");
   rbm = constructSlzem();
   if(rbm->numRBM() == 0) rbm = 0;
 }

 // for freqency sweep: need M, Muc, C, Cuc
 if(sinfo.doFreqSweep && sinfo.nFreqSweepRHS > 1) {
   //---- UH ----
   if(sinfo.freqSweepMethod == SolverInfo::PadeLanczos ||
      sinfo.freqSweepMethod == SolverInfo::GalProjection ||
      sinfo.freqSweepMethod == SolverInfo::KrylovGalProjection ||
      sinfo.freqSweepMethod == SolverInfo::QRGalProjection) {
     if (allOps.K)
       delete allOps.K;
     if (allOps.M) delete allOps.M;
     if (allOps.Muc) delete allOps.Muc;
     if (allOps.C_deriv) if (allOps.C_deriv[0]) delete allOps.C_deriv[0];

     if(sinfo.isCoupled) allOps.K = constructNBSparseMatrix<Scalar>(); // unsymmetric
     else allOps.K = constructDBSparseMatrix<Scalar>(); // symmetric
   }
   //---- UH ----
   if(sinfo.isCoupled) allOps.M = constructNBSparseMatrix<Scalar>();  // unsymmetric
   else allOps.M = constructDBSparseMatrix<Scalar>();  // symmetric
   allOps.Muc = constructCuCSparse<Scalar>();
   bool isDamped = (sinfo.alphaDamp != 0.0) || (sinfo.betaDamp != 0.0);
   if(isDamped || (numSommer > 0)) {
     int numC_deriv;
     if((numSommer > 0) && ((sommerfeldType == 2) || (sommerfeldType == 4)))
       numC_deriv = sinfo.nFreqSweepRHS - 1;
     else
       numC_deriv = 1;
     allOps.C_deriv = new GenSparseMatrix<Scalar> * [sinfo.nFreqSweepRHS - 1];
     for(int n=0; n<numC_deriv; ++n)
       allOps.C_deriv[n] = constructDBSparseMatrix<Scalar>();
     for(int n=numC_deriv; n<sinfo.nFreqSweepRHS - 1; ++n)
       allOps.C_deriv[n] = 0;
     if(c_dsa->size() > 0 && (c_dsa->size() - dsa->size()) != 0) {
       allOps.Cuc_deriv = new GenSparseMatrix<Scalar> * [sinfo.nFreqSweepRHS - 1];
       for(int n=0; n<numC_deriv; ++n)
         allOps.Cuc_deriv[n] = constructCuCSparse<Scalar>();
       for(int n=numC_deriv; n<sinfo.nFreqSweepRHS - 1; ++n)
         allOps.Cuc_deriv[n] = 0;
     }
   }
 }

 // ... Build stiffness matrix K and Kuc, etc...
 buildOps<Scalar>(allOps, 1.0, 0.0, 0.0, rbm, kelArray, (FullSquareMatrix *) NULL, factorize);
}

template<class Scalar>
void
Domain::makeStaticOpsAndSolver(AllOps<Scalar> &allOps, double Kcoef, double Mcoef,
                 double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
                 Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray)
{
  switch(sinfo.subtype) {
    default:
    case 0:
      spm = constructSkyMatrix<Scalar>(c_dsa,rbm);
      makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
      systemSolver  = (GenSkyMatrix<Scalar>*) spm;
      break;
    case 1:
      spm = constructBLKSparseMatrix<Scalar>(c_dsa, rbm);
      spm->zeroAll();
      makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
      systemSolver   = (GenBLKSparseMatrix<Scalar>*) spm;
      break;
    case 2:
      if(matrixTimers) matrixTimers->constructTime -= getTime();
      spm = constructSGISparseMatrix<Scalar>(rbm);
      if(matrixTimers) matrixTimers->constructTime += getTime();
      makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
      systemSolver   = (GenSGISparseMatrix<Scalar>*) spm;
      break;
    case 3:
#ifdef NO_COMPLEX
      spm = constructSGISkyMatrix(rbm);
      makeSparseOps<double>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
      systemSolver   = (SGISky*) spm;
#else
      fprintf(stderr,"ERROR: templated SGISkyMatrix class is not implemeted \n");
#endif
      break;
#ifdef USE_EIGEN3
    case 4:
      spm = constructEiSparseMatrix<Scalar>(c_dsa);
      makeSparseOps<Scalar>(allOps, Kcoef, Mcoef, Ccoef, spm, kelArray, melArray);
      systemSolver  = (GenEiSparseMatrix<Scalar>*) spm;
      break;
#endif
    case 5:
      makeFrontalOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,rbm,kelArray,melArray);
      systemSolver = allOps.sysSolver;
      break;
#ifdef EIGEN_SUPERLU_SUPPORT
    case 7:
      spm = constructEiSparseMatrix<Scalar>(c_dsa, nodeToNode, false);
      makeSparseOps<Scalar>(allOps, Kcoef, Mcoef, Ccoef, spm, kelArray, melArray);
      systemSolver  = (GenEiSparseMatrix<Scalar>*) spm;
      break;
#endif
#ifdef USE_SPOOLES
    case 8:
      spm = constructSpooles<Scalar>(c_dsa, rbm);
      makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
      systemSolver   = (GenSpoolesSolver<Scalar>*) spm;
      break;
#endif
#ifdef USE_MUMPS
    case 9:
#ifdef DISTRIBUTED
      spm = constructMumps<Scalar>(c_dsa, rbm, new FSCommunicator(structCom));
#else
      spm = constructMumps<Scalar>(c_dsa, rbm);
#endif
      makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
      systemSolver   = (GenMumpsSolver<Scalar>*) spm;
      break;
#endif
    case 10:
      //filePrint(stderr," ... Diagonal Solver is Selected    ...\n");
      spm = new GenDiagMatrix<Scalar>(c_dsa); // XML NEED TO DEAL WITH RBMS
      makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
      systemSolver   = (GenDiagMatrix<Scalar>*) spm;
      break;
    case 11:
      //filePrint(stderr," ... POD-GN Solver is Selected      ...\n");
      {
        Rom::GenGaussNewtonSolver<Scalar> * solver = constructGaussNewtonSolver<Scalar>();
        spm = solver;
        spm->zeroAll();
        makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
        systemSolver = solver;
      }
      break;
    case 12:
      //filePrint(stderr," ... POD-Galerkin Solver is Selected...\n");
      {
        Rom::GenGalerkinProjectionSolver<Scalar> * solver = constructGalerkinProjectionSolver<Scalar>();
        spm = solver;
        spm->zeroAll();
        makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
        systemSolver = solver;
      }
      break;
    case 13:
      //filePrint(stderr," ... Gappy-POD Solver is Selected   ...\n");
      {
        Rom::GenGappyProjectionSolver<Scalar> * solver = constructGappyProjectionSolver<Scalar>();
        spm = solver;
        spm->zeroAll();
        makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray);
        systemSolver = solver;
      }
      break;
  }
}

template<class Scalar>
void
Domain::makeDynamicOpsAndSolver(AllOps<Scalar> &allOps, double Kcoef, double Mcoef,
                 double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
                 Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray)
{
  switch(sinfo.iterSubtype) {
    case 2:
      filePrint(stderr," ... Node Based Sparse Matrix       ...\n");
      spm = constructNBSparseMatrix<Scalar>();
      break;
    default:
    case 3:
      filePrint(stderr," ... Dof Based Sparse Matrix        ...\n");
      spm = constructDBSparseMatrix<Scalar>();
      break;
#ifdef USE_EIGEN3
    case 4:
      filePrint(stderr," ... Eigen 3 Sparse Matrix          ...\n");
      spm = constructEiSparseMatrix<Scalar>();
      break;
#endif
  }
  switch(sinfo.precond) {
    case 1:
      filePrint(stderr," ... Diagonal Preconditioner        ...\n");
      GenDiagMatrix<Scalar> *diag = new GenDiagMatrix<Scalar>(c_dsa);
      allOps.prec = (GenSolver<Scalar> *) diag;
      allOps.spp = (GenSparseMatrix<Scalar> *) diag;
      break;
  }
  makeSparseOps<Scalar>(allOps, Kcoef, Mcoef, Ccoef, spm, kelArray, melArray);
  if(allOps.prec) allOps.prec->factor();
  if(sinfo.inpc) { systemSolver = 0; return; }
  switch(sinfo.iterType) {
    default:
    case 0: {
      filePrint(stderr," ... CG Solver is Selected          ...\n");
      if (sinfo.precond==3) {
        GenBLKSparseMatrix<Scalar> *prec_solver = constructBLKSparseMatrix<Scalar>(c_dsa, rbm);
        prec_solver->zeroAll();
        AllOps<Scalar> allOps_tmp;
        makeSparseOps<Scalar>(allOps_tmp,Kcoef,Mcoef,Ccoef,prec_solver,kelArray,melArray);
        prec_solver->factor();
        spm->setMeanSolver(prec_solver);
      }

      GenPCGSolver<Scalar,GenVector<Scalar>, GenSparseMatrix<Scalar> >
        *pcgSolver = constructPCGSolver<Scalar>(spm, rbm);
      systemSolver = pcgSolver;
      break;
    }
    case 4: {
      filePrint(stderr," ... Bi-CG Solver is Selected       ...\n");
      GenBCGSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar> >
        *bcgSolver = new GenBCGSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar> >(sinfo.maxit, sinfo.tol, spm, allOps.prec);
      systemSolver = bcgSolver;
      break;
    }
    case 5: {
      filePrint(stderr," ... CR Solver is Selected          ...\n");
      GenCRSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar> >
        *crSolver = new GenCRSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar> >(sinfo.maxit, sinfo.tol, spm, allOps.prec);
      systemSolver = crSolver;
      break;
    }
    case 1 : {
      filePrint(stderr," ... GMRES Solver is Selected       ...\n");
      GmresSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar>, GenSolver<Scalar> >
        *gmresSolver = new  GmresSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar>, GenSolver<Scalar>, GenSolver<Scalar> >
         (sinfo.maxit, sinfo.tol, spm, &GenSparseMatrix<Scalar>::matvec, allOps.prec, &GenSolver<Scalar>::apply, NULL, &GenSolver<Scalar>::solve, (FSCommunicator *) NULL);
      if(sinfo.maxvecsize > 0) gmresSolver->maxortho = sinfo.maxvecsize;
      gmresSolver->verbose = verboseFlag;
      gmresSolver->printNumber = sinfo.fetiInfo.printNumber;
      systemSolver = gmresSolver;
      break;
    }
  }
}


template<class Scalar>
void
Domain::addGravityForce(GenVector<Scalar> &force)
{
  // ... ADD ELEMENT MASS CONTRIBUTION ...
  Vector elementGravityForce(maxNumDOFs);
  int gravflg;
  for(int iele = 0; iele < numele; ++iele) {
    if(packedEset[iele]->getProperty() == 0) continue; // phantom element
       
    if(geoSource->consistentQFlag() && !(sinfo.isDynam() && packedEset[iele]->getMassType() == 0)) 
      gravflg = 2;                       // 2: consistent (for dynamics, consistent gravity should not be used if element only has a lumped mass matrix)
    else gravflg = geoSource->fixedEndM; // 1: lumped with fixed-end moments
                                         // 0: lumped without fixed-end moments

    elementGravityForce.zero();
    packedEset[iele]->getGravityForce(nodes, gravityAcceleration, elementGravityForce, gravflg);

    for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(cn >= 0)
        force[cn] += elementGravityForce[idof];
    }
  }

  // ... ADD DISCRETE MASS TO GRAVITY FORCE ...
  DMassData *current = firstDiMass;
  while(current != 0) {
    int thisDof = c_dsa->locate(current->node, (1 << current->dof));
    if(thisDof >= 0 ) {
      if(current->dof == 0)
        force[thisDof] += (current->diMass)*gravityAcceleration[0];
      if(current->dof == 1)
        force[thisDof] += (current->diMass)*gravityAcceleration[1];
      if(current->dof == 2)
        force[thisDof] += (current->diMass)*gravityAcceleration[2];
    }
    current = current->next;
  }
}

template<class Scalar>
void
Domain::addPressureForce(GenVector<Scalar> &force, double lambda)
{
  Vector elementPressureForce(maxNumDOFs);
  int cflg = 1; // NOW WE ALWAYS USE CONSISTENT PRESSURE

  if(pressureFlag()) {
    for(int iele = 0; iele < numele; ++iele) {
      // If there is a zero pressure defined, skip it.
      if(packedEset[iele]->getPressure() == 0.0) continue;

      // Otherwise, compute element pressure force
      elementPressureForce.zero();
      packedEset[iele]->computePressureForce(nodes, elementPressureForce, (GeomState *) 0, cflg);

      // Assemble element pressure forces into domain force vector
      for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
        if(cn >= 0)
          force[cn] += lambda*elementPressureForce[idof];
      }
    }
  }

  for(int iele = 0; iele < numNeum; ++iele) {
    
    if(dynamic_cast<QuadPressureBC*>(neum[iele]) == NULL && 
       dynamic_cast<TrianglePressureBC*>(neum[iele]) == NULL) continue;

    // Compute structural element distributed Neumann force
    elementPressureForce.zero();
    neum[iele]->neumVector(nodes, elementPressureForce);

    // Assemble element force vector into domain force vector
    int *dofs = neum[iele]->dofs(*c_dsa);
    for(int idof = 0; idof < neum[iele]->numDofs(); ++idof) {
      if(dofs[idof] >= 0)
        force[dofs[idof]] += lambda*elementPressureForce[idof];
    }
    delete [] dofs;
  }

}

template<class Scalar>
void
Domain::addAtddnbForce(GenVector<Scalar> &force, double lambda)
{
  Vector elementAtddnbForce(maxNumDOFs);

  for(int iele = 0; iele < numNeum; ++iele) {
    // Compute acoustic element distributed Neumann force
    elementAtddnbForce.zero();
    neum[iele]->neumVector(nodes, elementAtddnbForce);
    elementAtddnbForce *= sinfo.ATDDNBVal;
    
    // Assemble element force vector into domain force vector
    int *dofs = neum[iele]->dofs(*c_dsa);
    for(int idof = 0; idof < neum[iele]->numDofs(); ++idof) {
      if(dofs[idof] >= 0)
        force[dofs[idof]] += lambda*elementAtddnbForce[idof];
    }
    delete [] dofs;
  }
}

template<class Scalar>
void
Domain::addAtdrobForce(GenVector<Scalar> &force, double lambda)
{
  Vector elementAtdrobForce(maxNumDOFs);

  for(int iele = 0; iele < numScatter; ++iele) {
    // Compute acoustic element Robin distributed boundary condition
    elementAtdrobForce.zero();
    scatter[iele]->neumVector(nodes, elementAtdrobForce);
    elementAtdrobForce *= (sinfo.ATDROBVal/sinfo.ATDROBalpha);

    // Assemble element force vector into domain force vector
    int *dofs = scatter[iele]->dofs(*c_dsa);
    for(int idof = 0; idof < scatter[iele]->numDofs(); ++idof) {
      if(dofs[idof] >= 0)
        force[dofs[idof]] += lambda*elementAtdrobForce[idof];
    }
    delete [] dofs;
  }
}

template<class Scalar>
void
Domain::addThermalForce(GenVector<Scalar> &force)
{
  if(!temprcvd) initNodalTemperatures(); // XXXX to be moved
  Vector elementTemp(maxNumNodes);
  Vector elementThermalForce(maxNumDOFs);
  if(!elemToNode) elemToNode = new Connectivity(&packedEset);

  for(int iele = 0; iele < numele;  ++iele) {
    // By convention phantom elements do not have thermal load
    if(packedEset[iele]->getProperty() == 0) continue;

    // Extract the element nodal temperatures from temprcvd and/or element property ambient temperature
    for(int inod = 0; inod < elemToNode->num(iele); ++inod) {
      double t = temprcvd[(*elemToNode)[iele][inod]];
      elementTemp[inod] = (t == defaultTemp) ? packedEset[iele]->getProperty()->Ta : t;
    }

    // Compute element thermal force in the global coordinates
    elementThermalForce.zero();
    packedEset[iele]->getThermalForce(nodes, elementTemp, elementThermalForce, 0);

    // Assemble element thermal forces into the force vector
    for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
      int dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(dofNum >= 0)
        force[dofNum] += elementThermalForce[idof];
    }
  }
}

template<class Scalar>
void
Domain::addMpcRhs(GenVector<Scalar> &force)
{
  Vector elementForce(maxNumDOFs);

  for(int iele = 0; iele < numele; ++iele) {
    // If there is element is not a MpcElement, skip it.
    if(!packedEset[iele]->isMpcElement()) continue; // this also works for superelements

    // Otherwise, compute element force due to mpc rhs
    packedEset[iele]->computePressureForce(nodes, elementForce, (GeomState *) 0, 0);

    // Assemble element pressure forces into domain force vector
    for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(cn >= 0)
        force[cn] += elementForce[idof];
    }
  }
}

template<class Scalar>
void
Domain::scaleDisp(Scalar *u)
{
  // PJSA: 9-22-06 this is for multi-point pade with coupled fluid-structure
  for(int inode = 0; inode < numnodes; ++inode) {
    int cdofs[6];
    c_dsa->number(inode, DofSet::XYZdisp | DofSet::XYZrot, cdofs);
    for(int jdof = 0; jdof<6; ++jdof)
      if(cdofs[jdof] >= 0) u[cdofs[jdof]] *= coupledScaling;
  }
}


template<class Scalar>
void
Domain::scaleInvDisp(Scalar *u)
{
  // RT: 9-10-09 this is for multi-point projection with coupled fluid-structure
  for(int inode = 0; inode < numnodes; ++inode) {
    int cdofs[6];
    c_dsa->number(inode, DofSet::XYZdisp | DofSet::XYZrot, cdofs);
    for(int jdof = 0; jdof<6; ++jdof)
      if(cdofs[jdof] >= 0) u[cdofs[jdof]] /= coupledScaling;
  }
}


template<class Scalar>
int Domain::mergeDistributedDisp(Scalar (*xyz)[11], Scalar *u, Scalar *bcx)//DofSet::max_known_nonL_dof
{
  // PJSA 9-22-06 u is already scaled
  int inode, nodeI;
  int realNode = -1;
  for (inode = 0; inode < numnodes; ++inode){

    if(nodeToElem)
      if(nodeToElem->num(inode) <= 0) continue;
    realNode++;
    nodeI = (outFlag) ? realNode : inode;

    int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);
    int xLoc1 =   dsa->locate(inode, DofSet::Xdisp);

    if (xLoc >= 0)
      xyz[nodeI][0] = u[xLoc];          // free
    else if (xLoc1 >= 0 && bcx)
      xyz[nodeI][0] = bcx[xLoc1];       // constrained
    else
      xyz[nodeI][0] = 0.0;

    int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);
    int yLoc1 =   dsa->locate(inode, DofSet::Ydisp);

    if (yLoc >= 0)
      xyz[nodeI][1] = u[yLoc];
    else if (yLoc1 >= 0 && bcx)
      xyz[nodeI][1] = bcx[yLoc1];
    else
      xyz[nodeI][1] = 0.0;

    int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);
    int zLoc1 =   dsa->locate(inode, DofSet::Zdisp);

    if (zLoc >= 0)
      xyz[nodeI][2] = u[zLoc];
    else if (zLoc1 >= 0 && bcx)
      xyz[nodeI][2] = bcx[zLoc1];
    else
      xyz[nodeI][2] = 0.0;

    int xRot  = c_dsa->locate(inode, DofSet::Xrot);
    int xRot1 =   dsa->locate(inode, DofSet::Xrot);

    if (xRot >= 0)
      xyz[nodeI][3] = u[xRot];
    else if(xRot1 >= 0 && bcx)
      xyz[nodeI][3] = bcx[xRot1];
    else
      xyz[nodeI][3] = 0.0;

    int yRot  = c_dsa->locate(inode, DofSet::Yrot);
    int yRot1 =   dsa->locate(inode, DofSet::Yrot);

    if (yRot >= 0)
      xyz[nodeI][4] = u[yRot];
    else if (yRot1 >= 0 && bcx)
      xyz[nodeI][4] = bcx[yRot1];
    else
      xyz[nodeI][4] = 0.0;

    int zRot  = c_dsa->locate(inode, DofSet::Zrot);
    int zRot1 =   dsa->locate(inode, DofSet::Zrot);

    if (zRot >= 0)
      xyz[nodeI][5] = u[zRot];
    else if (zRot1 >= 0 && bcx)
      xyz[nodeI][5] = bcx[zRot1];
    else
      xyz[nodeI][5] = 0.0;

    int xTemp  = c_dsa->locate(inode, DofSet::Temp);
    int xTemp1 =   dsa->locate(inode, DofSet::Temp);

    if (xTemp >= 0)
      xyz[nodeI][6] = u[xTemp];
    else if (xTemp1 >= 0 && bcx)
      xyz[nodeI][6] = bcx[xTemp1];
    else
      xyz[nodeI][6] = 0.0;

    int xHelm  = c_dsa->locate(inode, DofSet::Helm);
    int xHelm1 =   dsa->locate(inode, DofSet::Helm);

    if (xHelm >= 0)
      xyz[nodeI][7] = u[xHelm];
    else if (xHelm1 >= 0 && bcx)
      xyz[nodeI][7] = bcx[xHelm1];
    else
      xyz[nodeI][7] = 0.0;

    int xPot  = c_dsa->locate(inode, DofSet::Potential);
    int xPot1 =   dsa->locate(inode, DofSet::Potential);

    if (xPot >= 0)
      xyz[nodeI][10] = u[xPot];
    else if (xPot1 >= 0 && bcx)
      xyz[nodeI][10] = bcx[xPot1];
    else
      xyz[nodeI][10] = 0.0;
  }

  return ++realNode;
}


template<class Scalar>
void
Domain::buildRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *kuc)
{
  if(! dynamic_cast<GenSubDomain<Scalar>*> (this))
    checkSommerTypeBC(this);
  // ... ZERO FORCE VECTOR INITIALLY
  force.zero();

  // ... COMPUTE EXTERNAL FORCE FROM REAL NEUMAN BC
  int i;
  for(i=0; i < numNeuman; ++i) {
    int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
    if(dof < 0) continue;
    ScalarTypes::addScalar(force[dof], nbc[i].val);
  }

  // ... COMPUTE EXTERNAL FORCE FROM COMPLEX NEUMAN BC
  for(i=0; i < numComplexNeuman; ++i) {
    int dof  = c_dsa->locate(cnbc[i].nnum, (1 << cnbc[i].dofnum));
    if(dof < 0) continue;
    ScalarTypes::addScalar(force[dof], cnbc[i].reval, cnbc[i].imval);
  }

  if (implicitFlag) {
    int i, iele;
    double *direction = getWaveDirection();
    ComplexVector elementNeumanScatterForce(this->maxNumDOFs,0.0);
    int* edofs = (int*) dbg_alloca(this->maxNumDOFs*4*sizeof(int));
    for(iele=0; iele<numNeum; ++iele) {
      double kappa = neum[iele]->el->getProperty()->kappaHelm; // PJSA 1-15-2008
      neum[iele]->dofs(*this->dsa,edofs);
      neum[iele]->neumVector(this->nodes,elementNeumanScatterForce,
                             kappa,direction[0],direction[1],direction[2],
                             pointSourceFlag);
      for(i=0;i<neum[iele]->numDofs();i++) {
        int cn = c_dsa->getRCN(edofs[i]);
        if(cn >= 0) {
          ScalarTypes::addScalar(force[cn], elementNeumanScatterForce[i].real(),
                                 elementNeumanScatterForce[i].imag());
        }
      }
    }

    ComplexVector elementWetInterfaceScatterForce(this->maxNumDOFs,0.0);
    int numWetDof;
    for(iele=0; iele<numWet; ++iele) {
      HelmElement *he = dynamic_cast<HelmElement *>(wet[iele]->el);
      if (he==0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numSolidDofs();
        wet[iele]->solidDofs(*this->dsa,edofs);
      } else if (he!=0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numDofs();
        wet[iele]->dofs(*this->dsa,edofs);
      } else {
        numWetDof = wet[iele]->numWetDofs();
        wet[iele]->wetDofs(*this->dsa,edofs);
      }

      complex<double> kappaw =
        sqrt(geoSource->shiftVal())/wet[iele]->soundSpeed;
      wet[iele]->wetInterfaceVector(this->nodes,
                        elementWetInterfaceScatterForce, real(kappaw),
                        direction[0], direction[1], direction[2],0,
                        pointSourceFlag);
      for(i=0;i<numWetDof;i++) {
        int cn = c_dsa->getRCN(edofs[i]);
        if(cn >= 0)
          ScalarTypes::addScalar(force[cn],
            elementWetInterfaceScatterForce[i].real(),
            elementWetInterfaceScatterForce[i].imag());
      }
    }
  }

  // ... ADD GRAVITY FORCES
  if(gravityFlag()) addGravityForce<Scalar>(force);

  // ... ADD THERMAL FORCES
  if(thermalFlag() && !sinfo.isNonLin()) addThermalForce<Scalar>(force);

  // ... ADD PRESSURE LOAD
  if(!sinfo.isNonLin()) addPressureForce<Scalar>(force);

  // ... ADD LMPC RHS
  if(lmpc.max_size() && !sinfo.isNonLin()) addMpcRhs<Scalar>(force);

  // scale RHS force for coupled domains
  if(sinfo.isCoupled) {
    int cdofs[6];  DofSet structdofs = DofSet::XYZdisp | DofSet::XYZrot;
    for(i=0; i<numnodes; ++i) {
      c_dsa->number(i, structdofs, cdofs);
      for(int j=0; j<6; ++j)
        if(cdofs[j] > -1) force[cdofs[j]] *= cscale_factor;
    }
  }

  // COMPUTE NON-HOMOGENEOUS FORCE CONTRIBUTION
  // IN THE CASE OF NONLINEAR, NON-HOMOGENEOUS (PRESCRIBED) FORCES
  // ARE TAKEN CARE OF USING THE GEOMSTATE CLASS, NOT BY
  // MODIFYING THE RHS VECTOR
  if(kuc && !sinfo.isNonLin()) {
    GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);

    // CONSTRUCT NON-HOMONGENOUS DIRICHLET BC VECTOR (PRESCRIBED)
    for(i=0; i<numDirichlet; ++i) {
      int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
      if(dof < 0) continue;
      dof = c_dsa->invRCN(dof);
      if(dof >= 0) {
        if(sinfo.isCoupled && dbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof], dbc[i].val/coupledScaling); else // PJSA 1-9-08
        ScalarTypes::initScalar(Vc[dof], dbc[i].val);
      }
    }

    // CONSTRUCT NON-HOMONGENOUS COMPLEX DIRICHLET BC VECTOR
    ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
    for(i=0; i<numComplexDirichlet; ++i) {
      int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
      if(dof2 < 0) continue;
      dof2 = c_dsa->invRCN(dof2);
      if(dof2 >= 0) {
        if(sinfo.isCoupled && cdbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval/coupledScaling, cdbcMRHS[i].imval/coupledScaling); else // PJSA 1-9-08
        ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
      }
    }

    kuc->multSubtract(Vc, force);
  }
}

template<class Scalar>
void
Domain::computeReactionForce(GenVector<Scalar> &fc, GenVector<Scalar> &Vu,
                             GenSparseMatrix<Scalar> *_kuc, GenSparseMatrix<Scalar> *_kcc)
{
  // COMPUTE NON-HOMOGENEOUS FORCE CONTRIBUTION
  GenCuCSparse<Scalar> *kuc = dynamic_cast<GenCuCSparse<Scalar> *>(_kuc);
  if(kuc) kuc->transposeMultNew(Vu.data(), fc.data()); // fc = kuc^T * Vu
  else fc.zero();

  // COMPUTE NON-HOMOGENEOUS FORCE CONTRIBUTION
  GenCuCSparse<Scalar> *kcc = dynamic_cast<GenCuCSparse<Scalar> *>(_kcc);
  if(kcc) {
    GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);

    // CONSTRUCT NON-HOMONGENOUS DIRICHLET BC VECTOR (PRESCRIBED)
    for(int i=0; i<numDirichlet; ++i) {
      int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
      if(dof < 0) continue;
      dof = c_dsa->invRCN(dof);
      if(dof >= 0) {
        if(sinfo.isCoupled && dbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof], dbc[i].val/coupledScaling); else // PJSA 1-9-08
        ScalarTypes::initScalar(Vc[dof], dbc[i].val);
      }
    }

    // CONSTRUCT NON-HOMONGENOUS COMPLEX DIRICHLET BC VECTOR
    ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
    for(int i=0; i<numComplexDirichlet; ++i) {
      int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
      if(dof2 < 0) continue;
      dof2 = c_dsa->invRCN(dof2);
      if(dof2 >= 0) {
        if(sinfo.isCoupled && cdbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval/coupledScaling, cdbcMRHS[i].imval/coupledScaling); else // PJSA 1-9-08
        ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
      }
    }

    kcc->multAddNew(Vc.data(), fc.data());
  }
}

template<class Scalar>
void
Domain::buildRHSForce(GenVector<Scalar> &force, GenVector<Scalar> &tmp,
                      GenSparseMatrix<Scalar> *kuc, GenSparseMatrix<Scalar> *muc, 
                      GenSparseMatrix<Scalar> **cuc_deriv, double omega, double delta_omega, GeomState *gs)
{
  if(! dynamic_cast<GenSubDomain<Scalar>*> (this))
    checkSommerTypeBC(this);
  // ... ZERO FORCE VECTOR INITIALLY
  force.zero();

  // ... COMPUTE EXTERNAL FORCE FROM REAL NEUMAN BC
  int i;
  for(i=0; i < numNeuman; ++i) {
    int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
    if(dof < 0) continue;
    //if(nbc[i].dofnum == 7) ScalarTypes::addScalar(force[dof], nbc[i].val/fluidDensity); else // PJSA 1-9-2008 temp fix (homogenous fluid only) 
    //                                                                                         // required for backward-compatibility with old versions (before acoustic equations divided by rho)
    ScalarTypes::addScalar(force[dof], nbc[i].val); 
  }

  // ... COMPUTE EXTERNAL FORCE FROM COMPLEX NEUMAN BC
  for(i=0; i < numComplexNeuman; ++i) {
    int dof  = c_dsa->locate(cnbc[i].nnum, (1 << cnbc[i].dofnum));
    if(dof < 0) continue;
    //if(cnbc[i].dofnum == 7) ScalarTypes::addScalar(force[dof], cnbc[i].reval/fluidDensity, cnbc[i].imval/fluidDensity); else // PJSA 1-9-2008
    ScalarTypes::addScalar(force[dof], cnbc[i].reval, cnbc[i].imval); 
  }

  // PJSA: new FETI-H stuff
  if (implicitFlag) {
    int i, iele;
    double *direction = getWaveDirection();
// RT
    ComplexVector elementNeumanScatterForce(this->maxNumDOFs,0.0);
    int* edofs = (int*) dbg_alloca(this->maxNumDOFs*sizeof(int));
    for(iele=0; iele<numNeum; ++iele) {
      double kappa = omega/ScalarTypes::Real(neum[iele]->soundSpeed);
    //double kappa = neum[iele]->el->getProperty()->kappaHelm; // PJSA 1-15-2008
      neum[iele]->dofs(*this->dsa,edofs);
      neum[iele]->neumVector(this->nodes,elementNeumanScatterForce,
                             kappa,direction[0],direction[1],direction[2],
                             pointSourceFlag);
      for(i=0;i<neum[iele]->numDofs();i++) { 
       int cn = c_dsa->getRCN(edofs[i]);
       if(cn >= 0) {
        ScalarTypes::addScalar(force[cn], elementNeumanScatterForce[i].real(),
                               elementNeumanScatterForce[i].imag());
       }
      }
    }
  }

  if (implicitFlag) {

    double *direction = getWaveDirection();

    int iele;
    ComplexVector elementWetInterfaceScatterForce(this->maxNumDOFs,0.0);
    int* edofs = (int*) alloca(this->maxNumDOFs*4*sizeof(int));
    int numWetDof;
    for(iele=0; iele<numWet; ++iele) {
      HelmElement *he = dynamic_cast<HelmElement *>(wet[iele]->el);
      if (he==0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numSolidDofs();
        wet[iele]->solidDofs(*this->dsa,edofs);
      } else if (he!=0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numDofs();
        wet[iele]->dofs(*this->dsa,edofs);
      } else {
        numWetDof = wet[iele]->numWetDofs();
        wet[iele]->wetDofs(*this->dsa,edofs);
      }

      complex<double> kappaw =
        omega/wet[iele]->soundSpeed;
      wet[iele]->wetInterfaceVector(this->nodes,
                        elementWetInterfaceScatterForce, real(kappaw),
                        direction[0], direction[1], direction[2],0,
                        pointSourceFlag);
      int i; 
      for(i=0;i<numWetDof;i++) {
       int cn = c_dsa->getRCN(edofs[i]);
       if(cn >= 0)
        ScalarTypes::addScalar(force[cn],
          elementWetInterfaceScatterForce[i].real(),
          elementWetInterfaceScatterForce[i].imag());
      }
    }
  }

  // ... ADD GRAVITY FORCES
  if(gravityFlag()) addGravityForce<Scalar>(force);

  // ... ADD THERMAL FORCES
  if(thermalFlag() && !sinfo.isNonLin()) addThermalForce<Scalar>(force);

  // ... ADD PRESSURE LOAD
  if(!sinfo.isNonLin()) addPressureForce<Scalar>(force);

  GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);

  // CONSTRUCT NON-HOMONGENOUS DIRICHLET BC VECTOR (PRESCRIBED)
  for(i=0; i<numDirichlet; ++i) {
    int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
    if(dof < 0) continue;
    dof = c_dsa->invRCN(dof);
    if(dof >= 0) {
      if(sinfo.isCoupled && dbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof], dbc[i].val/coupledScaling); else // PJSA 1-9-08
      ScalarTypes::initScalar(Vc[dof], dbc[i].val);
    }
  }

  // CONSTRUCT NON-HOMONGENOUS COMPLEX DIRICHLET BC VECTOR
  ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
  for(i=0; i<numComplexDirichlet; ++i) {
    int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
    if(dof2 < 0) continue;
    dof2 = c_dsa->invRCN(dof2);
    if(dof2 >= 0) {
      if(sinfo.isCoupled && cdbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval/coupledScaling, cdbcMRHS[i].imval/coupledScaling); else // PJSA 1-9-08
      ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
    }
  }

  // scale RHS force for coupled domains
  if(sinfo.isCoupled) {
    int cdofs[6];  DofSet structdofs = DofSet::XYZdisp | DofSet::XYZrot;
    for(i=0; i<numnodes; ++i) {
      c_dsa->number(i, structdofs, cdofs);
      for(int j=0; j<6; ++j) 
        if(cdofs[j] > -1) force[cdofs[j]] *= cscale_factor;
    }
  }

  // COMPUTE NON-HOMOGENEOUS FORCE CONTRIBUTION
  // IN THE CASE OF NONLINEAR, NON-HOMOGENEOUS (PRESCRIBED) FORCES
  // ARE TAKEN CARE OF USING THE GEOMSTATE CLASS, NOT BY
  // MODIFYING THE RHS VECTOR
  if(probType() != SolverInfo::NonLinStatic &&
     probType() != SolverInfo::NonLinDynam  &&
     probType() != SolverInfo::ArcLength) {
    if(kuc) {
      kuc->multSubtract(Vc, force);
    }
    if (muc) {
      tmp.zero();
      muc->multSubtract(Vc, tmp);
      tmp *= (omega-delta_omega)*(omega-delta_omega)-omega*omega;
      force += tmp;
    }
    if(cuc_deriv) if (cuc_deriv[0]) {
      tmp.zero();
      cuc_deriv[0]->multSubtract(Vc, tmp);
      Scalar c;
      ScalarTypes::initScalar(c,0.0,delta_omega);
      tmp *= c;
      force += tmp;
    }
  }
}

template<class Scalar>
void
Domain::buildFreqSweepRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *muc,
                               GenSparseMatrix<Scalar> **cuc_deriv, int iRHS, double omega)
{
  if(iRHS < 1) return; // this shouldn't happen

  // PJSA 5-15-06
  if (implicitFlag) {
    double *direction = getWaveDirection();
    ComplexVector elementNeumanScatterForce(this->maxNumDOFs,0.0);
    for(int iele=0; iele<numNeum; ++iele) {
      int *dofs = neum[iele]->dofs(*this->dsa);
      double kappa = neum[iele]->el->getProperty()->kappaHelm; // PJSA 1-15-2008
      neum[iele]->neumVectorDeriv(this->nodes,elementNeumanScatterForce,
                                  kappa,direction[0],direction[1],direction[2],iRHS, pointSourceFlag);
      for(int i=0;i<neum[iele]->numDofs();i++) {
       int cn = c_dsa->getRCN(dofs[i]);
       if(cn >= 0) {
        ScalarTypes::addScalar(force[cn],
             elementNeumanScatterForce[i].real() /
               pow(real(neum[iele]->soundSpeed),iRHS),
             elementNeumanScatterForce[i].imag() /
               pow(real(neum[iele]->soundSpeed),iRHS));
       }
      }
    }
  }

/*
  if (implicitFlag) {

    double *direction = getWaveDirection();
// RT - check that this is the heterogenous scattering
    int imat;
    double rho1, rho2;
    complex<double> kappa1(0.0,0.0), kappa2(0.0,0.0);
    for(imat=0;imat<geoSource->getNumProps();imat++) {
      StructProp *sp = &((geoSource->getStructProps())[imat]);
      if (sp!=0)
      if (sp->kappaHelm!=0.0) {
        if (kappa1==complex<double>(0.0,0.0)) {
          kappa1 = complex<double>(sp->kappaHelm,0.0);
          rho1 = sp->rho;
        } else {
          if (kappa1!=complex<double>(sp->kappaHelm,sp->kappaHelmImag)) {
            if (kappa2==complex<double>(0.0,0.0)) {
              kappa2 = complex<double>(sp->kappaHelm,sp->kappaHelmImag);
              rho2 = sp->rho;
            }
          }
        }
      }
    }

    complex<double> diri[3][3] ={ {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
    complex<double> coefi[3] = {0.0,0.0,0.0};
    complex<double> kappai[3] = {0.0,0.0,0.0};

    if (real(kappa2)!=0) {
fprintf(stderr,"Complex kappa2 not supported yet!\n");
// Assumptions - fluid 1 is in z<0, fluid 2 is in z>0
//             - incident wave is of form (dx,dy,dz), dz>0
//             - incident wave vector is normalized
//             - kappa1 is real
//             - |kappa2| < kappa1

       double alpha = acos(direction[2]);

       double ll = sqrt(direction[0]*direction[0]+direction[1]*direction[1]);
       double dx = 0.0;
       double dy = 0.0;
       if (ll!=0.0) {
         dx = direction[0]/ll;
         dy = direction[1]/ll;
       }

       diri[0][0] = kappa1*complex<double>(0.0,direction[0]);
       diri[0][1] = kappa1*complex<double>(0.0,direction[1]);
       diri[0][2] = kappa1*complex<double>(0.0,direction[2]);
       diri[1][0] = kappa1*complex<double>(0.0,direction[0]);
       diri[1][1] = kappa1*complex<double>(0.0,direction[1]);
       diri[1][2] = kappa1*complex<double>(0.0,-direction[2]);
       diri[2][0] = kappa1*complex<double>(0.0,sin(alpha)*dx);
       diri[2][1] = kappa1*complex<double>(0.0,sin(alpha)*dy);
       diri[2][2] = sqrt(sin(alpha)*sin(alpha)*kappa1*kappa1-kappa2*kappa2);
       if (real(diri[2][2])>1e-15) diri[2][2] = -diri[2][2];

       coefi[0] = complex<double>(1.0,0.0);
       coefi[1] = (rho2*kappa1/(rho1*kappa2)-
                  diri[2][2]/(kappa2*complex<double>(0.0,direction[2]))) /
                  (rho2*kappa1/(rho1*kappa2)+
                  diri[2][2]/(kappa2*complex<double>(0.0,direction[2])));
       coefi[2] = 2.0 /
                  (rho2*kappa1/(rho1*kappa2)+
                  diri[2][2]/(kappa2*complex<double>(0.0,direction[2])));
       coefi[2] *= rho2/rho1*kappa1/kappa2;
    } else {
       diri[2][0] = kappa1*complex<double>(0.0,direction[0]);
       diri[2][1] = kappa1*complex<double>(0.0,direction[1]);
       diri[2][2] = kappa1*complex<double>(0.0,direction[2]);
       coefi[2] = cscale_factor*complex<double>(1.0,0.0);
       kappai[2] = kappa1;
    }

    int iele;
    ComplexVector elementWetInterfaceScatterForce(this->maxNumDOFs,0.0);
    int* edofs = (int*) alloca(this->maxNumDOFs*4*sizeof(int));
    int numWetDof;
    for(iele=0; iele<numWet; ++iele) {
      HelmElement *he = dynamic_cast<HelmElement *>(wet[iele]->el);
      if (he==0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numSolidDofs();
        wet[iele]->solidDofs(*this->dsa,edofs);
      } else if (he!=0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numDofs();
        wet[iele]->dofs(*this->dsa,edofs);
      } else {
        numWetDof = wet[iele]->numWetDofs();
        wet[iele]->wetDofs(*this->dsa,edofs);
      }

      wet[iele]->wetInterfaceVectorDeriv(this->nodes,
         elementWetInterfaceScatterForce, diri,coefi,kappai,iRHS);

      int i;
      for(i=0;i<numWetDof;i++) {
        ScalarTypes::addScalar(force[edofs[i]],
          elementWetInterfaceScatterForce[i].real()/pow(domain->fluidCelerity,iRHS),
          elementWetInterfaceScatterForce[i].imag()/pow(domain->fluidCelerity,iRHS));
      }
    }
  }
*/

  if (implicitFlag) {

    double *direction = getWaveDirection();

    int iele;
    ComplexVector elementWetInterfaceScatterForce(this->maxNumDOFs,0.0);
    int* edofs = (int*) alloca(this->maxNumDOFs*4*sizeof(int));
    int numWetDof;
    for(iele=0; iele<numWet; ++iele) {
      HelmElement *he = dynamic_cast<HelmElement *>(wet[iele]->el);
      if (he==0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numSolidDofs();
        wet[iele]->solidDofs(*this->dsa,edofs);
      } else if (he!=0 && wet[iele]->el2==0) {
        numWetDof = wet[iele]->numDofs();
        wet[iele]->dofs(*this->dsa,edofs);
      } else {
        numWetDof = wet[iele]->numWetDofs();
        wet[iele]->wetDofs(*this->dsa,edofs);
      }

      complex<double> kappaw =
        sqrt(geoSource->shiftVal())/wet[iele]->soundSpeed;
      wet[iele]->wetInterfaceVector(this->nodes,
                        elementWetInterfaceScatterForce, real(kappaw),
                        direction[0], direction[1], direction[2],iRHS,
                        pointSourceFlag);
      if (he==0 && wet[iele]->el2==0) {
        for(int i=0;i<numWetDof;i++)
          elementWetInterfaceScatterForce[i] *= cscale_factor;
      } else if (he!=0 && wet[iele]->el2==0) {
      } else {
        for(int i=0;i<numWetDof/4;i++) {
          elementWetInterfaceScatterForce[4*i+0] *= cscale_factor;
          elementWetInterfaceScatterForce[4*i+1] *= cscale_factor;
          elementWetInterfaceScatterForce[4*i+2] *= cscale_factor;
        }
      }

      int i;
      for(i=0;i<numWetDof;i++) {
       int cn = c_dsa->getRCN(edofs[i]);
       if(cn >= 0)
        ScalarTypes::addScalar(force[cn],
          elementWetInterfaceScatterForce[i].real()/
           pow(real(wet[iele]->soundSpeed),iRHS),
          elementWetInterfaceScatterForce[i].imag()/
           pow(real(wet[iele]->soundSpeed),iRHS));
      }
    }
  }

  int i;
  GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);
  GenVector<Scalar> Vc_tmp(numDirichlet+numComplexDirichlet);

  // CONSTRUCT NON-HOMONGENOUS DIRICHLET BC VECTOR (PRESCRIBED)
  for(i=0; i<numDirichlet; ++i) {
    int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
    if(dof < 0) continue;
    dof = c_dsa->invRCN(dof);
    if(dof >= 0) ScalarTypes::initScalar(Vc[dof], dbc[i].val);
  }

  // CONSTRUCT NON-HOMONGENOUS COMPLEX DIRICHLET BC VECTOR
  ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
  for(i=0; i<numComplexDirichlet; ++i) {
    int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
    if(dof2 < 0) continue;
    dof2 = c_dsa->invRCN(dof2);
    if(dof2 >= 0) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
  }

  if(muc && iRHS <= 2) {
    if(iRHS==1) Vc_tmp.linC(Vc, -2.0*omega);
    else Vc_tmp.linC(Vc,-2.0);
    muc->multSubtract(Vc_tmp, force);
  }

  if(cuc_deriv && cuc_deriv[iRHS-1]) cuc_deriv[iRHS-1]->multSubtract(Vc, force);
}

template<class Scalar>
void
Domain::assembleATDROB(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops,double Kcoef)
{
//  fprintf(stderr, "ATDROB ok\n");
//  checkSommerTypeBC(this);
  double *v = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
  if(numScatter > 0) {
    int i;
    double temp = sinfo.ATDROBbeta/sinfo.ATDROBalpha;
    for(i=0; i<numScatter; i++) {
     int *dofs = scatter[i]->dofs(*dsa);
     FullSquareMatrix ms = scatter[i]->sommerMatrix(nodes,v);//sommerMatrix is negative definite...
     FullSquareMatrix mm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
     ms.multiply(mm,temp);
     updateMatrices<Scalar>(ops,K,dofs,&mm,0,Kcoef);
    }
  }
}

template<class Scalar>
void
Domain::assembleSommer(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops)
{
 checkSommerTypeBC(this); // TODO check

 if(numSommer > 0) {
   if(sommer[0]->dim() == 3) {
     if(sommerfeldType==1) sommerfeldType = 3;
     if(sommerfeldType==2) sommerfeldType = 4;
   }

   if(sommerfeldType == 1 )
     fprintf(stderr, " ... 1st order Bayliss-Turkel       ...\n");
   else if(sommerfeldType == 2 )
     fprintf(stderr, " ... 2nd order Bayliss-Turkel       ...\n");
   else if(sommerfeldType == 3 )
     fprintf(stderr, " ... 1st order 3D Bayliss-Turkel    ...\n");
   else if(sommerfeldType == 4 )
     fprintf(stderr, " ... 2nd order 3D Bayliss-Turkel    ...\n");

   if(sommerfeldType == 1 || sommerfeldType == 2)
     getCurvatures(this);
   if(curvatureFlag != 1) {
     if(sommerfeldType == 3 || sommerfeldType == 4)
       getCurvatures3D(this);
   }

   double *v = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
   double *vbt = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
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
     double ss = real(sommer[i]->el->getProperty()->soundSpeed);

     ms.multiply(mm,kappa);

     // "Zero-order" term
     updateMatrices(ops,K,dofs,0,&mm);

     double psi; // curvature of the boundary
     double HH, KK;
     // 1st order Bayliss-Turkel boundary condition
     if(sommerfeldType == 1 ) {
       psi = curvatures[i];
       HH = psi/2.0;
       ms.multiply(mm,-HH);
       updateMatrices(ops,K,dofs,&mm,0);
     }

     // 2nd order Bayliss-Turkel boundary condition
     else if(sommerfeldType == 2) {
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
       updateMatrices(ops,K,dofs,&mm,0);
       // mm = ms*c1 + ks*c3;
       ms.multiply(mm,c1);
       ks.multiply(mm1,c3);
       mm += mm1;
       updateMatrices(ops,K,dofs,0,&mm);
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
       updateMatrices(ops,K,dofs,&mm,0);
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
         updateMatrices(ops,K,dofs,&mm,0);
         ms.multiply(mm,c1);
         updateMatrices(ops,K,dofs,0,&mm);
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
         for (iDof=0; iDof<ms.dim(); iDof++) {
           for (jDof=0; jDof<ms.dim(); jDof++) {
             ksRe[iDof][jDof] = real(bt2Matrix[kDof]);
             ksIm[iDof][jDof] = imag(bt2Matrix[kDof]);
             kDof++;
           }
         }

         updateMatrices(ops,K,dofs,&ksRe,&ksIm);
       }

       double curv_e[3];
       double curv_f[3];
       double curv_g[3];
       double tau1[3][3];
       double tau2[3][3];
       if (curvatureFlag != 1) {
         int n1 = (*somElemToNode)[i][0];
         int n2 = (*somElemToNode)[i][1];
         int n3 = (*somElemToNode)[i][2];
         int nn[3];
         nn[0] = nodeToSommerNodeMap[n1];
         nn[1] = nodeToSommerNodeMap[n2];
         nn[2] = nodeToSommerNodeMap[n3];
         double *curv_normal[3];
         int iiNode;
         for (iiNode=0; iiNode<3; iiNode++) {
           int curvNode = nn[iiNode];
           curv_e[iiNode]      = curvatures_e[curvNode];
           curv_f[iiNode]      = curvatures_f[curvNode];
           curv_g[iiNode]      = curvatures_g[curvNode];
           curv_normal[iiNode] = curvatures_normal[nn[iiNode]];
           getTau(curv_normal[iiNode], tau1[iiNode], tau2[iiNode]);

           // This is for experiments with exact curvatures for ellipsoid
           if (curvatureFlag == 2) {
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
         for (iiNode=0; iiNode<3; iiNode++) {
           double nrmal[3];
           Node nd = nodes.getNode(nds[iiNode]);
           double l = sqrt(nd.x*nd.x+nd.y*nd.y+nd.z*nd.z);
           nrmal[0] = nd.x/l;
           nrmal[1] = nd.y/l;
           nrmal[2] = nd.z/l;
           getTau(nrmal, tau1[iiNode], tau2[iiNode]);
         }
       }
       if (curvatureFlag != 1) {
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
       if (curvatureFlag == 2) {
         sommer[i]->ellipsoidBT2(nodes, curvatureConst2, curvatureConst1, kappa, bt2Matrix);
       }

       int kDof=0;
       int iDof, jDof;
       for (iDof=0; iDof<ms.dim(); iDof++) {
         for (jDof=0; jDof<ms.dim(); jDof++) {
           ksRe[iDof][jDof] = real(-1.0/(2.0*ii*kappa)*bt2Matrix[kDof]);
           ksIm[iDof][jDof] = imag(-1.0/(2.0*ii*kappa)*bt2Matrix[kDof]);
           kDof++;
         }
       }

       updateMatrices(ops,K,dofs,&ksRe,&ksIm);
     }

     if(solInfo().doFreqSweep) {
       switch(sommerfeldType) {
         case 0:
         case 1:
         case 3:
           updateDampingMatrices(ops, dofs, 0, &ms, ss, 1); // zero- and first-order sommerfeld
           break;
         case 4:
           computeSommerDerivatives(HH, KK, curvatureFlag, dofs, ms, bt2nMatrix, kappa, ss, ops); // 3D second-order sommerfeld
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
Domain::computeSommerDerivatives(double HH, double KK, int curvatureFlag, int *dofs, FullSquareMatrix &ms,
                                 DComplex **bt2nMatrix, double kappa, double ss, AllOps<Scalar> *ops)
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
      updateDampingMatrices(ops,dofs,&mm,0,ss,n);
      ms.multiply(mm,imag(cm)); // mm = imaginary part of cm*ms
      updateDampingMatrices(ops,dofs,0,&mm,ss,n);
    }
    else cerr << " *** WARNING: 3D 2nd order Sommerfeld with curvatureFlag 2 is not supported for frequency sweep \n";
  }

/* PJSA DEBUG    
  GenFullSquareMatrix<DComplex> bt2(ms.dim(),(DComplex*)dbg_alloca(ms.dim()*ms.dim()*sizeof(DComplex)));
  int kDof=0;
  for(int iDof=0; iDof<ms.dim(); iDof++) {
    for(int jDof=0; jDof<ms.dim(); jDof++) {
      bt2[iDof][jDof] = bt2Matrix[kDof];
      kDof++;
    }
  }
  FullSquareMatrix ksRe(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
  FullSquareMatrix ksIm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));

  GenFullSquareMatrix<DComplex> bt2n(bt2,1.0);
  for(int n=1; n<=N; ++n) {
    GenFullSquareMatrix<DComplex> copy(bt2n, 1.0);
    copy.multiply(bt2,bt2n);  // bt2n = bt2n*bt2
    DComplex cm = -pow(1.0/(2.0*ii*kappa),n+1)*pow(-2.0*ii,n)*double(DFactorial(n));
    for(int iDof=0; iDof<ms.dim(); iDof++) {
      for(int jDof=0; jDof<ms.dim(); jDof++) {
        DComplex cmbt = cm*bt2n[iDof][jDof];
        ksRe[iDof][jDof] = real(cmbt);
        ksIm[iDof][jDof] = imag(cmbt);
      }
    }
    updateDampingMatrices(ops,dofs,&ksRe,&ksIm,ss,n);
  }
*/
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
    updateDampingMatrices(ops,dofs,&ksRe,&ksIm,ss,n);
  }
}

template<class Scalar>
void
Domain::updateMatrices(AllOps<Scalar> *ops, GenSparseMatrix<Scalar> *Z, int *dofs,
                       FullSquareMatrix *reEl, FullSquareMatrix *imEl, double Kcoef)
{
  if((sinfo.isATDARB()) || (sinfo.ATDROBalpha != 0.0)) {
    if(reEl) {
      FullSquareMatrix temp(reEl->dim(),(double*)dbg_alloca(reEl->dim()*reEl->dim()*sizeof(double)));
      reEl->multiply(temp, Kcoef);
      if(Z) Z->add(temp, dofs);
      if(ops && ops->spp) ops->spp->add(temp, dofs);
      if(ops && ops->Kuc) ops->Kuc->add(temp, dofs);
      if(ops && ops->Kcc) ops->Kcc->add(temp, dofs);
    }
    if(imEl) {
      FullSquareMatrix temp(imEl->dim(),(double*)dbg_alloca(imEl->dim()*imEl->dim()*sizeof(double)));
      imEl->multiply(temp, Kcoef);
      if(Z) Z->addImaginary(*imEl, dofs);
      if(ops && ops->spp) ops->spp->addImaginary(temp, dofs);
      if(ops && ops->Kuc) ops->Kuc->addImaginary(temp, dofs);
      if(ops && ops->Kcc) ops->Kcc->addImaginary(temp, dofs);
    }
  }
  else {
    if(reEl) {
      if(Z) Z->add(*reEl, dofs);
      if(ops && ops->spp) ops->spp->add(*reEl, dofs);
      if(ops && ops->Kuc) ops->Kuc->add(*reEl, dofs);
      if(ops && ops->Kcc) ops->Kcc->add(*reEl, dofs);
    }
    if(imEl) {
      if(Z) Z->addImaginary(*imEl, dofs);
      if(ops && ops->spp) ops->spp->addImaginary(*imEl, dofs);
      if(ops && ops->Kuc) ops->Kuc->addImaginary(*imEl, dofs);
      if(ops && ops->Kcc) ops->Kcc->addImaginary(*imEl, dofs);
    }
  }
}

template<class Scalar>
void
Domain::updateDampingMatrices(AllOps<Scalar> *ops, int *dofs, FullSquareMatrix *reEl,
                              FullSquareMatrix *imEl, double ss, int n)
{
  if(reEl) {
    *reEl /= pow(ss,n);
    if(ops && ops->C_deriv && ops->C_deriv[n-1]) ops->C_deriv[n-1]->add(*reEl, dofs);
    if(ops && ops->Cuc_deriv && ops->Cuc_deriv[n-1]) ops->Cuc_deriv[n-1]->add(*reEl, dofs);
  }
  if(imEl) {
    *imEl /= pow(ss,n);
    if(ops && ops->C_deriv && ops->C_deriv[n-1]) ops->C_deriv[n-1]->addImaginary(*imEl, dofs);
    if(ops && ops->Cuc_deriv && ops->Cuc_deriv[n-1]) ops->Cuc_deriv[n-1]->addImaginary(*imEl, dofs);
  }
}

//-------------------------------------------------------------------------------------
template<class Scalar>
int Domain::processDispTypeOutputs(OutputInfo &oinfo, Scalar (*glDisp)[11], int numNodes,
                                   int i, double time, double freq, int printFlag)  {

  int success = 0;
  double tag;
  switch (oinfo.type)  {

    case OutputInfo::Displacement:  
      tag = time;
      success = 1;
    case OutputInfo::EigenPair:  {
      if (success == 0)
        tag = freq;
/* this isn't necessary
      Scalar (*data)[3] = new Scalar[numNodes][3];
      for (int jj = 0; jj < numNodes; jj++)  {
        data[jj][0] = glDisp[jj][0]; 
        data[jj][1] = glDisp[jj][1]; 
        data[jj][2] = glDisp[jj][2]; 
      }
      if (oinfo.nodeNumber == -1)  // all nodes
        geoSource->outputNodeVectors(i, data, numNodes, tag);
      else    // one node
        geoSource->outputNodeVectors(i, &(data[oinfo.nodeNumber]), 1, tag);
      delete [] data;
*/
      if (oinfo.nodeNumber == -1)  // all nodes
        geoSource->outputNodeVectors(i, glDisp, numNodes, tag);
      else    // one node
        geoSource->outputNodeVectors(i, &(glDisp[oinfo.nodeNumber]), 1, tag);
      success = 1;
    }
      break;
    case OutputInfo::Disp6DOF:
      tag = time;
      success = 1;
    case OutputInfo::EigenPair6:  {
     if (success == 0)
        tag = freq;
/* this isn't necessary
      Scalar (*data)[6] = new Scalar[numNodes][6];
      for (int jj = 0; jj < numNodes; jj++)  
        for (int kk = 0; kk < 6; kk++) 
          data[jj][kk] = glDisp[jj][kk];
      
      if (oinfo.nodeNumber == -1)
        geoSource->outputNodeVectors6(i, data, numNodes, tag);
      else
        geoSource->outputNodeVectors6(i, &(data[oinfo.nodeNumber]), 1, tag);
      delete [] data;
*/
      if (oinfo.nodeNumber == -1)
        geoSource->outputNodeVectors6(i, glDisp, numNodes, tag);
      else
        geoSource->outputNodeVectors6(i, &(glDisp[oinfo.nodeNumber]), 1, tag);
      success = 1;
    }
      break;
    default: 
      break;
  }

  Scalar *globVal = 0;
  if (success == 0)  {
    globVal = new Scalar[numNodes];
    int dof = -1;

    switch (oinfo.type)  {

      case OutputInfo::DispX:
        if(dof==-1) dof = 0;
      case OutputInfo::DispY:
        if(dof==-1) dof = 1;
      case OutputInfo::DispZ:
        if(dof==-1) dof = 2;
      case OutputInfo::RotX:
        if(dof==-1) dof = 3;
      case OutputInfo::RotY:
        if(dof==-1) dof = 4;
      case OutputInfo::RotZ:
        if(dof==-1) dof = 5;
      case OutputInfo::Temperature:
        if(dof==-1) dof = 6;

        for (int iNode=0; iNode<numNodes; ++iNode)
          globVal[iNode] = glDisp[iNode][dof];
        success = 1;
      case OutputInfo::DispMod:
        if (success == 0)  {
          for (int iNode = 0; iNode < numNodes; ++iNode)
            globVal[iNode] = ScalarTypes::sqrt(glDisp[iNode][0]*glDisp[iNode][0] +
                                               glDisp[iNode][1]*glDisp[iNode][1] +
                                               glDisp[iNode][2]*glDisp[iNode][2]);
          success = 1;
        }
      case OutputInfo::RotMod:
        if (success == 0)  {
          for (int iNode = 0; iNode < numNodes; ++iNode)
            globVal[iNode] = ScalarTypes::sqrt(glDisp[iNode][3]*glDisp[iNode][3] +
                                               glDisp[iNode][4]*glDisp[iNode][4] +
                                               glDisp[iNode][5]*glDisp[iNode][5]);
          success = 1;
        }
      case OutputInfo::TotMod:
        if (success == 0)  {
          for (int iNode = 0; iNode < numNodes; ++iNode)
            globVal[iNode] = ScalarTypes::sqrt(glDisp[iNode][0]*glDisp[iNode][0] +
                                               glDisp[iNode][1]*glDisp[iNode][1] +
                                               glDisp[iNode][2]*glDisp[iNode][2] +
                                               glDisp[iNode][3]*glDisp[iNode][3] +
                                               glDisp[iNode][4]*glDisp[iNode][4] +
                                               glDisp[iNode][5]*glDisp[iNode][5]);
          success = 1;
        }
        if (oinfo.nodeNumber == -1)
          geoSource->outputNodeScalars(i, globVal, numNodes, time);
        else
          geoSource->outputNodeScalars(i, &(globVal[oinfo.nodeNumber]), 1, time);
        break;
      case OutputInfo::AcousticPressure:
        tag = time;
        if (dof==-1) dof = 7;
      case OutputInfo::EigenPressure:
      case OutputInfo::HelmholtzModes:
      case OutputInfo::Helmholtz:
        if (dof==-1) { dof = 7; tag = freq; }
     
        for (int iNode=0; iNode<numNodes; ++iNode) 
          globVal[iNode] = glDisp[iNode][dof];
      
        success = 1;
      case OutputInfo::EigenSlosh:
        if (success == 0)  {  
          if (dof == -1) dof = 10;
          tag = freq;
          for (int iNode=0; iNode<numNodes; ++iNode)
            globVal[iNode] = glDisp[iNode][dof];
          success = 1;
        }
        if(oinfo.nodeNumber == -1)
          geoSource->outputNodeScalars(i, globVal, numNodes, tag);
        else
          geoSource->outputNodeScalars(i, &(globVal[oinfo.nodeNumber]), 1, tag);
        break;
      default:
        break;
    }
    delete [] globVal;
  }

  return success;
}

template<class Scalar>
int Domain::processOutput(OutputInfo::Type &type, GenVector<Scalar> &d_n, Scalar *bcx, int i,
                          double time, double freq, int printFlag)  {

  int success = 1;
  switch(type) {
    case OutputInfo::SloshDisplacement:
      getSloshDispAll(d_n, bcx, i, freq);
      break;
    case OutputInfo::SloshDispX:
      getSloshDisp(d_n, bcx, i, SLDX, freq);
      break;
    case OutputInfo::SloshDispY:
      getSloshDisp(d_n, bcx, i, SLDY, freq);
      break;
    case OutputInfo::SloshDispZ:
      getSloshDisp(d_n, bcx, i, SLDZ, freq);
      break;
    case OutputInfo::YModulus:
      getElementAttr(i,YOUNG);
      break;
    case OutputInfo::MDensity:
      getElementAttr(i,MDENS);
      break;
    case OutputInfo::Thicknes:
      getElementAttr(i,THICK);
      break;
    case OutputInfo::StressXX:
      getStressStrain(d_n,bcx,i,SXX, time);
      break;
    case OutputInfo::StressYY:
      getStressStrain(d_n,bcx,i,SYY, time);
      break;
    case OutputInfo::StressZZ:
      getStressStrain(d_n,bcx,i,SZZ, time);
      break;
    case OutputInfo::StressXY:
      getStressStrain(d_n,bcx,i,SXY, time);
      break;
    case OutputInfo::StressYZ:
      getStressStrain(d_n,bcx,i,SYZ, time);
      break;
    case OutputInfo::StressXZ:
      getStressStrain(d_n,bcx,i,SXZ, time);
      break;
    case OutputInfo::StrainXX:
      getStressStrain(d_n,bcx,i,EXX, time);
      break;
    case OutputInfo::StrainYY:
      getStressStrain(d_n,bcx,i,EYY, time);
      break;
    case OutputInfo::StrainZZ:
      getStressStrain(d_n,bcx,i,EZZ, time);
      break;
    case OutputInfo::StrainXY:
      getStressStrain(d_n,bcx,i,EXY, time);
      break;
    case OutputInfo::StrainYZ:
      getStressStrain(d_n,bcx,i,EYZ, time);
      break;
    case OutputInfo::StrainXZ:
      getStressStrain(d_n,bcx,i,EXZ, time);
      break;
    case OutputInfo::StressVM:
      getStressStrain(d_n,bcx,i,VON, time);
      break;
    case OutputInfo::Damage:
      getStressStrain(d_n,bcx,i,DAMAGE, time);
      break;
    case OutputInfo::EquivalentPlasticStrain:
      getStressStrain(d_n, bcx, i, EQPLSTRN, time);
      break;
    case OutputInfo::StressPR1:
      getPrincipalStress(d_n,bcx,i,PSTRESS1,time);
      break;
    case OutputInfo::StressPR2:
      getPrincipalStress(d_n,bcx,i,PSTRESS2,time);
      break;
    case OutputInfo::StressPR3:
      getPrincipalStress(d_n,bcx,i,PSTRESS3,time);
      break;
    case OutputInfo::StrainPR1:
      getPrincipalStress(d_n,bcx,i,PSTRAIN1,time);
      break;
    case OutputInfo::StrainPR2:
      getPrincipalStress(d_n,bcx,i,PSTRAIN2,time);
      break;
    case OutputInfo::StrainPR3:
      getPrincipalStress(d_n,bcx,i,PSTRAIN3,time);
      break;
    case OutputInfo::InXForce:
      getElementForces(d_n, bcx, i, INX, time);
      break;
    case OutputInfo::InYForce:
      getElementForces(d_n, bcx, i, INY, time);
      break;
    case OutputInfo::InZForce:
      getElementForces(d_n, bcx, i, INZ, time);
      break;
    case OutputInfo::AXMoment:
      getElementForces(d_n, bcx, i, AXM, time);
      break;
    case OutputInfo::AYMoment:
      getElementForces(d_n, bcx, i, AYM, time);
      break;
    case OutputInfo::AZMoment:
      getElementForces(d_n, bcx, i, AZM, time);
      break;
    case OutputInfo::StrainVM:
      getStressStrain(d_n,bcx,i,STRAINVON, time);
      break;
    case OutputInfo::HeatFlXX:
      getHeatFlux(d_n, bcx, i, HFLX);
      break;
    case OutputInfo::HeatFlXY:
      getHeatFlux(d_n, bcx, i, HFLY);
      break;
    case OutputInfo::HeatFlXZ:
      getHeatFlux(d_n, bcx, i, HFLZ);
      break;
    case OutputInfo::GrdTempX:
      getHeatFlux(d_n, bcx, i, GRTX);
      break;
    case OutputInfo::GrdTempY:
      getHeatFlux(d_n, bcx, i, GRTY);
      break;
    case OutputInfo::GrdTempZ:
      getHeatFlux(d_n, bcx, i, GRTZ);
      break;
    case OutputInfo::HeatFlX:
      getTrussHeatFlux(d_n, bcx, i, HFLX);
      break;
    case OutputInfo::GrdTemp:
      getTrussHeatFlux(d_n, bcx, i, GRTX);
      break;
    default:
      success = 0;
      break;
  }
  return success;
}

//-------------------------------------------------------------------------------------

// Templated Post-processing for direct solver statics, frequency response, helmholtz and eigen
template<class Scalar>
void Domain::postProcessing(GenVector<Scalar> &sol, Scalar *bcx, GenVector<Scalar> &force,
                            int ndflag, int index, double time, double eigV,
                            GenSparseMatrix<Scalar> *kuc, GenSparseMatrix<Scalar> *kcc) {

  if(outFlag && !nodeTable) makeNodeTable(outFlag);
  int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes
  double freq;
  if (domain->probType() == SolverInfo::Modal) freq = eigV;
  else freq = domain->getFrequencyOrWavenumber();

  if (geoSource->isShifted() || domain->probType() == SolverInfo::Modal) time = freq;

  Scalar *globVal = 0;
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  if(numOutInfo && firstOutput && ndflag==0)
    filePrint(stderr," ... Postprocessing                 ...\n");

  // organize displacements
  int numNodeLim = myMax(numNodes,numnodes); 
    
  Scalar (*xyz)[11] = new Scalar[numNodeLim][11];//DofSet::max_known_nonL_dof
  int i;
  for(i = 0; i < numNodeLim; ++i)
    for (int j = 0 ; j < 11 ; j++)
      xyz[i][j] = 0.0;
  mergeDistributedDisp<Scalar>(xyz, sol.data(), bcx);
  int numNodesOut = (outFlag) ? exactNumNodes : numNodes;

  // Open files and write file headers in first time step
  if(firstOutput) geoSource->openOutputFiles();

  // Define Energy Terms
  // Total Energy = Wext+Wela+Wkin+Wdmp
  Wext = 0.0;  // external energy
  Waero = 0.0;  // aerodynamic force energy
  Wdmp = 0.0;  // damping energy
  double Wela=0.0;  // elastic energy
  double Wkin=0.0;  // kinetic energy

  int dof;
  int iNode;
  for(i = 0; i < numOutInfo; ++i)  {
    if(oinfo[i].ndtype != ndflag) continue;
    if(ndflag !=0 && oinfo[i].type != OutputInfo::Disp6DOF && oinfo[i].type !=  OutputInfo::Displacement) continue;
    // if non-deterministic and NOT displacement
    if(oinfo[i].interval == 1 
        || oinfo[i].type == OutputInfo::Farfield
        || oinfo[i].type == OutputInfo::Kirchhoff) {
      dof = -1;
      int success = processDispTypeOutputs(oinfo[i], xyz, numNodesOut, i, time, freq);
      if (success) continue;
      success = processOutput(oinfo[i].type, sol, bcx, i, time, freq);
      if (success) continue;
      success = 1;
      switch(oinfo[i].type)  {
        case OutputInfo::Energies: {
          Wext = ScalarTypes::Real(force*sol);   // Wext = external energy
          Wela =   0.5 * Wext;   // Wela = elastic energy
          double error = Wext+Wela+Wkin+Wdmp;
          geoSource->outputEnergies(i, freq, Wext, Waero, Wela, Wkin, Wdmp, error);
          }
          break;
        case OutputInfo::Farfield: case OutputInfo::Kirchhoff:
          outputFFP(sol, i);
          break;
        case OutputInfo::ModeError:
          break;  // This is handled in Problems.d/DynamDescr.C
        // The remaining cases are not officially supported in manual
        case OutputInfo::ElemToNode:
          if(elemToNode) elemToNode->print(oinfo[i].filptr, oinfo[i].nodeNumber);
          break;
        case OutputInfo::NodeToElem:
          if(nodeToElem) nodeToElem->print(oinfo[i].filptr, oinfo[i].nodeNumber);
          break;
        case OutputInfo::NodeToNode:
          if(nodeToNode) nodeToNode->print(oinfo[i].filptr, oinfo[i].nodeNumber);
          break;
        case OutputInfo::Reactions: {
          GenVector<Scalar> fc(numDirichlet+numComplexDirichlet);
          computeReactionForce(fc, sol, kuc, kcc);
          Scalar (*rxyz)[3] = new Scalar[numNodeLim][3];
          DofSet dofs[3] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp };
          for(int inode = 0, realNode = -1; inode < numnodes; ++inode) {
            if(nodeToElem && nodeToElem->num(inode) <= 0) continue;
            realNode++;
            int nodeI = (outFlag) ? realNode : inode;
            for(int k = 0; k < 3; ++k) {
              int dof =   dsa->locate(inode, dofs[k].list());
              int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
              rxyz[nodeI][k] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
            }
          }
          geoSource->outputNodeVectors(i, rxyz, numNodesOut, time);
          delete [] rxyz;
          } break;
        case OutputInfo::Reactions6: {
          GenVector<Scalar> fc(numDirichlet+numComplexDirichlet);
          computeReactionForce(fc, sol, kuc, kcc);
          Scalar (*rxyz)[6] = new Scalar[numNodeLim][6];
          DofSet dofs[6] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
                             DofSet::Xrot, DofSet::Yrot, DofSet::Zrot };
          for(int inode = 0, realNode = -1; inode < numnodes; ++inode) {
            if(nodeToElem && nodeToElem->num(inode) <= 0) continue;
            realNode++;
            int nodeI = (outFlag) ? realNode : inode;
            for(int k = 0; k < 6; ++k) {
              int dof =   dsa->locate(inode, dofs[k].list());
              int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
              rxyz[nodeI][k] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
            }
          }
          geoSource->outputNodeVectors6(i, rxyz, numNodesOut, time);
          delete [] rxyz;
          } break;
        case OutputInfo::HeatReactions: {
          GenVector<Scalar> fc(numDirichlet+numComplexDirichlet);
          computeReactionForce(fc, sol, kuc, kcc);
          Scalar *rxyz = new Scalar[numNodeLim];
          DofSet dofs[1] = { DofSet::Temp };
          for(int inode = 0, realNode = -1; inode < numnodes; ++inode) {
            if(nodeToElem && nodeToElem->num(inode) <= 0) continue;
            realNode++;
            int nodeI = (outFlag) ? realNode : inode;
            for(int k = 0; k < 1; ++k) {
              int dof =   dsa->locate(inode, dofs[k].list());
              int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
              rxyz[nodeI] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
            }
          }
          geoSource->outputNodeScalars(i, rxyz, numNodesOut, time);
          delete [] rxyz;
          } break;
        default:
          success = 0;
          break;
      }
      if (success == 0)
        fprintf(stderr, " *** AS.WRN: output %d is not supported \n", i);
    }
    if(globVal)
      { delete [] globVal; globVal = 0; }
  }

 // --- Print Problem statistics to the screen -------------------------------
 if(firstOutput) {
   if (!domain->solInfo().doEigSweep) {
   //  printStatistics();

     // ... CALCULATE STRUCTURE MASS IF REQUESTED
     if(sinfo.massFlag)  {
       double mass = computeStructureMass();
       filePrint(stderr," ... Total System mass = %10.4f ...\n",mass);
       filePrint(stderr," --------------------------------------\n");
     }
   }

   firstOutput = false;
 }

 if (xyz) delete [] xyz;
}

// XXXX this can and should be merged with buildRHSForce
template <class Scalar>
void
Domain::computeConstantForce(GenVector<Scalar>& cnst_f, GenSparseMatrix<Scalar>* kuc)
{
  // This is called for linear and nonlinear dynamics 
  // cnst_f is independent of t

  if(! dynamic_cast<GenSubDomain<Scalar>*> (this)) checkSommerTypeBC(this);

  cnst_f.zero();

  // ... COMPUTE FORCE FROM GRAVITY
  if(domain->gravityFlag()) addGravityForce(cnst_f);

  // ... COMPUTE FORCE FROM DISCRETE NEUMANN BOUNDARY CONDITIONS
  // note #1 when MFTT is present then FORCES contribution is not constant
  // note #2 when HFTT is present the FLUX contribution is not constant
  for(int i = 0; i < numNeuman; ++i) {
    int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
    if(dof < 0) continue;
    switch(nbc[i].type) {
      case(BCond::Forces) : if(!domain->mftval) cnst_f[dof] += nbc[i].val; break;
      case(BCond::Flux) :   if(!domain->hftval) cnst_f[dof] += nbc[i].val; break;
      case(BCond::Actuators) : case(BCond::Usdf) : break;  // these are never constant
      default : cnst_f[dof] += nbc[i].val;
    }
  }

  // ... COMPUTE FORCE FROM ACOUSTIC DISTRIBUTED NEUMANN BOUNDARY CONDITIONS
  // note #1: when MFTT is present this term is not constant (see computeExtForce4)
  if(sinfo.ATDDNBVal != 0.0 && !domain->mftval) addAtddnbForce(cnst_f);

  // ... COMPUTE FORCE FROM ACOUSTIC ROBIN BOUNDARY CONDITIONS
  //  note #1: when MFTT is present this term is not constant (see computeExtForce4)
  if(sinfo.ATDROBalpha != 0.0 && !domain->mftval) addAtdrobForce(cnst_f);

  // ... COMPUTE FORCE FROM PRESSURE
  // note #1: when MFTT is present this term is not constant (see computeExtForce4)
  // note #2: for NONLINEAR problems this term is not constant (see getStiffAndForce)
  if(!domain->mftval && !sinfo.isNonLin()) addPressureForce(cnst_f);

  // ... ADD RHS FROM LMPCs for linear statics
  if(lmpc.max_size() && !sinfo.isNonLin() && !sinfo.isDynam()) addMpcRhs(cnst_f);

  // ... COMPUTE FORCE FROM TEMPERATURES
  // note #1: for THERMOE problems TEMPERATURES are ignored 
  // note #2: for NONLINEAR problems this term is not constant (see getStiffAndForce)
  if(sinfo.thermalLoadFlag && !(sinfo.thermoeFlag >= 0) && !sinfo.isNonLin()) addThermalForce(cnst_f);

  // ... COMPUTE FORCE FROM NON-HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
  // note #1: when USDD is present this is term is not constant (see computeExtForce4)
  // note #2  for nonlinear this term is not constant (see getStiffAndForce) 
  if(numDirichlet && !(claw && claw->numUserDisp) && !sinfo.isNonLin() && kuc) {
    Vector Vc(numDirichlet, 0.0);
    // construct the non-homogeneous dirichlet bc vector
    for(int i = 0; i < numDirichlet; ++i) {
      int dof = dsa->locate(dbc[i].nnum, (1 << dbc[i].dofnum));
      if(dof >= 0) {
        int dof2 = c_dsa->invRCN(dof);
        if(dof2 >= 0) Vc[dof2] = dbc[i].val;
      }
    }
    // compute the non-homogeneous force
    kuc->multSubtract(Vc, cnst_f);
  }
}

template <class Scalar>
void
Domain::computeExtForce(GenVector<Scalar>& f, double t, GenSparseMatrix<Scalar>* kuc)
{
  f.zero();

  // ... COMPUTE FORCE FROM DISCRETE NEUMANN BOUNDARY CONDITIONS
  // note #1 when MFTT is not present FORCES contribution is constant (see computeConstantForce)
  // note #2 when HFTT is not present FLUX contribution is constant (see computeConstantForce)
  double mfttFactor = (domain->mftval) ? domain->mftval->getVal(t) : 1.0; // MFTT time dependent force coefficient
  double hfttFactor = (domain->hftval) ? domain->hftval->getVal(t) : 1.0; // HFTT time dependent flux coefficient
  if(numNeuman && (domain->mftval || domain->hftval || (claw && (claw->numUserForce || claw->numActuator)))) {
    for(int i = 0; i < numNeuman; ++i) {
      int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
      if(dof < 0) continue;
      switch(nbc[i].type) {
        case(BCond::Forces) : if(domain->mftval) f[dof] += mfttFactor*nbc[i].val; break;
        case(BCond::Flux)   : if(domain->hftval) f[dof] += hfttFactor*nbc[i].val; break;
        case(BCond::Actuators) : case(BCond::Usdf) : f[dof] += nbc[i].val; break;
        default : /* all other cases are constant */ ;
      }
    }
  }

  // COMPUTE FORCE FROM ACOUSTIC DISTRIBUTED NEUMANN BOUNDARY CONDITIONS
  // note #1: when MFTT not present this term is constant (see computeConstantForce)
  if(sinfo.ATDDNBVal != 0.0 && domain->mftval) addAtddnbForce(f, mfttFactor);

  // COMPUTE FORCE FROM ACOUSTIC ROBIN BOUNDARY CONDITIONS
  // note #1: when MFTT not present this term is constant (see computeConstantForce)
  if(sinfo.ATDROBalpha != 0.0 && domain->mftval) addAtdrobForce(f, mfttFactor);

  // COMPUTE FORCE FROM PRESSURE
  // note #1: when MFTT not present this term is constant (see computeConstantForce)
  // note #2: for NONLINEAR problems this term is follower (see getStiffAndForce)
  if(domain->mftval && !sinfo.isNonLin()) addPressureForce(f, mfttFactor);

  // ... ADD RHS FROM LMPCs for linear dynamics
  if(lmpc.max_size() && !sinfo.isNonLin() && sinfo.isDynam()) addMpcRhs(f);

  // COMPUTE FORCE FROM THERMOE
  // note #2: for NONLINEAR problems this term is follower (see getStiffAndForce)
  if(sinfo.thermoeFlag >= 0 && !sinfo.isNonLin()) addThermalForce(f);

  // COMPUTE FORCE FROM NON-HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
  // note #1: when USDD is not present this term is constant (see computeConstantForce)
  // note #2: for nonlinear this term is follower (see getStiffAndForce)
  if(numDirichlet && (claw && claw->numUserDisp) && !sinfo.isNonLin() && kuc) {
    Vector Vc(numDirichlet, 0.0);
    // construct the non-homogeneous dirichlet bc vector
    for(int i = 0; i < numDirichlet; ++i) {
      int dof = dsa->locate(dbc[i].nnum, (1 << dbc[i].dofnum));
      if(dof < 0) continue;
      int dof2 = c_dsa->invRCN(dof);
      if(dof2 >= 0) Vc[dof2] = dbc[i].val;
    }
// see Domain::updateUsddInDbc
//    if(userDefineDisp) {
//      int numDisp = claw->numUserDisp;
//      for(int i = 0; i < numDisp; ++i) {
//        int dof = dsa->locate(claw->userDisp[i].nnum, (1 << claw->userDisp[i].dofnum));
//        if(dof < 0) continue;
//        int dof2 = c_dsa->invRCN(dof);
//        if(dof2 >= 0)
//          if(userMap) {
//            Vc[dof2] += userDefineDisp[userMap[i]];
//          }
//          else
//            Vc[dof2] += userDefineDisp[i];
//      }
//    }

    // compute the non-homogeneous force
    kuc->multSubtract(Vc, f);
  }
}

template <class Scalar>
void
Domain::computeExtForce4(GenVector<Scalar>& f, const GenVector<Scalar>& constantForce,
                         double t, GenSparseMatrix<Scalar>* kuc)
{
  // This is called for linear and nonlinear dynamics 
  // doesn't include follower forces

  computeExtForce(f, t, kuc);

  // ADD CONSTANT FORCE
  f += constantForce;
}
