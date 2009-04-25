#include <stdio.h>
#include <Driver.d/DecDomain.h>
#include <Feti.d/Feti.h>
#include <Paral.d/MDDynam.h>
#include <Utils.d/BlockAlloc.h>
#include <Paral.d/DomainGroupTask.h>
#include <Paral.d/Assembler.h>


template<>
double
GenDecDomain<double>::computeStabilityTimeStep(GenMDDynamMat<double>& dMat)
{
      double eigmax;
      double relTol    = 1.0e-3;
      double preeigmax = 0.0;
      int numdofs = internalInfo.len;
      int maxIte  = numdofs * 50;
      GenDistrVector<double> v(internalInfo);
      GenDistrVector<double> z(internalInfo);
// Starts from an arbitrary array.
      int i;
      for (i=0; i<numdofs; ++i)
        v.data()[i] = (double) (i+1) / (double) numdofs;

 FSCommPattern<double> *pat = new FSCommPattern<double>(communicator, cpuToSub, myCPU, FSCommPattern<double>::CopyOnSend);
 for(int i=0; i<numSub; ++i) subDomain[i]->setDofPlusCommSize(pat);
 pat->finalize();
 BasicAssembler *ba = new BasicAssembler(numSub, subDomain, pat);
 dMat.K->setAssembler(ba); 
 dMat.M->setAssembler(ba);
// Power iteration loop
      for (i=0; i<maxIte; ++i) {

        dMat.K->mult(v,z);
        dMat.M->multInvertDiag(z);



        //for (j=0; j< numdofs; ++j)
        // z[j] /= dMat.M->diag(j);
// Normalize
/*
        double zmax = z[0];
        for (j=1; j< numdofs; ++j)
          if (abs(z.data()[j])>zmax) zmax = abs(z[j]);
#ifdef DISTRIBUTED
  zmax = communicator->globalMax(zmax);
#endif
*/
        double zmax = z.infNorm();
        eigmax = zmax;
        //v = (1.0/zmax)*z;
        v.linC(z,1.0/zmax);
        if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;
        preeigmax = eigmax;
      }
      // compute stability maximum time step
      double sdt = 2.0 / sqrt(eigmax);
  delete ba;
  delete pat;
  dMat.K->setAssembler((BasicAssembler *)0);
  dMat.M->setAssembler((BasicAssembler *)0);
  return sdt;
}

/* TEMPLATED AND MOVED TO Driver.d/DecDomain.C
template<>
void
GenDecDomain<DComplex>::buildOps(MDDynamMat& res, double coeM, double coeC, double coeK, Rbm **rbm, FullSquareMatrix **kelArray, bool make_feti)
{
  fprintf(stderr, "WARNING: GenDecDomain<DComplex>::buildOps(MDDynamMat& res, ...) not implemented \n");
}

template<>
void
GenDecDomain<double>::buildOps(MDDynamMat &res, double coeM, double coeC, double coeK, Rbm **rbms, FullSquareMatrix **kelArray, bool make_feti)
{
 BasicAssembler *ba = 0; 

 bool isFeti = domain->solInfo().type == 2;
 FetiInfo *finfo = &domain->solInfo().getFetiInfo();

 int solvertype = finfo->solvertype;

 int isFeti2 =  (isFeti && finfo->version == FetiInfo::feti2) ? 1 : 0;
 int isCtcOrDualMpc = (numDualMpc) ? 1 : 0;

 GenDomainGroupTask<double> dgt(numSub, subDomain, coeM, coeC, coeK, rbms, kelArray,
                                domain->solInfo().alphaDamp, domain->solInfo().betaDamp, isFeti2,
                                solvertype, isCtcOrDualMpc);

 filePrint(stderr," ... Assemble Subdomain Matrices    ... \n");
 if(isFeti && (finfo->version == FetiInfo::fetidp) && (finfo->augment == FetiInfo::Gs)) {
   // this is for sending and receiving the number of coarse grid modes
   FSCommPattern<int> *sPat = new FSCommPattern<int>(communicator, cpuToSub, 0, FSCommPattern<int>::CopyOnSend);
   for(int i=0; i<numSub; ++i) subDomain[i]->setCommSize(sPat,1);
   sPat->finalize();
   execParal(numSub, &dgt, &GenDomainGroupTask<double>::runFor1, make_feti, sPat);
   sPat->exchange(); 
   execParal(numSub, &dgt, &GenDomainGroupTask<double>::runFor2, make_feti, sPat);
   delete sPat;
 }
 else {
   execParal(numSub, &dgt, &GenDomainGroupTask<double>::runFor, make_feti);
 }

 if(domain->solInfo().inpc) { 
   FSCommPattern<double> *pat = new FSCommPattern<double>(communicator, cpuToSub, myCPU, FSCommPattern<double>::CopyOnSend);
   for(int i=0; i<numSub; ++i) subDomain[i]->setDofPlusCommSize(pat);
   pat->finalize();
   ba = new BasicAssembler(numSub, subDomain, pat);
 }
 res.K   = new SubDOp(numSub, dgt.K, ba);
 res.Kuc = new SubDOp(numSub, dgt.Kuc);

 if(dgt.C[0]) {
   res.C = new SubDOp(numSub, dgt.C);
   res.Cuc = new SubDOp(numSub, dgt.Cuc);
 }
 else {
   res.C   = 0; delete [] dgt.C;
   res.Cuc = 0; delete [] dgt.Cuc;
 }
 res.M   = new SubDOp(numSub, dgt.M);
 res.Muc = new SubDOp(numSub, dgt.Muc);

 if(isFeti) {
   if(make_feti) res.dynMat = getDynamicFetiSolver(dgt);
 } else
   res.dynMat = getDiagSolver(numSub, dgt.sd, dgt.dynMats);
}

template<>
void
GenDecDomain<DComplex>::rebuildOps(MDDynamMat &res, double coeM, double coeC, double coeK)
{
  cerr << "GenDecDomain<DComplex>::rebuildOps is not implemented\n";
}

template<>
void
GenDecDomain<double>::rebuildOps(MDDynamMat &res, double coeM, double coeC, double coeK)
{
 res.dynMat->reconstruct(); // do anything that needs to be done before zeroing and assembling the matrices

 execParal4R(numSub, this, &GenDecDomain<double>::subRebuildOps, res, coeM, coeC, coeK);

 res.dynMat->refactor(); // do anything that needs to be done after zeroing and assembling the matrices
}

template<>
void
GenDecDomain<double>::subRebuildOps(int iSub, MDDynamMat &res, double coeM, double coeC, double coeK)
{
  AllOps<double> allOps;

  if(res.K)  allOps.K = (*res.K)[iSub];
  if(res.C)  allOps.C = (*res.C)[iSub];
  if(res.Cuc)  allOps.Cuc = (*res.Cuc)[iSub];
  if(res.M)  allOps.M = (*res.M)[iSub]; 
  if(res.Muc)  allOps.Muc = (*res.Muc)[iSub];
  if(res.Mcc)  allOps.Mcc = (*res.Mcc)[iSub];
  if(res.Kuc)  allOps.Kuc = (*res.Kuc)[iSub]; 
 
  allOps.zero();

  GenMultiSparse<double> allMats(subDomain[iSub]->KrrSparse, subDomain[iSub]->KiiSparse, subDomain[iSub]->Kbb,
                                 subDomain[iSub]->Kib, subDomain[iSub]->Krc, subDomain[iSub]->Kcc);
  allMats.zeroAll();

  subDomain[iSub]->makeSparseOps<double>(allOps, coeK, coeM, coeC, &allMats);
}
*/

template<>
void
GenDecDomain<double>::buildLocalFFP(int iSub, GenDistrVector<double> *u,
                                    double **ffp, int *numSample, double (*dir)[3])
{
  fprintf(stderr, "WARNING: GenDecDomain<double>::buildLocalFFP not implemented \n");
}

template<>
void
GenDecDomain<DComplex>::buildLocalFFP(int iSub, GenDistrVector<DComplex> *u,
                                      DComplex **ffp, int *numSample, double (*dir)[3])
{
 subDomain[iSub]->ffp(subDomain[iSub], *numSample, ffp[iSub], dir, u->subData(iSub));
}

template<>
void
GenDecDomain<double>::buildFFP(GenDistrVector<double> &u, FILE *fffp)
{
 fprintf(stderr, "WARNING: GenDecDomain<double>::buildFFP not implemented \n"); 
}

template<>
void
GenDecDomain<DComplex>::buildFFP(GenDistrVector<DComplex> &u, FILE *fffp)
{
 if(domain->numFFPDirections == 0) {
   int i,j;
   int nsint = MAX(2, domain->nffp);
   // Dimension of the problem
   int dim= domain->scatter[0]->dim();
   DComplex ffpCoef;
   if(dim==3) ffpCoef = exp(DComplex(0.0,M_PI/4.0))/sqrt(8.0*M_PI*geoSource->kappa());
   else ffpCoef = DComplex(0.25/M_PI, 0.0);

   double (*ffpDir)[3];
   int numSamples;
   if(dim==2) {
     ffpDir = new double[nsint][3];
     for(i=0;i<nsint;i++) {
       ffpDir[i][0] = cos(2*M_PI*double(i)/double(nsint));
       ffpDir[i][1] = sin(2*M_PI*double(i)/double(nsint));
       ffpDir[i][2] = 0.0;
     }
     numSamples = nsint;
   } 
   else {
     int numTheta = nsint/2+1;
     ffpDir = new double[nsint*numTheta][3];
     for(i=0; i<numTheta; ++i) {
       double theta = M_PI*(-0.5+((double) i)/(numTheta-1.0));
       for(j=0;j<nsint;j++) {
         ffpDir[i*nsint+j][0] = cos(theta)*cos(2*j*M_PI/double(nsint));
         ffpDir[i*nsint+j][1] = cos(theta)*sin(2*j*M_PI/double(nsint));
         ffpDir[i*nsint+j][2] = sin(theta);
       }
     }
     numSamples = nsint*numTheta;
   }

   DComplex **localFFP = new DComplex * [numSub];
   for(i=0; i<numSub; ++i) {
     localFFP[i] = new DComplex[numSamples];
     for(j=0; j<numSamples; ++j) {
       localFFP[i][j] = DComplex(0.0,0.0);
     }
   }

   execParal(numSub, this, &GenDecDomain<DComplex>::buildLocalFFP, &u, localFFP, &numSamples, ffpDir);

   for(i=1; i<numSub; ++i) 
     for(j=0; j<numSamples; ++j)
       localFFP[0][j] += localFFP[i][j];
   for(j=0; j<numSamples; ++j)
     localFFP[0][j] *= ffpCoef;
#ifdef DISTRIBUTED
   communicator->globalSum(numSamples, localFFP[0]);
   if(communicator->cpuNum() == 0) {
#endif
   if(dim==2) {
     for(j=0; j<nsint; ++j) {
       fprintf(fffp,"%e  %e  %.10e  %.10e\n",
               2*j*M_PI/double(nsint),0.0,
               ScalarTypes::Real(localFFP[0][j]), ScalarTypes::Imag(localFFP[0][j]));
     }
   }
   else {
     int numTheta = nsint/2+1;
     for(i=0; i<numTheta; ++i) {
       for(j=0; j<nsint; ++j) {
//         double y = 10.0 * log(2.0*M_PI*abs(localFFP[0][numSample])*
//                    abs(localFFP[0][numSample]))/log(10.0);
         fprintf(fffp,"%e  %e  %.10e  %.10e\n",
                 2*j*M_PI/double(nsint),M_PI*(-0.5+((double)i)/(numTheta-1.0)),
                 ScalarTypes::Real(localFFP[0][i*nsint+j]),
                 ScalarTypes::Imag(localFFP[0][i*nsint+j]));
       }
     }
   }
#ifdef DISTRIBUTED
   }
#endif
   for(i=0; i<numSub; ++i) delete [] localFFP[i];
   delete [] localFFP;
 } 
 else {
   // RT: new style input/output
   int i,j;
   int numSamples = domain->numFFPDirections;
   DComplex **localFFP = new DComplex * [numSub];
   for(i=0; i<numSub; ++i) {
     localFFP[i] = new DComplex[numSamples];
     for(j=0; j<numSamples; ++j) {
       localFFP[i][j] = DComplex(0.0,0.0);
     }
   }

   double (*ffpDir)[3];
   ffpDir = new double[numSamples][3];
   for(i=0;i<numSamples;i++) {
     ffpDir[i][0] = domain->ffpDirections[i*3+0];
     ffpDir[i][1] = domain->ffpDirections[i*3+1];
     ffpDir[i][2] = domain->ffpDirections[i*3+2];
   }

   execParal(numSub, this, &GenDecDomain<DComplex>::buildLocalFFP, &u, localFFP,
           &numSamples, ffpDir);

   for(i=1; i<numSub; ++i)
     for(j=0; j<numSamples; ++j)
       localFFP[0][j] += localFFP[i][j];

#ifdef DISTRIBUTED
   communicator->globalSum(numSamples, localFFP[0]);
   if(communicator->cpuNum() == 0) {
#endif

   for(i=0;i<domain->numFFPDirections;i++) {
     fprintf(fffp,"%e %e %e   %e %e\n",
             domain->ffpDirections[i*3+0],
             domain->ffpDirections[i*3+1],
             domain->ffpDirections[i*3+2],
             real(localFFP[0][i]),imag(localFFP[0][i]));
   }
#ifdef DISTRIBUTED
   }
#endif

   for(i=0; i<numSub; ++i) delete [] localFFP[i];
   delete [] localFFP;
   delete [] ffpDir;
 }
}

