#include <iostream>
#include <Driver.d/SysState.h>
#include <Paral.d/MDDynam.h>
#include <Threads.d/Paral.h>
#include <Driver.d/Dynam.h>
#include <Paral.d/MDOp.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/Vector.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Timers.d/GetTime.h>
#include <Control.d/ControlInterface.h>
#include <Threads.d/PHelper.h>
#include <Paral.d/GenMS.h>
#include <Solvers.d/SolverFactory.h>
#ifdef DISTRIBUTED
#include <Utils.d/DistHelper.h>
#endif

class IntFullM;

template<class Scalar>
GenDomainGroupTask<Scalar>::GenDomainGroupTask(int _nsub, GenSubDomain<Scalar> **_sd, double _cm, 
                                               double _cc, double _ck, Rbm **_rbms, FullSquareMatrix **_kelArray,
                                               double _alpha, double _beta, int _numSommer, int _solvertype,
                                               FSCommunicator *_com, FullSquareMatrix **_melArray, FullSquareMatrix **_celArray)
{
  nsub = _nsub;
  sd = _sd;
  dynMats = new GenSolver<Scalar> *[nsub];
  spMats  = new GenSparseMatrix<Scalar> *[nsub];
  M    = new GenSparseMatrix<Scalar> *[nsub];
  Muc  = new GenSparseMatrix<Scalar> *[nsub];
  Mcc  = new GenSparseMatrix<Scalar> *[nsub];
  C    = new GenSparseMatrix<Scalar> *[nsub];
  Cuc  = new GenSparseMatrix<Scalar> *[nsub];
  Ccc  = new GenSparseMatrix<Scalar> *[nsub];
  C_deriv    = new GenSparseMatrix<Scalar> **[nsub];
  Cuc_deriv    = new GenSparseMatrix<Scalar> **[nsub];
  K    = new GenSparseMatrix<Scalar> *[nsub];
  if(domain->solInfo().solvercntl->precond) {
    spp = new GenSparseMatrix<Scalar> *[nsub];
    sps = new GenSolver<Scalar> *[nsub];
  }
  else { 
    spp = 0;
    sps = 0;
  }
  rbms = _rbms;
  kelArray = _kelArray;
  melArray = _melArray;
  celArray = _celArray;
  Kuc  = new GenSparseMatrix<Scalar> *[nsub];
  coeM    = _cm;
  coeC    = _cc;
  coeK    = _ck;
  numSommer = _numSommer;
  alpha   = _alpha; 
  beta    = _beta;
  solvertype = _solvertype;
  com = _com;
  makeC = (alpha != 0.0 || beta != 0.0 || (numSommer > 0) || domain->getElementSet().hasDamping());
// RT - 053013 - to enable multiple impedance section, build C_deriv whenever C
//  makeC_deriv = (makeC && domain->solInfo().doFreqSweep && domain->solInfo().getSweepParams()->nFreqSweepRHS > 1);
  makeC_deriv = (makeC && domain->solInfo().doFreqSweep );
}

template<class Scalar>
GenDomainGroupTask<Scalar>::~GenDomainGroupTask()
{
  // delete [] dynMats;
  //delete [] spMats;
  //delete [] rbms;
  // don't delete K,Kuc,C,Cuc,M,Muc
}

template<class Scalar> 
void
GenDomainGroupTask<Scalar>::runFor(int isub, bool make_feti) 
{
  DofSetArray     *dsa = sd[isub]->getDSA();
  ConstrainedDSA *cdsa = sd[isub]->getCDSA();

  K[isub] = 0;
  Kuc[isub] = 0;
  M[isub] = 0;
  Muc[isub] = 0;
  Mcc[isub] = 0;
  C[isub] = 0;
  Cuc[isub] = 0;
  Ccc[isub] = 0;
  C_deriv[isub] = 0;
  Cuc_deriv[isub] = 0;
  if(spp) spp[isub] = 0;
  if(sps) sps[isub] = 0;

  if((cdsa->size() - dsa->size()) != 0)
    Kuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();

  if(domain->solInfo().isDynam() || domain->solInfo().doFreqSweep || domain->solInfo().probType == SolverInfo::Modal || domain->solInfo().probType == SolverInfo::PodRomOffline) {

    // XML Need to introduce Mcc
    if(domain->solInfo().isCoupled)
      K[isub] = sd[isub]->template constructNBSparseMatrix<Scalar>(); // unsymmetric
    else
      K[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();

    if(domain->solInfo().newmarkBeta == 0.0) { // explict dynamics
      if(solvertype != 10) {
        int numN = sd[isub]->numNodes();
        // import a diagonal connectivity for the mass matrix
        Connectivity *connForMass = new Connectivity(numN);
        connForMass = connForMass->modify();
        M[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>(cdsa, connForMass);
        delete connForMass;
      }
      else {
        // if only spectral elements are used, the mass matrix is purely diagonal for implicit and explicit
        M[isub] = new GenDiagMatrix<Scalar>(cdsa);
      }
    }
    else {
      if(domain->solInfo().isCoupled)
        M[isub] = sd[isub]->template constructNBSparseMatrix<Scalar>();
      else
        M[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();
    }

    if((cdsa->size() - dsa->size()) != 0) {
      Muc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
      Mcc[isub] = sd[isub]->template constructCCSparse<Scalar>();
    }

    if(makeC) {
      C[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();

      if((cdsa->size() - dsa->size()) != 0) {
        Cuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
        Ccc[isub] = sd[isub]->template constructCCSparse<Scalar>();
      }

      if(makeC_deriv) {
        int numC_deriv, numRHS;
        numRHS = domain->solInfo().getSweepParams()->nFreqSweepRHS;
        if((numSommer > 0) && ((sd[isub]->sommerfeldType == 2) || (sd[isub]->sommerfeldType == 4)))
          numC_deriv = numRHS - 1;
        else
          numC_deriv = 1;
        C_deriv[isub] = new GenSparseMatrix<Scalar> * [numRHS - 1];
        for(int n = 0; n < numC_deriv; ++n) {
          C_deriv[isub][n] = sd[isub]->template constructDBSparseMatrix<Scalar>();
        }
        for(int n = numC_deriv; n < numRHS - 1; ++n)
          C_deriv[isub][n] = 0;
        if(cdsa->size() > 0 && (cdsa->size() - dsa->size()) != 0) {
          Cuc_deriv[isub] = new GenSparseMatrix<Scalar> * [numRHS - 1];
          for(int n = 0; n < numC_deriv; ++n) {
            Cuc_deriv[isub][n] = sd[isub]->template constructCuCSparse<Scalar>();
          }
          for(int n = numC_deriv; n < numRHS - 1; ++n)
            Cuc_deriv[isub][n] = 0;
        }
      }
    }
    else if(domain->solInfo().ATDARBFlag != -2.0) { // for acoustic damping
      if(solvertype != 10) {
       // build the modified damping matrix for implicit and explicit
       C[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>(cdsa, sd[isub]->getNodeToNode_sommer());
     }
     else
       C[isub] = new GenDiagMatrix<Scalar>(cdsa);
       if((cdsa->size() - dsa->size()) != 0)
         Cuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
    }
  }

  // builds the datastructures for Kii, Kib, Kbb
  if(domain->solInfo().solvercntl->type == 2 && make_feti) { // FETI
    if(sd[isub]->numMPCs() > 0)
      sd[isub]->makeKbbMpc();
    else {
      sd[isub]->makeKbb(sd[isub]->getCCDSA());
    }
  }

  GenMultiSparse<Scalar> *allMats = 0;
  if(make_feti) {
    if(domain->solInfo().solvercntl->type == 2 || domain->solInfo().solvercntl->type == 3) {
      // construct local solver for the subdomain
      dynMats[isub] = GenSolverFactory<Scalar>::getFactory()->createSolver(sd[isub]->getNodeToNode(), sd[isub]->getDSA(), sd[isub]->getCCDSA(),
                                                                           *domain->solInfo().solvercntl->fetiInfo.local_cntl, spMats[isub],
                                                                           (Rbm*) NULL, spp[isub], sps[isub]);
    }
    else if(domain->solInfo().solvercntl->type == 1) {
      dynMats[isub] = 0;
      switch(domain->solInfo().solvercntl->iterSubtype) {
        case 2 :
          spMats[isub] = sd[isub]->template constructNBSparseMatrix<Scalar>();
          break;
        default:
        case 3 :
          spMats[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();
          break;
#ifdef USE_EIGEN3
        case 4:
          spMats[isub] = sd[isub]->template constructEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >();
          break;
#endif
      }
      if(domain->solInfo().solvercntl->precond == 1) {
        GenDiagMatrix<Scalar> *dm = new GenDiagMatrix<Scalar>(sd[isub]->getCDSA());
        spp[isub] = (GenSparseMatrix<Scalar>*) dm;
        sps[isub] = (GenSolver<Scalar>*) dm;
      }
    }

    if(domain->solInfo().solvercntl->type == 2 && domain->solInfo().getFetiInfo().version == FetiInfo::fetidp) {
      sd[isub]->constructKcc();
      sd[isub]->constructKrc();
    }

    if(domain->solInfo().solvercntl->type == 2) {
      if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted)
        allMats = new GenMultiSparse<Scalar>(spMats[isub], sd[isub]->Krc, sd[isub]->Kcc);
      else 
        allMats = new GenMultiSparse<Scalar>(spMats[isub], sd[isub]->KiiSparse, sd[isub]->Kbb,
                                             sd[isub]->Kib, sd[isub]->Krc, sd[isub]->Kcc);
    }
  }

  AllOps<Scalar> allOps;

  if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted)
    allOps.K = new GenMultiSparse<Scalar>(K[isub], sd[isub]->KiiSparse, sd[isub]->Kbb, sd[isub]->Kib);
  else
    allOps.K = K[isub];
  allOps.C = C[isub]; 
  allOps.Cuc = Cuc[isub];
  allOps.Ccc = Ccc[isub];
  allOps.M = M[isub];
  allOps.Muc = Muc[isub];
  allOps.Mcc = Mcc[isub];
  allOps.Kuc = Kuc[isub];
  allOps.C_deriv = C_deriv[isub];
  allOps.Cuc_deriv = Cuc_deriv[isub];
  allOps.spp = (spp) ? spp[isub] : 0;
  FullSquareMatrix *subKelArray = (kelArray) ? kelArray[isub] : 0;
  FullSquareMatrix *subMelArray = (melArray) ? melArray[isub] : 0;
  FullSquareMatrix *subCelArray = (celArray) ? celArray[isub] : 0;
  if(domain->solInfo().solvercntl->type == 2)
    sd[isub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, allMats, subKelArray, subMelArray, subCelArray);
  else
    sd[isub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, spMats[isub], subKelArray, subMelArray, subCelArray);

  if(allMats) delete allMats;
}
