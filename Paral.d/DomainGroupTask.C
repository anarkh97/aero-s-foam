#include <iostream>
#include <Driver.d/DynamProbType.h>
#include <Paral.d/MDDynam.h>
#include <Threads.d/Paral.h>
#include <Driver.d/Dynam.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Paral.d/MDOp.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/Vector.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/SGISparseMatrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/Skyline.d/SGISky.h>
#include <Timers.d/GetTime.h>
#include <Control.d/ControlInterface.h>
#include <Threads.d/PHelper.h>
#include <Paral.d/GenMS.h>

#ifdef DISTRIBUTED
#include <Utils.d/DistHelper.h>
#endif

class IntFullM;

template<class Scalar>
GenDomainGroupTask<Scalar>::GenDomainGroupTask(int _nsub, GenSubDomain<Scalar> **_sd, double _cm, 
                                               double _cc, double _ck, Rbm **_rbms, FullSquareMatrix **_kelArray,
                                               double _alpha, double _beta, int _numSommer, int _solvertype,
                                               FSCommunicator *_com)
{
  nsub = _nsub;
  sd = _sd;
  dynMats = new GenSolver<Scalar> *[nsub];
  spMats  = new GenSparseMatrix<Scalar> *[nsub];
  M    = new GenSparseMatrix<Scalar> *[nsub];
  Muc  = new GenSparseMatrix<Scalar> *[nsub];
  // XML Mcc  = new GenSparseMatrix<Scalar> *[nsub];
  C    = new GenSparseMatrix<Scalar> *[nsub];
  Cuc  = new GenSparseMatrix<Scalar> *[nsub];
// RT
  C_deriv    = new GenSparseMatrix<Scalar> **[nsub];
  Cuc_deriv    = new GenSparseMatrix<Scalar> **[nsub];
// RT end
  K    = new GenSparseMatrix<Scalar> *[nsub];
  //rbms = new Rbm *[nsub];
  rbms = _rbms;
  kelArray = _kelArray;
  Kuc  = new GenSparseMatrix<Scalar> *[nsub];
  coeM    = _cm;
  coeC    = _cc;
  coeK    = _ck;
  numSommer = _numSommer;
  alpha   = _alpha; 
  beta    = _beta;
  solvertype = _solvertype;
  com = _com;
}

template<class Scalar>
GenDomainGroupTask<Scalar>::~GenDomainGroupTask()
{
  // delete [] dynMats;
  delete [] spMats;
  //delete [] rbms;
  // don't delete K,Kuc,C,Cuc,M,Muc
}

template<class Scalar> 
void
GenDomainGroupTask<Scalar>::runFor(int isub, bool make_feti) 
{
 DofSetArray     *dsa = sd[isub]->getDSA();
 ConstrainedDSA *cdsa = sd[isub]->getCDSA();

 // access to implicit/explicit - JFD for Acoustic Time Domain
 double betaNewmark=sd[isub]->solInfo().newmarkBeta;
 Connectivity *connForMass = 0;
 if (solvertype!=10) {//note: solvertype=10 is diagonal matrices for explicit
   if (betaNewmark==0.0) {
     int numN=sd[isub]->numNodes();
     // import a diagonal connectivity for the mass matrix
     connForMass = new Connectivity(numN);
     //connForMass->print();
     connForMass = connForMass->modify();
   }
 }

 // XML Need to introduce Mcc
 K[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();

 if(domain->solInfo().inpc) M[isub] = 0; 
 else {  
   if (betaNewmark == 0.0)//explict simulation
     if (solvertype != 10)
       M[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>(sd[isub]->getCDSA(),connForMass);
     else
       // if only spectral elements are used, the mass matrix is purely diagonal for implicit and explicit
       M[isub] = new GenDiagMatrix<Scalar>(sd[isub]->getCDSA());
   else {
     if (domain->solInfo().isCoupled)
       M[isub] = sd[isub]->template constructNBSparseMatrix<Scalar>();
     else
       M[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();
   }
 }

 if (solvertype != 10)
   delete connForMass;
 
 if((cdsa->size() - dsa->size()) == 0) {  // PJSA
   Kuc[isub] = 0;
   Muc[isub] = 0;
 } 
 else {
   Kuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
   if(!domain->solInfo().inpc) Muc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
   else Muc[isub] = 0;
 }

 if(domain->solInfo().inpc) { C[isub] = 0; Cuc[isub] = 0; }
 else 
 {
   if(alpha != 0.0 || beta != 0.0 || (numSommer > 0)) { // for rayleigh damping
     C[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();

     if((cdsa->size() - dsa->size()) == 0) // PJSA
       Cuc[isub] = 0;
     else 
       Cuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
     int numC_deriv, numRHS;
     numC_deriv = 1; numRHS = 1+1;
     C_deriv[isub] = new GenSparseMatrix<Scalar> * [numRHS - 1];
     for(int n=0; n<numC_deriv; ++n) {
       C_deriv[isub][n] = sd[isub]->template constructDBSparseMatrix<Scalar>();
     }
     for(int n=numC_deriv; n<numRHS - 1; ++n)
       C_deriv[isub][n] = 0;
     if((cdsa->size() - dsa->size()) != 0) {
       Cuc_deriv[isub] = new GenSparseMatrix<Scalar> * [numRHS - 1];
       for(int n=0; n<numC_deriv; ++n) {
         Cuc_deriv[isub][n] = sd[isub]->template constructCuCSparse<Scalar>();
       }
       for(int n=numC_deriv; n<numRHS - 1; ++n)
         Cuc_deriv[isub][n] = 0;
     } else {
       Cuc_deriv[isub] = 0;
     }
   }
   else if (domain->solInfo().ATDARBFlag != -2.0) {// for acoustic damping
     if (solvertype != 10) {
       //build the modified damping matrix for implicit and explicit
       C[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>(sd[isub]->getCDSA(),sd[isub]->getNodeToNode_sommer());
     }
     else
       C[isub] = new GenDiagMatrix<Scalar>(sd[isub]->getCDSA());
       if((cdsa->size() - dsa->size()) == 0) // PJSA
         Cuc[isub] = 0;
       else
         Cuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
   }
   else {
     C[isub] = 0;
     Cuc[isub] = 0;
     C_deriv[isub] = 0;
     Cuc_deriv[isub] = 0;
   }
 }

 // builds the datastructures for Kii, Kib, Kbb
 if(domain->solInfo().type == 2 && make_feti) { // FETI
   if(sd[isub]->numMPCs() > 0)
     sd[isub]->makeKbbMpc();  
   else {
     sd[isub]->makeKbb(sd[isub]->getCCDSA());
   }
 }

 GenMultiSparse<Scalar> *allMats;
 if(make_feti) {
   if(domain->solInfo().type != 0) {
     switch(solvertype) {
       default :
         cerr << " ... WARNING: in DomainGroupTask::runFor solvertype " << solvertype << " is  not implemented. Skyline solver is used instead ...\n";
       case 0 : {
         GenSkyMatrix<Scalar> *skmat = sd[isub]->template constructSkyMatrix<Scalar>(sd[isub]->getCCDSA(), 0);
         dynMats[isub] = skmat;
         spMats[isub] = skmat;
       } break;
       case 1 : {
         GenBLKSparseMatrix<Scalar> *bsmat = sd[isub]->template constructBLKSparseMatrix<Scalar>(sd[isub]->getCCDSA(), 0);
         bsmat->zeroAll();
         dynMats[isub] = bsmat;
         spMats[isub]  = bsmat;
       } break;
       case 2 : {
         GenSGISparseMatrix<Scalar> *sgimat = sd[isub]->template constructSGISparseMatrix<Scalar>(isub, 0);
         dynMats[isub] = sgimat;
         spMats[isub]  = sgimat;
       } break;
       case 3 : {
         SGISky *sgisky = sd[isub]->constructSGISkyMatrix(0);
         dynMats[isub] = dynamic_cast<GenSolver<Scalar> *>(sgisky);
         spMats[isub] = dynamic_cast<GenSparseMatrix<Scalar> *>(sgisky);
       } break;
#ifdef USE_SPOOLES
       case 8 : {
         GenSpoolesSolver<Scalar> *ssmat = sd[isub]->template constructSpooles<Scalar>(sd[isub]->getCCDSA());
         dynMats[isub] = ssmat;
         spMats[isub] = ssmat;
       } break;
#endif
#ifdef USE_MUMPS
       case 9 : {
         GenMumpsSolver<Scalar> *msmat = sd[isub]->template constructMumps<Scalar>(sd[isub]->getCCDSA());
         dynMats[isub] = msmat;
         spMats[isub] = msmat;
       } break;
#endif
       case 10 : {
         GenDiagMatrix<Scalar> *spm = new GenDiagMatrix<Scalar>(sd[isub]->getCCDSA()); 
         dynMats[isub] = spm;
         spMats[isub] = spm;
       } break;
     }
   }

   if(domain->solInfo().type == 2 && domain->solInfo().getFetiInfo().version == FetiInfo::fetidp) {
     sd[isub]->constructKcc();
     sd[isub]->constructKrc();
   }

   if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted)
     allMats = new GenMultiSparse<Scalar>(spMats[isub], sd[isub]->Krc, sd[isub]->Kcc);
   else 
     allMats = new GenMultiSparse<Scalar>(spMats[isub], sd[isub]->KiiSparse, sd[isub]->Kbb,
                                          sd[isub]->Kib, sd[isub]->Krc, sd[isub]->Kcc);
 } 
 else {
   allMats = 0;
 }

 AllOps<Scalar> allOps;

 if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted)
   allOps.K = new GenMultiSparse<Scalar>(K[isub], sd[isub]->KiiSparse, sd[isub]->Kbb, sd[isub]->Kib);
 else
   allOps.K = K[isub];
 allOps.C = C[isub]; 
 allOps.Cuc = Cuc[isub];
 allOps.M = M[isub];
 allOps.Muc = Muc[isub];
 allOps.Kuc = Kuc[isub];
 allOps.C_deriv = C_deriv[isub];
 allOps.Cuc_deriv = Cuc_deriv[isub];
 FullSquareMatrix *subKelArray = (kelArray) ? kelArray[isub] : 0;
 sd[isub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, allMats, subKelArray);

/* XXXX should be done in problem descriptor, eg. Paral.d/MDStatic.C
  if(make_feti && domain->solInfo().type == 2 && sd[isub]->solInfo().getFetiInfo().version != FetiInfo::fetidp)  // PJSA DEBUG
    rbms[isub] = sd[isub]->constructRbm(false);
*/

 if(allMats) delete allMats;
}
