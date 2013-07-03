#include <Problems.d/TempDescr.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Driver.d/Domain.h>
#include <Hetero.d/FlExchange.h>
#include <Driver.d/Dynam.h>

#include <Math.d/FullMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/State.h>
#include <Solvers.d/Rbm.h>

typedef FSFullMatrix FullMatrix;

void
SingleDomainTemp::tempprojector_prep(Rbm *rbms, SparseMatrix *M)
{
 if (!numR) return;

 int ndof = M->dim();
 
 fprintf(stderr," ... Building the HZEM Projector    ...\n");

 //fprintf(stderr," ... Number of HZEM     =   %d       ...\n",numR);

 Rmem = new double[ndof];
 for(int n=0; n<ndof; ++n) Rmem[n] = 1.;
 
 StackFSFullMatrix Rt(numR, ndof, Rmem);

 double *MRmem = new double[ndof];
 StackFSFullMatrix MR(numR, ndof, MRmem);

 M->mult(Rmem, MRmem);

 FullMatrix MRt = MR.transpose();

 FullMatrix RtMR(numR,numR);
 Rt.mult(MRt,RtMR);

 FullMatrix RtMRinverse = RtMR.invert();

 X = new FullMatrix(ndof,numR);
 MRt.mult(RtMRinverse,(*X));

}

void
SingleDomainTemp::temptrProject(Vector &f)
{
 if (!numR) return;

 int ndof = f.size();

 double *yMem = (double *) dbg_alloca(numR*sizeof(double));
 double *zMem = (double *) dbg_alloca(ndof*sizeof(double));
 StackVector y(ndof,yMem);
 StackVector z(ndof,zMem);

 StackFSFullMatrix Rt(numR, ndof, Rmem);

 // y = Rt*f
 Rt.mult(f,y);

 // z = X*y
 (*X).mult(y,z);

 // f = f - z;
 f.linC(1.0, f, -1.0, z);
}

void
SingleDomainTemp::tempProject(Vector &v)
{
 if (!numR) return;

 int ndof = v.size();

 double *yMem = (double *) dbg_alloca(numR*sizeof(double));
 double *zMem = (double *) dbg_alloca(ndof*sizeof(double));
 StackVector y(ndof,yMem);
 StackVector z(ndof,zMem);

 StackFSFullMatrix Rt(numR, ndof, Rmem);

 // y = Xt*v
 (*X).trMult(v,y);


 // z = R*y
 Rt.trMult(y,z);


 // v = v - z;
 v.linAdd(-1.0, z);
}

int
SingleDomainTemp::solVecInfo()
{
 return domain->numUncon();
}

int
SingleDomainTemp::dbcVecInfo()
{
 return domain->numdof();  // Watch  out!!!!!
}

void
SingleDomainTemp::getTempTimes(double &dtemp, double &tmax)
{
 dtemp = domain->solInfo().getTimeStep();
 tmax = domain->solInfo().tmax;
}

void
SingleDomainTemp::getTempNeum(double &epsiln)
{
 epsiln = domain->solInfo().alphaTemp;
}

int
SingleDomainTemp::getTimeIntegration()
{
 return domain->solInfo().timeIntegration;
}

int
SingleDomainTemp::getAeroheatFlag()
{
 return domain->solInfo().aeroheatFlag;
}

int
SingleDomainTemp::getThermohFlag()
{
 return domain->solInfo().thermohFlag;
}


int
SingleDomainTemp::getHzemFlag()
{
 return domain->solInfo().hzemFilterFlag;
}

int
SingleDomainTemp::getZEMFlag()
{
 return domain->solInfo().hzemFlag;
}

void
SingleDomainTemp::getSteadyStateParam(int &steadyFlag, int &steadyMin,
                                         int &steadyMax, double &steadyTol)
{
 steadyFlag  = domain->solInfo().steadyFlag;
 steadyMin   = domain->solInfo().steadyMin;
 steadyMax   = domain->solInfo().steadyMax;
 steadyTol   = domain->solInfo().steadyTol;
}

void
SingleDomainTemp::getQuasiStaticParameters(double &maxVel, double &qsbeta)
{
 maxVel  = domain->solInfo().qsMaxvel;
 qsbeta  = domain->solInfo().qsBeta;
}

void
SingleDomainTemp::tempInitState(TempState<Vector> &inState)
{
  domain->initTempVector(inState.getDisp(), inState.getVeloc(), inState.getPrevVeloc());
}

void
SingleDomainTemp::computeExtForce(Vector &ext_f, double t, int tIndex, Vector &prev_f)
{

 domain->computeExtForce(ext_f, t, tIndex, kuc, prev_f);

 // apply projector ONLY for Dynamics, hzemFilterFlag is set to zero
 // in the beginning for quasistatic

   int useFilter = domain->solInfo().hzemFilterFlag;
   if (useFilter) temptrProject(ext_f);

}

void
SingleDomainTemp::preProcess()
{

 // Makes renumbering, connectivities and dofsets

 domain->preProcessing();  //Defined in Static.C and Domain.h

 int ndof = domain->numdof();

 bc  = new int[ndof];
 bcx = new double[ndof];

 domain->make_bc(bc,bcx);

 domain->make_constrainedDSA(bc);

 delete[]bc;

 domain->makeAllDOFs();
}


DynamMat
SingleDomainTemp::buildOps(double coeM, double coeC, double coeK)
{
 AllOps<double> allOps;
 DynamMat dMat;
 Rbm *rbm = 0;
 double gamma = domain->solInfo().alphaTemp;
 int algType = domain->solInfo().timeIntegration;
 int useFilter = domain->solInfo().hzemFilterFlag;
 int useHzem   = domain->solInfo().hzemFlag;

 // construct HZEM
 if(useHzem || useFilter) {
    rbm = domain->constructHzem();
    numR = rbm->numRBM();
    if(numR == 0) rbm = 0;
 }

 // construct K matrix
 if(algType == 0 && (gamma == 0.0 || (domain->solInfo().iacc_switch && !geoSource->getCheckFileInfo()->lastRestartFile)))
   allOps.K   = domain->constructDBSparseMatrix<double>();

 // construct Kuc matrix
 if(true) // this should be a test for non homogeneous dirichlet boundary conditions
   allOps.Kuc = domain->constructCuCSparse<double>();

 // construct M matrix
 if(algType == 0 && (gamma != 0.0 || domain->solInfo().modeDecompFlag || (useFilter && numR > 0))) {
   if(geoSource->getMRatio() == 0.0) // lumped M matrix
     allOps.M = new DiagMatrix(domain->getCDSA());
   else
     allOps.M = domain->constructDBSparseMatrix<double>();
 }

 // construct solver to compute correctly the initial temperature gradient v^0 = M^{-1}(f^0 - K*u^0)
 if(algType == 0 && (gamma != 0.0 && domain->solInfo().iacc_switch && !geoSource->getCheckFileInfo()->lastRestartFile)) {
   if(geoSource->getMRatio() == 0.0) { // lumped M matrix
     DiagMatrix *m = new DiagMatrix(domain->getCDSA());
     allOps.Msolver = m;
     dMat.Msolver = m;
   }
   else {
     BLKSparseMatrix *m = domain->constructBLKSparseMatrix<double>(domain->getCDSA());
     m->zeroAll();
     allOps.Msolver = m;
     dMat.Msolver = m;
   }
 }

 // assemble operators (K,Kuc,M), also construct/assemble/factor Ax=b solver where A = coek*K+coeM*M 
 if(!useHzem || (domain->solInfo().timeIntegration != 1)) // only use for quasistatics
   domain->buildOps<double>(allOps, coeK, coeM, 0.0);
 else
   domain->buildOps<double>(allOps, coeK, coeM, 0.0, rbm);

 if(useFilter) {
   fprintf(stderr," ... HZEMFilter Requested           ...\n");
   tempprojector_prep(rbm, allOps.M);
 }

 // PJSA 5-19-2008 Modal decomposition preprocessing
 int decompFlag = domain->solInfo().modeDecompFlag;
 if(decompFlag) {
   fprintf(stderr," ... Modal decomposition requested ...\n");
   modeDecompPreProcess(allOps.M);
 }

 dMat.K      = allOps.K;
 dMat.M      = allOps.M;
 kuc         = allOps.Kuc;
 dMat.dynMat = allOps.sysSolver;
 if(dMat.Msolver) dMat.Msolver->factor();

 return dMat;
}

void
SingleDomainTemp::getInternalForce(Vector& d,Vector& f)
{
  f.zero();
  FullSquareMatrix *kelArray;
  kelArray = 0;

  domain->getKtimesU(d,bcx,f,1.0,kelArray);
}

void
SingleDomainTemp::aeroHeatPreProcess(Vector& d_n, Vector& v_n, Vector& v_p)
{
  domain->aeroHeatPreProcess(d_n, v_n, v_p, bcx);
}

void
SingleDomainTemp::thermohPreProcess(Vector& d_n, Vector& v_n, Vector& v_p)
{
  domain->thermohPreProcess(d_n, v_n, v_p, bcx);
}

void
SingleDomainTemp::getInitialTime(int &initTimeIndex, double &initTime)
{
 initTimeIndex = domain->solInfo().initialTimeIndex;
 initTime      = domain->solInfo().initialTime;
}

SDTempDynamPostProcessor *
SingleDomainTemp::getPostProcessor()
{
 return new SDTempDynamPostProcessor(domain,bcx);
}

void
SDTempDynamPostProcessor
::tempdynamOutput(int tIndex, DynamMat& dMat, Vector& ext_f, TempState<Vector> &state)
{
 domain->tempdynamOutput(tIndex,bcx,dMat,ext_f,state.getDisp(),state.getVeloc(), state.getPrevVeloc());
}

int
SingleDomainTemp::getModeDecompFlag()
{
 return domain->solInfo().modeDecompFlag;
}

void
SingleDomainTemp::modeDecompPreProcess(SparseMatrix *M)
{
// Read Eigenvectors from file EIGENMODES
// ======================================
    int eigsize;

    BinFileHandler modefile("EIGENMODES" ,"r");
    modefile.read(&maxmode, 1);
    fprintf(stderr,"Number of Modes = %d\n", maxmode);

    modefile.read(&eigsize, 1);
    fprintf(stderr,"Size of EigenVector = %d\n", eigsize);

    eigmodes = new double*[maxmode];
    int i;
    for (i=0; i<maxmode; ++i)
     eigmodes[i] = new double[eigsize];

// Check if the problem sizes are identical

    int numdof = solVecInfo();

    if (eigsize != numdof)
      fprintf(stderr, "... ERROR !! EigenVector and Problem sizes differ... ");
      fflush(stderr);

    for (i=0; i<maxmode; ++i)
      modefile.read(eigmodes[i], eigsize);

// Compute Transpose(Phi_i)*M once and for all
// ===========================================

   fprintf(stderr, "Preparing Transpose(Phi_i)*M\n");

   tPhiM =  new double*[maxmode];
     for (i=0; i<maxmode; ++i)
       tPhiM[i] = new double[eigsize];

   for (i = 0 ; i<maxmode; ++i)
     M->mult(eigmodes[i], tPhiM[i]);  // taking advantage of symmetry of M and computing
                                      //   M*Phi_i instead of transpose(Phi_i)*M
}

void
SingleDomainTemp::modeDecomp(double t, int tIndex, Vector& d_n)
{

// Compute Alpha and error only if their output file is specified,
// otherwise, it wouldn't make sense

   // Compute alfa_i=PhiDiag_i*d_n

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int i, j, k;

  double *alfa = 0;
  for (i=0; i < numOutInfo; i++) {
    if(oinfo[i].interval != 0 && tIndex % oinfo[i].interval == 0) {

      int w = oinfo[i].width;
      int p = oinfo[i].precision;

      switch(oinfo[i].type) {
        case OutputInfo::ModeAlpha: {

          if(alfa == 0) {
            alfa = new double[maxmode];
            for(k = 0; k < maxmode; ++k) alfa[k] = 0.0;
            for(j = 0; j < maxmode; ++j)
              for(k = 0; k < d_n.size(); ++k)
                 alfa[j] += tPhiM[j][k]*d_n[k];
          }

          fprintf(oinfo[i].filptr, "%e  ", t);
          for(j = 0; j < maxmode; ++j)
            fprintf(oinfo[i].filptr, "% *.*E ", w, p, alfa[j]);
          fprintf(oinfo[i].filptr, "\n");

          fflush(oinfo[i].filptr);
        } break;

        case OutputInfo::ModeError: {

          if(alfa == 0) {
            alfa = new double[maxmode];
            for(k = 0; k < maxmode; ++k) alfa[k] = 0.0;
            for(j = 0; j < maxmode; ++j)
              for(k = 0; k < d_n.size(); ++k)
                 alfa[j] += tPhiM[j][k]*d_n[k];
          }
          double sumerror = 0;
          double normerror = 0;
          double sumdisp = 0;
          double normdisp = 0;
          int ersize = d_n.size();
          double *sumalfa = new double[ersize];
          double *error = new double[ersize];
          for(k = 0; k < ersize; ++k) sumalfa[k] = 0;

          for(k = 0; k < ersize; ++k)
          for(j = 0; j < maxmode; ++j)
            sumalfa[k] += alfa[j]*eigmodes[j][k];

          for(j=0; j < ersize; ++j) {
            error[j] = d_n[j]-sumalfa[j];
            sumerror += error[j]*error[j];
            sumdisp += d_n[j]*d_n[j];
          }

          normdisp = sqrt(sumdisp);
          if (normdisp == 0.0)  normerror = 0.0;
          else normerror = sqrt(sumerror)/normdisp;

          fprintf(oinfo[i].filptr, "%e % *.*E\n", t, w, p, normerror);
          fflush(oinfo[i].filptr);
          delete [] sumalfa;
          delete [] error;
        } break;
  
        default:
          break;
      }
    }
  }
  if(alfa) delete [] alfa;
}

int
SingleDomainTemp::cmdComHeat(int cmdFlag)
{
  return domain->getFileExchanger()->cmdComHeat(cmdFlag);
}
