#include <Problems.d/ModalBase.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/ModeData.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/DynamProbType.h>

extern ModeData modeData;

void DiagonalMatrix::setDiag(double val){
/*PRE: *d is allocated to double[neq]
 POST: set all entires of d to val
*/
  for(int i = 0; i < neq; ++i)
    d[i] = val;
}

//------------------------------------------------------------------------------

void DiagonalMatrix::mult(Vector &v, Vector &Av){
/*PRE: the data in Av have been allocated
 POST: return in Av, the matrix-vector product (*this)*v
*/
  for(int i = 0; i < neq; ++i)
    Av[i] = d[i] * v[i];
}

//------------------------------------------------------------------------------

void DiagonalMatrix::invertDiag(){
/*PRE: d is populated, preferably with non-zero entries
 POST: inversion of the diagonal done in place, ie entries in d are over-written
*/
  for(int i = 0; i < neq; ++i)
    d[i] = 1.0 / d[i];
}

//------------------------------------------------------------------------------

void DiagonalMatrix::reSolve(Vector &rhs){
/*PRE: d is populated with the inverse of A
 POST: returns solution to Ax = rhs; rhs is over-written with x
*/
  for(int i = 0; i < neq; ++i)
    rhs[i] *= d[i];
}

//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

void ModalBase::preProcessBase(){
/*PRE: none
 POST: initialize and populate various data members
*/
  numModes = numRBM = numFlex = 0;
  domain->preProcessing();

  const int numDofs = domain->numdof();

  int *bc = (int *) dbg_alloca(sizeof(int)*numDofs);
  bcx = new double[numDofs];
  vcx = new double[numDofs];

  int i;
  for(i = 0; i < numDofs; ++i) vcx[i] = 0.0;

  domain->make_bc(bc, bcx);
//  domain->make_constrainedDSA(bc);
  domain->make_constrainedDSA(1);
  domain->makeAllDOFs();

  numModes = modeData.numModes;  // this is assigment is temporary; it is
                                 //   overrulled in populateFlexModes
  fullTmpF.setData(new double[numDofs], numDofs);
  fullTmpGrav.setData(new double[numDofs], numDofs);
  fullAeroF.setData(new double[numDofs], numDofs);
  fullDsp.setData(new double[numDofs], numDofs);
  fullVel.setData(new double[numDofs], numDofs);
  fullAcc.setData(new double[numDofs], numDofs);
  fullPrevVel.setData(new double[numDofs], numDofs);

  fullTmpF.zeroAll();
  fullTmpGrav.zeroAll();
  fullAeroF.zeroAll();
  fullDsp.zeroAll();
  fullVel.zeroAll();
  fullAcc.zeroAll();
  fullPrevVel.zeroAll();

  prevFrc = new PrevFrc(numDofs);
  prevFrcBackup = new PrevFrc(numDofs);

  numConstr = domain->nDirichlet();
  cDofIdx = new int[numConstr];
  for(i = 0; i < numConstr; ++i){
    cDofIdx[i] = domain->getCDSA()->locate(domain->getDBC()[i].nnum,
      1 << domain->getDBC()[i].dofnum);
  }
}

//------------------------------------------------------------------------------

void ModalBase::populateRBModes(){

  // compute the rigid body modes
//  Rbm *rbm = domain->constructAllRbm();
  Rbm *rbm = domain->constructRbm();
  numRBM   = rbm->numRBM();
  modesRB  = new Vector[numRBM];

  rbm->getRBMs(modesRB);
  rbm->getxyzRot(0, cg);      // cg temporarily stores the point about
                              //   which rotation RBMs were calculated
}

//------------------------------------------------------------------------------

void ModalBase::populateFlexModes(double scale, bool readAll){
/*PRE: preProcessBase has been called
       scale has default value of 1.0; readAll has default value of 0
       the modes in modeData are : mode^T.M.mode = identity
 POST: populate modesFl with data from modeData multiplied by scale
       populate freqs with the circular frequency of each flexible mode
 NOTE: if readAll, then also include in modesFL those modes with zero frequency
*/
  const int numNodes   = modeData.numNodes;
  const int numDofs  = domain->numdof();

  // count the number of flexible modes
  numFlex = 0;
  int iMode;
  for(iMode = 0; iMode < numModes; ++iMode){
    if( readAll || (modeData.frequencies[iMode] != 0.0) )
      ++numFlex;
  }
  modesFl = new Vector[numFlex];
  freqs   = new double[numFlex];

  int iModeFl = 0, iNode, dof;
  for(iMode = 0; iMode < numModes; ++iMode){
    if( readAll || (modeData.frequencies[iMode] != 0.0) ){
      modesFl[iModeFl].setData(new double[numDofs], numDofs);
      modesFl[iModeFl].zeroAll();
      freqs[iModeFl] = modeData.frequencies[iMode] * 8 * atan(1.);

      for(iNode = 0; iNode < numNodes; ++iNode){

        dof = domain->getCDSA()->locate(iNode, DofSet::Xdisp);
        if(dof >= 0)   modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][0];

        dof = domain->getCDSA()->locate(iNode, DofSet::Ydisp);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][1];

        dof = domain->getCDSA()->locate(iNode, DofSet::Zdisp);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][2];

        dof = domain->getCDSA()->locate(iNode, DofSet::Xrot);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][3];

        dof = domain->getCDSA()->locate(iNode, DofSet::Yrot);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][4];

        dof = domain->getCDSA()->locate(iNode, DofSet::Zrot);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][5];
      }
      ++iModeFl;
    }
  }
  
}

//------------------------------------------------------------------------------

void ModalBase::initStateBase(Vector& dsp, Vector& vel,
  Vector& acc, Vector& vel_p, int idxOffset){

  dsp.zero();
  vel.zero();
  acc.zero();
  vel_p.zero();

  ControlInfo *cinfo = geoSource->getCheckFileInfo();
  SolverInfo &sinfo = domain->solInfo();

  if (cinfo->lastRestartFile) {
     fprintf(stderr, " ... Restarting From a Previous Run ...\n");
     int fn = open(cinfo->lastRestartFile,O_RDONLY );
     if(fn >= 0) {
       int vsize, restartTIndex;
       double restartT;
       int readSize = read(fn, &vsize, sizeof(int));
       if(vsize != dsp.size() || readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       readSize = read(fn, &restartTIndex, sizeof(int));
       if(readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTimeIndex = 0;
       else
         sinfo.initialTimeIndex = restartTIndex;
       readSize = read(fn, &restartT, sizeof(double));
       if(readSize != sizeof(double))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTime = 0.0;
       else
         sinfo.initialTime = restartT;
       readSize = read(fn, dsp.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, vel.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, vel_p.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double))) {
         fprintf(stderr," *** WARNING: Older version of restart file"
                        " -- Missing velocity field is set to zero\n");
         vel_p.zero();
       }

       close(fn);
     } 
     else {
       fprintf(stderr, " *** ERROR: Restart file could not be opened\n");
     }
  }
  else  {

    int j;
    for(j = 0; j <  domain->numInitDisp(); ++j)
      dsp[domain->getInitDisp()[j].nnum + idxOffset] += domain->getInitDisp()[j].val;

    for(j = 0; j <  domain->numInitVelocity(); ++j)
      vel[domain->getInitVelocity()[j].nnum + idxOffset] += domain->getInitVelocity()[j].val;
  }
}

//------------------------------------------------------------------------------

void ModalBase::outputModal(SysState<Vector>& state, Vector& extF, int tIndex){
/*PRE:
 POST:
*/
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  SolverInfo &sinfo = domain->solInfo();
  int i, w, p, iMode;
  double time = tIndex * domain->solInfo().getTimeStep();
  Vector &mDsp = state.getDisp();
  Vector &mVel = state.getVeloc();
  Vector &mVel_p = state.getPrevVeloc();

  if (sinfo.nRestart > 0)
    domain->writeRestartFile(time, tIndex, mDsp, mVel, mVel_p);

  for(i = 0; i < numOutInfo; ++i){
    if((oinfo[i].interval != 0) && (tIndex % oinfo[i].interval == 0)){

      w = oinfo[i].width;
      p = oinfo[i].precision;

      switch(oinfo[i].type){

        case OutputInfo::ModalDsp:
//          fprintf(oinfo[i].filptr, "# modal displacements\n");
          fprintf(oinfo[i].filptr, " % *.*E ", w, p, time);
          for(iMode = 0; iMode < mDsp.size(); ++iMode){
            fprintf(oinfo[i].filptr, "  % *.*E", w, p, mDsp[iMode]);
          }
          fprintf(oinfo[i].filptr, "\n");
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::ModalExF:
//          fprintf(oinfo[i].filptr, "# modal external forces\n");
          fprintf(oinfo[i].filptr, " % *.*E ", w, p, time);
          for(iMode = 0; iMode < extF.size(); ++iMode){
            fprintf(oinfo[i].filptr, "  % *.*E", w, p, extF[iMode]);
          }
          fprintf(oinfo[i].filptr, "\n");
          fflush(oinfo[i].filptr);
          break;

        default:
          break;
      }
    }
  }
}

//------------------------------------------------------------------------------
