#include <cmath>
#include <climits> //--- UH

#include <Driver.d/Domain.h>
#include <Problems.d/NLModalDescr.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/Dynam.h>
#include <Driver.d/DynamProbType.h>

void NLModalOpSolver::reSolve(Vector &rhs){
/*PRE:
 POST:
*/
#ifdef VICTORY
  fflush(stderr);
  rhs.print("OpSolver.reSolve rhs", "rhs");
  fflush(stderr);
#endif

//  exit(0);
  mat.reSolve(rhs.data());

#ifdef VICTORY
  rhs.print("OpSolver.reSolve solution", "sol");
  fflush(stderr);
  fprintf(stderr, "end OpSolver resolve ---------------------------------------\n");
#endif
}


//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

NLModalDescr::NLModalDescr(Domain *d) : ModalBase(d){}

//------------------------------------------------------------------------------

void NLModalDescr::rotateVector(Vector& globalF, Vector& rotF, double glR[3][3]){
/*PRE: globalF is full-sized and is in inertial coordinates
 POST: apply the transpose of glR to globalF and return in rotF
*/
  int i, xdof, ydof, zdof;
  for(i = 0; i < domain->numNode(); ++i){
    xdof = domain->getCDSA()->locate(i, DofSet::Xdisp);
    ydof = domain->getCDSA()->locate(i, DofSet::Ydisp);
    zdof = domain->getCDSA()->locate(i, DofSet::Zdisp);

    rotF[xdof] = glR[0][0]*globalF[xdof] + glR[1][0]*globalF[ydof]
      + glR[2][0]*globalF[zdof];
    rotF[ydof] = glR[0][1]*globalF[xdof] + glR[1][1]*globalF[ydof]
      + glR[2][1]*globalF[zdof];
    rotF[zdof] = glR[0][2]*globalF[xdof] + glR[1][2]*globalF[ydof]
      + glR[2][2]*globalF[zdof];

    xdof = domain->getCDSA()->locate(i, DofSet::Xrot);
    ydof = domain->getCDSA()->locate(i, DofSet::Yrot);
    zdof = domain->getCDSA()->locate(i, DofSet::Zrot);

    rotF[xdof] = glR[0][0]*globalF[xdof] + glR[1][0]*globalF[ydof]
      + glR[2][0]*globalF[zdof];
    rotF[ydof] = glR[0][1]*globalF[xdof] + glR[1][1]*globalF[ydof]
      + glR[2][1]*globalF[zdof];
    rotF[zdof] = glR[0][2]*globalF[xdof] + glR[1][2]*globalF[ydof]
      + glR[2][2]*globalF[zdof];
  }
}

//------------------------------------------------------------------------------

void NLModalDescr::projectForceNL(Vector& fullF,
  Vector& modalF, ModalGeomState* mgs){
/*PRE: fullF is the force vector obtained from Domain::computeExtForce
 POST: return in modalF the force conjugate to the non-lin modal formulation
 SIDE: fullFrot, a private member of *this, is modified by call to rotateVector
 NOTE: assumes all RBMs are present
 TODO: use mgs.q and modesPr in projection onto rigid body rotations
*/
  int i;
  for(i = 0; i < 3; ++i)
    modalF[i] = modesRB[i] * fullF;

  rotateVector(fullF, fullFrot, mgs->glR);

  for(i = 3; i < 6; ++i)
    modalF[i] = modesRB[i] * fullFrot;

  for(i = 0; i < numFlex; ++i)
    modalF[i+numRBM] = modesFl[i] * fullFrot;

  for(i = 0; i < numConstr; ++i)
    modalF[i+numRBM+numFlex] = 0.0;
} 

//------------------------------------------------------------------------------

void NLModalDescr::filterForce(Vector& fullF){
/*PRE:
 POST: fullF filtered appropriately
*/
//  modalF[1] = 0.0;  // hard-coded test case to remove force on y-translation
  fprintf(stderr, "rbmFilters: %d %d %d %d %d %d\n",
    domain->solInfo().rbmFilters[0], domain->solInfo().rbmFilters[1],
    domain->solInfo().rbmFilters[2], domain->solInfo().rbmFilters[3],
    domain->solInfo().rbmFilters[4], domain->solInfo().rbmFilters[5]);

  Vector tmpF(fullF.size(), 0.0);
  // double tmp1[numFilters], tmp2[numFilters];
  double *tmp1 = (double *) dbg_alloca(sizeof(double)*numFilters);
  double *tmp2 = (double *) dbg_alloca(sizeof(double)*numFilters);
  int i, j;

  // store in tmp1 the product Psi^T.fullF
  for(i = 0; i < numFilters; ++i)
    tmp1[i] = modesRB[rbmIdx[i]] * fullF;

  // store in tmp2 the product invPsiTMPsi.tmp1
  for(i = 0; i < numFilters; ++i){
    tmp2[i] = 0;
    for(j = 0; j < numFilters; ++j)
      tmp2[i] += (*invPsiTMPsi)[i][j] * tmp1[j];
  }

  // store in tmpF the product MPsi.tmp2
  for(i = 0; i < fullF.size(); ++i){
    for(j = 0; j < numFilters; ++j){
      tmpF[i] += (*MPsi)[i][j] * tmp2[j];
    }
  }
  fullF -= tmpF;
}


//------------------------------------------------------------------------------

void NLModalDescr::expandFlex(Vector& modalV, Vector& fullV){
/*PRE: modalV has size numFlex; fullV has size domain->numUncon()
 POST: project modalV into fullV using modes modesFl
       if there are zero flexible modes, then return fullV populated with zeros
*/
  fullV.zeroAll();
  for(int b = 0; b < numFlex; ++b)
    fullV.linAdd(modalV[b], modesFl[b]);
}

//------------------------------------------------------------------------------

void NLModalDescr::calcTmps(ModalGeomState &mgs){
/*PRE: called by getStiffAndForce
 POST: populate the temporary variables Bq, Bqq, Bhq, Bho, Btq, Bqo, c22a, c22b
*/
  Vector &vel = mgs.vel;
  int b, g, i, j;
  for(g = 0; g < numFlex; ++g){
    for(i = 0; i < 3; ++i){

      // populate Bho
      Bho[g][i] = 0.0;
      for(j = 0; j < 3; ++j){
        Bho[g][i] += Bh[g][i][j] * vel[j+3];

        // populate Bq
        Bq[g][i][j] = 0.0;
        for(b = 0; b < numFlex; ++b){
          Bq[g][i][j] += Bp[b][g][i][j] * mgs.q[b];
        }
      }

      // populate Btq
      Btq[g][i] = 0.0;
      for(b = 0; b < numFlex; ++b){
        Btq[g][i] += Bt[g][b][i] * mgs.q[b];

        // populate Bo
        Bo[b][g][i] = 0.0;
        for(j = 0; j < 3; ++j){
          Bo[b][g][i] += Bp[b][g][i][j] * vel[j+3];
        }
      }
    }

    for(b = 0; b < numFlex; ++b){
      Bto[b][g] = 0.0;
      for(i = 0; i < 3; ++i){
        Bto[b][g] += Bt[b][g][i] * vel[i+3];
      }
    }
  }

  for(i = 0; i < 3; ++i){

    c22b[i] = 0.0;
    q_modesPr[i].zeroAll();

    for(j = 0; j < 3; ++j){

      // populate Bqq and Bhq
      Bqq[i][j] = 0.0;
      Bhq[i][j] = 0.0;
      for(b = 0; b < numFlex; ++b){
        Bqq[i][j] += Bq[b][i][j] * mgs.q[b];
        Bhq[i][j] += Bh[b][i][j] * mgs.q[b];
      }
      // begin populating c22b
      c22b[i] += ( J[i][j] + Bqq[i][j] + Bhq[i][j] ) * vel[j+3];
      // populate c22a
      c22a[i][j] = J[i][j] + Bqq[i][j] + 2*Bhq[i][j];
    }

    // finish populating c22b; populate q_modesPr
    for(g = 0; g < numFlex; ++g){
      c22b[i] += Btq[g][i] * vel[numRBM+g];

      q_modesPr[i].linAdd(mgs.q[g], modesPr[g][i]);

      // populate Bqo
      Bqo[g][i] = 0.0;
      for(j = 0; j < 3; ++j){
        Bqo[g][i] += Bq[g][i][j] * vel[j+3];
      }
    }
  }
  evalConstraintsAndDerivatives(mgs);
}

//------------------------------------------------------------------------------

void NLModalDescr::processLastOutput() {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

//------------------------------------------------------------------------------

void NLModalDescr::preProcess(){
/*PRE: none
 POST: various, see comments throughout this member function
*/
  preProcessBase();
  fullFrot.setData(new double[fullTmpF.size()], fullTmpF.size());
  mass = domain->computeStructureMass();
  populateRBModes();
  populateFlexModes(sqrt(mass));
  numModes = numRBM + numFlex;

  solver = new NLModalOpSolver;
  solver->mat.setNewSize(solVecInfo(), solVecInfo(), 0.0);
  tolerance = domain->solInfo().getNLInfo().tolRes;

  // populate modesPr
  modesPr = new Vector*[numFlex];
  const int numNodes = domain->numNode();
  const int numDofs  = domain->numdof();

  int b, iNode, xdof, ydof, zdof; //, xrot, yrot, zrot;
  for(b = 0; b < numFlex; ++b){

    modesPr[b] = new Vector[3];
    modesPr[b][0].setData(new double[numDofs], numDofs);
    modesPr[b][1].setData(new double[numDofs], numDofs);
    modesPr[b][2].setData(new double[numDofs], numDofs);

    modesPr[b][0].zeroAll(); modesPr[b][1].zeroAll(); modesPr[b][2].zeroAll();

    for(iNode = 0; iNode < numNodes; ++iNode){

      xdof = domain->getCDSA()->locate(iNode, DofSet::Xdisp);
      ydof = domain->getCDSA()->locate(iNode, DofSet::Ydisp);
      zdof = domain->getCDSA()->locate(iNode, DofSet::Zdisp);
/*
      xrot = domain->getCDSA()->locate(iNode, DofSet::Xrot);
      yrot = domain->getCDSA()->locate(iNode, DofSet::Yrot);
      zrot = domain->getCDSA()->locate(iNode, DofSet::Zrot);
*/
      if(xdof >= 0 && ydof >= 0){
        modesPr[b][2][xdof] = -modesFl[b][ydof];
        modesPr[b][2][ydof] =  modesFl[b][xdof];
      }
      if(ydof >= 0 && zdof >= 0){
        modesPr[b][0][ydof] = -modesFl[b][zdof];
        modesPr[b][0][zdof] =  modesFl[b][ydof];
      }
      if(zdof >= 0 && xdof >= 0){
        modesPr[b][1][zdof] = -modesFl[b][xdof];
        modesPr[b][1][xdof] =  modesFl[b][zdof];
      }
    }
  }

  // build the finite element mass matrix
  AllOps<double> allOps;
  allOps.M = domain->constructDBSparseMatrix<double>();
  domain->makeSparseOps<double>(allOps, 0.0, 1.0, 0.0); // NULL, NULL);

  GenSparseMatrix<double> *massMat = allOps.M;

// testing stuff ************************************************
/*
  ((DBSparseMatrix*) massMat)->print();

  allOps.K = domain->constructDBSparseMatrix<double>();
  domain->makeSparseOps<double>(allOps, 1.0, 0.0, 0.0);

  GenSparseMatrix<double> *stiffnessMat = allOps.K;
  fprintf(stderr, "\nStiffness matrix:\n");
  ((DBSparseMatrix*) stiffnessMat)->print();

  exit(0);
*/
/*
  Vector tmp(12), M_tmp(12);
  for(b = 0; b < 12; ++b){
    tmp.zero();
    tmp[b] = 1.0;
    massMat->mult(tmp, M_tmp);
    for(iNode = 0; iNode < 12; ++iNode)
      fprintf(stderr, "%28.20e", M_tmp[iNode]);
    fprintf(stderr, "\n");
  }
  exit(0);
*/
//***************************************************************

  // locate the numeric center of gravity
  Vector M_phir[3];
  M_phir[0].setData(new double[numDofs], numDofs);
  M_phir[1].setData(new double[numDofs], numDofs);
  M_phir[2].setData(new double[numDofs], numDofs);

  int i, j;
  for(i = 0; i < 3; ++i){
    massMat->mult(modesRB[i+3], M_phir[i]); // NOTE: modesRB[i+3] assumes
  }                                          // all RBMs are present
  double dxyz[3][3];
  for(i = 0; i < 3; ++i)
    for(j = 0; j < 3; ++j)
      dxyz[i][j] = -1.0/mass * (modesRB[i] * M_phir[j]);

  cg[0] = cg[0] + (dxyz[2][1] - dxyz[1][2]) / 2.0;
  cg[1] = cg[1] + (dxyz[0][2] - dxyz[2][0]) / 2.0;
  cg[2] = cg[2] + (dxyz[1][0] - dxyz[0][1]) / 2.0;

  // re-center origin at cg
  for(iNode = 0; iNode < numNodes; ++iNode){
    domain->getNodes().getNode(iNode).x -= cg[0];
    domain->getNodes().getNode(iNode).y -= cg[1];
    domain->getNodes().getNode(iNode).z -= cg[2];
  }

  // adjust rotational modes : the point about which they rotate is the cg
  // NOTE: the following three lines assume all RBMs are present
  modesRB[3].linAdd(dxyz[1][0], modesRB[1], dxyz[2][0], modesRB[2]);
  modesRB[4].linAdd(dxyz[0][1], modesRB[0], dxyz[2][1], modesRB[2]);
  modesRB[5].linAdd(dxyz[0][2], modesRB[0], dxyz[1][2], modesRB[1]);

  // calc the mass moment of inertia tensor
  fprintf(stderr, "mass moment of inertia\n");
  for(j = 0; j < 3; ++j){
    massMat->mult(modesRB[j+3], M_phir[j]);
    for(i = 0; i < 3; ++i){
      J[i][j] = modesRB[i+3] * M_phir[j]; // NOTE: assumes all RBMs are present
      fprintf(stderr, "  %20.12e", J[i][j]);
    }
    fprintf(stderr, "\n");
  }

  // initialize the coupling coefficients -- Bp, Bh and Bt -- initialize temporary
  //   variables and create temporary var M_phip
  Bp = new double***[numFlex];    // Bp[numFlex][numFlex][3][3]
  Bh = new double**[numFlex];     // Bh[numFlex][3][3]
  Bt = new double**[numFlex];     // Bt[numFlex][numFlex][3]

  Bq  = new double**[numFlex];
  Bo  = new double**[numFlex];
  Bho = new double*[numFlex];
  Btq = new double*[numFlex];
  Bqo = new double*[numFlex];
  Bto = new double*[numFlex];

  q_modesPr[0].setData(new double[numDofs], numDofs);
  q_modesPr[1].setData(new double[numDofs], numDofs);
  q_modesPr[2].setData(new double[numDofs], numDofs);

  Vector **M_phip = new Vector*[numFlex];

  int g;
  for(b = 0; b < numFlex; ++b){
    Bp[b] = new double**[numFlex];
    Bh[b] = new double*[3];
    Bt[b] = new double*[numFlex];

    Bq[b]  = new double*[3];
    Bo[b]  = new double*[numFlex];
    Bho[b] = new double[3];
    Btq[b] = new double[3];
    Bqo[b] = new double[3];
    Bto[b] = new double[numFlex];

    for(g = 0; g < numFlex; ++g){
      Bp[b][g] = new double*[3];
      Bt[b][g] = new double[3];
      Bo[b][g] = new double[3];

      for(i = 0; i < 3; ++i)
        Bp[b][g][i] = new double[3];
    }

    M_phip[b] = new Vector[3];

    for(j = 0; j < 3; ++j){
      Bh[b][j] = new double[3];
      Bq[b][j] = new double[3];
      M_phip[b][j].setData(new double[numDofs], numDofs);
      massMat->mult(modesPr[b][j], M_phip[b][j]);
    }
  }

  // populate the coupling coefficients
  for(b = 0; b < numFlex; ++b){

    for(g = 0; g < numFlex; ++g){

      for(j = 0; j < 3; ++j){

        for(i = 0; i < 3; ++i){
          Bp[b][g][i][j] = 0.5 * ( (modesPr[b][i] * M_phip[g][j])
            + (modesPr[g][i] * M_phip[b][j]) );
        }
        Bt[b][g][j] = 0.5* ( (modesFl[b] * M_phip[g][j])
          - (modesFl[g] * M_phip[b][j]) );
      }
    }
    for(i = 0; i < 3; ++i){
      for(j = 0; j < 3; ++j){
        Bh[b][i][j] = 0.5 * ( (modesRB[i+3] * M_phip[b][j]) // NOTE: assumes all
          + (modesRB[j+3] * M_phip[b][i]) );                //   RBMs are present
      }
    }
  }

  // initialize dxhat, dflexhat, flexhat and constr
  dxhat    = new double*[numConstr];    // dxhat[constr#][rotation#]
  dflexhat = new double**[numConstr];   // dflexhat[constr#][flexmode#][rotation#]
  flexhat  = new double*[numConstr];    // flexhat[constr#][felxmode#]
  constr   = new double[numConstr];     // constr[constr#]

  dxhatnp1    = new double*[numConstr];
  dflexhatnp1 = new double**[numConstr];
  flexhatnp1  = new double*[numConstr];

// new stuff 040902 ------------------------------
  dhdp  = new double*[numConstr];   // dhdp is constant in time, so maybe populate it in preProcess
  dhdth = new double*[numConstr];
  dhdq  = new double*[numConstr];
  d2h   = new double**[numConstr];  // d2h[constr#][flexmode#][rotation#]

  dhdp_np1  = new double*[numConstr];
  dhdth_np1 = new double*[numConstr];
  dhdq_np1  = new double*[numConstr];
// -----------------------------------------------

  for(i = 0; i < numConstr; ++i){
    dxhat[i] = new double[3];
    dflexhat[i] = new double*[numFlex];
    flexhat[i] = new double[numFlex];

    dxhatnp1[i] = new double[3];
    dflexhatnp1[i] = new double*[numFlex];
    flexhatnp1[i] = new double[numFlex];


// new stuff 040902 ------------------------------
    dhdp[i]  = new double[3];
    dhdth[i] = new double[3];
    dhdq[i]  = new double[numFlex];
    d2h[i]   = new double*[numFlex];

    dhdp_np1[i]  = new double[3];
    dhdth_np1[i] = new double[3];
    dhdq_np1[i]  = new double[numFlex];
// -----------------------------------------------


    for(b = 0; b < numFlex; ++b){
      dflexhat[i][b] = new double[3];
      dflexhatnp1[i][b] = new double[3];
      d2h[i][b] = new double[3];
    }
  }
/*
  fprintf(stderr, "total mass: %23.16e\n", mass);
  char msg[100];
  for(b = 0; b < 6; ++b){
    sprintf(msg, "modeRB %d", b);
    modesRB[b].print(msg);
  }
  fprintf(stderr, "\n");
  for(b = 0; b < numFlex; ++b){
    sprintf(msg, "modeFl %d", b);
    modesFl[b].print(msg);
  }

  exit(0);
*/
  buildRBMFilter((DBSparseMatrix*)massMat);
}

//------------------------------------------------------------------------------

void NLModalDescr::buildRBMFilter(DBSparseMatrix* massMat){
/*PRE: domain.sinfo.rbmFilters[6] is populated with 0s and 1s only; 0 indicates
         that rbm is to be left untouched; 1, the mode is to be removed
 POST: populate MPsi and invPsiTMPsi, Psi being the matrix collecting the rigid
         body modes to be filtered
*/
  int i, j, k;
  int dim = fullTmpF.size();

  // count the number of rigid body modes to be filtered
  numFilters = domain->solInfo().rbmFilters[0] + domain->solInfo().rbmFilters[1]
    + domain->solInfo().rbmFilters[2] + domain->solInfo().rbmFilters[3]
    + domain->solInfo().rbmFilters[4] + domain->solInfo().rbmFilters[5];

  if(numFilters == 0){
    invPsiTMPsi = 0;
    MPsi = 0;
    return;
  }
  invPsiTMPsi = new FullM(numFilters);
  MPsi = new FullM(dim, numFilters);

  // store the indices of the rbms to be filtered
  rbmIdx = new int[numFilters];
  j = 0;
  for(k = 0; k < 6; ++k){
    if(domain->solInfo().rbmFilters[k]){
      rbmIdx[j] = k;
      ++j;
    }
  }
/*
  for(k = 0; k < numFilters; ++k)
    fprintf(stderr, "rbmIdx[%d]: %d\n", k, rbmIdx[k]);
  exit(0);
*/
  Vector tmpVec(dim);
  for(k = 0; k < numFilters; ++k){
    massMat->mult(modesRB[rbmIdx[k]], tmpVec);
    for(j = 0; j < dim; ++j)
      (*MPsi)[j][k] = tmpVec[j];
  }

  for(i = 0; i < numFilters; ++i){
    for(j = 0; j < numFilters; ++j){
      (*invPsiTMPsi)[i][j] = 0.0;
      for(k = 0; k < dim; ++k){
        (*invPsiTMPsi)[i][j]+= modesRB[rbmIdx[i]][k] * (*MPsi)[k][j];
      }
    }
  }

  fprintf(stderr, "total mass: %e\n", mass);
  MPsi->print("MPsi");

  invPsiTMPsi->print("inPsiTMPsi pre-factor");
  invPsiTMPsi->factor();
  invPsiTMPsi->print("inPsiTMPsi post-factor");
}

//------------------------------------------------------------------------------

void NLModalDescr::computeTimeInfo(){
/*PRE:
 POST:
*/
  tFinal = domain->solInfo().tmax;
  dt     = domain->solInfo().getTimeStep();
  delta  = 0.5 * dt;

  mcoef = 1.0;
  ccoef = delta;
  kcoef = delta*delta;
  hcoef = 1.0;//kcoef;

  maxStep = (int) (tFinal/dt);
  remainingTime = tFinal - maxStep*dt;
}

//------------------------------------------------------------------------------

void NLModalDescr::getConstForce(Vector &constF){
/*PRE: none
 POST: populate constF with zeros; leave for now
 NOTE: constF has size numConstr longer than fullTmpGrav 
*/
  domain->computeConstantForce(constF);
/*
  fullTmpGrav.zero();
  constF.zero();

  if( domain->gravityFlag()  ) domain->buildGravityForce<double>(fullTmpGrav);
  if( domain->pressureFlag() ) domain->buildPressureForce<double>(fullTmpGrav);
*/

//  projectForceNL(fullTmpGrav, constF);
}

//------------------------------------------------------------------------------

void NLModalDescr::getExternalForce(Vector &extF, Vector &gravF, int tIndex,
  double time, ModalGeomState* mgs, Vector &elemIntF, Vector& aeroF){
/*PRE: nothing special
 POST: non-linear, modalized force in extF
 NOTE: arguements to this function are of size numRBM+numFlex+numConstr
         while args for domain->computeExtForce are size domain->numdof
 TODO: write filter force member function and use it in this memeber function
*/
  initIterState(*mgs);
  domain->computeExtForce4<double>(fullTmpF, fullTmpGrav, time);

  if(MPsi && invPsiTMPsi){ filterForce(fullTmpF); }
  projectForceNL(fullTmpF, extF, mgs);
}

//------------------------------------------------------------------------------

int NLModalDescr::getInitState(Vector &dsp, Vector &vel,
  Vector &acc, Vector &vel_p){
/*PRE: arguements to this function are size returned by solVecInfo()
 POST: populate the arguments of this function with the initial state
*/
  initStateBase(dsp, vel, acc, vel_p, 6);
//  vel[0] = 7.5;
//  vel[5] = 5.;
  return 0;
}

//------------------------------------------------------------------------------
void NLModalDescr::readRestartFile(Vector &dsp, Vector &vel, Vector &acc,
  Vector &vel_p, ModalGeomState &mgs){
/*PRE: arguements vel and acc are populated
 POST: populate mgs with vel and acc
 NOTE: this member function does not actually read the restart file
       when this function is called from the driver, dsp, vel and acc have been
         populated by getInitState; this function is simply a way of
         populateding mgs with the specified initial values
       initial displacements are populated by updatePrescribedDisplacements
*/
  mgs.vel = vel;
  mgs.acc = acc;
}

//------------------------------------------------------------------------------

void NLModalDescr::updatePrescribedDisplacement(ModalGeomState *mgs){
/*PRE: 
 POST: populate *mgs with initial displacements and initialize data members of
         *mgs that store the displacements at time step n+1
*/
  int j;
  for(j = 0; j <  domain->numInitDisp(); ++j)
    mgs->q[domain->getInitDisp()[j].nnum] += domain->getInitDisp()[j].val;
}

//------------------------------------------------------------------------------

void NLModalDescr::initIterState(ModalGeomState &mgs){
/*PRE: nothing special
 POST: populate mgs with the state appropriate to begin the Newton iterations
 NOST: mgs.vel are left unchanged
*/
#ifdef VICTORY
  mgs.printState(" ^^^^^ pre initIterState ^^^^^^^^^^^^^^^");
#endif

  mgs.acc.copy(0.0);
  mgs.lam.copy(0.0);

  mgs.glTnp1[0] = mgs.glT[0] + 2*delta*mgs.vel[0];
  mgs.glTnp1[1] = mgs.glT[1] + 2*delta*mgs.vel[1];
  mgs.glTnp1[2] = mgs.glT[2] + 2*delta*mgs.vel[2];

  mgs.glT[0] += delta*mgs.vel[0];
  mgs.glT[1] += delta*mgs.vel[1];
  mgs.glT[2] += delta*mgs.vel[2];

  double dtheta[3] = {2*delta*mgs.vel[3], 2*delta*mgs.vel[4], 2*delta*mgs.vel[5]};
  double dR[3][3];

  vec_to_mat(dtheta, dR);
  mat_mult_mat(dR, mgs.glR, mgs.glRnp1, 0);

  dtheta[0] = delta*mgs.vel[3];
  dtheta[1] = delta*mgs.vel[4];
  dtheta[2] = delta*mgs.vel[5];

  inc_rottensor(dtheta, mgs.glR);

  for(int i = 0; i < numFlex; ++i){
    mgs.qnp1[i] = mgs.q[i] + 2*delta*mgs.vel[numRBM+i];
    mgs.q[i] += delta*mgs.vel[numRBM+i];
  }

#ifdef VICTORY
  mgs.printState(" ^^^^^ post initIterState ^^^^^^^^^^^^^^^");
#endif
}

//------------------------------------------------------------------------------

double NLModalDescr::getStiffAndForce(ModalGeomState &mgs,
  Vector &res, double midtime, ModalGeomState*){
/*PRE:
 POST:
*/
  Vector &vel = mgs.vel;
  Vector &acc = mgs.acc;

  calcTmps(mgs);

//  projectForceNL(fullTmpF, res, &mgs);

//************************************************
//  res.print("res from getStiffAndForce", "res");
//************************************************

  int b, i, s;
  // populate blocks 11 and 22
  for(s = 0; s < 3; ++s){

    solver->mat[s][s] = -mcoef * mass;

    for(i = 0; i < 3; ++i){
      solver->mat[i+3][s+3] = mcoef * (-J[i][s] - Bqq[i][s] - 2*Bhq[i][s]);
      for(b = 0; b < numFlex; ++b){
        solver->mat[i+3][s+3] +=
          ccoef * ( -2 * vel[numRBM+b] * (Bq[b][i][s] + Bh[b][i][s]) );
      }
    }

    solver->mat[3][s+3] += ccoef * ( c22a[1][s]*vel[5] - c22a[2][s]*vel[4] );
    solver->mat[4][s+3] += ccoef * ( c22a[2][s]*vel[3] - c22a[0][s]*vel[5] );
    solver->mat[5][s+3] += ccoef * ( c22a[0][s]*vel[4] - c22a[1][s]*vel[3] );

    res[s]   = -res[s];
    res[s+3] = -res[s+3];
  }

  solver->mat[3][4] += -ccoef * c22b[2];  solver->mat[3][5] +=  ccoef * c22b[1];
  solver->mat[4][3] +=  ccoef * c22b[2];  solver->mat[4][5] += -ccoef * c22b[0];
  solver->mat[5][3] += -ccoef * c22b[1];  solver->mat[5][4] +=  ccoef * c22b[0];

  // populate block 33
  int g, h, z;
  for(z = 0; z < numFlex; ++z){

    b = numRBM + z;
    for(h = 0; h < numFlex; ++h){

      g = numRBM + h;
      solver->mat[b][g] = ccoef * 2 * Bto[h][z];

      for(i = 0; i < 3; ++i){
        solver->mat[b][g] += kcoef*(-Bt[z][h][i]*acc[i+3] + Bo[z][h][i]*vel[i+3]);
      }
    }
    solver->mat[b][b] -= mass * ( mcoef + kcoef*freqs[z]*freqs[z]);
//   for derivTest
//    solver->mat[b][b] -= mass * kcoef*freqs[z]*freqs[z];
  }

//************************************************
//  solver->mat.print("pre block 33");
//************************************************

  // populate blocks 23 and 32
  double tmp, k23c[3];           // tmp variables; k23c is the cross term in
  for(h = 0; h < numFlex; ++h){  //   the 23 block of the K matrix
                                 //   tmp is 2(Bq + Bh).omega
    g = numRBM + h;
    for(s = 0; s < 3; ++s){

      solver->mat[s+3][g] = solver->mat[g][s+3] = -mcoef * Btq[h][s];
//      for derivTest
//      solver->mat[s+3][g] = solver->mat[g][s+3] = 0.0;
      k23c[s] = 0.0;
      for(i = 0; i < 3; ++i){
        tmp = 2 * (Bq[h][s][i] + Bh[h][s][i]) * vel[i+3];
        solver->mat[s+3][g] -= ccoef * tmp;
        solver->mat[g][s+3] += ccoef * tmp;

        solver->mat[s+3][g] += -2 * kcoef * (Bq[h][s][i] + Bh[h][s][i]) * acc[i+3];

        k23c[s] += tmp;
      }

      for(b = 0; b < numFlex; ++b){
        z = numRBM + b;
        solver->mat[g][s+3] += 2 * ccoef * Bt[b][h][s] * vel[z];
        solver->mat[s+3][g] += kcoef * ( -2 * Bo[h][b][s] * vel[z]
          - Bt[b][h][s] * acc[z] );

        k23c[s] += Bt[b][h][s] * vel[z];
      }
      // note to self: following line is in question
      solver->mat[s+3][g] += kcoef * (modesPr[h][s] * fullFrot);
    }

    solver->mat[3][g] += (ccoef*Btq[h][1] + kcoef*k23c[1]) * vel[5]
      - (ccoef*Btq[h][2] + kcoef*k23c[2]) * vel[4];

    solver->mat[4][g] += (ccoef*Btq[h][2] + kcoef*k23c[2]) * vel[3]
      - (ccoef*Btq[h][0] + kcoef*k23c[0]) * vel[5];

    solver->mat[5][g] += (ccoef*Btq[h][0] + kcoef*k23c[0]) * vel[4]
      - (ccoef*Btq[h][1] + kcoef*k23c[1]) * vel[3];

    res[g] = -res[g] + mass * mgs.q[h] * freqs[h] * freqs[h];
  }

  // terms due to constraints
  int r, ar, be;
  for(r = 0; r < numConstr; ++r){
      ar = r + numRBM + numFlex;
    for(b = 0; b < numFlex; ++b){
      be = b + numRBM;
      for(s = 0; s < 3; ++s){
        // blocks 23 and 32
      	solver->mat[s+3][be] += hcoef * kcoef * d2h[r][b][s] * mgs.lam[r];
        solver->mat[be][s+3] += hcoef * kcoef * d2h[r][b][s] * mgs.lam[r];
      }
      // blocks 34 and 43
      solver->mat[be][ar] =   hcoef * kcoef * dhdq[r][b];
      solver->mat[ar][be] = 2*hcoef * kcoef * dhdq_np1[r][b];
    }
    for(i = 0; i < 3; ++i){
      // blocks 14 and 41
      solver->mat[i][ar] =   hcoef * kcoef * dhdp[r][i];
      solver->mat[ar][i] = 2*hcoef * kcoef * dhdp_np1[r][i];

      // blocks 24 and 42
      solver->mat[i+3][ar] = hcoef * kcoef * dhdth[r][i];
      solver->mat[ar][i+3] = 2*hcoef * kcoef * dhdth_np1[r][i];

      // blocks 24 and 42 cont
/*      for(b = 0; b < numFlex; ++b){
      	solver->mat[i+3][ar] +=   hcoef * kcoef * mgs.q[b] * dflexhat[r][b][i];
      	solver->mat[ar][i+3] += 2*hcoef * kcoef * mgs.qnp1[b] * dflexhatnp1[r][b][i];
      }
*/
    }
    // block 44, reset to zero because factorization overwrites the entries
    for(s = 0; s < numConstr; ++s){
      solver->mat[ar][s + numRBM + numFlex] = 0.0;
    }
  }

#ifdef VICTORY
  fflush(stderr);
  solver->mat.print("from getStiffAndForce Solver matrix:");
  fflush(stderr);
  fprintf(stderr, "end getStiffAndForce ------------------------------\n");
#endif
  return  res.norm();

}

//------------------------------------------------------------------------------

void NLModalDescr::evalRHS(Vector &res, Vector &rhs, ModalGeomState &mgs){
/*PRE:
 POST:
*/
//************************************************
//  mgs.printRotation("rotation tensor from evalRHS");
//************************************************

  Vector &vel = mgs.vel;
  Vector &acc = mgs.acc;

  int b, i, r;
  for(i = 0; i < 3; ++i){

    // translational components of rhs; cont below
    rhs[i] = kcoef * ( res[i] + mass*acc[i] );

    // rotational components of rhs; cont below
    rhs[i+3] = kcoef * res[i+3];
    for(r = 0; r < 3; ++r){
      rhs[i+3] += kcoef *( J[i][r] + Bqq[i][r] + 2*Bhq[i][r] ) * acc[r+3];
    }
    for(b = 0; b < numFlex; ++b){
      rhs[i+3] += kcoef * ( Btq[b][i]*acc[numRBM+b]
        + 2*vel[numRBM+b]*(Bqo[b][i] + Bho[b][i])
        - mgs.q[b] * (modesPr[b][i] * fullFrot) );  // note to self: this line in question
    }
  }

  // rotational components continued; cont below
//  fprintf(stderr, "evalRHS c22b: %20.12e %20.12e %20.12e\n", c22b[0], c22b[1], c22b[2]);
  rhs[3] -= kcoef * ( c22b[1]*vel[5] - c22b[2]*vel[4] );
  rhs[4] -= kcoef * ( c22b[2]*vel[3] - c22b[0]*vel[5] );
  rhs[5] -= kcoef * ( c22b[0]*vel[4] - c22b[1]*vel[3] );

  int g, z;
  // deformational components of rhs; cont below
  for(b = 0; b < numFlex; ++b){

    z = numRBM + b;

    rhs[z] = kcoef * ( res[z] + mass*acc[z] );
//    next line for derivTest
//    rhs[z] = kcoef * mass * ( mgs.q[b]*freqs[b]*freqs[b] + acc[z] );
    for(i = 0; i < 3; ++i){
      rhs[z] += kcoef * ( Btq[b][i]*acc[i+3] - (Bqo[b][i] + Bho[b][i])*vel[i+3] );
    }
    for(g = 0; g < numFlex; ++g){
      rhs[z] -= kcoef * 2*Bto[g][b]*vel[numRBM+g];
    }
  }

  // constraint terms of rhs
  //double xdof, ydof, zdof;
  for(r = 0; r < numConstr; ++r){

    for(i = 0; i < 3; ++i){

      // translation cont
      rhs[i] -= hcoef * kcoef * dhdp[r][i] * mgs.lam[r];

      // rotation cont
      rhs[i+3] -= hcoef * kcoef * dhdth[r][i] * mgs.lam[r];
//      for(b = 0; b < numFlex; ++b)
//      	rhs[i+3] -= hcoef * kcoef * mgs.qnp1[b] * dflexhat[r][b][i] * mgs.lam[r];
    }
    // deformation cont
    for(b = 0; b < numFlex; ++b){
      rhs[b+numRBM] -= hcoef * kcoef * dhdq[r][b] * mgs.lam[r];
    }
    rhs[r + numRBM + numFlex] = -hcoef * kcoef * constr[r];
  }
}

//------------------------------------------------------------------------------

void NLModalDescr::formRHSpredictor(Vector &res, Vector &rhs,
  ModalGeomState &mgs){
/*PRE:
 POST:
*/
  calcTmps(mgs);
  evalRHS(res, rhs, mgs);

}

//------------------------------------------------------------------------------

double NLModalDescr::formRHScorrector(Vector& res, Vector& rhs,
  ModalGeomState& mgs){
/*PRE:
 POST:
*/
  evalRHS(res, rhs, mgs);
  return rhs.norm();
}

//------------------------------------------------------------------------------

int NLModalDescr::checkConvergence(int iter, double normRes, Vector &residual,
  Vector &dvec, double time){
/*PRE:
 POST:
*/
  if(iter == 0){ 
    firstRes = normRes;
#ifdef VICTORY_CONVERGE
    fprintf(stderr, "checkConvergence time: %e ----------------------------------\n", time);
#endif
//    firstRes = dvec.norm();
  }

  double relRes = (firstRes == 0) ? 0.0 : normRes/firstRes;
//  double relRes = (firstRes == 0) ? 0.0 : dvec.norm()/firstRes;

  int converged = 0;
  if(relRes <= tolerance)
    converged = 1;
  else if(relRes > 10000)   // check for divergence
    converged = -1;

#ifdef VICTORY_CONVERGE
  fprintf(stderr, "  Iteration # %d\n", iter);
  fprintf(stderr, "  normRes = %20.12e  relRes = %20.12e  firstRes = %20.12e  tol = %e\n",
    normRes, relRes, firstRes, tolerance);
  fprintf(stderr, "end checkConvergence ------------------------------------------\n");
#endif

  return converged;
}

//------------------------------------------------------------------------------

void NLModalDescr::dynamOutput(ModalGeomState* mgs, Vector& vel, Vector& vel_p,
  double time, int tIndex, Vector& extF, Vector &aeroF, Vector &acc){
/*PRE:
 POST:
 NOTE: valid only for displacements, leave for now
*/
  expandFlex(mgs->q, fullDsp);

  double x0, y0, z0, tmp[3], rottensor[3][3], result[3][3], theta[3];
  int xdof, ydof, zdof, iNode;

  for(iNode = 0; iNode < domain->numNode(); ++iNode){
    x0 = domain->getNodes().getNode(iNode).x;
    y0 = domain->getNodes().getNode(iNode).y;
    z0 = domain->getNodes().getNode(iNode).z;

    xdof = domain->getCDSA()->locate(iNode, DofSet::Xdisp);
    ydof = domain->getCDSA()->locate(iNode, DofSet::Ydisp);
    zdof = domain->getCDSA()->locate(iNode, DofSet::Zdisp);

    tmp[0] = x0 + fullDsp[xdof];
    tmp[1] = y0 + fullDsp[ydof];
    tmp[2] = z0 + fullDsp[zdof];

    fullDsp[xdof] = mgs->glR[0][0]*tmp[0] + mgs->glR[0][1]*tmp[1]
      + mgs->glR[0][2]*tmp[2] - x0
      + mgs->glT[0]*modesRB[0][xdof] + mgs->glT[1]*modesRB[1][xdof] + mgs->glT[2]*modesRB[2][xdof];

    fullDsp[ydof] = mgs->glR[1][0]*tmp[0] + mgs->glR[1][1]*tmp[1]
      + mgs->glR[1][2]*tmp[2] - y0
      + mgs->glT[0]*modesRB[0][ydof] + mgs->glT[1]*modesRB[1][ydof] + mgs->glT[2]*modesRB[2][ydof];

    fullDsp[zdof] = mgs->glR[2][0]*tmp[0] + mgs->glR[2][1]*tmp[1]
      + mgs->glR[2][2]*tmp[2] - z0
      + mgs->glT[0]*modesRB[0][zdof] + mgs->glT[1]*modesRB[1][zdof] + mgs->glT[2]*modesRB[2][zdof];

    xdof = domain->getCDSA()->locate(iNode, DofSet::Xrot);
    ydof = domain->getCDSA()->locate(iNode, DofSet::Yrot);
    zdof = domain->getCDSA()->locate(iNode, DofSet::Zrot);
    tmp[0] = fullDsp[xdof];  tmp[1] = fullDsp[ydof];  tmp[2] = fullDsp[zdof];

    vec_to_mat(tmp, rottensor);
    mat_mult_mat(mgs->glR, rottensor, result, 0);
    mat_to_vec(result, theta);

    fullDsp[xdof] = theta[0];
    fullDsp[ydof] = theta[1];
    fullDsp[zdof] = theta[2];
  }
  DynamMat dumDMat;
  domain->dynamOutput(tIndex+1, bcx, dumDMat, fullTmpF, fullAeroF,
    fullDsp, fullVel, fullAcc, fullPrevVel, vcx);

// testing
/*
  int i, j;
  printf(" --- step: %d from dynamOutput - glR: ---------\n", tIndex+1);
  for(i = 0; i < 3; ++i){
    for(j = 0; j < 3; ++j)
      printf("%20.12e", mgs->glR[i][j]);
    printf("\n");
  }
*/
// end testing

  SysState<Vector> state(mgs->q, mgs->vel, mgs->acc, mgs->qnp1);
  outputModal(state, extF, tIndex+1);
}

//------------------------------------------------------------------------------

void NLModalDescr::dRdTheta(double R[3][3], double dR[3][3][3]){
/*PRE: R[3][3] is the current rigid rotation tensor
 POST: return in dR the derivative dR/dtheta evaluated at the given rotation
       index order : dR[i][j]/dth[k] = dR[i][j][k]
*/
  int i, j;
  for(i = 0; i < 3; ++i){
    for(j = 0; j < 3; ++j){
//      dR[i][j][i] = 0.0;
      dR[j][i][i] = 0.0;
    }
/*
    dR[1][i][0] = -R[2][i];
    dR[2][i][0] =  R[1][i];

    dR[0][i][1] =  R[2][i];
    dR[2][i][1] = -R[0][i];

    dR[0][i][2] = -R[1][i];
    dR[1][i][2] =  R[0][i];
*/

    dR[i][1][0] =  R[i][2];
    dR[i][2][0] = -R[i][1];

    dR[i][0][1] = -R[i][2];
    dR[i][2][1] =  R[i][0];

    dR[i][0][2] =  R[i][1];
    dR[i][1][2] = -R[i][0];

  }
}

//------------------------------------------------------------------------------

void NLModalDescr::evalConstraintsAndDerivatives(ModalGeomState &mgs){
/*PRE:
 POST:
*/
  double dR[3][3][3], dRnp1[3][3][3], x0[3];
  double rVecNode[3], rVecTmp[3], rVecNode_np1[3];
  double Rtmp[3][3], glRRtmp[3][3], dRRtmp[3];
  int cdof, cnode;

  dRdTheta(mgs.glRnp1, dRnp1);
  dRdTheta(mgs.glR, dR);
  
  int g, b, i, j, k, m;
  for(g = 0; g < numConstr; ++g){

    cdof = domain->getDBC()[g].dofnum;
    cnode = domain->getDBC()[g].nnum;
    constr[g] = -domain->getDBC()[g].val;

    for(i = 0; i < 3; ++i){
      dhdp_np1[g][i] = modesRB[i][cDofIdx[g]];
      dhdp[g][i] = modesRB[i][cDofIdx[g]];
    }

    if(cdof > 2){	 // ie if gth constr does not correspond to a translation
      cdof -= 3;    
      rVecNode_np1[0] = rVecNode_np1[1] = rVecNode_np1[2] = 0.0;
      rVecNode[0] = rVecNode[1] = rVecNode[2] = 0.0;

      i = (cdof+2) % 3;
      j = (cdof+1) % 3;

      for(b = 0; b < numFlex; ++b){
        rVecTmp[0] = modesFl[b][domain->getCDSA()->locate(cnode, DofSet::Xrot)];
        rVecTmp[1] = modesFl[b][domain->getCDSA()->locate(cnode, DofSet::Yrot)];
        rVecTmp[2] = modesFl[b][domain->getCDSA()->locate(cnode, DofSet::Zrot)];
        vec_to_mat(rVecTmp, Rtmp);  // TODO store "magnitude" of rVecTmp for recovery later

        mat_mult_mat(mgs.glRnp1, Rtmp, glRRtmp, 0);
        mat_to_vec(glRRtmp, x0);
        dhdq_np1[g][b] = x0[cdof];  // TODO use "magnitude" of rVecTmp here; figure out how to use it

        mat_mult_mat(mgs.glR, Rtmp, glRRtmp, 0);
        mat_to_vec(glRRtmp, x0);
        dhdq[g][b] = x0[cdof];

        for(k = 0; k < 3; ++k){

          for(m = 0; m < 3; ++m)
            dRRtmp[m] = dR[m][0][k]*Rtmp[0][j] + dR[m][1][k]*Rtmp[1][j] + dR[m][2][k]*Rtmp[2][j];
          d2h[g][b][k] = glRRtmp[0][i]*dRRtmp[0] + glRRtmp[1][i]*dRRtmp[1] + glRRtmp[2][i]*dRRtmp[2];
/*
          d2h[g][b][k] = 0.0;
          for(i = 0; i < 3; ++i)
            for(j = 0; j < 3; ++j)
              d2h[g][b][k] += dR[i][j][k] * Rtmp[i][j];
*/
        }
        rVecNode_np1[0] += mgs.qnp1[b]*rVecTmp[0];  rVecNode[0] += mgs.q[b]*rVecTmp[0];
        rVecNode_np1[1] += mgs.qnp1[b]*rVecTmp[1];  rVecNode[1] += mgs.q[b]*rVecTmp[1];
        rVecNode_np1[2] += mgs.qnp1[b]*rVecTmp[2];  rVecNode[2] += mgs.q[b]*rVecTmp[2];
      }

      vec_to_mat(rVecNode_np1, Rtmp);
      mat_mult_mat(mgs.glRnp1, Rtmp, glRRtmp, 0);
      mat_to_vec(glRRtmp, x0);
      constr[g] += x0[cdof];

      for(k = 0; k < 3; ++k){
        for(m = 0; m < 3; ++m)
          dRRtmp[m] = dRnp1[m][0][k]*Rtmp[0][j] + dRnp1[m][1][k]*Rtmp[1][j] + dRnp1[m][2][k]*Rtmp[2][j];
        dhdth_np1[g][k] = glRRtmp[0][i]*dRRtmp[0] + glRRtmp[1][i]*dRRtmp[1] + glRRtmp[2][i]*dRRtmp[2];
      }

      vec_to_mat(rVecNode, Rtmp);
      mat_mult_mat(mgs.glR, Rtmp, glRRtmp, 0);
        
      for(k = 0; k < 3; ++k){
        for(m = 0; m < 3; ++m)
          dRRtmp[m] = dR[m][0][k]*Rtmp[0][j] + dR[m][1][k]*Rtmp[1][j] + dR[m][2][k]*Rtmp[2][j];
        dhdth[g][k] = glRRtmp[0][i]*dRRtmp[0] + glRRtmp[1][i]*dRRtmp[1] + glRRtmp[2][i]*dRRtmp[2];
      }
    }
    else{
      for(i = 0; i < 3; ++i)
        constr[g] += mgs.glTnp1[i] * modesRB[i][cDofIdx[g]];

      x0[0] = domain->getNodes().getNode(cnode).x;
      x0[1] = domain->getNodes().getNode(cnode).y;
      x0[2] = domain->getNodes().getNode(cnode).z;

      constr[g] += mgs.glRnp1[cdof][0] * x0[0] + mgs.glRnp1[cdof][1] * x0[1]
      	+ mgs.glRnp1[cdof][2] * x0[2] - x0[cdof];

      for(i = 0; i < 3; ++i){
        dhdth_np1[g][i] = dRnp1[cdof][0][i] * x0[0] + dRnp1[cdof][1][i] * x0[1] + dRnp1[cdof][2][i] * x0[2];
        dhdth[g][i] = dR[cdof][0][i] * x0[0] + dR[cdof][1][i] * x0[1] + dR[cdof][2][i] * x0[2];
      }

      for(b = 0; b < numFlex; ++b){
        x0[0] = modesFl[b][domain->getCDSA()->locate(cnode, DofSet::Xdisp)];
        x0[1] = modesFl[b][domain->getCDSA()->locate(cnode, DofSet::Ydisp)];
        x0[2] = modesFl[b][domain->getCDSA()->locate(cnode, DofSet::Zdisp)];

      	dhdq_np1[g][b] = mgs.glRnp1[cdof][0] * x0[0] + mgs.glRnp1[cdof][1] * x0[1] + mgs.glRnp1[cdof][2] * x0[2];
      	constr[g] += mgs.qnp1[b] * dhdq_np1[g][b];

        dhdq[g][b] = mgs.glR[cdof][0] * x0[0] + mgs.glR[cdof][1] * x0[1] + mgs.glR[cdof][2] * x0[2];

      	for(i = 0; i < 3; ++i){
          dhdth_np1[g][i] += mgs.qnp1[b] * (dRnp1[cdof][0][i] * x0[0] + dRnp1[cdof][1][i] * x0[1] + dRnp1[cdof][2][i] * x0[2]);
          d2h[g][b][i] = dR[cdof][0][i] * x0[0] + dR[cdof][1][i] * x0[1] + dR[cdof][2][i] * x0[2];
          dhdth[g][i] += mgs.q[b] * d2h[g][b][i];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

void NLModalDescr::test2(ModalGeomState *mgsref){
/*PRE:
 POST:
*/
  if(!mgsref)
    mgsref = new ModalGeomState(numRBM, numFlex, numConstr);

  const double pi = 4*atan(1.0);
  const double eps = 1.0e-6;

  double thref[3] = {0, 0, pi/6.}, thp[3], thm[3];
  vec_to_mat(thref, mgsref->glR);
  thp[0] = thm[0] = thref[0];
  thp[1] = thm[1] = thref[1];
  thp[2] = thm[2] = thref[2];
/*
  mgsref->vel[4] = 0.8;
  mgsref->vel[5] = 1.0;
  mgsref->vel[numRBM] = 0.1;
  mgsref->vel[numRBM+3] = 0.08;
*/  
  double glRp[3][3], glRm[3][3];
  int i, j, k, r;

  FullM matdiff;
  matdiff.setNewSize(solVecInfo(), solVecInfo(), 0.0);
  Vector restmp(solVecInfo(), 0.0), rhsp(solVecInfo()), rhsm(solVecInfo());

  mcoef = 0.0; ccoef = 0.0; kcoef = 1.0;
  getStiffAndForce(*mgsref, restmp);

  fflush(stderr);
  fprintf(stderr, "moment of inertia:\n");
  for(i = 0; i < 3; ++i){
    for(j = 0; j < 3; ++j)
      fprintf(stderr, "%18.10e", J[i][j]);
    fprintf(stderr, "\n");
  }
  fflush(stderr);
  solver->mat.print("Analytic");

  ModalGeomState mgsp(*mgsref);
  ModalGeomState mgsm(*mgsref);


  mcoef = 0.0; ccoef = 0.0; kcoef = 1.0;
  for(k = 0; k < 3; ++k){
    mgsp.glT[k] += eps;
    mgsm.glT[k] -= eps;

    calcTmps(mgsp);
    evalRHS(restmp, rhsp, mgsp);

    calcTmps(mgsm);
    evalRHS(restmp, rhsm, mgsm);

    for(r = 0; r < solVecInfo(); ++r)
      matdiff[r][k] = (-rhsp[r] + rhsm[r]) / (2*eps);

    mgsp.glT[k] = mgsref->glT[k];
    mgsm.glT[k] = mgsref->glT[k];
  }

  for(k = 0; k < 3; ++k){
    for(j = 0; j < 3; ++j){
      thp[j] =  eps * mgsref->glR[j][k];
      thm[j] = -eps * mgsref->glR[j][k];
    }
    vec_to_mat(thp, glRp);
    vec_to_mat(thm, glRm);

    for(j = 0; j < 3; ++j){
      for(r = 0; r < 3; ++r){
      	mgsp.glR[j][r] = glRp[j][0] * mgsref->glR[0][r] + glRp[j][1] * mgsref->glR[1][r] + glRp[j][2] * mgsref->glR[2][r];
      	mgsm.glR[j][r] = glRm[j][0] * mgsref->glR[0][r] + glRm[j][1] * mgsref->glR[1][r] + glRm[j][2] * mgsref->glR[2][r];
      }
    }
    calcTmps(mgsp);
    evalRHS(restmp, rhsp, mgsp);
    calcTmps(mgsm);
    evalRHS(restmp, rhsm, mgsm);

    for(r = 0; r < solVecInfo(); ++r)
      matdiff[r][k+3] = (-rhsp[r] + rhsm[r]) / (2*eps);
  }

  thp[0] = thm[0] = thref[0];
  thp[1] = thm[1] = thref[1];
  thp[2] = thm[2] = thref[2];
  vec_to_mat(thp, mgsp.glR);
  vec_to_mat(thm, mgsm.glR);

  for(k = 0; k < numFlex; ++k){

    mgsp.q[k] += eps;
    calcTmps(mgsp);
    evalRHS(restmp, rhsp, mgsp);

    mgsm.q[k] -= eps;
    calcTmps(mgsm);
    evalRHS(restmp, rhsm, mgsm);

    for(r = 0; r < solVecInfo(); ++r)
      matdiff[r][k+numRBM] = (-rhsp[r] + rhsm[r]) / (2*eps);

    mgsp.q[k] = mgsref->q[k];
    mgsm.q[k] = mgsref->q[k];
  }

  for(k = 0; k < numConstr; ++k){

    mgsp.lam[k] += eps;
    calcTmps(mgsp);
    evalRHS(restmp, rhsp, mgsp);

    mgsm.lam[k] -= eps;
    calcTmps(mgsm);
    evalRHS(restmp, rhsm, mgsm);
    for(r = 0; r < solVecInfo(); ++r)
      matdiff[r][k+numRBM+numFlex] = (-rhsp[r] + rhsm[r]) / (2*eps);

    mgsp.lam[k] = mgsm.lam[k] = mgsref->lam[k];

  }

  matdiff.print("Finite difference");
  
}

//------------------------------------------------------------------------------

void NLModalDescr::test(ModalGeomState *mgsref){

  if(!mgsref)
    mgsref = new ModalGeomState(numRBM, numFlex, numConstr);

  const double pi = 4*atan(1.0);
  const double eps = 1.0e-6;

  int i, j, k, r;

  double thref[3] = {1.6*pi, -pi/3.7, pi/6}, thp[3], thm[3];
  double glRp[3][3], glRm[3][3];

  fflush(stderr);
  fprintf(stderr, "moment of inertia:\n");
  for(i = 0; i < 3; ++i){
    for(j = 0; j < 3; ++j)
      fprintf(stderr, "%18.10e", J[i][j]);
    fprintf(stderr, "\n");
  }
  
  vec_to_mat(thref, mgsref->glR);
  fflush(stderr);
  fprintf(stderr, "thref: %18.10e%18.10e%18.10e\n", thref[0], thref[1], thref[2]);
  fflush(stderr);
  mgsref->printRotation("ref rotation:");
  fflush(stderr);
/*
  mgsref->glT[0] = 1.0;
  mgsref->glT[1] = 2.0;
  mgsref->glT[2] = 3.0;
*/
  mgsref->q[0] = 1.0;
//  mgsref->q[1] = 1.1;
//  mgsref->q[2] = 1.2;

  FullM matdiff, dhdiff, dhanalytic;
  matdiff.setNewSize(solVecInfo(), solVecInfo(), 0.0);
  dhdiff.setNewSize(numConstr, solVecInfo(), 0.0);
  dhanalytic.setNewSize(numConstr, solVecInfo(), 0.0);

  Vector restmp(solVecInfo(), 0.0), rhsp(solVecInfo()), rhsm(solVecInfo());
  calcTmps(*mgsref);

  for(k = 0; k < numConstr; ++k){
    for(i = 0; i < 3; ++i){
      dhanalytic[k][i] = modesRB[i][cDofIdx[k]];
      dhanalytic[k][i+3] = dxhat[k][i];
      for(j = 0; j < numFlex; ++j)
        dhanalytic[k][i+3] += mgsref->q[j] * dflexhat[k][j][i];
    }
    for(j = 0; j < numFlex; ++j){
      dhanalytic[k][numRBM+j] = flexhat[k][j];

//      if(domain->getDBC()[k].dofnum > 2)
//        dhanalytic[k][numRBM+j] += 1.0;
    }
  }

/*
  mcoef = 1.0; ccoef = 0.0; kcoef = 0.0;
  getStiffAndForce(*mgsref, restmp);

  solver->mat.print("Analytic");
  fflush(stderr);
*/
  ModalGeomState mgsp(*mgsref);
  ModalGeomState mgsm(*mgsref);

  kcoef = 1.0;

  for(k = 0; k < 3; ++k){
    for(j = 0; j < 3; ++j){
      mgsp.glT[j] = mgsref->glT[j];
      mgsm.glT[j] = mgsref->glT[j];
    }

    mgsp.glT[k] += eps;
    calcTmps(mgsp);
    for(i = 0; i < numConstr; ++i){ dhdiff[i][k] = constr[i]/(2*eps); }

    mgsm.glT[k] -= eps;
    calcTmps(mgsm);
    for(i = 0; i < numConstr; ++i){ dhdiff[i][k] -= constr[i]/(2*eps); }

  }

  for(i = 0; i < 3; ++i){
    mgsp.glT[i] = mgsref->glT[i];
    mgsm.glT[i] = mgsref->glT[i];
  }

  for(k = 0; k < 3; ++k){
    for(j = 0; j < 3; ++j){
      thp[j] =  eps * mgsref->glR[j][k];
      thm[j] = -eps * mgsref->glR[j][k];
    }
    vec_to_mat(thp, glRp);
    vec_to_mat(thm, glRm);

    for(j = 0; j < 3; ++j){
      for(r = 0; r < 3; ++r){
      	mgsp.glR[j][r] = glRp[j][0] * mgsref->glR[0][r] + glRp[j][1] * mgsref->glR[1][r] + glRp[j][2] * mgsref->glR[2][r];
      	mgsm.glR[j][r] = glRm[j][0] * mgsref->glR[0][r] + glRm[j][1] * mgsref->glR[1][r] + glRm[j][2] * mgsref->glR[2][r];
      }
    }

    calcTmps(mgsp);
    for(i = 0; i < numConstr; ++i){ dhdiff[i][3+k] = constr[i]/(2*eps); }
    calcTmps(mgsm);
    for(i = 0; i < numConstr; ++i){ dhdiff[i][3+k] -= constr[i]/(2*eps); }

  }

  for(i = 0; i < 3; ++i){
    thp[i] = thref[i];
    thm[i] = thref[i];
  }
  vec_to_mat(thref, mgsp.glR);
  vec_to_mat(thref, mgsm.glR);
  for(k = 0; k < numFlex; ++k){

    mgsp.q = mgsref->q;
    mgsm.q = mgsref->q;

    mgsp.q[k] += eps;
    calcTmps(mgsp);
    for(i = 0; i < numConstr; ++i){ dhdiff[i][numRBM+k] = constr[i]/(2*eps); }

    mgsm.q[k] -= eps;
    calcTmps(mgsm);
    for(i = 0; i < numConstr; ++i){ dhdiff[i][numRBM+k] -= constr[i]/(2*eps); }

  }

//  matdiff.print("\nfinite diff");
  dhanalytic.print("analytic derivative of the constraints");
  dhdiff.print("------------\nfinite difference");
      

}

//------------------------------------------------------------------------------

void NLModalDescr::printCoefs(){

  int b, g, i;
/*
  int j;
  fprintf(stderr, "Bp ------------------------------\n");
  for(b = 0; b < numFlex; ++b){
    for(i = 0; i < 3; ++i){
      for(g = 0; g < numFlex; ++g){
        for(j = 0; j < 3; ++j){
          fprintf(stderr, "%16.8e", Bp[b][g][i][j]);
        }
        fprintf(stderr, "    ");
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }


  fprintf(stderr, "Bh ------------------------------\n");
  for(b = 0; b < numFlex; ++b){
    for(i = 0; i < 3; ++i){
      for(j = 0; j < 3; ++j){
        fprintf(stderr, "%16.8e", Bh[b][i][j]);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "Bt ------------------------------\n");
  for(b = 0; b < numFlex; ++b){
    for(i = 0; i < numFlex; ++i){
      for(j = 0; j < 3; ++j){
        fprintf(stderr, "%16.8e", Bt[b][i][j]);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
*/
  fprintf(stderr, "dflexhat ----------------------\n");
  for(b = 0; b < numFlex; ++b){
    for(g = 0; g < numFlex; ++g){
      for(i = 0 ; i < 3; ++i){
        fprintf(stderr, "%16.8e", dflexhat[b][g][i]);
      }
      fprintf(stderr, "   ");
    }
    fprintf(stderr, "\n");
  }


}
