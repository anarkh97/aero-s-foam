#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

// New include files for Restart file
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#include <Corotational.d/Corotator.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/mathUtility.h>
#include <Math.d/matrix.h>
#include <Control.d/ControlInterface.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>

#include <Driver.d/GeoSource.h>
#include <Corotational.d/MatNLCorotator.h>
#ifdef USE_EIGEN3
#include <Element.d/Dimass.d/InertialForceFunction.h>
#endif
#include <algorithm>

void
Domain::getElemStiffAndForce(const GeomState &geomState, double time,
                             const GeomState *refState, const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
    const_cast<Corotator &>(elemCorot).getStiffAndForce(
        const_cast<GeomState *>(refState),
        const_cast<GeomState &>(geomState),
        nodes, elemStiff, elemForce, domain->solInfo().getTimeStep(), time);
}

void
Domain::getElemStiffAndForce(const GeomState &geomState, double time,
                             const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
  const_cast<Corotator &>(elemCorot).getStiffAndForce(
      const_cast<GeomState &>(geomState),
      nodes, elemStiff, elemForce, domain->solInfo().getTimeStep(), time);
}

void
Domain::getStiffAndForce(GeomState &geomState, Vector& elementForce,
                         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time,
                         GeomState *refState, Vector *reactions, FullSquareMatrix *mel)
/*******************************************************************
 *
 * Purpose :
 *
 *  Compute element tangential stiffness and element internal force
 *  and assemble element internal force into global internal force.
 *  Also compute follower external force contribution to both 
 *  tangential stiffness matrix and residual
 *
 * Input :
 *
 *  corotators : element corotator array
 *  geomState  : current node geometric state
 *  nodes      : undeformed nodal coordinate set
 *
 * Output :
 *
 *  residual   : residual vector = external force - internal force
 *  kel        : array of element tangential stiffness matrices
 *               in current configuration
 *
 *****************************************************************/

{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // mpc needs lambda for nonlinear statics

  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    // Get updated tangent stiffness matrix and element internal force
    if (const Corotator *elemCorot = corotators[iele]) {
      getElemStiffAndForce(geomState, pseudoTime, refState, *elemCorot, elementForce.data(), kel[iele]);
      if(domain->solInfo().galerkinPodRom && packedEset[iele]->hasRot()) {
        transformElemStiffAndForce(geomState, elementForce.data(), kel[iele], iele, true);
      }
    }
    // Compute k and internal force for an element with x translation (or temperature) dofs
    else if(solInfo().soltyp == 2) {
      kel[iele].zero();
      Vector temp(packedEset[iele]->numNodes());
      for(int i=0; i<packedEset[iele]->numNodes(); ++i) {
        temp[i] = geomState[packedEset[iele]->nodes()[i]].x;
      }
      kel[iele] = packedEset[iele]->stiffness(nodes, kel[iele].data());
      kel[iele].multiply(temp, elementForce, 1.0); // elementForce = kel*temp
    }
    // Assemble element internal force into residual force vector
    for(int idof = 0; idof < kel[iele].dim(); ++idof) {
      int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(uDofNum >= 0)
        residual[uDofNum] -= elementForce[idof];
      else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] += elementForce[idof];
      }
    }
  }

  getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, refState, reactions, true);

  if(sinfo.isDynam() && mel) getFictitiousForce(geomState, kel, residual, time, refState, reactions, mel, true);

  if(!solInfo().getNLInfo().unsymmetric && solInfo().newmarkBeta != 0)
    for(int iele = 0; iele < numele; ++iele)
      kel[iele].symmetrize();
}

void
Domain::getFollowerForce(GeomState &geomState, Vector& elementForce,
                         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time,
                         GeomState *refState, Vector *reactions, bool compute_tangents)
{
  if(domain->pressureFlag()) {
    double cflg = (sinfo.newmarkBeta == 0.0) ? 0.0 : 1.0;
    double loadFactor = (domain->mftval && sinfo.isDynam()) ? lambda*domain->mftval->getVal(std::max(time,0.0)) : lambda;
    double p0;
    for(int iele = 0; iele < numele;  ++iele) {
      // If there is a zero pressure defined, skip the element
      if((p0 = packedEset[iele]->getPressure()) == 0) continue;

      // Compute (linear) element pressure force in the local coordinates
      elementForce.zero();
      packedEset[iele]->setPressure(p0*loadFactor, domain->getMFTT(), sinfo.ConwepOnOff);
      packedEset[iele]->computePressureForce(nodes, elementForce, &geomState, 1, time);
      packedEset[iele]->setPressure(p0, domain->getMFTT(), sinfo.ConwepOnOff);

      // Include the "load stiffness matrix" in kel[iele]
      if(compute_tangents) {
        FullSquareMatrix elementLoadStiffnessMatrix(kel[iele].dim());
        elementLoadStiffnessMatrix.zero();
        corotators[iele]->getDExternalForceDu(geomState, nodes, elementLoadStiffnessMatrix,
                                              elementForce.data());
        for(int i=0; i<kel[iele].dim(); ++i)
          for(int j=0; j<kel[iele].dim(); ++j)
            kel[iele][i][j] += elementLoadStiffnessMatrix[i][j];
      }

      // Determine the elemental force for the corrotated system
      corotators[iele]->getExternalForce(geomState, nodes, elementForce.data());

      // Assemble element pressure forces into residual force vector
      for(int idof = 0; idof < packedEset[iele]->numDofs(); ++idof) {
        int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
        if(uDofNum >= 0)
          residual[uDofNum] += elementForce[idof];
        else if(reactions) {
          int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
          if(cDofNum >= 0)
            (*reactions)[cDofNum] -= elementForce[idof];
        }
      }
    }
  }

  // pressure using surfacetopo
  int* edofs = (int*) dbg_alloca(maxNumDOFs*sizeof(int));
  double mfttFactor = (domain->mftval && sinfo.isDynam()) ? domain->mftval->getVal(std::max(time,0.0)) : 1.0;
  SubDomain *subCast = (numNeum > 0) ? dynamic_cast<SubDomain*>(this) : NULL;
  for(int iele = 0; iele < numNeum; ++iele) {
    neum[iele]->dofs(*dsa, edofs);
    elementForce.zero();
    neum[iele]->neumVector(nodes, elementForce, 0, &geomState);

    // Include the "load stiffness matrix" in kel[iele]
    if(compute_tangents) {
      int jele = (subCast) ? subCast->globalToLocalElem(neum[iele]->getAdjElementIndex()) : neum[iele]->getAdjElementIndex();
      if(jele > -1) {
        FullSquareMatrix elementLoadStiffnessMatrix(neum[iele]->numDofs());
        elementLoadStiffnessMatrix.zero();
        neum[iele]->neumVectorJacobian(nodes, elementLoadStiffnessMatrix, 0, &geomState);
        int *eledofs = new int[neum[iele]->numDofs()];
        for(int j = 0; j < neum[iele]->numDofs(); ++j) {
          for(int k = 0; k < allDOFs->num(jele); ++k)
            if(edofs[j] == (*allDOFs)[jele][k]) { eledofs[j] = k; break; }
        }
        for(int i=0; i<neum[iele]->numDofs(); ++i)
          for(int j=0; j<neum[iele]->numDofs(); ++j)
            kel[jele][eledofs[i]][eledofs[j]] -= lambda*mfttFactor*elementLoadStiffnessMatrix[i][j];
        delete [] eledofs;
      }
    }

    for(int idof = 0; idof < neum[iele]->numDofs(); ++idof) {
      int uDofNum = c_dsa->getRCN(edofs[idof]);
      if(uDofNum >= 0)
        residual[uDofNum] += lambda*mfttFactor*elementForce[idof];
      else if(reactions) {
        int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
        if(cDofNum >= 0)
          (*reactions)[cDofNum] -= lambda*mfttFactor*elementForce[idof];
      }
    }
  }

  // treatment of nodal moments
  for(int i = 0; i < numNeuman; ++i) {
    if((nbc[i].type == BCond::Forces || nbc[i].type == BCond::Usdf || nbc[i].type == BCond::Actuators) 
       && (nbc[i].dofnum == 3 || nbc[i].dofnum == 4 || nbc[i].dofnum == 5)) {
      int dofs[3];
      dsa->number(nbc[i].nnum, DofSet::XYZrot, dofs);
      double m0[3] = { 0, 0, 0 }, m[3], r[3], rotvar[3][3];
      m0[nbc[i].dofnum-3] = lambda*mfttFactor*nbc[i].val;

      switch(nbc[i].mtype) {
        case BCond::Axial : // axial (constant) moment: m = m0
          for(int j=0; j<3; ++j) m[j] = m0[j];
          break;
        case BCond::Rotational : { // rotational moment: m = T^{-1}*m0
          mat_to_vec(geomState[nbc[i].nnum].R,r);
          pseudorot_var(r, rotvar);
          mat_mult_vec(rotvar,m0,m,1);
        } break;
        case BCond::Follower : // follower moment: m = R*m0
          mat_mult_vec(geomState[nbc[i].nnum].R,m0,m,0);
          break;
        default :
          std::cerr << " *** WARNING: selected moment type is not supported\n";
      }
      for(int j = 0; j < 3; ++j) {
        int uDofNum = c_dsa->getRCN(dofs[j]);
        if(uDofNum >= 0)
          residual[uDofNum] += m[j];
        else if(reactions) {
          int cDofNum = c_dsa->invRCN(dofs[j]);
          if(cDofNum >= 0)
            (*reactions)[cDofNum] -= m[j];
        }
      }
      // tangent stiffness contribution: 
      if(compute_tangents) {
        switch(nbc[i].mtype) {
          case BCond::Axial : { // axial (constant) moment
            double skewm0[3][3] = { {     0, -m0[2],  m0[1] },
                                    {  m0[2],     0, -m0[0] },
                                    { -m0[1],  m0[0],    0  } };
            for(int j=0; j<3; ++j)
              for(int k=0; k<3; ++k) 
                rotvar[j][k] = 0.5*skewm0[j][k];
          } break;
          case BCond::Rotational : { // rotational moment
            double scndvar[3][3];
            pseudorot_2var(r, m0, scndvar);
            for(int j=0; j<3; ++j)
              for(int k=0; k<3; ++k)
                rotvar[j][k] = 0.5*(scndvar[j][k] + scndvar[k][j]);
          } break;
          case BCond::Follower : { // follower moment
            double skewm[3][3] = { {     0, -m[2],  m[1] },
                                   {  m[2],     0, -m[0] },
                                   { -m[1],  m[0],    0  } };
            for(int j=0; j<3; ++j)
              for(int k=0; k<3; ++k)
                rotvar[j][k] = -0.5*skewm[j][k];
          } break;
        }
        for(int inode = 0; inode < nodeToElem->num(nbc[i].nnum); ++inode) { // loop over the elements attached to the node
                                                                            // at which the nodal moment is applied
          int iele = (*nodeToElem)[nbc[i].nnum][inode];
          int eledofs[3] = { -1, -1, -1 };
          for(int j = 0; j < 3; ++j) {
            for(int k = 0; k < allDOFs->num(iele); ++k)
              if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
          }
          if(eledofs[0] != -1 && eledofs[1] != -1 && eledofs[2] != -1) {
            // found an element with the 3 rotation dofs of node nbc[i].nnum so we can add the load stiffness
            // contribution of the nodal moment to the tangent stiffness matrix of this element
            for(int j = 0; j < 3; ++j)
              for(int k = 0; k < 3; ++k)
                kel[iele][eledofs[j]][eledofs[k]] -= rotvar[j][k];
            break;
          }
        }
      }
    }
  }

  // treatment of nodal follower forces
  for(int i = 0; i < numNeuman; ++i) {
    if((nbc[i].type == BCond::Forces || nbc[i].type == BCond::Usdf || nbc[i].type == BCond::Actuators) 
       && (nbc[i].dofnum == 0 || nbc[i].dofnum == 1 || nbc[i].dofnum == 2)
       && nbc[i].mtype == BCond::Follower) {
      int dofs[6];
      dsa->number(nbc[i].nnum, DofSet::XYZdisp | DofSet::XYZrot, dofs);
      double f0[3] = { 0, 0, 0 }, f[3] = { 0, 0, 0 }, r[3], rotvar[3][3];
      f0[nbc[i].dofnum] = lambda*mfttFactor*nbc[i].val;
      mat_mult_vec(geomState[nbc[i].nnum].R,f0,f,0); // f = R*f0
      for(int j = 0; j < 3; ++j) {
        int uDofNum = c_dsa->getRCN(dofs[j]);
        if(uDofNum >= 0)
          residual[uDofNum] += f[j];
        else if(reactions) {
          int cDofNum = c_dsa->invRCN(dofs[j]);
          if(cDofNum >= 0)
            (*reactions)[cDofNum] -= f[j];
        }
      }
      // tangent stiffness contribution: 
      if(compute_tangents) {
        double skewf[3][3] = { {     0, -f[2],  f[1] },
                               {  f[2],     0, -f[0] },
                               { -f[1],  f[0],    0  } };

        for(int inode = 0; inode < nodeToElem->num(nbc[i].nnum); ++inode) { // loop over the elements attached to the node
                                                                            // at which the nodal moment is applied
          int iele = (*nodeToElem)[nbc[i].nnum][inode];
          int eledofs[6] = { -1, -1, -1, -1, -1, -1 };
          for(int j = 0; j < 6; ++j) {
            for(int k = 0; k < allDOFs->num(iele); ++k)
              if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
          }
          if(eledofs[0] != -1 && eledofs[1] != -1 && eledofs[2] != -1 &&
             eledofs[3] != -1 && eledofs[4] != -1 && eledofs[5] != -1) {
            // found an element with the 6 translation & rotation dofs of node nbc[i].nnum so we can add the load stiffness
            // contribution of the nodal force to the tangent stiffness matrix of this element
            for(int j = 0; j < 3; ++j)
              for(int k = 0; k < 3; ++k)
                kel[iele][eledofs[j]][eledofs[3+k]] -= -skewf[j][k];
            break;
          }
        }
      }
    }
  }

  if(domain->thermalFlag()) {
    if(!temprcvd) initNodalTemperatures(); // XXXX to be moved
    Vector elementTemp(maxNumNodes);
    for(int iele = 0; iele < numele;  ++iele) {
      // By convention phantom elements do not have thermal load
      if(packedEset[iele]->getProperty() == 0) continue;

      // Extract the element nodal temperatures from temprcvd and/or element property ambient temperature
      for(int inod = 0; inod < elemToNode->num(iele); ++inod) {
        double t = temprcvd[(*elemToNode)[iele][inod]];
        elementTemp[inod] = (t == defaultTemp) ? packedEset[iele]->getProperty()->Ta : t;
      }

      // Compute element thermal force in the local coordinates
      elementForce.zero();
      packedEset[iele]->getThermalForce(nodes, elementTemp, elementForce, 1);
/*
      elementForce *= lambda;

      // Include the "load stiffness matrix" in kel[iele]
      if(compute_tangents)
        corotators[iele]->getDExternalForceDu(geomState, nodes, kel[iele],
                                              elementForce.data());
*/
      // Determine the elemental force for the corrotated system
      corotators[iele]->getExternalForce(geomState, nodes, elementForce.data());

      // Assemble element thermal forces into residual force vector
      for(int idof = 0; idof < packedEset[iele]->numDofs(); ++idof) {
        int uDofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
        if(uDofNum >= 0)
          residual[uDofNum] += elementForce[idof];
        else if(reactions) {
          int cDofNum = c_dsa->invRCN((*allDOFs)[iele][idof]);
          if(cDofNum >= 0)
            (*reactions)[cDofNum] -= elementForce[idof];
        }
      }
    }
  }
}

void
Domain::getWeightedStiffAndForceOnly(const std::map<int, double> &weights,
                                     GeomState &geomState, Vector& elementForce,
                                     Corotator **corotators, FullSquareMatrix *kel,
                                     Vector &residual, double lambda, double time,
                                     GeomState *refState, FullSquareMatrix *mel)
{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // MPC needs lambda for nonlinear statics
  
  for (std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iElem = it->first;

    // Get updated tangent stiffness matrix and element internal force
    if (const Corotator *elementCorot = corotators[iElem]) {
      elementForce.zero();

      FullSquareMatrix &elementStiff = kel[iElem];
      getElemStiffAndForce(geomState, pseudoTime, refState, *elementCorot, elementForce.data(), elementStiff);
      if(domain->solInfo().galerkinPodRom && packedEset[iElem]->hasRot()) {
        transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iElem, true);
      }
   
      // Apply lumping weight 
      const double lumpingWeight = it->second;
      elementForce *= lumpingWeight;
      elementStiff *= lumpingWeight;

      const int elemDofCount = elementStiff.dim();
      for(int iDof = 0; iDof < elemDofCount; ++iDof) {
        const int dofId = c_dsa->getRCN((*allDOFs)[iElem][iDof]);
        if (dofId >= 0) {
          residual[dofId] -= elementForce[iDof];
        }
      }
    }
  }

  getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, refState, NULL, true);

  if(sinfo.isDynam() && mel) getFictitiousForce(geomState, kel, residual, time, refState, NULL, mel, true);
}

void
Domain::applyResidualCorrection(GeomState &geomState, Corotator **corotators, Vector &residual, double rcoef)
{
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) {
      Vector residualCorrection(packedEset[iele]->numDofs());
      residualCorrection.zero();
      corotators[iele]->getResidualCorrection(geomState, residualCorrection.data());
    
      // Assemble element residualCorrection into global residual vector
      for(int idof = 0; idof < packedEset[iele]->numDofs(); ++idof) {
        int dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
        if(dofNum >= 0)
          residual[dofNum] += rcoef*residualCorrection[idof];
      }
    }
  }
}

void
Domain::updateStates(GeomState *refState, GeomState &geomState, Corotator **corotators)
{
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) corotators[iele]->updateStates(refState, geomState, nodes);
  }
}

void
Domain::initializeParameters(GeomState &geomState, Corotator **corotators)
{
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) corotators[iele]->initMultipliers(geomState);
  }

  SPropContainer &sProps = geoSource->getStructProps();
  SPropContainer::iterator it = sProps.begin();
  while(it != sProps.end()) {
    StructProp* p = &(it->second);
    p->penalty = p->initialPenalty;
    it++;
  }
  if(p) p->penalty = sinfo.penalty; 
}

void
Domain::updateParameters(GeomState &geomState, Corotator **corotators)
{
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) corotators[iele]->updateMultipliers(geomState);
  }

  SPropContainer &sProps = geoSource->getStructProps();
  SPropContainer::iterator it = sProps.begin();
  while(it != sProps.end()) {
    StructProp* p = &(it->second);
    p->penalty *= sinfo.penalty_beta;
    it++;
  }
  if(p) p->penalty *= sinfo.penalty_beta;
  // TODO this doesn't allow for elements that have a property that is not in the sProps container
  // or StructProps that don't use the global default penalty parameter from SolverInfo
}

double
Domain::getError(Corotator **corotators)
{
  // returns the largest constraint violation
  double err = 0;
  for(int iele = 0; iele < numele; ++iele) {
    if(corotators[iele]) err = std::max(err,corotators[iele]->getError());
  }
  return err;
}

// used in nonlinear statics
void
Domain::createKelArray(FullSquareMatrix *&kArray)
{
  // Allocate array of pointers to FullSquareMatrix to store the
  // element stiffness matrices
  kArray = new FullSquareMatrix[numele];

  // Allocate the correct size for each elements stiffness matrix
  int iele;
  for(iele = 0; iele<numele; ++iele) {
    kArray[iele].setSize(packedEset[iele]->numDofs());
    kArray[iele].zero(); 
  }
}

// used in nonlinear dynamics

void
Domain::createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray)
{

 // Allocate array of pointers to FullSquareMatrix to store
 // the element stiffness matrices and element mass matrices
 kArray = new FullSquareMatrix[numele];
 mArray = new FullSquareMatrix[numele];

 // Allocate the correct size for each element's stiffness & mass matrix
 int iele;
 for(iele = 0; iele<numele; ++iele) {
   int dimension = packedEset[iele]->numDofs();
   kArray[iele].setSize(dimension);
   kArray[iele].zero();
   mArray[iele].setSize(dimension);
 }

 // Form and store element mass matrices into an array
 for(iele=0; iele<numele; ++iele) {
   // note: only lumped mass matrix is supported currently for elements with rotation dofs in nonlinear dynamics
   //       (only the euler beam element is affected)
   double mratio = (packedEset[iele]->hasRot() && sinfo.isNonLin() && sinfo.isDynam()) ? 0 : geoSource->getMRatio();
   mArray[iele].copy(packedEset[iele]->massMatrix(nodes, mArray[iele].data(), mratio));
 }

 // zero rotational degrees of freedom within element mass matrices
 // for nonlinear implicit dynamics
 if(sinfo.zeroRot && sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0) {
   int *dofType = dsa->makeDofTypeArray();

   int i,j;
   for(iele=0; iele<numele; ++iele) {
     for(i=0; i<mArray[iele].dim(); ++i)
       for(j=0; j<mArray[iele].dim(); ++j)
          if( dofType[ (*allDOFs)[iele][i] ] == 1 ||
              dofType[ (*allDOFs)[iele][j] ] == 1) {
             mArray[iele][i][j] = 0.0;
          }
   }
 }

}

// used in nonlinear dynamics
void
Domain::createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray, FullSquareMatrix *&cArray)
{
 // Allocate array of pointers to FullSquareMatrix to store
 // the element stiffness matrices, element damping matrices and element mass matrices
 kArray = new FullSquareMatrix[numele];
 mArray = new FullSquareMatrix[numele];
 cArray = new FullSquareMatrix[numele];

 // Allocate the correct size for each elements stiffness & mass matrix
 int iele;
 for(iele = 0; iele<numele; ++iele) {
   int dimension = packedEset[iele]->numDofs();
   kArray[iele].setSize(dimension);
   kArray[iele].zero();
   mArray[iele].setSize(dimension);
   cArray[iele].setSize(dimension);
 }

 // Form and store element damping matrices and mass matrices into arrays
 for(iele=0; iele<numele; ++iele) {
   // note: only lumped mass matrix is supported currently for elements with rotation dofs in nonlinear dynamics
   //       (only the euler beam element is affected)
   double mratio = (packedEset[iele]->hasRot() && sinfo.isNonLin() && sinfo.isDynam()) ? 0 : geoSource->getMRatio();
   mArray[iele] = packedEset[iele]->massMatrix(nodes, mArray[iele].data(), mratio);
   cArray[iele] = packedEset[iele]->dampingMatrix(nodes, cArray[iele].data());
 }

 // add Rayleigh damping
 int i,j;
 double alpha, beta;
 double *karray = new double[maxNumDOFs*maxNumDOFs];
 FullSquareMatrix kel;
 for(iele=0; iele<numele; ++iele) {
   if(!packedEset[iele]->getProperty()) continue; // phantom
   kel = packedEset[iele]->stiffness(nodes, karray);
   if(packedEset[iele]->isConstraintElement()) cArray[iele].zero();
   else {
     alpha = (packedEset[iele]->isDamped()) ? packedEset[iele]->getProperty()->alphaDamp : sinfo.alphaDamp;
     beta  = (packedEset[iele]->isDamped()) ? packedEset[iele]->getProperty()->betaDamp : sinfo.betaDamp;
     for(i=0; i<cArray[iele].dim(); ++i)
       for(j=0; j<cArray[iele].dim(); ++j)
         cArray[iele][i][j] += alpha*mArray[iele][i][j] + beta*kel[i][j];
   }
 }
 delete [] karray;

 // zero rotational degrees of freedom within element mass matrices and damping matrices
 // for nonlinear implicit dynamics
 if(sinfo.zeroRot && sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0) {
   int *dofType = dsa->makeDofTypeArray();
   for(iele=0; iele<numele; ++iele) {
     for(i=0; i<mArray[iele].dim(); ++i)
       for(j=0; j<mArray[iele].dim(); ++j)
          if( dofType[ (*allDOFs)[iele][i] ] == 1 ||
              dofType[ (*allDOFs)[iele][j] ] == 1) {
             mArray[iele][i][j] = 0.0;
             cArray[iele][i][j] = 0.0;
         }
   } 
 }
}

bool
Domain::reactionsReqd(double time, int step)
{
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for(int iInfo = 0; iInfo < numOutInfo; ++iInfo) {
    if((oinfo[iInfo].type == OutputInfo::Reactions || oinfo[iInfo].type == OutputInfo::Reactions6)
       && (step%oinfo[iInfo].interval == 0 || time == 0.0)) return true;
  }
  return false;
}

void
Domain::postProcessing(GeomState *geomState, Vector& force, Vector &aeroForce,
                       double time, int step, double* velocity, double *vcx,
                       Corotator **allCorot, FullSquareMatrix *mel, double *acceleration,
                       double *acx, GeomState *refState, Vector *reactions)
{
  if(time == sinfo.initialTime) {
    geoSource->openOutputFiles();
    //printStatistics();
  }

  if( sinfo.nRestart > 0 && velocity !=0) {
    StackVector v_n(velocity, numUncon());
    writeRestartFile(time, step, v_n, geomState);
  }

  int numOutInfo = geoSource->getNumOutInfo();
  for(int iInfo = 0; iInfo < numOutInfo; ++iInfo)
  {
    postProcessingImpl(iInfo, geomState, force, aeroForce, time, step, velocity, vcx,
                       allCorot, mel, acceleration, acx, refState, reactions);
  }

}

void
Domain::postProcessingImpl(int iInfo, GeomState *geomState, Vector& force, Vector &aeroForce,
                           double time, int step, double* velocity, double *vcx,
                           Corotator **allCorot, FullSquareMatrix *mel, double *acceleration,
                           double *acx, GeomState *refState, Vector *reactions)
{
 if(outFlag && !nodeTable) makeNodeTable(outFlag);
 int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes

 //fprintf(stderr,"Running postPro in NLSTATIC\n");
 enum {SXX=0,SYY=1,SZZ=2,SXY= 3,SYZ= 4,SXZ= 5,VON=6,
       EXX=7,EYY=8,EZZ=9,EXY=10,EYZ=11,EXZ=12,STRAINVON=13,
       VONTOP=14,VONBOT=15};

 enum {INX,INY,INZ,AXM,AYM,AZM};

 enum {PSTRESS1=0,PSTRESS2=1,PSTRESS3=2,
       PSTRAIN1=3,PSTRAIN2=4,PSTRAIN3=5};

  int i, nodeI, realNode; 
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // Check output interval
  if(step%oinfo[iInfo].interval != 0 && time != 0.0) return;

  int w = oinfo[iInfo].width;
  int p = oinfo[iInfo].precision;

  int first_node, last_node;
  if (oinfo[iInfo].nodeNumber == -1) {
    first_node=0;
    last_node=numNodes;
  }
  else  {
    first_node=oinfo[iInfo].nodeNumber;
    last_node=first_node+1;
  }
  int nNodes = last_node - first_node;
  int nPrintNodes = (outFlag && oinfo[iInfo].nodeNumber == -1) ? nodes.nnz() : nNodes;

  switch(oinfo[iInfo].type) {
    case OutputInfo::Displacement:  {
      double (*data)[3] = new double[nPrintNodes][3];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        data[nodeI][0] = (nodes[iNode] && iNode<geomState->numNodes()) ? (*geomState)[iNode].x-nodes[iNode]->x : 0;
        data[nodeI][1] = (nodes[iNode] && iNode<geomState->numNodes()) ? (*geomState)[iNode].y-nodes[iNode]->y : 0;
        data[nodeI][2] = (nodes[iNode] && iNode<geomState->numNodes()) ? (*geomState)[iNode].z-nodes[iNode]->z : 0;
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Temperature:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        data[nodeI] = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].x : 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } 
      break;
    case OutputInfo::Disp6DOF:  {
      double (*data)[6] = new double[nPrintNodes][6];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes()) {
          if (nodes[iNode]) {
            data[nodeI][0] = (*geomState)[iNode].x - nodes[iNode]->x;
            data[nodeI][1] = (*geomState)[iNode].y - nodes[iNode]->y;
            data[nodeI][2] = (*geomState)[iNode].z - nodes[iNode]->z;
          } else {
            std::fill_n(&data[nodeI][0], 3, 0.0);
          }
          mat_to_vec((*geomState)[iNode].R, &data[nodeI][3]);
        } else {
          std::fill_n(&data[nodeI][0], 6, 0.0);
        }
      }
      geoSource->outputNodeVectors6(iInfo, data, nPrintNodes, time);
      delete [] data;
    } 
      break;
    case OutputInfo::RotationMatrix:  {
      double (*data)[9] = new double[nPrintNodes][9];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if (iNode < geomState->numNodes() && nodes[iNode]) {
          data[nodeI][0] = (*geomState)[iNode].R[0][0];
          data[nodeI][1] = (*geomState)[iNode].R[0][1];
          data[nodeI][2] = (*geomState)[iNode].R[0][2];
          data[nodeI][3] = (*geomState)[iNode].R[1][0];
          data[nodeI][4] = (*geomState)[iNode].R[1][1];
          data[nodeI][5] = (*geomState)[iNode].R[1][2];
          data[nodeI][6] = (*geomState)[iNode].R[2][0];
          data[nodeI][7] = (*geomState)[iNode].R[2][1];
          data[nodeI][8] = (*geomState)[iNode].R[2][2];
        } else {
          std::fill_n(&data[nodeI][0], 9, 0.0);
        }
      }
      geoSource->outputNodeVectors9(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Velocity6: { 
      if(!velocity) break;
      StackVector v_n(velocity, numUncon());
      double (*data)[6] = new double[nPrintNodes][6];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        getOrAddDofForPrint(false, v_n, vcx, first_node+iNode, data[nodeI], &DofSet::Xdisp,
                            data[nodeI]+1, &DofSet::Ydisp, data[nodeI]+2, &DofSet::Zdisp);
        getOrAddDofForPrint(false, v_n, vcx, first_node+iNode, data[nodeI]+3, &DofSet::Xrot,
                              data[nodeI]+4, &DofSet::Yrot, data[nodeI]+5, &DofSet::Zrot);
        if(oinfo[iInfo].angularouttype == OutputInfo::spatial) {
          // transform from convected to spatial angular velocity
          if (iNode < geomState->numNodes() && nodes[iNode]) {
            double V[3] = { data[nodeI][3], data[nodeI][4], data[nodeI][5] };
            mat_mult_vec((*geomState)[iNode].R,V,data[nodeI]+3,0); // v = R*V
          }
        }
        else if(oinfo[iInfo].angularouttype == OutputInfo::total) {
          // transform from convected angular velocity to time derivative of total rotation vector
          if (iNode < geomState->numNodes() && nodes[iNode]) {
#ifdef USE_EIGEN3
            Eigen::Vector3d V, Psi;
            Eigen::Map<Eigen::Vector3d> Psidot(data[nodeI]+3);
            Eigen::Matrix3d R, T;
            V << data[nodeI][3], data[nodeI][4], data[nodeI][5];
            R << (*geomState)[iNode].R[0][0], (*geomState)[iNode].R[0][1], (*geomState)[iNode].R[0][2],
                 (*geomState)[iNode].R[1][0], (*geomState)[iNode].R[1][1], (*geomState)[iNode].R[1][2],
                 (*geomState)[iNode].R[2][0], (*geomState)[iNode].R[2][1], (*geomState)[iNode].R[2][2];
            mat_to_vec(R, Psi);
            tangential_transf(Psi, T);
            Psidot = T.inverse()*V;
#else
            data[nodeI][3] = data[nodeI][4] = data[nodeI][5] = 0;
#endif
          }
        }
      }
      geoSource->outputNodeVectors6(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Velocity:  {
      if(!velocity) break;
      StackVector v_n(velocity, numUncon());
      double (*data)[3] = new double[nPrintNodes][3];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        getOrAddDofForPrint(false, v_n, vcx, first_node+iNode, data[nodeI], &DofSet::Xdisp,
                            data[nodeI]+1, &DofSet::Ydisp, data[nodeI]+2, &DofSet::Zdisp);
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    } 
      break;
    case OutputInfo::TemperatureFirstTimeDerivative: {
      if(!velocity) break;
      StackVector v_n(velocity, numUncon());
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode) {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        getOrAddDofForPrint(false, v_n, vcx, first_node+iNode, data+nodeI, &DofSet::Temp);
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } 
      break;
    case OutputInfo::Accel6:  {
      if(!acceleration) break;
      StackVector a_n(acceleration, numUncon());
      double (*data)[6] = new double[nPrintNodes][6];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        getOrAddDofForPrint(false, a_n, acx, first_node+iNode, data[nodeI],
                            &DofSet::Xdisp, data[nodeI]+1, &DofSet::Ydisp,
                            data[nodeI]+2, &DofSet::Zdisp);
        getOrAddDofForPrint(false, a_n, (double *) acx, first_node+iNode, data[nodeI]+3,
                            &DofSet::Xrot, data[nodeI]+4, &DofSet::Yrot, data[nodeI]+5,
                            &DofSet::Zrot);
        if(oinfo[iInfo].angularouttype == OutputInfo::spatial) {
          // transform from convected to spatial angular acceleration
          if (iNode < geomState->numNodes() && nodes[iNode]) {
            double A[3] = { data[nodeI][3], data[nodeI][4], data[nodeI][5] };
            mat_mult_vec((*geomState)[iNode].R,A,data[nodeI]+3,0); // a = R*A
          }
        }
        else if(oinfo[iInfo].angularouttype == OutputInfo::total) {
          // transform from convected angular acceleration to second time derivative of total rotation vector
          if (iNode < geomState->numNodes() && nodes[iNode]) {
#ifdef USE_EIGEN3
            Eigen::Vector3d A, V, Psi, Psidot;
            Eigen::Map<Eigen::Vector3d> Psiddot(data[nodeI]+3);
            Eigen::Matrix3d R, T, Tdot;
            A << data[nodeI][3], data[nodeI][4], data[nodeI][5];
            V << (*geomState)[iNode].v[3], (*geomState)[iNode].v[4], (*geomState)[iNode].v[5];
            R << (*geomState)[iNode].R[0][0], (*geomState)[iNode].R[0][1], (*geomState)[iNode].R[0][2],
                 (*geomState)[iNode].R[1][0], (*geomState)[iNode].R[1][1], (*geomState)[iNode].R[1][2],
                 (*geomState)[iNode].R[2][0], (*geomState)[iNode].R[2][1], (*geomState)[iNode].R[2][2];
            mat_to_vec(R, Psi);
            tangential_transf(Psi, T);
            Psidot = T.inverse()*V;
            tangential_transf_dot(Psi, Psidot, Tdot);
            Psiddot = T.inverse()*(A - Tdot*Psidot);
#else
            data[nodeI][3] = data[nodeI][4] = data[nodeI][5] = 0;
#endif
          }
        }
      }
      geoSource->outputNodeVectors6(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Acceleration: {
      if(!acceleration) break;
      StackVector a_n(acceleration, numUncon());
      double (*data)[3] = new double[nPrintNodes][3];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode) {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        getOrAddDofForPrint(false, a_n, acx, first_node+iNode, data[nodeI], &DofSet::Xdisp, data[nodeI]+1, &DofSet::Ydisp, data[nodeI]+2, &DofSet::Zdisp);
      }
      geoSource->outputNodeVectors(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispX:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        data[nodeI] = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].x - nodes[iNode]->x : 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispY:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        data[nodeI] = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].y - nodes[iNode]->y : 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispZ:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        data[nodeI] = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].z - nodes[iNode]->z : 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotX:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if(iNode < geomState->numNodes()) {
          double rot[3];
          mat_to_vec((*geomState)[iNode].R,rot);
          data[nodeI] = rot[0];
        }
        else data[nodeI] = 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotY:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if(iNode < geomState->numNodes()) {
          double rot[3];
          mat_to_vec((*geomState)[iNode].R,rot);
          data[nodeI] = rot[1];
        }
        else data[nodeI] = 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotZ:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if(iNode < geomState->numNodes()) {
          double rot[3];
          mat_to_vec((*geomState)[iNode].R,rot);
          data[i] = rot[2];
        } 
        else data[nodeI] = 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::DispMod:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        double x = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].x - nodes[iNode]->x : 0;
        double y = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].y - nodes[iNode]->y : 0;
        double z = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].z - nodes[iNode]->z : 0;
        data[nodeI] = sqrt(x*x+y*y+z*z);
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::RotMod:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        if(iNode < geomState->numNodes()) {
          double rot[3];
          mat_to_vec((*geomState)[iNode].R,rot);
          data[nodeI] = sqrt(rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
        }
        else data[nodeI] = 0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::TotMod:  {
      double *data = new double[nPrintNodes];
      for (i = 0, realNode = -1; i < nNodes; ++i) {
        int iNode = first_node+i;
        if(outFlag) { if(nodes[iNode] == 0) continue; nodeI = ++realNode; } else nodeI = i;
        double x = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].x - nodes[iNode]->x : 0;
        double y = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].y - nodes[iNode]->y : 0;
        double z = (nodes[iNode] && iNode < geomState->numNodes()) ? (*geomState)[iNode].z - nodes[iNode]->z : 0;
        double rot[3];
        if(iNode < geomState->numNodes()) 
          mat_to_vec((*geomState)[iNode].R,rot);
        else {
          rot[0] = rot[1] = rot[2] = 0;
        }
        data[nodeI] = sqrt(x*x+y*y+z*z+rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    }
      break;
    case OutputInfo::Rigid:
/*
         Vector rigid(maxNumDOFs);
         Vector rDisp(numUncon(),0.0);
         for(i=0; i<numele; ++i) {
           allCorot[i]->extractRigidMotion(geomState, nodes, rigid.data());
           int NodesPerElement = packedEset[iele]->numNodes();
           for(k=0; k<NodesPerElement; ++k) {
             rDisp[(*elemToNode)[iele][k]] += rigid[k];
             (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
           }
         }
*/
    case OutputInfo::StressXX:
      getStressStrain(*geomState, allCorot,  iInfo, SXX, time, refState);
      break;
    case OutputInfo::StressYY:
      getStressStrain(*geomState, allCorot,  iInfo, SYY, time, refState);
      break;
    case OutputInfo::StressZZ:
      getStressStrain(*geomState, allCorot,  iInfo, SZZ, time, refState);
      break;
    case OutputInfo::StressXY:
      getStressStrain(*geomState, allCorot,  iInfo, SXY, time, refState);
      break;
    case OutputInfo::StressYZ:
      getStressStrain(*geomState, allCorot,  iInfo, SYZ, time, refState);
      break;
    case OutputInfo::StressXZ:
      getStressStrain(*geomState, allCorot,  iInfo, SXZ, time, refState);
      break;
    case OutputInfo::StrainXX:
      getStressStrain(*geomState, allCorot,  iInfo, EXX, time, refState);
      break;
    case OutputInfo::StrainYY:
      getStressStrain(*geomState, allCorot,  iInfo, EYY, time, refState);
      break;
    case OutputInfo::StrainZZ:
      getStressStrain(*geomState, allCorot,  iInfo, EZZ, time, refState);
      break;
    case OutputInfo::StrainXY:
      getStressStrain(*geomState, allCorot,  iInfo, EXY, time, refState);
      break;
    case OutputInfo::StrainYZ:
      getStressStrain(*geomState, allCorot,  iInfo, EYZ, time, refState);
      break;
    case OutputInfo::StrainXZ:
      getStressStrain(*geomState, allCorot,  iInfo, EXZ, time, refState);
      break;
    case OutputInfo::StressVM:
      getStressStrain(*geomState, allCorot,  iInfo, VON, time, refState);
      break;
    case OutputInfo::StrainVM:
      getStressStrain(*geomState, allCorot,  iInfo, STRAINVON, time, refState);
      break;
    case OutputInfo::StressPR1:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS1, time);
      break;
    case OutputInfo::StressPR2:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS2, time);
      break;
    case OutputInfo::StressPR3:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS3, time);
      break;
    case OutputInfo::StrainPR1:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN1, time);
      break;
    case OutputInfo::StrainPR2:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN2, time);
      break;
    case OutputInfo::StrainPR3:
      getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN3, time);
      break;
    case OutputInfo::EquivalentPlasticStrain:
      getStressStrain(*geomState, allCorot,  iInfo, EQPLSTRN, time, refState);
      break;
    case OutputInfo::InXForce:
      getElementForces(*geomState, allCorot, iInfo, INX, time);
      break;
    case OutputInfo::InYForce:
      getElementForces(*geomState, allCorot, iInfo, INY, time);
      break;
    case OutputInfo::InZForce:
      getElementForces(*geomState, allCorot, iInfo, INZ, time);
      break;
    case OutputInfo::AXMoment:
      getElementForces(*geomState, allCorot, iInfo, AXM, time);
      break;
    case OutputInfo::AYMoment:
      getElementForces(*geomState, allCorot, iInfo, AYM, time);
      break;
    case OutputInfo::AZMoment:
      getElementForces(*geomState, allCorot, iInfo, AZM, time);
      break;
    case OutputInfo::Energies: {

      // Since this file is used for both Non-linear Statics
      // and Non-linear Dynamics
      // we need both included in this routine.

      // Build displacement Vector
      Vector sol(numUncon(), 0.0 );
      int i;
      for(i=0; i<numNodes; ++i) {
        int xloc  = c_dsa->locate(i, DofSet::Xdisp);
        if(xloc >= 0)
          sol[xloc]  = ( (*geomState)[i].x - nodes[i]->x);
        int yloc  = c_dsa->locate(i, DofSet::Ydisp);
        if(yloc >= 0)
          sol[yloc]  = ( (*geomState)[i].y - nodes[i]->y);
        int zloc  = c_dsa->locate(i, DofSet::Zdisp);
        if(zloc >= 0)
          sol[zloc]  = ( (*geomState)[i].z - nodes[i]->z);
        double rot[3];
        mat_to_vec((*geomState)[i].R,rot);
        int xrot  = c_dsa->locate(i, DofSet::Xrot);
        if(xrot >= 0)
          sol[xrot]  = rot[0];
        int yrot  = c_dsa->locate(i, DofSet::Yrot);
        if(yrot >= 0)
          sol[yrot]  = rot[1];
        int zrot  = c_dsa->locate(i, DofSet::Zrot);
        if(zrot >= 0)
          sol[zrot]  = rot[2];
      }

      // Non-Linear Dynamics
      if (sinfo.probType == SolverInfo::NonLinDynam) {
        double dW=0.0, dWaero = 0.0, Wela=0.0, Wkin=0.0;

        if(!previousExtForce) {
          Wext=0.0;
          Waero=0.0;
          Wdmp=0.0;
          pWela=0.0;
          pWkin=0.0;
          previousExtForce = new Vector(force);
          previousAeroForce = new Vector(aeroForce);
          dW = force*sol;
          dWaero = aeroForce*sol;
        } else {
          double c = solInfo().newmarkGamma;
/*
     NOTE: The formula for dW should be (1-c)*f^n + c*f^{n+1}.
           However, force stores f^{n+1/2}
           and previousExtForce stores f^{n-1/2}. For this reason,
           the formula below looksk
           different. Furthermore, we had to assume
           f^n = (f^{n-1/2} + f^{n+1/2})/2 to minimize
           code changes. On the other hand, f^{n+1/2} = (f^n + f^{n+1})/2
           holds without any assumption.
*/

          dW = (c*force + (1.0-c)*(*previousExtForce)) * (sol - (*previousDisp));
          dWaero = (c*aeroForce + (1.0-c)*(*previousAeroForce)) *(sol - (*previousDisp));

          if(step==sinfo.initialTimeIndex) { dW*=2.0; dWaero *= 2.0; }
        }

        Wext += dW;
        Waero += dWaero;

        Vector vtmp(velocity,numUncon());
        Vector tmpVec(numUncon(),0.0);
        int iele, idof, jdof, dofn1, dofn2;

        // Compute Kinetic Energy
        // Compute vtmp^t M vtmp
        for(iele = 0; iele < numele; ++iele) {
          for(idof = 0; idof <  mel[iele].dim(); ++idof) {
            dofn1 = c_dsa->getRCN((*allDOFs)[iele][idof]);
            if(dofn1 >= 0) {
              for(jdof = 0; jdof <  mel[iele].dim(); ++jdof) {
                dofn2 = c_dsa->getRCN((*allDOFs)[iele][jdof]);
                if(dofn2 >= 0)
                  tmpVec[dofn1] += vtmp[dofn2]*mel[iele][idof][jdof];
              }
            }
          }
        }
      Wkin = 0.5 * (vtmp * tmpVec);

      // Compute Internal Energy
      // This is done at the element level.
      double EleWela;
      for(iele = 0; iele < numele; ++iele) {
        EleWela = 0.0;
        EleWela = allCorot[iele]->getElementEnergy(*geomState,nodes);
        Wela += EleWela;
      }

      double error = (time==sinfo.initialTime) ? 0.0 : (Wela+Wkin)-(pWela+pWkin)-dW;
      geoSource->outputEnergies(iInfo, time, Wext, Waero, Wela, Wkin, 0.0, error);

      pWela=Wela;
      pWkin=Wkin;

      if(!previousDisp) {
        previousDisp = new Vector(sol);
      } else {
        (*previousExtForce) = force;
        (*previousAeroForce) = aeroForce;
        (*previousDisp)     = sol;
      }

      // Non-Linear Statics
    } else {
      double lambda = time;
      if (time == 0.0) {
        double deltaLambda = solInfo().getNLInfo().dlambda;
        double maxLambda = solInfo().getNLInfo().maxLambda;
        if(deltaLambda == maxLambda) lambda = 1.0;
      }
      Wext=0.0;
      Waero=0.0;
      Wdmp=0.0;
      pWela=0.0;
      pWkin=0.0;

      // Wkin = kinetic energy
      // Total Energy = Wext+Wela+Wkin
      double Wkin=0.0;

      // Wext = external energy
      Wext = lambda*force * sol;
      Waero = lambda*aeroForce * sol;

      // Compute Internal Energy
      // This is done at the element level.
      // Wela = elastic energy
      int iele;
      double Wela = 0.0;
      double EleWela;
        for(iele = 0; iele < numele; ++iele) {
          EleWela = 0.0;
          EleWela = allCorot[iele]->getElementEnergy(*geomState,nodes);
          Wela += EleWela;
        }

          double error = Wext+Wela+Wkin;
          geoSource->outputEnergies(iInfo,time,Wext, Waero, Wela,Wkin,0.0,error);
        }
    } break;
    case OutputInfo::AeroForce: break; // this is done in FlExchange.C
    case OutputInfo::AeroXForce:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
        data[nodeI]  = (xloc >= 0) ? aeroForce[xloc] : 0.0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroYForce:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
        data[nodeI]  = (yloc >= 0) ? aeroForce[yloc] : 0.0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroZForce:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
        data[nodeI] = (zloc >= 0) ? aeroForce[zloc] : 0.0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroXMom:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
        data[nodeI] = (xrot >= 0) ? aeroForce[xrot] : 0.0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroYMom:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
        data[nodeI] = (yrot >= 0) ? aeroForce[yrot] : 0.0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::AeroZMom:  {
      double *data = new double[nPrintNodes];
      for (int iNode = 0, realNode = -1; iNode < nNodes; ++iNode)  {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
        data[nodeI] = (zrot >= 0) ? aeroForce[zrot] : 0.0;
      }
      geoSource->outputNodeScalars(iInfo, data, nPrintNodes, time);
      delete [] data;
    } break;
    case OutputInfo::Reactions: {
      if(!reactions) break;
      double (*rxyz)[3] = new double[nPrintNodes][3];
      DofSet dofs[3] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp };
      for(int iNode = 0, realNodes = -1; iNode < nNodes; ++iNode) {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        for(int k = 0; k < 3; ++k) {
          int dof =   dsa->locate(first_node+iNode, dofs[k].list());
          int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
          rxyz[nodeI][k] = (cdof >= 0) ? (*reactions)[cdof] : 0;     // constrained
        }
      }
      geoSource->outputNodeVectors(iInfo, rxyz, nPrintNodes, time);
      delete [] rxyz;
    } break;
    case OutputInfo::Reactions6: {
      if(!reactions) break;
      double (*rxyz)[6] = new double[nPrintNodes][6];
      DofSet dofs[6] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
                         DofSet::Xrot, DofSet::Yrot, DofSet::Zrot };
      for(int iNode = 0, realNode = -1; iNode < nNodes; ++iNode) {
        if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
        for(int k = 0; k < 6; ++k) {
          int dof =   dsa->locate(first_node+iNode, dofs[k].list());
          int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
          rxyz[nodeI][k] = (cdof >= 0) ? (*reactions)[cdof] : 0;     // constrained
        }
      }
      geoSource->outputNodeVectors6(iInfo, rxyz, nPrintNodes, time);
      delete [] rxyz;
    } break;
     case OutputInfo::Statevector:
        break;
     case OutputInfo::Residual:
        break;
     case OutputInfo::Jacobian:
        break;
     case OutputInfo::RobData:
        break;
     case OutputInfo::SampleMesh:
        break;
    default:
      fprintf(stderr," *** WARNING: Output case %d not implemented for non-linear direct solver \n", iInfo);
      break;
  }

}

void
Domain::createCorotators(Corotator **allCorot)
{
 // Allocate memory for element stiffness matrix
 double *kmatrix = new double[maxNumDOFs*maxNumDOFs];

 // Loop over elements to get their Corotators
 for(int iele = 0; iele < numele; ++iele)
   allCorot[iele] = packedEset[iele]->getCorotator(nodes, kmatrix,
                                      sinfo.getNLInfo().fitAlgShell,
                                      sinfo.getNLInfo().fitAlgBeam);

 delete [] kmatrix;
}

void
Domain::getGeometricStiffness(GeomState &geomState, Vector& elementInternalForce,
                              Corotator **allCorot, FullSquareMatrix *&geomKelArray)
{

   // Get Geometric Stiffness

   int iele;
   for(iele = 0; iele < numele; ++iele) {

      // Get updated tangent stiffness matrix and element internal force
      allCorot[iele]->formGeometricStiffness(geomState, nodes,
                                             geomKelArray[iele],
                                             elementInternalForce.data());
      if(!solInfo().getNLInfo().unsymmetric) geomKelArray[iele].symmetrize();
   }

   if(domain->pressureFlag()) {
     for(iele = 0; iele < numele;  ++iele) {
       // If there is a zero pressure defined, skip the element
       if(packedEset[iele]->getPressure() == 0) continue;
 
       // Compute (linear) element pressure force in the local coordinates
       elementInternalForce.zero();
       packedEset[iele]->computePressureForce(nodes, elementInternalForce, &geomState, 1);
       elementInternalForce *= -1.0; // due to IDISP6 -1

       // Include the "load stiffness matrix" in kel[iele]
       allCorot[iele]->getDExternalForceDu(geomState, nodes, geomKelArray[iele],
                                           elementInternalForce.data());
     }
   }

}

void
Domain::computeGeometricPreStress(Corotator **&allCorot, GeomState *&geomState,
                                  FullSquareMatrix *&kelArray, StaticTimers *times,
                                  FullSquareMatrix *&geomKelArray, FullSquareMatrix *&melArray,
                                  bool melFlag)
{
   // ... ALLOCATE MEMORY FOR THE ARRAY OF COROTATORS
   times->corotatorTime -= getTime();
   allCorot = new Corotator *[numElements()];

   // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
   createCorotators(allCorot);
   times->corotatorTime += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Create Element Corotators %19.5f s\n", times->corotatorTime/1000.0);
#endif

   // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
   times->kelArrayTime -= getTime();
   if(melFlag) createKelArray(kelArray, melArray);
   else createKelArray(kelArray);
   times->kelArrayTime += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Create Element Stiffness Array %14.5f s\n", times->kelArrayTime/1000.0);
#endif

   times->timeGeom -= getTime();
   geomState = new GeomState( *getDSA(), *getCDSA(), getNodes(), &getElementSet() );
   times->timeGeom += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Build GeomState %29.5f s\n", times->timeGeom/1000.0);
#endif

   // geometric prestress for linear: update geomState with initial displacements and compute element stiffness matrices
   if(!sinfo.isNonLin() && (numInitDisp6() > 0) && (sinfo.gepsFlg == 1)) {
     times->timePresc -= getTime();
     geomState->updatePrescribedDisplacement(getInitDisp6(), numInitDisp6(), 1.0);
     times->timePresc += getTime();
#ifdef PRINT_NLTIMERS
     fprintf(stderr," ... Update Geometry to initial state %12.5f s\n", times->timePresc/1000.0);
#endif
   }

   times->buildStiffAndForce -= getTime();
   Vector elementInternalForce(maxNumDOF(), 0.0);
   Vector residual(numUncon(), 0.0);
   getStiffAndForce(*geomState, elementInternalForce, allCorot,
                    kelArray, residual, 1.0, 0.0, geomState, (Vector*) NULL, melArray);
   times->buildStiffAndForce += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Build Element Tangent Stiffness %13.5f s\n",
           times->buildStiffAndForce/1000.0);
#endif

   // Buckling analysis if requested
   if(geomKelArray) {
     Vector elementInternalForce(maxNumDOF(), 0.0);
     getGeometricStiffness(*geomState, elementInternalForce, allCorot, geomKelArray);
   }
}


void
Domain::getStressStrain(GeomState &geomState, Corotator **allCorot,
                        int fileNumber, int stressIndex, double time,
                        GeomState *refState)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS, WHERE
  // ... ONLY THE AXIAL STRAIN (EXX) AND AXIAL STRESS (SXX) ARE CALCULATED

  int avgnum  = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  int k;

  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;

  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;

  // ... WRITE CURRENT TIME VALUE
  //if(oinfo[fileNumber].nodeNumber == -1)
   // fprintf(oinfo[fileNumber].filptr,"  % *.*E\n",w,p,time);

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(stress == 0)
    stress = new Vector(numnodes,0.0);
  else 
    stress->zero();

  if(weight == 0)
    weight = new Vector(numnodes,0.0);
  else 
    weight->zero();

  if(elDisp == 0)
    elDisp = new Vector(maxNumDOFs,0.0);

  int iele;
  if((elstress == 0)||(elweight == 0)) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = myMax(maxNodesPerElement, NodesPerElement);
    }
    if(elstress == 0) elstress = new Vector(maxNodesPerElement, 0.0);
    if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
  }


  int flag;
  for(iele = 0; iele < numElements(); ++iele) {
    if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isMpcElement()) continue;
    elDisp->zero();
    elstress->zero();
    elweight->zero();

     // extract deformations from current Geometry State of structure
     allCorot[iele]->extractDeformations( geomState, nodes,
                                          elDisp->data(), flag);
//---------------------------------------------------------------------------

   int *nodeNumbers = new int[maxNumNodes];
   Vector elemNodeTemps(maxNumNodes);
   elemNodeTemps.zero();

   packedEset[iele]->nodes(nodeNumbers);

   int NodesPerElement = packedEset[iele]->numNodes();
   double *nodalTemperatures = 0;
   // Either get the nodal temperatures from the input file or
   // from the thermal model
   if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
   if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

   int iNode;
   if (sinfo.thermalLoadFlag || (sinfo.thermoeFlag>=0))
     for (iNode = 0; iNode < NodesPerElement; ++iNode) {
       if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
         elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
       else
         elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
     }

//----------------------------------------------------------------------------

     if (flag == 1) {
// USE LINEAR STRESS ROUTINE
// ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
     packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                   *elDisp, stressIndex, surface,
				   elemNodeTemps.data(), ylayer,
                                   zlayer, avgnum);

     } else if (flag == 2) {
// USE NON-LINEAR STRESS ROUTINE
     allCorot[iele]->getNLVonMises(*elstress, *elweight, geomState,
                                   refState, nodes, stressIndex, surface,
                                   elemNodeTemps.data(), ylayer, zlayer,
                                   avgnum);

     } else {
// NO STRESS RECOVERY
     }

// ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

//     int NodesPerElement = packedEset[iele]->numNodes();

     for(k=0; k<NodesPerElement; ++k) {
       (*stress)[(*elemToNode)[iele][k]] += (*elstress)[k];
       (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
     }

// ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
     if(avgnum == 0) {
       for(k=0; k<NodesPerElement; ++k)
         fprintf(oinfo[fileNumber].filptr," % *.*E",w,p,(*elstress)[k]);
       fprintf(oinfo[fileNumber].filptr,"\n");
     }
  }


// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

   if(avgnum == 1 || avgnum == 2) {

     if(oinfo[fileNumber].nodeNumber == -1) {
       int numNodes = geoSource->numNode();
       int numNodesOut = (outFlag) ? exactNumNodes : numNodes;
       double *data = new double[numNodesOut];
       for(k=0; k<numNodesOut; ++k) data[k] = 0;
       for(k=0; k<numNodes; ++k) {
          int l = (outFlag) ? nodeTable[k]-1 : k;
          if(l < 0) continue;
          if(k < numnodes && (*weight)[k] != 0)
            data[l] = (*stress)[k]/=(*weight)[k];
       }
       geoSource->outputNodeScalars(fileNumber, data, numNodesOut, time);
       delete [] data;
     } else {
       if((*weight)[oinfo[fileNumber].nodeNumber] == 0.0)
         fprintf(oinfo[fileNumber].filptr," %*.*E % *.*E\n",w,p,time,w,p,0.0);
       else
         fprintf(oinfo[fileNumber].filptr," %*.*E % *.*E\n",w,p,time,w,p,(*stress)[oinfo[fileNumber].nodeNumber]/=(*weight)[oinfo[fileNumber].nodeNumber]);
     }
   }
   fflush(oinfo[fileNumber].filptr);
}

void
Domain::getPrincipalStress(GeomState &geomState, Corotator **allCorot,
                         int fileNumber, int stressIndex, double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // set stress VS. strain for element subroutines
  int strInd;
  int strDir;
  if ((stressIndex==0)||(stressIndex==1)||(stressIndex==2)) {
    strInd = 0;
    strDir = stressIndex+1;
  }
  else if ((stressIndex==3)||(stressIndex==4)||(stressIndex==5)) {
    strInd = 1;
    strDir = stressIndex-2;
  }
  else {
    fprintf(stderr,"*** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS

  int avgnum  = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  int j,k;

  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;
  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;

  // ... OUTPUT NODE NUMBER
  int n = oinfo[fileNumber].nodeNumber;

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(p_stress == 0) p_stress = new FullM(numnodes,6);
  if(weight == 0) weight = new Vector(numnodes,0.0);
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  int iele;
  if((p_elstress == 0)||(elweight == 0)) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = myMax(maxNodesPerElement, NodesPerElement);
    }
    if(p_elstress == 0) p_elstress = new FullM(maxNodesPerElement,9);
    if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
  }

  // zero the vectors
  p_stress->zero();
  weight->zero();

  // ... WRITE CURRENT TIME VALUE
  if(n == -1) {
    fprintf(oinfo[fileNumber].filptr,"  % *.*E\n",w,p,time);
  }

  int flag;

  for(iele = 0; iele < numElements(); ++iele) {
    if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isMpcElement()) continue;
    elDisp->zero();
    p_elstress->zero();
    elweight->zero();

     // extract deformations from current Geometry State of structure
     allCorot[iele]->extractDeformations( geomState, nodes,
                                          elDisp->data(), flag);

     if (flag == 1) {
// USE LINEAR STRESS ROUTINE
// ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
     packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                    *elDisp, strInd, surface);

     } else if (flag == 2) {
// USE NON-LINEAR STRESS ROUTINE
     allCorot[iele]->getNLAllStress(*p_elstress, *elweight, geomState,
                                    nodes, strInd);

     } else {
// NO STRESS RECOVERY
     }

// ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

    int NodesPerElement = elemToNode->num(iele);

    for(k=0; k<NodesPerElement; ++k) {
      for(j=0; j<6; ++j) {
        (*p_stress)[(*elemToNode)[iele][k]][j] += (*p_elstress)[k][j];
      }
      (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
    }

// ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
//     THIS WRITES THE CHOSEN PRINCIPAL STRESS FOR EACH ELEMENT
    if(avgnum == 0) {
      geoSource->outputNodeScalars(fileNumber, (*p_elstress)[0]+(5+strDir), 1);
      //fprintf(oinfo[fileNumber].filptr," % *.*E\n", w,p,(*p_elstress)[0][5+strDir]);
    }
  }

// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

  if(avgnum == 1 || avgnum == 2) {

    if(n == -1) {
      for(k=0; k<numnodes; ++k) {
        if((*weight)[k] == 0.0) {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] = 0.0;
          }
        } else  {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] /= (*weight)[k];
          }
        }
      }
    } else {
      if((*weight)[n] == 0.0) {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] = 0.0;
        }
      } else {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] /= (*weight)[n];
        }
      }
    }

// ... CALCULATE PRINCIPALS AT EACH NODE

    double svec[6], pvec[3];
    if(n == -1) {
      int numNodes = geoSource->numNode();
      int numNodesOut = (outFlag) ? exactNumNodes : numNodes;
      double *globalPVec = new double[numnodes];
      for(k=0; k<numNodes; ++k) {
        int l = (outFlag) ? nodeTable[k]-1 : k;
        if(l < 0) continue;
        for (j=0; j<6; ++j) {
          svec[j] = (*p_stress)[k][j];
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        globalPVec[l] = pvec[strDir-1];
        //fprintf(oinfo[fileNumber].filptr," % *.*E\n",w,p,pvec[strDir-1]);
      }
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodesOut, time);
    }
    else {
      for (j=0; j<6; ++j) {
        svec[j] = (*p_stress)[n][j];
      }
// Convert Engineering to Tensor Strains
      if(strInd != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      geoSource->outputNodeScalars(fileNumber, pvec+strDir-1, 1);
    }
  }
  else {
    fflush(oinfo[fileNumber].filptr);
  }
}


// Nonlinear version of getElementForces
void
Domain::getElementForces( GeomState &geomState, Corotator **allCorot,
                          int fileNumber, int forceIndex, double time)
{

  // allocate integer array to store node numbers
  // make sure that maxNumNodes is computed
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();
  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag)  nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  // ... FORCES ARE CALCULATED FOR BARS AND BEAMS CURRENTLY

  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  // ... Watch out for this as the max number of nodes may increase
  int NodesPerElement = 2;

  FullM forces(numele, NodesPerElement);

  int i, iele;

  // note, we are reusing elstress here for an elemental force vector
  if(elstress == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = myMax(maxNodesPerElement, NodesPerElement);
    }
    elstress = new Vector(maxNodesPerElement, 0.0);
  }

  int flag;

  for(iele=0; iele<numele; ++iele) {

  packedEset[iele]->nodes(nodeNumbers);

// ... DETERMINE ELEMENT DISPLACEMENT VECTOR

    allCorot[iele]->extractDeformations(geomState, nodes,
                                        elDisp->data(), flag);

// ... CALCULATE INTERNAL FORCE VALUE FOR EACH ELEMENT
     if(sinfo.thermalLoadFlag || (sinfo.thermoeFlag >=0)) {
       int iNode;
       for(iNode=0; iNode<NodesPerElement; ++iNode) {
        if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
     }

    packedEset[iele]->getIntrnForce(*elstress,nodes,elDisp->data(),forceIndex,elemNodeTemps.data());

// ... COPY ELEMENT'S NODAL FORCES INTO A TOTAL FORCE MATRIX

    for(i=0; i<NodesPerElement; ++i)
      forces[iele][i] = (*elstress)[i];

  }

// ... PRINT THE ELEMENT FORCES TO A FILE
   geoSource->outputElemVectors(fileNumber, forces.data(), numele, time);
}

// Nonlinear restart file
void
Domain::writeRestartFile(double time, int timeIndex, Vector &v_n,
                         GeomState *geomState, const char *ext)
{
// either test for pointer or frequency > 0

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if((timeIndex % sinfo.nRestart == 0) || (time >= sinfo.tmax-0.1*domain->solInfo().getTimeStep())) {

   int fn;
   if(strlen(ext) != 0) {
     char *currentRestartFile = new char[strlen(cinfo->currentRestartFile)+strlen(ext)+1];
     strcpy(currentRestartFile, cinfo->currentRestartFile);
     strcat(currentRestartFile, ext);
     fn = open(currentRestartFile, O_WRONLY | O_CREAT, 0666);
     delete [] currentRestartFile;
   } else
   fn = open(cinfo->currentRestartFile, O_WRONLY | O_CREAT, 0666);
   if(fn >= 0) {
     int writeSize;
     writeSize = write(fn, &timeIndex, sizeof(int));
     if(writeSize != sizeof(int))
       fprintf(stderr," *** ERROR: Writing restart file time index\n");

     writeSize = write(fn, &time, sizeof(double));
     if(writeSize != sizeof(double))
       fprintf(stderr," *** ERROR: Writing restart file time\n");

      writeSize = write(fn, v_n.data(), v_n.size()*sizeof(double));
      if(int(writeSize) != int(v_n.size()*sizeof(double)))
        fprintf(stderr," *** ERROR: Writing restart file velocity\n");

     double *positions = (double *) dbg_alloca(sizeof(double)*3*numnodes);
     geomState->getPositions(positions);
     writeSize = write(fn, positions,numnodes*3*sizeof(double));
     if(int(writeSize) != int(numnodes*3*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file geometry_state\n");

     double *rotations = (double *) dbg_alloca(sizeof(double)*9*numnodes);
     geomState->getRotations(rotations);
     writeSize = write(fn, rotations,numnodes*9*sizeof(double));
     if(int(writeSize) != int(numnodes*9*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file geometry_state\n");

     // new method of storing the element states in the GeomState object
     int numElemStates = geomState->getTotalNumElemStates();
     double *elemStates = (double *) dbg_alloca(sizeof(double)*numElemStates);
     geomState->getElemStates(elemStates);
     writeSize = write(fn, elemStates, numElemStates*sizeof(double));
     if(int(writeSize) != int(numElemStates*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file geometry_state\n");

     // PJSA 9-17-2010 (note: this idea of the element storing the internal states is deprecated
     // and will eventually be removed
     int numEle = packedEset.last();
     for(int i = 0; i < numEle; ++i)
       packedEset[i]->writeHistory(fn);

     close(fn);
   } else {
      perror(" *** ERROR: Restart file could not be opened: ");
      exit(-1);
   }
 }
}

void
Domain::readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                        Vector &v_p, double *bcx, double *vcx,
                        GeomState &geomState, const char *ext)
{
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(cinfo->lastRestartFile) {
   int fn;
   if(strlen(ext) != 0) {
     char *lastRestartFile = new char[strlen(cinfo->lastRestartFile)+strlen(ext)+1];
     strcpy(lastRestartFile, cinfo->lastRestartFile);
     strcat(lastRestartFile, ext);
     fn = open(lastRestartFile, O_RDONLY );
     delete [] lastRestartFile;
   } else
   fn = open(cinfo->lastRestartFile,O_RDONLY );
   if(fn >= 0) {
     int restartTIndex;
     double restartT;
     int readSize = read(fn, &restartTIndex, sizeof(int));
     if(readSize != sizeof(int))
       fprintf(stderr," *** ERROR: Inconsistent restart file 1\n");
     sinfo.initialTimeIndex = restartTIndex;

     readSize = read(fn, &restartT, sizeof(double));
     if(readSize != sizeof(double))
       fprintf(stderr," *** ERROR: Inconsistent restart file 2\n");
     //fprintf(stderr,"Initial Time = %f\n\n",restartT);
     sinfo.initialTime = restartT;

     v_n.zero();
     readSize = read(fn, v_n.data(), sizeof(double)*v_n.size());
     if(int(readSize) != int(sizeof(double)*v_n.size()))
       fprintf(stderr," *** ERROR: Inconsistent restart file 2.5\n");

     double *positions = (double *)dbg_alloca(sizeof(double)*3*numnodes);
     readSize = read(fn, positions, numnodes*3*sizeof(double));
     if(int(readSize) != int(numnodes*3*sizeof(double)))
       fprintf(stderr," *** ERROR: Inconsistent restart file 3\n");
     geomState.setPositions(positions);

     double *rotations = (double *)dbg_alloca(sizeof(double)*9*numnodes);
     readSize = read(fn, rotations, numnodes*9*sizeof(double));
     if(int(readSize) != int(numnodes*9*sizeof(double)))
       fprintf(stderr," *** ERROR: Inconsistent restart file 4\n");
     geomState.setRotations(rotations);

     // new method of storing the element states in the GeomState object
     int numElemStates = geomState.getTotalNumElemStates();
     double *elemStates = (double *) dbg_alloca(sizeof(double)*numElemStates);
     readSize = read(fn, elemStates, numElemStates*sizeof(double));
     if(int(readSize) != int(numElemStates*sizeof(double)))
       fprintf(stderr," *** ERROR: Inconsistent restart file 5\n");
     geomState.setElemStates(elemStates);

     // PJSA 9-17-2010
     int numEle = packedEset.last();
     for(int i = 0; i < numEle; ++i) 
       packedEset[i]->readHistory(fn);

     close(fn);

     d_n.zero();
     a_n.zero();
     v_p.zero();
     for(int i = 0; i < numNodes(); ++i) {

       int xloc  = c_dsa->locate(i, DofSet::Xdisp );
       int xloc1 =   dsa->locate(i, DofSet::Xdisp );

       if(xloc >= 0)
         d_n[xloc]  = ( (geomState)[i].x - nodes[i]->x);
       else if (xloc1 >= 0)
         bcx[xloc1] = ( (geomState)[i].x - nodes[i]->x);

       int yloc  = c_dsa->locate(i, DofSet::Ydisp );
       int yloc1 =   dsa->locate(i, DofSet::Ydisp );

       if(yloc >= 0)
         d_n[yloc]  = ( (geomState)[i].y - nodes[i]->y);
       else if (yloc1 >= 0)
         bcx[yloc1] = ( (geomState)[i].y - nodes[i]->y);

       int zloc  = c_dsa->locate(i, DofSet::Zdisp);
       int zloc1 =   dsa->locate(i, DofSet::Zdisp);

       if(zloc >= 0)
         d_n[zloc]  = ( (geomState)[i].z - nodes[i]->z);
       else if (zloc1 >= 0)
         bcx[zloc1] = ( (geomState)[i].z - nodes[i]->z);
     }

     if(solInfo().aeroFlag >= 0 && solInfo().newmarkBeta != 0) // for nonlinear explicit aeroPreProcess is called in the main driver
       aeroPreProcess( d_n, v_n, a_n, v_p, bcx, vcx );

   } else {
      perror(" *** ERROR: Restart file could not be opened: ");
      //exit(-1);
   }

 }
 //else {
 //    fprintf(stderr, " ... No restart                     ...\n");
 //}


}

// nonlinear statics
void
Domain::computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
                             FullSquareMatrix *kel, double lambda, GeomState *refState)
{
  // TODO: include non-follower external forces on constrained dofs
  Vector elementInternalForce(maxNumDOF(), 0.0);
  Vector residual(numUncon(), 0.0);
  fc.zero();
  getInternalForce(*geomState, elementInternalForce, corotators, kel, residual, lambda, 0.0, refState, &fc);
}

// nonlinear dynamics
void
Domain::computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
                             FullSquareMatrix *kel, double time, GeomState *refState,
                             Vector &Vu, Vector &Au, double *vcx, double *acx,
                             SparseMatrix *_cuc, SparseMatrix *_ccc,
                             SparseMatrix *_muc, SparseMatrix *_mcc)
{
  // TODO: include non-follower external forces on constrained dofs
  Vector elementInternalForce(maxNumDOF(), 0.0);
  Vector residual(numUncon(), 0.0);
  fc.zero();

  getInternalForce(*geomState, elementInternalForce, corotators, kel, residual, 1.0, time, refState, &fc);

  CuCSparse *cuc = dynamic_cast<CuCSparse *>(_cuc);
  if(cuc) cuc->transposeMultAddNew(Vu.data(), fc.data()); // fc += Cuc^T * Vu

  CuCSparse *muc = dynamic_cast<CuCSparse *>(_muc);
  if(muc) muc->transposeMultAddNew(Au.data(), fc.data()); // fc += Muc^T * Vu

  CuCSparse *ccc = dynamic_cast<CuCSparse *>(_ccc);
  CuCSparse *mcc = dynamic_cast<CuCSparse *>(_mcc);
  Vector Dc(numDirichlet, 0.0);
  Vector Vc(numDirichlet, 0.0);
  Vector Ac(numDirichlet, 0.0);

  for(int i=0; i<numDirichlet; ++i) {
    int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
    if(dof < 0) continue;
    int cdof = c_dsa->invRCN(dof);
    if(cdof >= 0) {
      Vc[cdof] = vcx[dof];
      Ac[cdof] = (acx) ? acx[dof] : 0;
    }
  }

  if(ccc) ccc->multAddNew(Vc.data(), fc.data()); // fc += Ccc * Vc
  if(mcc) mcc->multAddNew(Ac.data(), fc.data()); // fc += Mcc * Ac

}

void
Domain::transformElemStiffAndForce(const GeomState &geomState, double *elementForce,
                                   FullSquareMatrix &kel, int iele, bool compute_tangents, FullSquareMatrix *mel)
{
#ifdef USE_EIGEN3
  // Convert from eulerian spatial to total lagrangian or updated lagrangian spatial
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> >
    G(elementForce, packedEset[iele]->numDofs());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>,Eigen::RowMajor>
    H(kel.data(),packedEset[iele]->numDofs(),packedEset[iele]->numDofs());
  int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
  int *nodes = packedEset[iele]->nodes();
  for(int k = 0; k < numNodes; ++k) {
    Eigen::Matrix3d R;
    R << geomState[nodes[k]].R[0][0], geomState[nodes[k]].R[0][1], geomState[nodes[k]].R[0][2],
         geomState[nodes[k]].R[1][0], geomState[nodes[k]].R[1][1], geomState[nodes[k]].R[1][2],
         geomState[nodes[k]].R[2][0], geomState[nodes[k]].R[2][1], geomState[nodes[k]].R[2][2];
    Eigen::Vector3d Psi;
    mat_to_vec(R, Psi);
    Eigen::Matrix3d T;
    tangential_transf(Psi, T);

    Eigen::Vector3d V = G.segment<3>(6*k+3);
    if(sinfo.newmarkBeta == 0) {
      if(domain->solInfo().galerkinPodRom && mel) {
        Eigen::Matrix3d M;
        M << (*mel)[6*k+3][6*k+3], (*mel)[6*k+3][6*k+4], (*mel)[6*k+3][6*k+5],
             (*mel)[6*k+4][6*k+3], (*mel)[6*k+4][6*k+4], (*mel)[6*k+4][6*k+5],
             (*mel)[6*k+5][6*k+3], (*mel)[6*k+5][6*k+4], (*mel)[6*k+5][6*k+5];

        G.segment<3>(6*k+3) = M*(T.inverse()*M.inverse()*T.transpose().inverse())*T*V;
      }
      else {
        G.segment<3>(6*k+3) = T.transpose()*V;
      }
    }
    else
      G.segment<3>(6*k+3) = T*V;
    if(compute_tangents) {
      Eigen::Matrix3d C1;
      directional_deriv1(Psi, V, C1);

      for(int l=0; l<2*numNodes; ++l) {
        H.block<3,3>(6*k+3,3*l) = (T*H.block<3,3>(6*k+3,3*l)).eval();
        H.block<3,3>(3*l,6*k+3) = (H.block<3,3>(3*l,6*k+3)*T.transpose()).eval();
      }
      for(int l=0; l<packedEset[iele]->numInternalNodes(); ++l) {
        H.block<3,1>(6*k+3,6*numNodes+l) = (T*H.block<3,1>(6*k+3,6*numNodes+l)).eval();
        H.block<1,3>(6*numNodes+l,6*k+3) = (H.block<1,3>(6*numNodes+l,6*k+3)*T.transpose()).eval();
      }
      H.block<3,3>(6*k+3,6*k+3) += 0.5*(C1 + C1.transpose());
    }
  }
  delete [] nodes;
#else
  cerr << "USE_EIGEN3 is not defined here in Domain::transformElemStiffAndForce\n";
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
void
Domain::transformNodalMoment(const GeomState &geomState, Eigen::Vector3d &G,
                             Eigen::Matrix3d &H, int inode, bool compute_tangents)
{
  // transform from eulerian spatial description to total lagrangian or updated lagrangian spatial description

  Eigen::Matrix3d R;
  R << geomState[inode].R[0][0], geomState[inode].R[0][1], geomState[inode].R[0][2],
       geomState[inode].R[1][0], geomState[inode].R[1][1], geomState[inode].R[1][2],
       geomState[inode].R[2][0], geomState[inode].R[2][1], geomState[inode].R[2][2];
  Eigen::Vector3d Psi;
  mat_to_vec(R, Psi);

  Eigen::Matrix3d T;
  tangential_transf(Psi, T);

  Eigen::Vector3d V = G;
  G = T*V;
  if(compute_tangents) {
    Eigen::Matrix3d C1;
    directional_deriv1(Psi, V, C1);

    H = (T*H*T.transpose()).eval() + 0.5*(C1 + C1.transpose());
  }
}
#endif
