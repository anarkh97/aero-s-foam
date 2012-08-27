#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Corotational.d/Corotator.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/mathUtility.h>
#include <Math.d/matrix.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>

#include <Driver.d/GeoSource.h>
#include <Corotational.d/MatNLCorotator.h>

#include <algorithm>

void
Domain::getElemInternalForce(const GeomState &geomState, double time,
                             const GeomState *refState, const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
    const_cast<Corotator &>(elemCorot).getInternalForce(
        const_cast<GeomState *>(refState),
        const_cast<GeomState &>(geomState),
        nodes, elemStiff, elemForce, sinfo.getTimeStep(), time);
}

void
Domain::getElemInternalForce(const GeomState &geomState, double time,
                             const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
  const_cast<Corotator &>(elemCorot).getInternalForce(
      const_cast<GeomState &>(geomState),
      nodes, elemStiff, elemForce, sinfo.getTimeStep(), time);
}

void
Domain::getInternalForce(GeomState &geomState, Vector& elementForce,
                         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time,
                         GeomState *refState, Vector *reactions)
/*******************************************************************
 *
 * Purpose :
 *
 *  Compute element internal force
 *  and assemble element internal force into global internal force.
 *  Also compute follower external force contribution to
 *  residual
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
 *
 *****************************************************************/

{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // mpc needs lambda for nonlinear statics

  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    // Get updated tangent stiffness matrix and element internal force
    if (const Corotator *elemCorot = corotators[iele]) {
      getElemInternalForce(geomState, pseudoTime, refState, *elemCorot, elementForce.data(), kel[iele]);
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

  if(domain->pressureFlag()) {
    double cflg = (sinfo.newmarkBeta == 0.0) ? 0.0 : 1.0;
    double loadFactor = (domain->mftval && sinfo.isDynam()) ? lambda*domain->mftval->getVal(std::max(time,0.0)) : lambda;
    double p0;
    for(int iele = 0; iele < numele;  ++iele) {
      // If there is a zero pressure defined, skip the element
      if((p0 = packedEset[iele]->getPressure()) == 0) continue;

      // Compute (linear) element pressure force in the local coordinates
      elementForce.zero();
      packedEset[iele]->setPressure(p0*loadFactor);
      packedEset[iele]->computePressureForce(nodes, elementForce, &geomState, 1);
      packedEset[iele]->setPressure(p0);

      // Determine the elemental force for the corrotated system
      corotators[iele]->getExternalForce(geomState, nodes, elementForce.data());

      // Assemble element pressure forces into residual force vector
      for(int idof = 0; idof < kel[iele].dim(); ++idof) {
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
  for(int iele = 0; iele < numNeum; ++iele) {
    neum[iele]->dofs(*dsa, edofs);
    elementForce.zero();
    neum[iele]->neumVector(nodes, elementForce, 0, &geomState);
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

  // In order to make the nodal moments non-follower we need to make a correction...
  for(int i = 0; i < numNeuman; ++i) {
    if((nbc[i].type == BCond::Forces || nbc[i].type == BCond::Usdf || nbc[i].type == BCond::Actuators) 
       && (nbc[i].dofnum == 3 || nbc[i].dofnum == 4 || nbc[i].dofnum == 5)) {
      int dofs[3];
      dsa->number(nbc[i].nnum, DofSet::XYZrot, dofs);
      double m0[3] = { 0, 0, 0 }, m[3], r[3], rotvar[3][3];
      m0[nbc[i].dofnum-3] = lambda*mfttFactor*nbc[i].val;
      mat_to_vec(geomState[nbc[i].nnum].R,r);
      pseudorot_var(r, rotvar);
      mat_mult_vec(rotvar,m0,m,1);
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
*/
      // Determine the elemental force for the corrotated system
      corotators[iele]->getExternalForce(geomState, nodes, elementForce.data());

      // Assemble element thermal forces into residual force vector
      for(int idof = 0; idof < kel[iele].dim(); ++idof) {
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
Domain::getWeightedInternalForceOnly(const std::map<int, double> &weights,
                                     GeomState &geomState, Vector& elementForce,
                                     Corotator **corotators, FullSquareMatrix *kel,
                                     Vector &residual, double lambda, double time,
                                     GeomState *refState)
{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // MPC needs lambda for nonlinear statics
  
  for (std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iElem = it->first;

    // Get updated tangent stiffness matrix and element internal force
    if (const Corotator *elementCorot = corotators[iElem]) {
      elementForce.zero();

      FullSquareMatrix &elementStiff = kel[iElem];
      getElemInternalForce(geomState, pseudoTime, refState, *elementCorot, elementForce.data(), elementStiff);
   
      // Apply lumping weight 
      const double lumpingWeight = it->second;
      elementForce *= lumpingWeight;

      const int elemDofCount = elementStiff.dim();
      for(int iDof = 0; iDof < elemDofCount; ++iDof) {
        const int dofId = c_dsa->getRCN((*allDOFs)[iElem][iDof]);
        if (dofId >= 0) {
          residual[dofId] -= elementForce[iDof];
        }
      }
    }
  }
}

