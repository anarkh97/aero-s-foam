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

  getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, refState, reactions);
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

  getFollowerForce(geomState, elementForce, corotators, (FullSquareMatrix *) NULL, residual, lambda, time, refState, NULL);
}
