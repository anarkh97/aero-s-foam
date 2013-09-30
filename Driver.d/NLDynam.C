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
#ifdef USE_EIGEN3
#include <Element.d/Function.d/InertialForce.d/InertialForceFunction.h>
#include <Element.d/Function.d/InertialForce.d/InertialForceFunctionExp.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#endif
#include <algorithm>

void
Domain::getElemInternalForce(const GeomState &geomState, double time,
                             const GeomState *refState, const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
    const_cast<Corotator &>(elemCorot).getInternalForce(
        const_cast<GeomState *>(refState),
        const_cast<GeomState &>(geomState),
        nodes, elemStiff, elemForce, domain->solInfo().getTimeStep(), time);
}

void
Domain::getElemInternalForce(const GeomState &geomState, double time,
                             const Corotator &elemCorot,
                             double *elemForce, FullSquareMatrix &elemStiff) {
  const_cast<Corotator &>(elemCorot).getInternalForce(
      const_cast<GeomState &>(geomState),
      nodes, elemStiff, elemForce, domain->solInfo().getTimeStep(), time);
}

void
Domain::getInternalForce(GeomState &geomState, Vector& elementForce,
                         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time,
                         GeomState *refState, Vector *reactions, FullSquareMatrix *mel)
/*******************************************************************
 *
 * Purpose :
 *
 *  Compute element internal force
 *  and assemble element internal force into global internal force.
 *  Also compute follower external force contribution to
 *  residual.
 *  Also for dynamics, compute fictitious (inertial) force
 *  contribution to residual.
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
 *                                 - fictitious force
 *
 *****************************************************************/

{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // mpc needs lambda for nonlinear statics
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  if(elemAdj.empty()) makeElementAdjacencyLists();

  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    // Get updated tangent stiffness matrix and element internal force
    if(corotators[iele] && !solInfo().getNLInfo().linearelastic) {
      getElemInternalForce(geomState, pseudoTime, refState, *corotators[iele], elementForce.data(), kel[iele]);
    }
    // Or, get linear elastic element internal force
    else {
      Vector disp(packedEset[iele]->numDofs());
      getElementDisp(iele, geomState, disp);
      kel[iele].multiply(disp, elementForce, 1.0);
    }

    // Add configuration-dependent external forces and their element stiffness contributions
    getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
                         (corotators[iele]), kel[iele], lambda, time, false, conwep);

    if(domain->solInfo().galerkinPodRom && packedEset[iele]->hasRot() && !solInfo().getNLInfo().linearelastic) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), kel[iele], iele, false);
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

  if(sinfo.isDynam() && mel) getFictitiousForce(geomState, elementForce, kel, residual, time, refState, reactions, mel, false);
}

void
Domain::getWeightedInternalForceOnly(const std::map<int, double> &weights,
                                     GeomState &geomState, Vector& elementForce,
                                     Corotator **corotators, FullSquareMatrix *kel,
                                     Vector &residual, double lambda, double time,
                                     GeomState *refState, FullSquareMatrix *mel)
{
  const double pseudoTime = sinfo.isDynam() ? time : lambda; // MPC needs lambda for nonlinear statics
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  if(elemAdj.empty()) makeElementAdjacencyLists();
  
  for(std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iele = it->first;

    elementForce.zero();
    FullSquareMatrix &elementStiff = kel[iele];
    
    // Get updated tangent stiffness matrix and element internal force
    if(const Corotator *elementCorot = corotators[iele]) {
      getElemInternalForce(geomState, pseudoTime, refState, *elementCorot, elementForce.data(), elementStiff);
    }

    // Add configuration-dependent external forces and their element stiffness contributions
    if(domain->solInfo().reduceFollower) {
      getElemFollowerForce(iele, geomState, elementForce.data(), elementForce.size(),
                           (corotators[iele]), elementStiff, lambda, time, false, conwep);
    }

    if(packedEset[iele]->hasRot()) {
      // Transform element stiffness and force to solve for the increment in the total rotation vector
      transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iele, false);
    }
   
    // Apply lumping weight 
    const double lumpingWeight = it->second;
    elementForce *= lumpingWeight;

    // Assemble element internal force into residual force vector
    const int elemDofCount = elementStiff.dim();
    for(int iDof = 0; iDof < elemDofCount; ++iDof) {
      const int dofId = c_dsa->getRCN((*allDOFs)[iele][iDof]);
      if (dofId >= 0) {
        residual[dofId] -= elementForce[iDof];
      }
    }
  }

  if(!domain->solInfo().reduceFollower) {
    getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, refState, NULL, false);
  }

  if(sinfo.isDynam() && mel) getWeightedFictitiousForceOnly(weights, geomState, elementForce, kel, residual, time, refState, NULL, mel, false);
}

void
Domain::getFictitiousForce(GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel, Vector &residual,
                           double time, GeomState *refState, Vector *reactions, FullSquareMatrix *mel,
                           bool compute_tangents)
{
  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    getElemFictitiousForce(iele, geomState, elementForce.data(), kel[iele], time, refState, mel[iele], compute_tangents);

    // Assemble element force into residual force vector
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
}

void
Domain::getElemFictitiousForce(int iele, GeomState &geomState, double *_f, FullSquareMatrix &kel,
                               double time, GeomState *refState, FullSquareMatrix &mel,
                               bool compute_tangents)
{
#ifdef USE_EIGEN3
  // add the correction to the residual and tangent stiffness due to the inertial effects of 
  // element iele with rotation dofs. Currently implemented for lumped mass matrix only
  // (or more specifically, mass matrices with decoupled rotational and translational diagonal blocks),
  // mass-proportional damping only, and elements with 6 dofs per node.
  if(packedEset[iele]->hasRot() && !packedEset[iele]->isSpring()) {
    int *nodes = packedEset[iele]->nodes();
    int numNodes = packedEset[iele]->numNodes() - packedEset[iele]->numInternalNodes();
    for(int i=0; i<numNodes; ++i) {

      Eigen::Matrix3d M;
      M << mel[6*i+3][6*i+3], mel[6*i+3][6*i+4], mel[6*i+3][6*i+5],
           mel[6*i+4][6*i+3], mel[6*i+4][6*i+4], mel[6*i+4][6*i+5],
           mel[6*i+5][6*i+3], mel[6*i+5][6*i+4], mel[6*i+5][6*i+5];
      if((M.array() == 0).all()) continue;

      Eigen::Vector3d f;
      Eigen::Matrix3d K;
      getNodeFictitiousForce(nodes[i], geomState, time, refState, M, f, K, compute_tangents);

      for(int j = 0; j < 3; ++j) {
        _f[6*i+3+j] += f[j];
        if(compute_tangents) {
          for(int k = 0; k < 3; ++k)
            kel[6*i+3+j][6*i+3+k] += K(j,k);
        }
      }
    }
    delete [] nodes;
  }

  if(iele >= elemAdj.size()) return;

  // treatment of discrete inertias adjacent to the element
  for(std::vector<std::pair<DMassData*,std::vector<int> > >::iterator it = elemAdj[iele].dimass.begin(); it != elemAdj[iele].dimass.end(); ++it) {

    DMassData *current = it->first;
    std::vector<int> &eledofs = it->second;

    int idof = current->dof;
    int jdof = (current->jdof > -1) ? current->jdof : idof;

    Eigen::Matrix3d M = Eigen::Matrix3d::Zero(); 
    M(idof-3,jdof-3) = current->diMass;
    if(idof != jdof) M(jdof-3,idof-3) = current->diMass;

    Eigen::Vector3d f;
    Eigen::Matrix3d K;
    getNodeFictitiousForce(current->node, geomState, time, refState, M, f, K, compute_tangents);

    for(int j = 0; j < 3; ++j) {
      _f[eledofs[j]] += f[j];
      if(compute_tangents) {
        for(int k = 0; k < 3; ++k)
          kel[eledofs[j]][eledofs[k]] += K(j,k);
      }
    }
  }
#endif
}

#ifdef USE_EIGEN3
void
Domain::getNodeFictitiousForce(int inode, GeomState &geomState, double time, GeomState *refState, Eigen::Matrix3d &M,
                               Eigen::Vector3d &f, Eigen::Matrix3d &K, bool compute_tangents)
{
  double &beta = sinfo.newmarkBeta,
         &gamma = sinfo.newmarkGamma,
         &alphaf = sinfo.newmarkAlphaF,
         &alpham = sinfo.newmarkAlphaM,
          dt = domain->solInfo().getTimeStep();

  Eigen::Matrix3d R, T, Tdot, R_n;
  R << geomState[inode].R[0][0], geomState[inode].R[0][1], geomState[inode].R[0][2],
       geomState[inode].R[1][0], geomState[inode].R[1][1], geomState[inode].R[1][2],
       geomState[inode].R[2][0], geomState[inode].R[2][1], geomState[inode].R[2][2];
  if(refState) {
    R_n << (*refState)[inode].R[0][0],(*refState)[inode].R[0][1], (*refState)[inode].R[0][2],
           (*refState)[inode].R[1][0],(*refState)[inode].R[1][1], (*refState)[inode].R[1][2],
           (*refState)[inode].R[2][0],(*refState)[inode].R[2][1], (*refState)[inode].R[2][2];
  }

  Eigen::Vector3d Psi, Psidot, V, A, Psi_n, V_n, A_n;
  Psi << geomState[inode].theta[0], geomState[inode].theta[1], geomState[inode].theta[2];
  V << geomState[inode].v[3], geomState[inode].v[4], geomState[inode].v[5];
  A << geomState[inode].a[3], geomState[inode].a[4], geomState[inode].a[5];
  if(refState) {
    Psi_n << (*refState)[inode].theta[0], (*refState)[inode].theta[1], (*refState)[inode].theta[2];
    V_n << (*refState)[inode].v[3], (*refState)[inode].v[4], (*refState)[inode].v[5];
    A_n << (*refState)[inode].a[3], (*refState)[inode].a[4], (*refState)[inode].a[5];
  }

  if(beta == 0) { // compute the fictitious force for explicit central difference
    // V is either the convected angular velocity at t^{n+1/2} for FOM or ROM model II or model III,
    //   or the convected angular velocity at current snapshot after projection for explicit ROM "training"
    if(domain->solInfo().galerkinPodRom || domain->solInfo().samplingPodRom) {
      tangential_transf(Psi, T);
      Psidot = T.inverse()*V;
      tangential_transf_dot(Psi, Psidot, Tdot);
      f = T.transpose()*(M*Tdot*Psidot + V.cross(M*V));
    }
    else {
      f = R*V.cross(M*V);
    }

    // TODO: compute tangents for explict (for critical timestep estimate)
    if(compute_tangents) K.setZero();
  }
  else { // compute the fictitious force for implicit generalized-alpha
    if(domain->solInfo().samplingPodRom) {
      // V and A are the convected angular velocity and acceleration at current snapshot after projection
      tangential_transf(Psi, T);
      f = (T.transpose() - Eigen::Matrix3d::Identity())*M*(A+sinfo.alphaDamp*V) + T.transpose()*V.cross(M*V);
    }
    else {
      Eigen::Vector3d f0;
      if(time == 0) {
        Psi.setZero();
        f0.setZero();
        compute_tangents = false;
      }
      else if(domain->solInfo().galerkinPodRom) {
        Eigen::Vector3d incd = Psi - Psi_n;
        // compute the total angular velocity at t^{n+1-alphaf}
        V = gamma/(dt*beta)*incd + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
        // compute the total angular acceleration at t^{n+1-alpham}
        A = (1-alpham)/(dt*dt*beta*(1-alphaf))*incd - (1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;
        f0 = M*A;
      }
      else {
        Eigen::Vector3d incd;
        Eigen::Matrix3d dR = R_n.transpose()*R;
        mat_to_vec(dR, incd);
        // compute the convected angular velocity at t^{n+1-alphaf}
        V = gamma/(dt*beta)*incd + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
        // compute the convected angular acceleration at t^{n+1-alpham}
        A = (1-alpham)/(dt*dt*beta*(1-alphaf))*incd - (1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;
      }

      // compute the fictitious force and the correction to the inertial+viscous force computed in probDesc->formRHScorrector which is (M*A + C*V)
      if(domain->solInfo().galerkinPodRom) {
        // the correct inertia+viscous force is T*R*(M*A + C*V + V.cross(M*V)). note T*R = T.transpose()
        tangential_transf(Psi, T);
        tangential_transf_dot(Psi, V, Tdot);
        f = T.transpose()*( M*(T*A + Tdot*V) + (T*V).cross(M*T*V) ) - f0;
      }
      else {
        // the correct inertia+viscous force is R*(M*A + C*V + V.cross(M*V))
        f = (R - Eigen::Matrix3d::Identity())*M*(A+sinfo.alphaDamp*V) + R*V.cross(M*V);
      }

      if(compute_tangents) { // tangent stiffness contribution of the fictitious force and correct linearization of rotary inertia+viscous force

        Eigen::Matrix<double,3,1> q;
        if(domain->solInfo().galerkinPodRom) {
          Eigen::Array<double,24,1> dconst;
          Eigen::Array<int,0,1> iconst;
          dconst << M(0,0), M(0,1), M(0,2), M(1,0), M(1,1), M(1,2), M(2,0), M(2,1), M(2,2),
                    A_n[0], A_n[1], A_n[2],
                    V_n[0], V_n[1], V_n[2],
                    Psi_n[0], Psi_n[1], Psi_n[2],
                    beta, gamma, alphaf, alpham, dt, sinfo.alphaDamp;

          // evaluate the jacobian of the inertial+viscous force
          Simo::Jacobian<double,InertialForceFunctionExp> dFdq(dconst,iconst);
          q << geomState[inode].theta[0], geomState[inode].theta[1], geomState[inode].theta[2];
          K = dFdq(q, time);
        }
        else {
          Eigen::Array<double,39,1> dconst;
          Eigen::Array<int,0,1> iconst;
          dconst << M(0,0), M(0,1), M(0,2), M(1,0), M(1,1), M(1,2), M(2,0), M(2,1), M(2,2),
                    A_n[0], A_n[1], A_n[2],
                    V_n[0], V_n[1], V_n[2],
                    R_n(0,0), R_n(0,1), R_n(0,2), R_n(1,0), R_n(1,1), R_n(1,2), R_n(2,0), R_n(2,1), R_n(2,2),
                    R(0,0), R(0,1), R(0,2), R(1,0), R(1,1), R(1,2), R(2,0), R(2,1), R(2,2),
                    beta, gamma, alphaf, alpham, dt, sinfo.alphaDamp;

          // evaluate the jacobian of the inertial+viscous force
          Simo::Jacobian<double,InertialForceFunction> dFdq(dconst,iconst);
          q = Eigen::Vector3d::Zero();
          K = dFdq(q, time);
        }

        // subtract the part which is added to the dynamic tangent stiffness in probDesc->reBuild
        K -= ((1-alpham)/((1-alphaf)*(dt*dt*beta)) + gamma/(dt*beta)*sinfo.alphaDamp)*M;
      }
    }
  }
}
#endif

void
Domain::getWeightedFictitiousForceOnly(const std::map<int, double> &weights, GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel,
                                       Vector &residual, double time, GeomState *refState, Vector *reactions,
                                       FullSquareMatrix *mel, bool compute_tangents)
{
  for (std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iele = it->first;
    const double lumpingWeight = it->second;

    elementForce.zero();
    FullSquareMatrix kel2(kel[iele].dim());
    kel2.zero();

    getElemFictitiousForce(iele, geomState, elementForce.data(), kel2, time, refState, mel[iele], compute_tangents);

    elementForce *= lumpingWeight;

    if(compute_tangents) {
      kel2 *= lumpingWeight;
      kel[iele] += kel2;
    }

    // Assemble element force into residual force vector
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
}
