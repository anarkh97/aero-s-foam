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
#include <Element.d/Dimass.d/InertialForceFunction.h>
#include <Element.d/Dimass.d/InertialForceFunctionExp.h>
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
      if(domain->solInfo().galerkinPodRom && packedEset[iele]->hasRot()) {
        transformElemStiffAndForce(geomState, elementForce.data(), kel[iele], iele, false);
      }
    }
    // Compute k and internal force for an element with x translation (or temperature) dofs
    else if(solInfo().soltyp == 2) {
      kel[iele].zero();
      Vector temp(packedEset[iele]->numNodes());
      int *nn = packedEset[iele]->nodes();
      for(int i=0; i<packedEset[iele]->numNodes(); ++i) {
        temp[i] = geomState[nn[i]].x;
      }
      kel[iele] = packedEset[iele]->stiffness(nodes, kel[iele].data());
      kel[iele].multiply(temp, elementForce, 1.0); // elementForce = kel*temp
      delete [] nn;
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

  getFollowerForce(geomState, elementForce, corotators, kel, residual, lambda, time, refState, reactions, false);

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
  
  for (std::map<int, double>::const_iterator it = weights.begin(), it_end = weights.end(); it != it_end; ++it) {
    const int iElem = it->first;
    
    // Get updated tangent stiffness matrix and element internal force
    if (const Corotator *elementCorot = corotators[iElem]) {
      elementForce.zero();

      FullSquareMatrix &elementStiff = kel[iElem];
      getElemInternalForce(geomState, pseudoTime, refState, *elementCorot, elementForce.data(), elementStiff);
      if (domain->solInfo().galerkinPodRom && packedEset[iElem]->hasRot()) {
        transformElemStiffAndForce(geomState, elementForce.data(), elementStiff, iElem, false);
      }
   
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

  getFollowerForce(geomState, elementForce, corotators, (FullSquareMatrix *) NULL, residual, lambda, time, refState, NULL, false);

  if(sinfo.isDynam() && mel) getWeightedFictitiousForceOnly(weights, geomState, elementForce, kel, residual, time, refState, NULL, mel, false);
}

void
Domain::getFictitiousForce(GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel, Vector &residual,
                           double time, GeomState *refState, Vector *reactions, FullSquareMatrix *mel,
                           bool compute_tangents)
{
  // TODO: consider case of t=0
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

  getDMassFictitiousForce(geomState, kel, residual, time, refState, reactions, compute_tangents);
}

void
Domain::getElemFictitiousForce(int iele, GeomState &geomState, double *_f, FullSquareMatrix &kel,
                               double time, GeomState *refState, FullSquareMatrix &mel,
                               bool compute_tangents)
{
#ifdef USE_EIGEN3
  double &beta = sinfo.newmarkBeta,
         &gamma = sinfo.newmarkGamma,
         &alphaf = sinfo.newmarkAlphaF,
         &alpham = sinfo.newmarkAlphaM,
          dt = domain->solInfo().getTimeStep();

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

      Eigen::Matrix3d R, K, T, Tdot;
      R << geomState[nodes[i]].R[0][0], geomState[nodes[i]].R[0][1], geomState[nodes[i]].R[0][2],
           geomState[nodes[i]].R[1][0], geomState[nodes[i]].R[1][1], geomState[nodes[i]].R[1][2],
           geomState[nodes[i]].R[2][0], geomState[nodes[i]].R[2][1], geomState[nodes[i]].R[2][2];

      Eigen::Vector3d f, Psi, V, A;

      if(beta == 0) { // compute the fictitious force for explicit central difference        
        // V is either the convected angular velocity at t^{n+1/2} for FOM or ROM model II or model III,
        //   or the convected angular velocity at current snapshot after projection for explicit ROM "training"
        V << geomState[nodes[i]].v[3], geomState[nodes[i]].v[4], geomState[nodes[i]].v[5];
        if(domain->solInfo().galerkinPodRom) { // ROM
          mat_to_vec(R, Psi);
          tangential_transf(Psi, T);
          Eigen::Vector3d Psidot;
          Psidot = T.inverse()*V;
          Eigen::Matrix3d Tdot;
          tangential_transf_dot(Psi, Psidot, Tdot);
          f = T.transpose()*(M*Tdot*Psidot + V.cross(M*V));
        }
        else { // FOM
          f = R*V.cross(M*V);
        }

        // TODO: compute tangents for explict (for critical timestep estimate)
        if(compute_tangents) K.setZero();
      }
      else { // compute the fictitious force for implicit generalized-alpha
        if(domain->solInfo().samplingPodRom) { // implicit ROM "training"
          // V and A are the convected angular velocity and acceleration at current snapshot after projection
          V << geomState[nodes[i]].v[3], geomState[nodes[i]].v[4], geomState[nodes[i]].v[5];
          A << geomState[nodes[i]].a[3], geomState[nodes[i]].a[4], geomState[nodes[i]].a[5];
          //mat_to_vec(R, Psi);
          Psi << geomState[nodes[i]].theta[0], geomState[nodes[i]].theta[1], geomState[nodes[i]].theta[2];
          tangential_transf(Psi, T);
          f = (T.transpose() - Eigen::Matrix3d::Identity())*M*(A+sinfo.alphaDamp*V) + T.transpose()*V.cross(M*V);
        }
        else {
          Eigen::Vector3d V_n, A_n, Psi_n, f0; // convected angular velocity and acceleration at t^{n}
          Eigen::Matrix3d R_n;
          if(time == domain->solInfo().initialTime) {
            V << geomState[nodes[i]].v[3], geomState[nodes[i]].v[4], geomState[nodes[i]].v[5];
            A << geomState[nodes[i]].a[3], geomState[nodes[i]].a[4], geomState[nodes[i]].a[5];
            Psi.setZero();
            f0.setZero();
            compute_tangents = false;
          }
          else if(domain->solInfo().galerkinPodRom) {
            V_n << (*refState)[nodes[i]].v[3], (*refState)[nodes[i]].v[4], (*refState)[nodes[i]].v[5];
            A_n << (*refState)[nodes[i]].a[3], (*refState)[nodes[i]].a[4], (*refState)[nodes[i]].a[5];
            Psi_n << (*refState)[nodes[i]].theta[0], (*refState)[nodes[i]].theta[1], (*refState)[nodes[i]].theta[2];
            R_n << (*refState)[nodes[i]].R[0][0], (*refState)[nodes[i]].R[0][1], (*refState)[nodes[i]].R[0][2],
                   (*refState)[nodes[i]].R[1][0], (*refState)[nodes[i]].R[1][1], (*refState)[nodes[i]].R[1][2],
                   (*refState)[nodes[i]].R[2][0], (*refState)[nodes[i]].R[2][1], (*refState)[nodes[i]].R[2][2];

            Psi << geomState[nodes[i]].theta[0], geomState[nodes[i]].theta[1], geomState[nodes[i]].theta[2];

            Eigen::Vector3d incd = Psi - Psi_n;
            // compute the total angular velocity at t^{n+1-alphaf}
            V = gamma/(dt*beta)*incd + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
            // compute the total angular acceleration at t^{n+1-alpham}
            A = (1-alpham)/(dt*dt*beta*(1-alphaf))*incd - (1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;
            f0 = M*A;
/*
            if(Psi.squaredNorm() != 0 && Psi.norm() > M_PI) {
              Eigen::Matrix3d B = complement_transf(Psi);
              Eigen::Matrix3d Bdot = complement_transf_dot(Psi,V);
              A = (B*A + Bdot*V).eval();
              V = (B*V).eval();
              Psi = complement_rot_vec(Psi);
            }
*/
          }
          else {
            V_n << (*refState)[nodes[i]].v[3], (*refState)[nodes[i]].v[4], (*refState)[nodes[i]].v[5];
            A_n << (*refState)[nodes[i]].a[3], (*refState)[nodes[i]].a[4], (*refState)[nodes[i]].a[5];
            R_n << (*refState)[nodes[i]].R[0][0], (*refState)[nodes[i]].R[0][1], (*refState)[nodes[i]].R[0][2],
                   (*refState)[nodes[i]].R[1][0], (*refState)[nodes[i]].R[1][1], (*refState)[nodes[i]].R[1][2],
                   (*refState)[nodes[i]].R[2][0], (*refState)[nodes[i]].R[2][1], (*refState)[nodes[i]].R[2][2];

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
            //mat_to_vec(R, Psi); // XXX
            tangential_transf(Psi, T);
            //f = (T.transpose() - Eigen::Matrix3d::Identity())*M*(A+sinfo.alphaDamp*V) + T.transpose()*V.cross(M*V);
            tangential_transf_dot(Psi, V, Tdot);
            //f = T.transpose()*M*(T*A + Tdot*V) + V.cross(T.inverse()*M*T*V) - M*A;
            f = T.transpose()*( M*(T*A + Tdot*V) + (T*V).cross(M*T*V) ) - f0;
          }
          else {
            // the correct inertia+viscous force is R*(M*A + C*V + V.cross(M*V))
            f = (R - Eigen::Matrix3d::Identity())*M*(A+sinfo.alphaDamp*V) + R*V.cross(M*V);
          }

          if(compute_tangents) { // tangent stiffness contribution of the fictitious force and correct linearization of rotary inertia+viscous force

            if(domain->solInfo().galerkinPodRom) {
              Eigen::Array<double,24,1> dconst;
              Eigen::Array<int,0,1> iconst;
              dconst << M(0,0), M(0,1), M(0,2), M(1,0), M(1,1), M(1,2), M(2,0), M(2,1), M(2,2),
                        A_n[0], A_n[1], A_n[2],
                        V_n[0], V_n[1], V_n[2],
                        Psi_n[0], Psi_n[1], Psi_n[2], 
                        beta, gamma, alphaf, alpham, dt, sinfo.alphaDamp;

              // evaluate the jacobian of the inertial+viscous force
              VectorValuedFunctionJacobian<double,InertialForceFunctionExp> dFdq(dconst,iconst,time);
              Eigen::Matrix<double,9,1> jacF;
              Eigen::Matrix<double,3,1> q;
              q << geomState[nodes[i]].theta[0], geomState[nodes[i]].theta[1], geomState[nodes[i]].theta[2];
              dFdq(q, jacF);
              //if(q.squaredNorm() > 0 && q.norm() > M_PI) std::cerr << "#1 jacF = " << jacF.transpose() << std::endl;

              for(int j = 0; j < 3; ++j)
                for(int k = 0; k < 3; ++k)
                  K(j,k) = jacF[j+k*3];
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
              VectorValuedFunctionJacobian<double,InertialForceFunction> dFdq(dconst,iconst,time);
              Eigen::Matrix<double,9,1> jacF;
              Eigen::Vector3d q = Eigen::Vector3d::Zero();
              dFdq(q, jacF);

              for(int j = 0; j < 3; ++j)
                for(int k = 0; k < 3; ++k)
                  K(j,k) = jacF[j+k*3];
            }
/*
            if(domain->solInfo().galerkinPodRom) {
              // transform K (note: Psi and T have already been computed...)
              Eigen::Vector3d F = R*(M*(A + sinfo.alphaDamp*V) + V.cross(M*V));
              Eigen::Matrix3d C1;
              directional_deriv1(Psi, F, C1);
              K = (T*K*T.transpose()).eval() + 0.5*(C1 + C1.transpose());
            }
*/

            // subtract the part which is added to the dynamic tangent stiffness in probDesc->reBuild
            K -= ((1-alpham)/((1-alphaf)*(dt*dt*beta)) + gamma/(dt*beta)*sinfo.alphaDamp)*M;
          }
        }
      }

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
#endif
}


void
Domain::getDMassFictitiousForce(GeomState &geomState, FullSquareMatrix *kel, Vector &residual, double time,
                                GeomState *refState, Vector *reactions, bool compute_tangents)
{
#ifdef USE_EIGEN3
  double &beta = sinfo.newmarkBeta,
         &gamma = sinfo.newmarkGamma,
         &alphaf = sinfo.newmarkAlphaF,
         &alpham = sinfo.newmarkAlphaM,
          dt = domain->solInfo().getTimeStep();

  // treatment of discrete inertias
  if(firstDiMass != NULL) {
    DMassData *current = firstDiMass;
    while(current != 0) {
      int idof = current->dof;
      int jdof = (current->jdof > -1) ? current->jdof : idof;
      if((idof == 3 || idof == 4 || idof == 5) && (jdof == 3 || jdof == 4 || jdof == 5)) {

        Eigen::Matrix3d M = Eigen::Matrix3d::Zero(); 
        M(idof-3,jdof-3) = current->diMass;
        if(idof != jdof) M(jdof-3,idof-3) = current->diMass;

        Eigen::Matrix3d R, K, T, Tdot;
        R << geomState[current->node].R[0][0], geomState[current->node].R[0][1], geomState[current->node].R[0][2],
             geomState[current->node].R[1][0], geomState[current->node].R[1][1], geomState[current->node].R[1][2],
             geomState[current->node].R[2][0], geomState[current->node].R[2][1], geomState[current->node].R[2][2];

        Eigen::Vector3d f, Psi, V, A;
        int dofs[3];
        dsa->number(current->node, DofSet::XYZrot, dofs);

        if(beta == 0) { // compute the fictitious force for explicit central difference
        // V is either the convected angular velocity at t^{n+1/2} for FOM or ROM model II or model III,
        //   or the convected angular velocity at current snapshot after projection for explicit ROM "training"
          V << geomState[current->node].v[3], geomState[current->node].v[4], geomState[current->node].v[5];
          if(domain->solInfo().galerkinPodRom) {
            mat_to_vec(R, Psi);
            tangential_transf(Psi, T);
            Eigen::Vector3d Psidot;
            Psidot = T.inverse()*V;
            Eigen::Matrix3d Tdot;
            tangential_transf_dot(Psi, Psidot, Tdot);
            f = T.transpose()*(M*Tdot*Psidot + V.cross(M*V));
          }
          else { // FOM
            f = R*V.cross(M*V);
          }

          // TODO: compute tangents for explict (for critical timestep estimate)
         if(compute_tangents) K.setZero();
        }
        else { // compute the fictitious force for implicit generalized-alpha

          Eigen::Vector3d V_n, A_n, Psi_n, f0; // convected angular velocity and acceleration at t^{n}
          Eigen::Matrix3d R_n;
          if(time == domain->solInfo().initialTime) {
            V << geomState[current->node].v[3], geomState[current->node].v[4], geomState[current->node].v[5];
            A << geomState[current->node].a[3], geomState[current->node].a[4], geomState[current->node].a[5];
            Psi.setZero();
            f0.setZero();
            compute_tangents = false;
          }
          else if(domain->solInfo().galerkinPodRom) {
            V_n << (*refState)[current->node].v[3], (*refState)[current->node].v[4], (*refState)[current->node].v[5];
            A_n << (*refState)[current->node].a[3], (*refState)[current->node].a[4], (*refState)[current->node].a[5];
            Psi_n << (*refState)[current->node].theta[0], (*refState)[current->node].theta[1], (*refState)[current->node].theta[2];
            R_n << (*refState)[current->node].R[0][0],(*refState)[current->node].R[0][1], (*refState)[current->node].R[0][2],
                   (*refState)[current->node].R[1][0],(*refState)[current->node].R[1][1], (*refState)[current->node].R[1][2],
                   (*refState)[current->node].R[2][0],(*refState)[current->node].R[2][1], (*refState)[current->node].R[2][2];

            Psi << geomState[current->node].theta[0], geomState[current->node].theta[1], geomState[current->node].theta[2];

            Eigen::Vector3d incd = Psi - Psi_n;
            // compute the total angular velocity at t^{n+1-alphaf}
            V = gamma/(dt*beta)*incd + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
            // compute the total angular acceleration at t^{n+1-alpham}
            A = (1-alpham)/(dt*dt*beta*(1-alphaf))*incd - (1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;
            f0 = M*A;
/*
            if(Psi.squaredNorm() != 0 && Psi.norm() > M_PI) {
              Eigen::Matrix3d B = complement_transf(Psi);
              Eigen::Matrix3d Bdot = complement_transf_dot(Psi,V);
              A = (B*A + Bdot*V).eval();
              V = (B*V).eval();
              Psi = complement_rot_vec(Psi);
            }
*/
          }
          else {
            V_n << (*refState)[current->node].v[3], (*refState)[current->node].v[4], (*refState)[current->node].v[5];
            A_n << (*refState)[current->node].a[3], (*refState)[current->node].a[4], (*refState)[current->node].a[5];
            R_n << (*refState)[current->node].R[0][0],(*refState)[current->node].R[0][1], (*refState)[current->node].R[0][2],
                   (*refState)[current->node].R[1][0],(*refState)[current->node].R[1][1], (*refState)[current->node].R[1][2],
                   (*refState)[current->node].R[2][0],(*refState)[current->node].R[2][1], (*refState)[current->node].R[2][2];

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
            //mat_to_vec(R, Psi); // XXX
            tangential_transf(Psi, T);
            tangential_transf_dot(Psi, V, Tdot);
            //f = (T.transpose() - Eigen::Matrix3d::Identity())*M*(A+sinfo.alphaDamp*V) + T.transpose()*V.cross(M*V);
            //f = T.transpose()*M*(T*A + Tdot*V) + V.cross(T.inverse()*M*T*V) - M*A;
            f = T.transpose()* ( M*(T*A + Tdot*V) + (T*V).cross(M*T*V) ) - f0;
          }
          else {
            // the correct inertia+viscous force is R*(M*A + C*V + V.cross(M*V))
            f = (R - Eigen::Matrix3d::Identity())*M*(A+sinfo.alphaDamp*V) + R*V.cross(M*V);
          }

          if(compute_tangents) { // tangent stiffness contribution of the fictitious force and correct linearization of rotary inertia+viscous force

            if(domain->solInfo().galerkinPodRom) {
              Eigen::Array<double,24,1> dconst;
              Eigen::Array<int,0,1> iconst;
              dconst << M(0,0), M(0,1), M(0,2), M(1,0), M(1,1), M(1,2), M(2,0), M(2,1), M(2,2),
                        A_n[0], A_n[1], A_n[2],
                        V_n[0], V_n[1], V_n[2],
                        Psi_n[0], Psi_n[1], Psi_n[2],
                        beta, gamma, alphaf, alpham, dt, sinfo.alphaDamp;

              // evaluate the jacobian of the inertial+viscous force
              VectorValuedFunctionJacobian<double,InertialForceFunctionExp> dFdq(dconst,iconst,time);
              Eigen::Matrix<double,9,1> jacF;
              Eigen::Matrix<double,3,1> q;
              q << geomState[current->node].theta[0], geomState[current->node].theta[1], geomState[current->node].theta[2];
              dFdq(q, jacF);
              //if(q.squaredNorm() > 0 && q.norm() > M_PI) std::cerr << "#2 jacF = " << jacF.transpose() << std::endl;

              for(int i = 0; i < 3; ++i)
                for(int j = 0; j < 3; ++j)
                  K(i,j) = jacF[i+j*3];
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
              VectorValuedFunctionJacobian<double,InertialForceFunction> dFdq(dconst,iconst,time);
              Eigen::Matrix<double,9,1> jacF;
              Eigen::Vector3d q = Eigen::Vector3d::Zero();
              dFdq(q, jacF);

              for(int i = 0; i < 3; ++i)
                for(int j = 0; j < 3; ++j)
                  K(i,j) = jacF[i+j*3];
            }
/*
            if(domain->solInfo().galerkinPodRom) {
              // transform K (note: Psi and T have already been computed...)
              Eigen::Vector3d F = R*(M*(A + sinfo.alphaDamp*V) + V.cross(M*V));
              Eigen::Matrix3d C1;
              directional_deriv1(Psi, F, C1);
              K = (T*K*T.transpose()).eval() + 0.5*(C1 + C1.transpose());
            }
*/
            // subtract the part which is added to the dynamic tangent stiffness in probDesc->reBuild
            K -= ((1-alpham)/((1-alphaf)*(dt*dt*beta)) + gamma/(dt*beta)*sinfo.alphaDamp)*M;

            for(int inode = 0; inode < nodeToElem->num(current->node); ++inode) { // loop over the elements attached to the node
                                                                                  // at which the discrete mass is located
              int iele = (*nodeToElem)[current->node][inode];
              int eledofs[3] = { -1, -1, -1 };
              for(int j = 0; j < 3; ++j) {
                for(int k = 0; k < allDOFs->num(iele); ++k)
                  if(dofs[j] == (*allDOFs)[iele][k]) { eledofs[j] = k; break; }
              }
              if(eledofs[0] != -1 && eledofs[1] != -1 && eledofs[2] != -1) {
                // found an element with the 3 rotation dofs of current->node so we can add the inertial stiffness
                // contribution of the discrete mass to the tangent stiffness matrix of this element
                for(int j = 0; j < 3; ++j)
                  for(int k = 0; k < 3; ++k)
                    kel[iele][eledofs[j]][eledofs[k]] += K(j,k);
                break;
              }
            }
          }
        }
  
        // assemble f into the global residual
        for(int j = 0; j < 3; ++j) {
          int uDofNum = c_dsa->getRCN(dofs[j]);
          if(uDofNum >= 0)
            residual[uDofNum] -= f[j];
          else if(reactions) {
            int cDofNum = c_dsa->invRCN(dofs[j]);
            if(cDofNum >= 0)
              (*reactions)[cDofNum] += f[j];
          }
        }
      }
      current = current->next;
    }
  }
#endif
}

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

  // XXX treatment of discrete inertias
  getDMassFictitiousForce(geomState, kel, residual, time, refState, reactions, compute_tangents);
}
