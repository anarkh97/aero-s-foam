#include <Element.d/MpcElement.d/MpcElement.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/Domain.h>
#include <Driver.d/EFrameData.h>

extern Domain * domain;

MpcElement::MpcElement(int _nNodes, DofSet nodalDofs, int* _nn)
 : nNodes(_nNodes), LMPCons(0, 0.0)
{
  // this constructor is for a constraint involving the same DofSet on each node
  nn = new int[nNodes+1];
  for(int i = 0; i < nNodes; ++i)
    nn[i] = _nn[i];
  nn[nNodes] = -1;
  addTerms(nodalDofs);
}

void
MpcElement::addTerms(DofSet nodalDofs)
{
  terms.clear();
  nterms = 0;
  for(int i = 0; i < nNodes; ++i) {
    for(int j = 0; j < nodalDofs.count(); ++j)
      for(int k = 0; k < 11; ++k) 
        if(nodalDofs.contains(1 << k)) {
          LMPCTerm t(nn[i], k, 0.0);
          addterm(&t);
        }
  }
}

MpcElement::MpcElement(int _nNodes, DofSet *nodalDofs, int* _nn)
 : nNodes(_nNodes), LMPCons(0, 0.0)
{
  // this constructor is for a constraint involving a different DofSet on each nodes
  nn = new int[nNodes+1];
  for(int i = 0; i < nNodes; ++i)
    nn[i] = _nn[i];
  nn[nNodes] = -1;
  addTerms(nodalDofs);
}

void
MpcElement::addTerms(DofSet *nodalDofs)
{
  terms.clear();
  nterms = 0;
  for(int i = 0; i < nNodes; ++i) {
    for(int j = 0; j < nodalDofs[i].count(); ++j)
      for(int k = 0; k < 11; ++k)
        if(nodalDofs[i].contains(1 << k)) {
          LMPCTerm t(nn[i], k, 0.0);
          addterm(&t);
        }
  }
}

MpcElement::MpcElement(LMPCons *mpc, bool nlflag)
 : LMPCons(*mpc)
{
  // count the number of nodes. also fill node numbers array
  nn = new int[mpc->nterms + 1];
  nNodes = 0;
  for(int i = 0; i < mpc->nterms; i++) {
    int nnum = (mpc->terms)[i].nnum;
    bool found = false;
    for(int j = 0; j < nNodes; ++j) {
      if(nnum == nn[j]) { found = true; break; }
    }
    if(!found) nn[nNodes++] = nnum;
  }
  nn[nNodes] = -1;

  if(nlflag) {
    // store the original coefficients of the rotation dofs terms
    for(int i = 0; i < nterms; ++i) {
      if(terms[i].dofnum == 3 || terms[i].dofnum == 4 || terms[i].dofnum == 5) {
        int nodeNumber = terms[i].nnum;
        std::map<int,std::vector<double> >::iterator it = rotation_coefs.find(nodeNumber);
        if(it == rotation_coefs.end()) {
          std::vector<double> v(3); v[0] = v[1] = v[2] = 0;
          v[terms[i].dofnum-3] = terms[i].coef.r_value;
          rotation_coefs.insert(it, std::pair<int,std::vector<double> >(nodeNumber,v));
          std::vector<int> k(3); k[0] = k[1] = k[2] = -1;
          k[terms[i].dofnum-3] = i;
          rotation_indices[nodeNumber] = k;
        }
        else {
          it->second[terms[i].dofnum-3] = terms[i].coef.r_value;
          rotation_indices[it->first][terms[i].dofnum-3] = i;
        }
      }
    }
    // insert any missing terms (all three rotation dofs should be present even if original coefs are zero)
    for(std::map<int,std::vector<int> >::iterator it = rotation_indices.begin(); it != rotation_indices.end(); ++it) {
      for(int j=0; j<3; ++j)
        if(it->second[j] == -1) {
          LMPCTerm t(it->first, j+3, 0.0);
          addterm(&t);
          it->second[j] = nterms-1;
        }
    }
  }
}

MpcElement::~MpcElement()
{
  delete [] nn;
}

int
MpcElement::getNumMPCs()
{
  if(prop->penalty == 0) return 1; // lagrange multiplier and elimination 
  else return 0; // penalty and augmented lagrangian
}

LMPCons** 
MpcElement::getMPCs()
{
  LMPCons **mpcs = new LMPCons*[1];
  mpcs[0] = new LMPCons(*this);
  return mpcs;
}

int
MpcElement::numInternalNodes()
{
  // method of multipliers and augmented lagrangian have an internal node to store the lagrange multiplier
  return (prop->lagrangeMult) ? 1 : 0;
}

void
MpcElement::setInternalNodes(int* in)
{
  if(prop->lagrangeMult)
    nn[nNodes] = in[0];
}

int
MpcElement::numNodes()
{
  return (prop->lagrangeMult) ? nNodes+1 : nNodes;
}

void
MpcElement::renum(int* table)
{
  for(int i = 0; i < numNodes(); ++i)
    if(nn[i] > -1)
      nn[i] = table[nn[i]];
  for(int i = 0; i < nterms; ++i)
    terms[i].nnum = table[terms[i].nnum];

  std::map<int,std::vector<int> > rotation_indices_renum;
  for(std::map<int,std::vector<int> >::iterator it = rotation_indices.begin(); it != rotation_indices.end(); ++it) {
    rotation_indices_renum[table[it->first]] = it->second;
  }
  rotation_indices = rotation_indices_renum;

  std::map<int,std::vector<double> > rotation_coefs_renum;
  for(std::map<int,std::vector<double> >::iterator it = rotation_coefs.begin(); it != rotation_coefs.end(); ++it) {
    rotation_coefs_renum[table[it->first]] = it->second;
  }
  rotation_coefs = rotation_coefs_renum;
}

void
MpcElement::renum(EleRenumMap& table)
{
  for(int i = 0; i < numNodes(); ++i)
    if(nn[i] > -1)
      nn[i] = table[nn[i]];
  for(int i = 0; i < nterms; ++i)
    terms[i].nnum = table[terms[i].nnum];

  std::map<int,std::vector<int> > rotation_indices_renum;
  for(std::map<int,std::vector<int> >::iterator it = rotation_indices.begin(); it != rotation_indices.end(); ++it) {
    rotation_indices_renum[table[it->first]] = it->second;
  }
  rotation_indices = rotation_indices_renum;

  std::map<int,std::vector<double> > rotation_coefs_renum;
  for(std::map<int,std::vector<double> >::iterator it = rotation_coefs.begin(); it != rotation_coefs.end(); ++it) {
    rotation_coefs_renum[table[it->first]] = it->second;
  }
  rotation_coefs = rotation_coefs_renum;
}

int*
MpcElement::nodes(int* p)
{
  if(p == 0) p = new int[numNodes()];
  for(int i = 0; i < numNodes(); ++i) p[i] = nn[i];
  return p;
}

int
MpcElement::numDofs()
{
  // note: augmented lagrangian does NOT have a degree of freedom for the lagrange multiplier
  return (prop->lagrangeMult && prop->penalty == 0) ? nterms+1 : nterms;
}

int *
MpcElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];
  for(int i = 0; i < nterms; i++)
    dsa.number(terms[i].nnum, 1 << terms[i].dofnum, p+i);
  if(prop->lagrangeMult && prop->penalty == 0) {
    if(type == 0 && prop->relop == 0)
      dsa.number(nn[nNodes], DofSet::LagrangeE, p+nterms);
    else
      dsa.number(nn[nNodes], DofSet::LagrangeI, p+nterms);
  }
  return p;
}

void
MpcElement::markDofs(DofSetArray &dsa)
{
  for(int i = 0; i < nterms; i++)
    dsa.mark(terms[i].nnum, 1 << terms[i].dofnum);
  if(prop->lagrangeMult && prop->penalty == 0) {
    if(type == 0 && prop->relop == 0)
      dsa.mark(nn[nNodes], DofSet::LagrangeE);
    else
      dsa.mark(nn[nNodes], DofSet::LagrangeI);
  }
}

FullSquareMatrix
MpcElement::stiffness(CoordSet& c0, double* karray, int)
{
  // this function computes the constraint stiffness matrix for linear statics and dynamics
  // see comments in getStiffAndForce (nonlinear version)
  FullSquareMatrix ret(numDofs(), karray);
  ret.zero();
  if(prop->lagrangeMult || prop->penalty != 0) {
    double lambda = 0;
    if(prop->penalty != 0 && (type == 1 && -rhs.r_value <= -lambda/prop->penalty)) {
      //if(prop->lagrangeMult) {
      //  ret[nterms][nterms] = -1/prop->penalty;
      //}
    }
    else {
      if(prop->penalty != 0) lambda += prop->penalty*(-rhs.r_value);
      GeomState c1(c0);
      FullSquareMatrix H(nterms);
      if(lambda != 0) getHessian(&c1, c1, c0, H, 0); else H.zero();
      for(int i = 0; i < nterms; ++i) {
        for(int j = 0; j < nterms; ++j) {
          ret[i][j] = lambda*H[i][j];
          if(prop->penalty != 0) ret[i][j] += prop->penalty*terms[i].coef.r_value*terms[j].coef.r_value;
        }
        if(prop->lagrangeMult && prop->penalty == 0) ret[i][nterms] = ret[nterms][i] = terms[i].coef.r_value;
      }
    }
  }

  return ret;
}

void
MpcElement::getGravityForce(CoordSet&, double*, Vector& f, int, GeomState*)
{
  f.zero();
}

Corotator*
MpcElement::getCorotator(CoordSet&, double*, int, int)
{
  return this;
}

void
MpcElement::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* f, double delt, double t)
{
  getStiffAndForce(NULL, c1, c0, Ktan, f, delt, t);
}

void
MpcElement::getStiffAndForce(GeomState* refState, GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* f, double delt, double t)
{
  Ktan.zero();
  for(int i = 0; i < numDofs(); ++i) f[i] = 0.0;
  if(getSource() != mpc::ContactSurfaces)
    update(refState, c1, c0, t); // update rhs and coefficients to the value and gradient the constraint function, respectively 
  if(t == 0 && delt > 0 && type == 0 && prop->penalty == 0)
    rhs.r_value = this->getAccelerationConstraintRhs(refState, c1, c0, t); // for dynamics we compute initial acceleration at t=0 for equality constraints

  // general augmented lagrangian implementation from RT Rockafellar "Lagrange multipliers and optimality" Siam Review 1993, eq 6.7
  // NOTES:
  //  1. penalty method is the particular case with prop->lagrangeMult is set to false
  //  2. multipliers method is the particular case with prop->penalty set to 0
  //  3. for direct elimination set prop->lagrangeMult to false and prop->penalty to 0

  if(prop->lagrangeMult || prop->penalty != 0) {
    double lambda = (prop->lagrangeMult) ? c1[nn[nNodes]].x : 0; // y is the lagrange multiplier (if used)
    if(prop->penalty != 0 && (type == 1 && -rhs.r_value <= -lambda/prop->penalty)) { //
      //if(prop->lagrangeMult) {
      //  Ktan[nterms][nterms] = -1/prop->penalty;
      //  f[nterms] = -lambda/prop->penalty;
      //}
    }
    else {
      if(prop->penalty != 0) lambda += prop->penalty*(-rhs.r_value);
      FullSquareMatrix H(nterms);
      if(lambda != 0) getHessian(refState, c1, c0, H, t); else H.zero();
      for(int i = 0; i < nterms; ++i) {
        for(int j = 0; j < nterms; ++j) {
          Ktan[i][j] = lambda*H[i][j];
          if(prop->penalty != 0) Ktan[i][j] += prop->penalty*terms[i].coef.r_value*terms[j].coef.r_value;
        }
        if(prop->lagrangeMult && prop->penalty == 0) Ktan[i][nterms] = Ktan[nterms][i] = terms[i].coef.r_value;
        if(!(type == 1 && prop->lagrangeMult && prop->penalty == 0)) f[i] = lambda*terms[i].coef.r_value;
                                                   // note: for inequalities we solve for lambda^{k} at every SQP iteration
                                                   // but for equalities we solve for the increment (lambda^{k}-lambda^{k-1})

      }
      if(prop->lagrangeMult && prop->penalty == 0) f[nterms] = -rhs.r_value;
    }
  }
}

void
MpcElement::getInternalForce(GeomState *refState, GeomState& c1, CoordSet& c0, FullSquareMatrix&, double* f, double delt, double t)
{
  for(int i = 0; i < numDofs(); ++i) f[i] = 0.0;
  if(getSource() != mpc::ContactSurfaces)
    update(refState, c1, c0, t); // update rhs and coefficients to the value and gradient the constraint function, respectively 
  if(t == 0 && delt > 0 && type == 0 && prop->penalty == 0)
    rhs.r_value = this->getAccelerationConstraintRhs(refState, c1, c0, t); // for dynamics we compute initial acceleration at t=0 for equality constraints

  // general augmented lagrangian implementation from RT Rockafellar "Lagrange multipliers and optimality" Siam Review 1993, eq 6.7
  // NOTES:
  //  1. penalty method is the particular case with prop->lagrangeMult is set to false
  //  2. multipliers method is the particular case with prop->penalty set to 0
  //  3. for direct elimination set prop->lagrangeMult to false and prop->penalty to 0

  if(prop->lagrangeMult || prop->penalty != 0) {
    double lambda = (prop->lagrangeMult) ? c1[nn[nNodes]].x : 0; // y is the lagrange multiplier (if used)
    if(prop->penalty != 0 && (type == 1 && -rhs.r_value <= -lambda/prop->penalty)) { // note: if -rhs == -lambda/penalty
                                                                                     // then the derivative is indeterminate
      //if(prop->lagrangeMult) {
      //  f[nterms] = -lambda/prop->penalty;
      //}
    }
    else {
      if(prop->penalty != 0) lambda += prop->penalty*(-rhs.r_value);
      for(int i = 0; i < nterms; ++i) {
        if(!(type == 1 && prop->lagrangeMult && prop->penalty == 0)) f[i] = lambda*terms[i].coef.r_value; // for inequalities we solve for lambda^{k}
                                                   // at every SQP iteration but for equalities we solve for the increment (lambda^{k}-lambda^{k-1})

      }
      if(prop->lagrangeMult && prop->penalty == 0) f[nterms] = -rhs.r_value;
    }
  }
}

void 
MpcElement::update(GeomState* refState, GeomState& c1, CoordSet& c0, double t) 
{ 
  // this is for a linear constraint in nonlinear analysis. Nonlinear constraints must overload this function
  rhs = original_rhs;

  // data structures to store rotation vector and updates coefficients for nodes with rotational lmpc terms
  std::map<int, std::vector<double> > theta, updated_coefs;

  // initialize the rotation dofs from the current state
  std::map<int,std::vector<double> >::iterator it;
  for(it = rotation_coefs.begin(); it != rotation_coefs.end(); ++it) {
    int nodeNumber = it->first;
    std::vector<double> v(3);
    v[0] = c1[nodeNumber].theta[0];
    v[1] = c1[nodeNumber].theta[1];
    v[2] = c1[nodeNumber].theta[2];
    if(NFrameData *cd = c0.dofFrame(nodeNumber)) cd->transformVector3(&v[0]); // transform from basic frame to DOF_FRM
    theta[nodeNumber] = v;
    double rotvar[3][3];
    pseudorot_var(&v[0], rotvar);
    mat_mult_vec(rotvar, &it->second[0], &v[0], 1);
    updated_coefs[nodeNumber] = v;
  }

  for(int i = 0; i < nterms; ++i) {
    double u;
    NFrameData *cd = c0.dofFrame(terms[i].nnum);
    switch(terms[i].dofnum) {
      case 0 : 
        if(!cd) {
          u = c1[terms[i].nnum].x-c0[terms[i].nnum]->x; 
          rhs.r_value -= terms[i].coef.r_value*u;
        }
        else {
          double u[3] = { c1[terms[i].nnum].x-c0[terms[i].nnum]->x,
                          c1[terms[i].nnum].y-c0[terms[i].nnum]->y,
                          c1[terms[i].nnum].z-c0[terms[i].nnum]->z };
          cd->transformVector3(u);
          rhs.r_value -= terms[i].coef.r_value*u[0];
        }
        break;
      case 1 : 
        if(!cd) {
          u = c1[terms[i].nnum].y-c0[terms[i].nnum]->y;
          rhs.r_value -= terms[i].coef.r_value*u;
        }
        else {
          double u[3] = { c1[terms[i].nnum].x-c0[terms[i].nnum]->x,
                          c1[terms[i].nnum].y-c0[terms[i].nnum]->y,
                          c1[terms[i].nnum].z-c0[terms[i].nnum]->z };
          cd->transformVector3(u);
          rhs.r_value -= terms[i].coef.r_value*u[1];
        }
        break;
      case 2 : 
        if(!cd) {
          u = c1[terms[i].nnum].z-c0[terms[i].nnum]->z;
          rhs.r_value -= terms[i].coef.r_value*u;
        }
        else {
          double u[3] = { c1[terms[i].nnum].x-c0[terms[i].nnum]->x,
                          c1[terms[i].nnum].y-c0[terms[i].nnum]->y,
                          c1[terms[i].nnum].z-c0[terms[i].nnum]->z };
          cd->transformVector3(u);
          rhs.r_value -= terms[i].coef.r_value*u[2];
        }
        break;
      case 3 : case 4 : case 5 : {
        u = theta[terms[i].nnum][terms[i].dofnum-3];
        rhs.r_value -= rotation_coefs[terms[i].nnum][terms[i].dofnum-3]*u;
        terms[i].coef.r_value = updated_coefs[terms[i].nnum][terms[i].dofnum-3];
      } break;
    }
  }
}

void 
MpcElement::getHessian(GeomState*, GeomState& c1, CoordSet& c0, FullSquareMatrix& _H, double t) 
{
#ifdef USE_EIGEN3
  if(getSource() == mpc::ContactSurfaces && H.size() > 0) {
    for(int i=0; i<H.rows(); ++i)
      for(int j=0; j<H.cols(); ++j)
        _H[i][j] = H(i,j);
  }
  else
#endif
  {
    _H.zero();

    // compute second variation of the rotational dof terms
    double v[3], rotvar[3][3], scndvar[3][3];
    std::map<int,std::vector<double> >::iterator it;
    for(it = rotation_coefs.begin(); it != rotation_coefs.end(); ++it) {
      int nodeNumber = it->first;
      std::vector<int> dofs = rotation_indices[nodeNumber];
      v[0] = c1[nodeNumber].theta[0];
      v[1] = c1[nodeNumber].theta[1];
      v[2] = c1[nodeNumber].theta[2];
      if(NFrameData *cd = c0.dofFrame(nodeNumber)) cd->transformVector3(&v[0]); // transform from basic frame to DOF_FRM
      pseudorot_2var(v, &it->second[0], scndvar);
      for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
          _H[dofs[i]][dofs[j]] = 0.5*(scndvar[i][j]+scndvar[j][i]);
    }
  }
}

double
MpcElement::getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double)
{
  return 0;
}

double
MpcElement::getAccelerationConstraintRhs(GeomState* refState, GeomState& gState, CoordSet& cs, double t)
{
  // compute rhs = -(G(q)*qdot)_q * qdot assuming other terms are zero. Overload function if this assumption is not correct
  // XXX consider nodal frames (careful though, this function is used by both LMPCs and constraint function elements)
  FullSquareMatrix H(nterms);
  getHessian(refState, gState, cs, H, t);
  Vector v(nterms);
  // fill v with initial velocity 
  for(int i = 0; i < nterms; ++i) {
   if(terms[i].dofnum == 3 || terms[i].dofnum == 4 || terms[i].dofnum == 5) {
      // compute spatial angular velocity
      double omega[3];
      mat_mult_vec(gState[terms[i].nnum].R, &gState[terms[i].nnum].v[3], omega);
      v[i] = omega[terms[i].dofnum-3];
    }
    else {
      v[i] = gState[terms[i].nnum].v[terms[i].dofnum];
    }
  }

  Vector Hv(nterms,0.0);
  for(int i = 0; i < nterms; ++i)
    for(int j = 0; j < nterms; ++j)
      Hv[i] += H[i][j]*v[j];

  return -(v*Hv);
}

void
MpcElement::getResidualCorrection(GeomState& c1, double* r)
{
  // note #1: r = [0; -f] + dr = [-G^t*lambda; -(pos_part<f>-neg_part<lambda>)]
  // therefore dr = [-G^t*lambda; -(-neg_part<f>-neg_part<lambda>)]
  //              = [-G^t*lambda; -pos_part<rhs>+neg_part<lambda>]
  // note #2: rhs = -f (-ve value of the constraint function f <= 0)
  // 
  // r = [-G^t*lambda; pos_part<f>+neg_part<lambda>]
  if(prop->lagrangeMult && prop->penalty == 0 && type == 1) {
    double lambda = c1[nn[nNodes]].x;
    for(int i = 0; i < nterms; ++i)
      r[i] -= lambda*terms[i].coef.r_value;
    if(rhs.r_value > 0) r[nterms] -= rhs.r_value;
    /*when using exact QP solver works then lambda is always dual-feasible
    if(lambda < 0) r[nterms] += lambda; */
  }
}

double
MpcElement::getElementEnergy(GeomState& c1, CoordSet& c0)
{
  // TODO
  return 0;
}

void
MpcElement::computePressureForce(CoordSet&, Vector& f, GeomState*, int, double)
{
  // this function computes the constraint force vector for linear statics and dynamics
  // see comments in getStiffAndForce (nonlinear version)
  f.zero();
  if(prop->lagrangeMult || prop->penalty != 0) {
    double lambda = 0;
    if(prop->penalty != 0 && (type == 1 && -rhs.r_value <= -lambda/prop->penalty)) {
/*
      if(prop->lagrangeMult) {
        f[nterms] = lambda/prop->penalty;
      }
*/
    }
    else {
      if(prop->penalty != 0) lambda += prop->penalty*(-rhs.r_value);
      for(int i = 0; i < nterms; ++i) {
        f[i] = -lambda*terms[i].coef.r_value;
      }
      if(prop->lagrangeMult && prop->penalty == 0) f[nterms] = rhs.r_value;
    }
  }
}

void
MpcElement::getNLVonMises(Vector& stress, Vector& weight, GeomState&, CoordSet&, int)
{
  stress.zero();
  weight.zero();
}

void
MpcElement::getNLAllStress(FullM& stress, Vector& weight, GeomState&, CoordSet&, int)
{
  stress.zero();
  weight.zero();
}

void
MpcElement::initMultipliers(GeomState& c1)
{
  // augmented lagrangian
  if(prop->lagrangeMult && prop->penalty != 0) {
    c1[nn[nNodes]].x = 0;
  }
}

void
MpcElement::updateMultipliers(GeomState& c1)
{
  // augmented lagrangian
  if(prop->lagrangeMult && prop->penalty != 0) {
    c1[nn[nNodes]].x += prop->penalty*(-rhs.r_value);
    if(type == 1 && c1[nn[nNodes]].x < 0) c1[nn[nNodes]].x = 0;
  }
}

double
MpcElement::getError()
{
  return (type == 1) ? std::max(0.0,-rhs.r_value) : std::abs(-rhs.r_value);
}
