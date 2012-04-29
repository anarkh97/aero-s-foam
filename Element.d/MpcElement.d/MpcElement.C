#include <Element.d/MpcElement.d/MpcElement.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>


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

MpcElement::MpcElement(LMPCons *mpc)
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
}

MpcElement::~MpcElement()
{
  delete [] nn;
}

int
MpcElement::getNumMPCs()
{
  return 1;
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
  return (prop->lagrangeMult) ? nterms+1 : nterms;
}

int *
MpcElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];
  for(int i = 0; i < nterms; i++)
    dsa.number(terms[i].nnum, 1 << terms[i].dofnum, p+i);
  if(prop->lagrangeMult) {
    if(type == 0)
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
  if(prop->lagrangeMult) {
    if(type == 0)
      dsa.mark(nn[nNodes], DofSet::LagrangeE);
    else
      dsa.mark(nn[nNodes], DofSet::LagrangeI);
  }
}

FullSquareMatrix
MpcElement::stiffness(CoordSet& c0, double* karray, int)
{
  // this function computes the constraint stiffness matrix for linear statics and dynamics
  // see comments in ::getStiffAndForce (nonlinear version)
  FullSquareMatrix ret(numDofs(), karray);
  ret.zero();
  if(prop->lagrangeMult || prop->penalty != 0) {
    double lambda = 0;
    if(prop->penalty != 0.0 && (type == 1 && -rhs.r_value <= lambda/prop->penalty)) {
      if(prop->lagrangeMult) {
        ret[nterms][nterms] = -1/prop->penalty;
      }
    }
    else {
      if(prop->penalty != 0) lambda += prop->penalty*(-rhs.r_value);
      GeomState c1(c0);
      FullSquareMatrix H(nterms);
      getHessian(c1, c0, H);
      for(int i = 0; i < nterms; ++i) {
        for(int j = 0; j < nterms; ++j) {
          ret[i][j] = lambda*H[i][j];
          if(prop->penalty != 0) ret[i][j] += prop->penalty*terms[i].coef.r_value*terms[j].coef.r_value;
        }
        if(prop->lagrangeMult) ret[i][nterms] = ret[nterms][i] = terms[i].coef.r_value;
      }
    }
  }

  return ret;
}

Corotator*
MpcElement::getCorotator(CoordSet&, double*, int, int)
{
  return this;
}

void
MpcElement::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* f, double delt, double t)
{
  Ktan.zero();
  for(int i = 0; i < numDofs(); ++i) f[i] = 0.0;
  if(getSource() != mpc::ContactSurfaces)
    update(c1, c0, t); // update rhs and coefficients to the value and gradient the constraint function, respectively 
  if(t == 0 && delt > 0 && type == 0)
    rhs.r_value = getAccelerationConstraintRhs(&c1, c0, t); // for dynamics we compute initial acceleration at t=0 for equality constraints

  // general augmented lagrangian implementation from RT Rockafellar "Lagrange multipliers and optimality" Siam Review 1993, eq 6.7
  // NOTES:
  //  1. penalty method is the particular case with prop->lagrangeMult is set to false
  //  2. multipliers method is the particular case with prop->penalty set to 0
  //  3. for direct elimination set prop->lagrangeMult to false and prop->penalty to 0

  if(prop->lagrangeMult || prop->penalty != 0.0) {
    double lambda = (prop->lagrangeMult) ? c1[nn[nNodes]].x : 0; // y is the lagrange multiplier (if used)
    if(prop->penalty != 0.0 && (type == 1 && -rhs.r_value <= lambda/prop->penalty)) { //
      if(prop->lagrangeMult) {
        Ktan[nterms][nterms] = -1/prop->penalty;
        f[nterms] = -lambda/prop->penalty;
      }
    }
    else {
      if(prop->penalty != 0) lambda += prop->penalty*(-rhs.r_value);
      FullSquareMatrix H(nterms);
      getHessian(c1, c0, H);
      for(int i = 0; i < nterms; ++i) {
        for(int j = 0; j < nterms; ++j) {
          Ktan[i][j] = lambda*H[i][j];
          if(prop->penalty != 0) Ktan[i][j] += prop->penalty*terms[i].coef.r_value*terms[j].coef.r_value;
        }
        if(prop->lagrangeMult) Ktan[i][nterms] = Ktan[nterms][i] = terms[i].coef.r_value;
        if(!(type == 1 && prop->lagrangeMult)) f[i] = lambda*terms[i].coef.r_value; // for inequalities we solve for lambda^{k} at every SQP iteration
                                                   // but for equalities we solve for the increment (lambda^{k}-lambda^{k-1})

      }
      if(prop->lagrangeMult) f[nterms] = -rhs.r_value;
    }
  }
}

void 
MpcElement::update(GeomState& c1, CoordSet& c0, double t) 
{ 
  // this is for a linear constraint. Nonlinear constraints must overload this function
  rhs = original_rhs;

  for(int i = 0; i < nterms; ++i) {
    double u;
    switch(terms[i].dofnum) {
      case 0 : u = c1[terms[i].nnum].x-c0[terms[i].nnum]->x; break;
      case 1 : u = c1[terms[i].nnum].y-c0[terms[i].nnum]->y; break;
      case 2 : u = c1[terms[i].nnum].z-c0[terms[i].nnum]->z; break;
      case 3 : case 4 : case 5 : {
        double theta[3];
        mat_to_vec(c1[terms[i].nnum].R, theta);
        u = theta[terms[i].dofnum-3];
      } break;
    }
    rhs.r_value -= terms[i].coef.r_value*u;
  }
}

void 
MpcElement::getHessian(GeomState&, CoordSet&, FullSquareMatrix& _H) 
{
  //std::cerr << "here in GenMpcElement<Scalar>::getHessian\n";
  //std::cerr << "eigenvalues of H: " << this->H.eigenvalues().transpose() << std::endl;
#ifdef USE_EIGEN3
  if(getSource() == mpc::ContactSurfaces && H.size() > 0) {
    //std::cerr << "H = \n" << this->H << std::endl;
    for(int i=0; i<H.rows(); ++i)
      for(int j=0; j<H.cols(); ++j)
        _H[i][j] = H(i,j);
  }
  else
#endif
    _H.zero();
}

double
MpcElement::getVelocityConstraintRhs(GeomState *gState, CoordSet& cs, double t)
{
  return 0;
}

double
MpcElement::getAccelerationConstraintRhs(GeomState *gState, CoordSet& cs, double t)
{
  // compute rhs = -(G(q)*qdot)_q * qdot assuming other terms are zero. Overload function if this assumption is not correct
  FullSquareMatrix H(nterms);
  getHessian(*gState, cs, H);
  Vector v(nterms);
  // fill v with initial velocity 
  for(int i = 0; i < nterms; ++i)
    v[i] = (*gState)[terms[i].nnum].v[terms[i].dofnum];

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
  if(prop->lagrangeMult && type == 1) {
    double lambda = c1[nn[nNodes]].x;
    for(int i = 0; i < nterms; ++i)
      r[i] -= lambda*terms[i].coef.r_value;
    if(rhs.r_value > 0) r[nterms] -= rhs.r_value;
    /*when using exact QP solver works then lambda is always dual-feasible
    if(lambda < 0) r[nterms] += lambda; */
  }
}

void
MpcElement::computePressureForce(CoordSet&, Vector& f, GeomState*, int, double)
{
  // this function computes the constraint force vector for linear statics and dynamics
  // see comments in ::getStiffAndForce (nonlinear version)
  f.zero();
  if(prop->lagrangeMult || prop->penalty != 0.0) {
    double lambda = 0;
    if(prop->penalty != 0.0 && (type == 1 && -rhs.r_value <= lambda/prop->penalty)) {
      if(prop->lagrangeMult) {
        f[nterms] = lambda/prop->penalty;
      }
    }
    else {
      if(prop->penalty != 0) lambda += prop->penalty*(-rhs.r_value);
      for(int i = 0; i < nterms; ++i) {
        f[i] = -lambda*terms[i].coef.r_value;
      }
      if(prop->lagrangeMult) f[nterms] = rhs.r_value;
    }
  }
}

void
MpcElement::getNLVonMises(Vector& stress, Vector& weight,
                          GeomState &, CoordSet &, int)
{
  stress.zero();
  weight.zero();
}

