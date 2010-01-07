#include <Element.d/MpcElement.d/ConstraintElement.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>

ConstraintElement::ConstraintElement(int _nNodes, DofSet _nodalDofs, int* _nn)
 : nNodes(_nNodes), nDofs(_nNodes*_nodalDofs.count()), mode(1), LMPCons(0, 0.0)
{
  // this constructor is for a constraint involving the same DofSet on each node
  nDofsPerNode = new int[nNodes];
  nodalDofs = new DofSet[nNodes];
  nn = new int[nNodes+1];
  for(int i = 0; i < nNodes; ++i) {
    nDofsPerNode[i] = _nodalDofs.count();
    nodalDofs[i] = _nodalDofs;
    nn[i] = _nn[i];
    for(int j = 0; j < nodalDofs[i].count(); ++j)
      for(int k = 0; k < 11; ++k)
        if(nodalDofs[i].contains(1 << k)) addterm(new LMPCTerm(nn[i], k, 0.0));
  }
}

ConstraintElement::ConstraintElement(LMPCons *mpc)
 : LMPCons(*mpc)
{
}

ConstraintElement::~ConstraintElement()
{
  delete [] nn;
}

int
ConstraintElement::getNumMPCs()
{
  return 1;
}

LMPCons** 
ConstraintElement::getMPCs()
{
  mode = 0; // XXXX
  LMPCons **mpcs = new LMPCons*[1];
  mpcs[0] = this;
  return mpcs;
}

int
ConstraintElement::numInternalNodes()
{
  return (mode == 0) ? 0 : 1;
}

void
ConstraintElement::setInternalNodes(int* in)
{
  if(mode != 0) nn[nNodes] = in[0];
}

int
ConstraintElement::numNodes()
{
  return (mode == 0) ? nNodes : nNodes+1;
}

void
ConstraintElement::renum(int* table)
{
  for(int i = 0; i < nNodes; ++i) nn[i] = table[nn[i]];
  if(mode != 0) {
    nn[nNodes] = table[nn[nNodes]];
    for(int j = 0; j < nterms; ++j)
      terms[j].nnum = table[terms[j].nnum];
  }
}

int*
ConstraintElement::nodes(int* p)
{
  if(p == 0) p = new int[numNodes()];
  for(int i = 0; i < nNodes; ++i) p[i] = nn[i];
  if(mode != 0) p[nNodes] = nn[nNodes];
  return p;
}

int
ConstraintElement::numDofs()
{
  return (mode == 0) ? nDofs : nDofs + 1;
}

int *
ConstraintElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];

  int k=0;
  for(int i = 0; i < nNodes; ++i) {
    dsa.number(nn[i], nodalDofs[i], p+k);
    k += nDofsPerNode[i];
  }
  if(mode != 0)
    dsa.number(nn[nNodes], DofSet::Lagrange1, p+k);
  return p;
}

void
ConstraintElement::markDofs(DofSetArray &dsa)
{
  for(int i = 0; i < nNodes; ++i)
    dsa.mark(nn+i, 1, nodalDofs[i].list());
  if(mode != 0)
    dsa.mark(nn+nNodes, 1, DofSet::Lagrange1);
}

FullSquareMatrix
ConstraintElement::stiffness(CoordSet&, double* karray, int)
{
  FullSquareMatrix ret(numDofs(), karray);
  ret.zero();
  if(mode != 0) {
    int j, k, l;
    for(j = 0; j < nterms; ++j) {
      l = 0;
      for(k = 0; k < nNodes; ++k) { if (terms[j].nnum == nn[k]) break; l += nDofsPerNode[k]; }
      l += nodalDofs[k].locate(1 << terms[j].dofnum);
      ret[l][nDofs] = ret[nDofs][l] = terms[j].coef.r_value;
    }
  }
  //cerr << "here in ConstraintElement::stiffness\n";
  //ret.print();
  return ret;
}

Corotator*
ConstraintElement::getCorotator(CoordSet&, double*, int, int)
{
  return this;
}

void
ConstraintElement::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* f)
{
  Ktan.zero();
  for(int i = 0; i < numDofs(); ++i) f[i] = 0.0;

  if(mode != 0) {
    updateLMPCs(c1, c0); // update the rhs and coefficients of each mpc to the value and gradient of it's constraint function, respectively 
    int j, k, l;
    FullSquareMatrix H(nDofs);
    getHessian(c1, c0, H); // compute the hessian of the ith constraint function and store in H
    double lambda = c1[nn[nNodes]].x;
    for(j = 0; j < nDofs; ++j)
      for(k = 0; k < nDofs; ++k) 
        Ktan[j][k] += lambda*H[j][k];
    for(j = 0; j < nterms; ++j) {
      l = 0;
      for(k = 0; k < nNodes; ++k) { if (terms[j].nnum == nn[k]) break; l += nDofsPerNode[k]; }
      l += nodalDofs[k].locate(1 << terms[j].dofnum);
      Ktan[l][nDofs] = Ktan[nDofs][l] = terms[j].coef.r_value;
      f[l] += lambda*terms[j].coef.r_value;
    }
    f[nDofs] = rhs.r_value; // value of the ith constraint function
  }
  //cerr << "here in ConstraintElement::getStiffAndForce\n";
  //Ktan.print();
}

void 
ConstraintElement::updateLMPCs(GeomState&, CoordSet&) 
{ 
  rhs.r_value = 0.0;
}

void 
ConstraintElement::getHessian(GeomState&, CoordSet&, FullSquareMatrix& H) 
{ 
  H.zero(); 
}

