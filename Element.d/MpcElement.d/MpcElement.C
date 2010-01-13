#include <Element.d/MpcElement.d/MpcElement.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>

MpcElement::MpcElement(int _nNodes, DofSet nodalDofs, int* _nn)
 : nNodes(_nNodes), LMPCons(0, 0.0)
{
  // this constructor is for a constraint involving the same DofSet on each node
  nn = new int[nNodes+1];
  for(int i = 0; i < nNodes; ++i) {
    nn[i] = _nn[i];
    for(int j = 0; j < nodalDofs.count(); ++j)
      for(int k = 0; k < 11; ++k)
        if(nodalDofs.contains(1 << k)) addterm(new LMPCTerm(nn[i], k, 0.0));
  }
}

MpcElement::MpcElement(int _nNodes, DofSet *nodalDofs, int* _nn)
 : nNodes(_nNodes), LMPCons(0, 0.0)
{
  // this constructor is for a constraint involving a different DofSet on each nodes
  nn = new int[nNodes+1];
  for(int i = 0; i < nNodes; ++i) {
    nn[i] = _nn[i];
    for(int j = 0; j < nodalDofs[i].count(); ++j)
      for(int k = 0; k < 11; ++k)
        if(nodalDofs[i].contains(1 << k)) addterm(new LMPCTerm(nn[i], k, 0.0));
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
  mpcs[0] = this;
  return mpcs;
}

int
MpcElement::numInternalNodes()
{
  return 1;
}

void
MpcElement::setInternalNodes(int* in)
{
  nn[nNodes] = in[0];
}

int
MpcElement::numNodes()
{
  return nNodes+1;
}

void
MpcElement::renum(int* table)
{
  for(int i = 0; i < nNodes+1; ++i) nn[i] = table[nn[i]];
  for(int i = 0; i < nterms; ++i) terms[i].nnum = table[terms[i].nnum];
}

int*
MpcElement::nodes(int* p)
{
  if(p == 0) p = new int[nNodes+1];
  for(int i = 0; i < nNodes+1; ++i) p[i] = nn[i];
  return p;
}

int
MpcElement::numDofs()
{
  return nterms+1;
}

int *
MpcElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[nterms+1];
  for(int i = 0; i < nterms; i++) dsa.number(terms[i].nnum, 1 << terms[i].dofnum, p+i);
  dsa.number(nn[nNodes], DofSet::Lagrange, p+nterms);
  return p;
}

void
MpcElement::markDofs(DofSetArray &dsa)
{
  for(int i = 0; i < nterms; i++) dsa.mark(terms[i].nnum, 1 << terms[i].dofnum);
  dsa.mark(nn[nNodes], DofSet::Lagrange);
}

FullSquareMatrix
MpcElement::stiffness(CoordSet&, double* karray, int)
{
  FullSquareMatrix ret(nterms+1, karray);
  ret.zero();
  for(int i = 0; i < nterms; ++i) ret[i][nterms] = ret[nterms][i] = terms[i].coef.r_value;
  return ret;

}

Corotator*
MpcElement::getCorotator(CoordSet&, double*, int, int)
{
  return this;
}

void
MpcElement::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* f)
{
  Ktan.zero();
  for(int i = 0; i < nterms+1; ++i) f[i] = 0.0;

  update(c1, c0); // update rhs and coefficients to the value and gradient the constraint function, respectively 
  FullSquareMatrix H(nterms);
  getHessian(c1, c0, H); // H is the hessian of the constraint function
  double lambda = c1[nn[nNodes]].x;
  for(int i = 0; i < nterms; ++i)
    for(int j = 0; j < nterms; ++j) 
      Ktan[i][j] += lambda*H[i][j];
  for(int i = 0; i < nterms; ++i) {
    Ktan[i][nterms] = Ktan[nterms][i] = terms[i].coef.r_value;
    f[i] += lambda*terms[i].coef.r_value;
  }
  f[nterms] = rhs.r_value; // value of the constraint function
}

void 
MpcElement::update(GeomState&, CoordSet&) 
{ 
  rhs.r_value = 0.0;
}

void 
MpcElement::getHessian(GeomState&, CoordSet&, FullSquareMatrix& H) 
{ 
  H.zero(); 
}

