#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <Element.d/MpcElement.d/MpcElement.h>
#include <Driver.d/Domain.h>
#include <Corotational.d/Corotator.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

MpcElement::MpcElement(LMPCons* _mpc)
{
  mpc = _mpc;
  nn = new int[mpc->nterms + 1];

  // count the number of nodes. also fill node numbers array
  nnodes = 0;
  int i,j;
  for(i = 0; i < mpc->nterms; i++) {
    int nnum = (mpc->terms)[i].nnum;
    bool found = false;
    for(j=0; j<nnodes; ++j) {
      if(nnum == nn[j]) { found = true; break; } 
    }
    if(!found) nn[nnodes++] = nnum;
  }
  ndofs = mpc->nterms + 1;
  nn[nnodes] = -1; // internal node, number is set later
  // cerr << " ndofs = " << ndofs << ", nnodes = " << nnodes << ", nodes: "; for(i=0; i<nnodes; ++i) cerr << nn[i] << " "; cerr << endl;

  renumTable = 0;
  prop = new StructProp();
}

void
MpcElement::renum(int *table)
{
  int i;
  renumTable = new int[mpc->nterms+1];
  for(i=0; i<nnodes; ++i) nn[i] = table[nn[i]];
  for(i = 0; i < mpc->nterms; i++) {
    int nnum = (mpc->terms)[i].nnum;
    renumTable[i] = table[nnum];
  }
  renumTable[mpc->nterms] = table[nn[nnodes-1]];
}

FullSquareMatrix
MpcElement::stiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(ndofs,k);
  ret.zero();
  
  int i;
  for(i = 0; i < mpc->nterms; i++) {
    if((mpc->terms)[i].isComplex) {
      double val = (mpc->terms)[i].coef.c_value.real();
      ret[i][ndofs-1] = ret[ndofs-1][i] = val;
    }
    else {
      double val = (mpc->terms)[i].coef.r_value;
      ret[i][ndofs-1] = ret[ndofs-1][i] = val;
    }
  }
  //ret *= 1.0e3; // XXXX
  return ret;
}

FullSquareMatrix
MpcElement::imagStiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(ndofs,k);
  ret.zero();

  int i;
  for(i = 0; i < mpc->nterms; i++) {
    if((mpc->terms)[i].isComplex) {
      double val = (mpc->terms)[i].coef.c_value.imag();
      ret[i][ndofs-1] = ret[ndofs-1][i] = val;
    }
  }
  //cerr << "MpcElement imaginary stiffness matrix = \n"; ret.print();
  return ret;
}

int *
MpcElement::nodes(int *p)
{
  if(p == 0) p = new int[nnodes];
  int i;
  for(i=0; i<nnodes; ++i) p[i] = nn[i];
 
  // cerr << "nnodes = " << nnodes << ", MpcElement nodes: "; for(i=0; i<nnodes; ++i) cerr << p[i] << " "; cerr << endl; 
  return p;
}

int *
MpcElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[ndofs];

  int i;
  for(i = 0; i < mpc->nterms; i++) {
    int nnum = (mpc->terms)[i].nnum;
    if(renumTable) nnum = renumTable[i];
    int dofnum = (mpc->terms)[i].dofnum;
    dsa.number(nnum, (1 << dofnum), p+i);
  }
  dsa.number(nn[nnodes-1],DofSet::Lagrange1, p+ndofs-1);

  // cerr << "MpcElement dofs: "; for(i=0; i<ndofs; ++i) cerr << p[i] << " "; cerr << endl;

  return p;
}

void
MpcElement::markDofs(DofSetArray &dsa)
{
  int i;
  for(i = 0; i < mpc->nterms; i++) {
    int nnum = (mpc->terms)[i].nnum;
    if(renumTable) nnum = renumTable[i];
    int dofnum = (mpc->terms)[i].dofnum;
    dsa.mark(nnum, (1 << dofnum));
  }
  dsa.mark(nn[nnodes-1],DofSet::Lagrange1);
}

void 
MpcElement::getStiffAndForce(GeomState &gState, CoordSet &cs, FullSquareMatrix &Ktan, double *f)
{
  NodeState lagrangeNode = gState[nn[nnodes-1]];
  Ktan.zero();

  int i;
  if(!mpc->isComplex) {
    double error = -mpc->rhs.r_value;
    for(i = 0; i < mpc->nterms; i++) {
      double val = (mpc->terms)[i].coef.r_value;
      Ktan[i][ndofs-1] = Ktan[ndofs-1][i] = val;
      f[i] = val*lagrangeNode.x;
      int nnum = (mpc->terms)[i].nnum;
      int dofnum = (mpc->terms)[i].dofnum;
      if(renumTable) nnum = renumTable[i];
      Node &node1 = cs.getNode(nnum);
      NodeState ns1 = gState[nnum];
      error += val * ns1.diff(node1, dofnum);
    }
    f[ndofs-1] = error;
  }
  else cerr << " *** WARNING: getStiffAndForce not implemented for complex MpcElement \n";
}
