#include <Element.d/MpcElement.d/RigidMpcElement.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

RigidMpcElement::RigidMpcElement(int _nNodes, int _nDofsPerNode, int _nodalDofs, int _nMpcs, int* _nn)
 : nNodes(_nNodes), nDofsPerNode(_nDofsPerNode), nodalDofs(_nodalDofs), nMpcs(_nMpcs), mode(1)
{
  nn = new int[nNodes+nMpcs];
  nn_orig = new int[nNodes];
  for(int i = 0; i < nNodes; ++i) 
    nn[i] = nn_orig[i] = _nn[i];
  mpcs = new LMPCons*[nMpcs];
}

RigidMpcElement::~RigidMpcElement()
{
  delete [] nn;
  delete [] nn_orig;
  if(mpcs) delete [] mpcs;
}

void
RigidMpcElement::setProp(StructProp*)
{
  /* do nothing */
}

int
RigidMpcElement::getNumMPCs()
{
  return nMpcs;
}

LMPCons** 
RigidMpcElement::getMPCs()
{
  mode = 0; // XXXX
  return mpcs;
}

int
RigidMpcElement::numInternalNodes()
{
  return (mode == 0) ? 0 : nMpcs;
}

void
RigidMpcElement::setInternalNodes(int* in)
{
  if(mode != 0) {
    for(int i = 0; i < nMpcs; ++i) 
      nn[nNodes+i] = in[i];
  }
}

int
RigidMpcElement::numNodes()
{
  return (mode == 0) ? nNodes : nNodes+nMpcs;
}

void
RigidMpcElement::renum(int* table)
{
  for(int i = 0; i < nNodes; ++i) nn[i] = table[nn[i]];
  if(mode != 0) {
    for(int i = 0; i < nMpcs; ++i) 
      nn[nNodes+i] = table[nn[nNodes+i]];
  }
}

int*
RigidMpcElement::nodes(int* p)
{
  if(p == 0) p = new int[numNodes()];
  for(int i = 0; i < nNodes; ++i) p[i] = nn[i];
  if(mode != 0) {
    for(int i = 0; i < nMpcs; ++i) p[nNodes+i] = nn[nNodes+i];
  }
  return p;
}

int
RigidMpcElement::numDofs()
{
  return (mode == 0) ? nNodes*nDofsPerNode : nNodes*nDofsPerNode + nMpcs;
}

int *
RigidMpcElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];

  for(int i = 0; i < nNodes; ++i) 
    dsa.number(nn[i], nodalDofs /*DofSet::XYZdisp | DofSet::XYZrot*/, p+(nDofsPerNode*i) );
  if(mode != 0) {
    for(int i = 0; i < nMpcs; ++i) 
      dsa.number(nn[nNodes+i], DofSet::Lagrange1, p+(nNodes*nDofsPerNode+i) );
  }
  return p;
}

void
RigidMpcElement::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, nNodes, nodalDofs /*DofSet::XYZdisp | DofSet::XYZrot*/);
  if(mode != 0) {
    dsa.mark(nn+nNodes, nMpcs, DofSet::Lagrange1);
  }
}

FullSquareMatrix
RigidMpcElement::stiffness(CoordSet&, double* k, int)
{
  FullSquareMatrix ret(numDofs(), k);
  ret.zero();
  if(mode != 0) {
    int offset = nNodes*nDofsPerNode;
    int i, j, k, l;
    for(i = 0; i < nMpcs; ++i) {
      for(j = 0; j < mpcs[i]->nterms; ++j) {
        for(k = 0; k < nNodes; ++k) if ( mpcs[i]->terms[j].nnum == nn_orig[k]) break;
        l = k*nDofsPerNode + mpcs[i]->terms[j].dofnum;
        ret[l][offset+i] = ret[offset+i][l] = mpcs[i]->terms[j].coef.r_value;
      }
    }
  }
  return ret;
}

Corotator*
RigidMpcElement::getCorotator(CoordSet&, double*, int, int)
{
  return this;
}

void
RigidMpcElement::getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix& Ktan, double* f)
{
  Ktan.zero();
  for(int i = 0; i < numDofs(); ++i) f[i] = 0.0;
}

bool 
RigidMpcElement::isRigidMpcElement(const DofSet& ds, bool forAllNodes)
{ 
  return  ds == DofSet::nullDofset || ds == nodalDofs;
}

