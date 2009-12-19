#include <Element.d/MpcElement.d/RigidMpcElement.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

RigidMpcElement::RigidMpcElement(int _nNodes, int _nDofsPerNode, int _nodalDofs, int _nMpcs, int* _nn)
 : nNodes(_nNodes), nDofsPerNode(_nDofsPerNode), nodalDofs(_nodalDofs), nMpcs(_nMpcs), mode(1)
{
  nn = new int[nNodes+nMpcs];
  for(int i = 0; i < nNodes; ++i) 
    nn[i] = _nn[i];
  mpcs = new LMPCons*[nMpcs];
  first = true;
}

RigidMpcElement::~RigidMpcElement()
{
  delete [] nn;
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
    for(int i = 0; i < nMpcs; ++i) {
      nn[nNodes+i] = table[nn[nNodes+i]];
      for(int j = 0; j < mpcs[i]->nterms; ++j)
        mpcs[i]->terms[j].nnum = table[mpcs[i]->terms[j].nnum];
    }
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
RigidMpcElement::stiffness(CoordSet&, double* karray, int)
{
  FullSquareMatrix ret(numDofs(), karray);
  ret.zero();
  if(mode != 0) {
    int i, j, k, l, m;
    m = nNodes*nDofsPerNode;
    for(i = 0; i < nMpcs; ++i) {
      for(j = 0; j < mpcs[i]->nterms; ++j) {
        for(k = 0; k < nNodes; ++k) if ( mpcs[i]->terms[j].nnum == nn[k]) break;
        l = k*nDofsPerNode + mpcs[i]->terms[j].dofnum;
        ret[l][m+i] = ret[m+i][l] = mpcs[i]->terms[j].coef.r_value;
      }
    }
  }
  //cerr << "here in RigidMpcElement::stiffness\n";
  //ret.print();
  return ret;
}

Corotator*
RigidMpcElement::getCorotator(CoordSet&, double*, int, int)
{
  return this;
}

void
RigidMpcElement::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan, double* f)
{
  Ktan.zero();
  for(int i = 0; i < numDofs(); ++i) f[i] = 0.0;

  if(mode != 0) {
    updateLMPCs(c1, c0); // update the rhs and coefficients of each mpc to the value and gradient of it's constraint function, respectively 
    int i, j, k, l, m;
    m = nNodes*nDofsPerNode;
    FullSquareMatrix H(m);
    for(i = 0; i < nMpcs; ++i) {
      getHessian(c1, c0, i, H); // compute the hessian of the ith constraint function and store in H
      double lambda = c1[nn[nNodes+i]].x;
      for(j = 0; j < m; ++j)
        for(k = 0; k < m; ++k) 
          Ktan[j][k] += lambda*H[j][k];
      for(j = 0; j < mpcs[i]->nterms; ++j) {
        for(k = 0; k < nNodes; ++k) if ( mpcs[i]->terms[j].nnum == nn[k]) break;
        l = k*nDofsPerNode + mpcs[i]->terms[j].dofnum;
        Ktan[l][m+i] = Ktan[m+i][l] = mpcs[i]->terms[j].coef.r_value;
        f[l] += lambda*mpcs[i]->terms[j].coef.r_value;
      }
      f[m+i] = mpcs[i]->rhs.r_value; // value of the ith constraint function
    }
  }
  //cerr << "here in RigidMpcElement::getStiffAndForce\n";
  //Ktan.print();
}

bool 
RigidMpcElement::isRigidMpcElement(const DofSet& ds, bool forAllNodes)
{ 
  return  ds == DofSet::nullDofset || ds == nodalDofs;
}

void 
RigidMpcElement::updateLMPCs(GeomState&, CoordSet&) 
{ 
  //cerr << "here in RigidMpcElement::updateLMPCs(GeomState&, CoordSet&)\n";
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i]->rhs.r_value = 0.0;
}

void 
RigidMpcElement::getHessian(GeomState&, CoordSet&, int, FullSquareMatrix& H) 
{ 
  H.zero(); 
}

