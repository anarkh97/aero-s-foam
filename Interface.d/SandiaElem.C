#include <Interface.d/SandiaElem.h>
StructProp refProp;

#include <Element.d/Truss.d/TwoNodeTruss.C>
StructProp sp;

SandiaElem::SandiaElem(int nn, int *nd, unsigned short int *dm)
{
  prop = &refProp;
  nnodes = nn;
  // nds    = nd;
  nds = new int[nnodes]; for(int i=0; i<nnodes; ++i) nds[i] = nd[i];
  flags  = dm;
  //proxy = new TwoNodeTruss(nd);
  int dummy[2] = { 0, 0 }; 
  proxy = new TwoNodeTruss(dummy); // PJSA: 12-6-05 nnodes can be < 2 !!!
  sp.E = 1e7;
  sp.nu = 0.3;
  proxy->setProp(&sp);
} 

int
SandiaElem::numNodes() 
{
  return nnodes;
}

bool 
SandiaElem::isSafe() 
{
  if(nnodes == 2)
    if((numDofs()/numNodes()) == 6)
      return true; // beam element
    else 
      return false;  // truss or spring element
  else 
    return true;
}

bool
SandiaElem::isRotMidSideNode(int iNode)
{
  // this returns true for midside nodes with rotational dofs (used in CornerMaker for rotCorners)
  if(nnodes == 6 && numDofs() == 36) { // 6 node tri Shell
    if(iNode < 3) return false;
    else return true;
  }
  else if((nnodes == 8) && (numDofs() == 48)) { // 8 node quad shell
    if(iNode < 4) return false;
    else return true;
  }
  else return false; 
}

int *
SandiaElem::nodes(int *nd)
{
  if(nd == 0) { nd = new int[nnodes]; }
  for(int i=0; i < nnodes; ++i) nd[i] = nds[i];
  return nd;
}

void
SandiaElem::markDofs(DofSetArray &dsa)
{
  for(int i=0; i < nnodes; ++i) 
    dsa.mark(nds[i], int(flags[nds[i]]));
}

void
SandiaElem::renum(int *) { /* not implemented */ }

FullSquareMatrix
SandiaElem::stiffness(CoordSet&cs, double *d, int f) 
{
  return proxy->stiffness(cs,d,f);
}

FullSquareMatrix
SandiaElem::massMatrix(CoordSet&, double *,int) 
{
  return FullSquareMatrix(1);
}

int*
SandiaElem::dofs(DofSetArray &ds, int *p)
{
  p = proxy->dofs(ds,p);
  return p;
}

int
SandiaElem::numDofs()
{
  int tot = 0;
  for(int i=0; i < nnodes; ++i)
    tot += DofSet(int(flags[nds[i]])).count();
  return tot;
}

