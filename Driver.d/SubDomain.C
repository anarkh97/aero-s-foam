#include <typeinfo>
#include <cstdio>
#ifdef SUN10
#include <typeinfo.h>
#endif 
#include <cmath>
#include <Utils.d/dbg_alloca.h>

#include <climits> //--- UH

#include <Element.d/Element.h>
#include <Utils.d/dofset.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/VectorSet.h>
#include <Math.d/matrix.h>
#include <Math.d/SparseSet.h>
#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Math.d/MpcSparse.h>
#include <Corotational.d/TemperatureState.h>
#include <Solvers.d/SolverFactory.h>

//#define DEBUG_MPC

extern int verboseFlag;
extern int isFeti3;
extern Domain * domain;
extern int salinasFlag;

extern "C" {
  void _FORTRAN(dgemm)(const char &, const char &, const int &,const int &,
                       const int &, const double &, double *, const int &,
                       double *, const int &, const double &, double *, const int &);

  void _FORTRAN(zgemm)(const char &, const char &, const int &,const int &,
                       const int &, const complex<double> &, complex<double> *, const int &,
                       complex<double> *, const int &, const complex<double> &, complex<double> *,
                       const int &);

  void _FORTRAN(dgemv)(const char &, const int &,const int &,
                       const double &, double *, const int &,
                       double *, const int &, const double &, double *, const int &);

  void _FORTRAN(zgemv)(const char &, const int &,const int &,
                       const complex<double> &, complex<double> *, const int &,
                       complex<double> *, const int &, const complex<double> &, complex<double> *,
                       const int &);
}

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(int sn, int lsn) :
  Domain(), // virtual base first
  BaseSub()
{
  initialize();
  subNumber      =  sn;
  localSubNumber = lsn;
}

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(Domain &dom, int sn, Connectivity &con, Connectivity &nds, int gn) :
  Domain(dom, con.num(gn), con[gn], nds.num(gn), nds[gn]),  // virtual base first
  BaseSub(dom, sn, con, nds, gn)
{
  initialize();
}

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(Domain &dom, int sn, int nNodes, int *nds, int nElems, int *elems, int gn) :
  Domain(dom, nElems, elems, nNodes, nds), //virtual base first
  BaseSub(dom, sn, nNodes, nds, nElems, elems, gn)
{
  initialize();
}

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(Domain &dom, int sn, CoordSet* _nodes, Elemset* _elems, int *glNodeNums, int *glElemNums, int gn) :
  Domain(dom, _elems, _nodes), //virtual base first
  BaseSub(dom, sn, _nodes, _elems, glNodeNums, glElemNums, gn)
{
  initialize();
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeAllDisp(Scalar (*xyz)[11], Scalar *u, Scalar (*xyz_loc)[11])
{
 // this is for either Helmholtz or other solver types
 // note: coupledScaling always has a default value of 1.0
 // xyz should be initialized to zero before being passed into this function
 // note: u is already scaled
 int inode, nodeI;
 for(inode = 0; inode < numnodes; ++inode){
   nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]]-1 : glNums[inode];

   int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);
   int xLoc1 =   dsa->locate(inode, DofSet::Xdisp);

   if(xLoc >= 0)
     xyz[nodeI][0] = u[xLoc];           // free
   else if(xLoc1 >= 0)
     xyz[nodeI][0] = Bcx(xLoc1);	// constrained

   int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);
   int yLoc1 =   dsa->locate(inode, DofSet::Ydisp);

   if(yLoc >= 0)
     xyz[nodeI][1] = u[yLoc];
   else if(yLoc1 >= 0)
     xyz[nodeI][1] = Bcx(yLoc1);

   int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);
   int zLoc1 =   dsa->locate(inode, DofSet::Zdisp);

   if(zLoc >= 0)
     xyz[nodeI][2] = u[zLoc];
   else if(zLoc1 >= 0)
     xyz[nodeI][2] = Bcx(zLoc1);

   int xRot  = c_dsa->locate(inode, DofSet::Xrot);
   int xRot1 =   dsa->locate(inode, DofSet::Xrot);

   if(xRot >= 0)
     xyz[nodeI][3] = u[xRot];
   else if(xRot1 >= 0)
     xyz[nodeI][3] = Bcx(xRot1);

   int yRot  = c_dsa->locate(inode, DofSet::Yrot);
   int yRot1 =   dsa->locate(inode, DofSet::Yrot);

   if(yRot >= 0)
     xyz[nodeI][4] = u[yRot];
   else if(yRot1 >= 0)
     xyz[nodeI][4] = Bcx(yRot1);

   int zRot  = c_dsa->locate(inode, DofSet::Zrot);
   int zRot1 =   dsa->locate(inode, DofSet::Zrot);

   if(zRot >= 0)
     xyz[nodeI][5] = u[zRot];
   else if(zRot1 >= 0)
     xyz[nodeI][5] = Bcx(zRot1);

   int xTemp  = c_dsa->locate(inode, DofSet::Temp);
   int xTemp1 =   dsa->locate(inode, DofSet::Temp);

   if(xTemp >= 0)
     xyz[nodeI][6] = u[xTemp];
   else if(xTemp1 >= 0)
     xyz[nodeI][6] = Bcx(xTemp1);

   int xHelm  = c_dsa->locate(inode, DofSet::Helm);
   int xHelm1 =   dsa->locate(inode, DofSet::Helm);

   if(xHelm >= 0)
     xyz[nodeI][7] = u[xHelm];
   else if(xHelm1 >= 0)
     xyz[nodeI][7] = Bcx(xHelm1);

   // transform displacements and rotations (if present) from DOF_FRM to basic coordinates
   // and keep a copy 
   if(!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
      && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
     if(xyz_loc) for(int j=0; j<11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
     bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
     transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
   }
 }
}


template<class Scalar>
void
GenSubDomain<Scalar>::forceContinuity(Scalar *u, Scalar (*xyz)[11])//DofSet::max_known_nonL_dof
{
 int inode;
 for(inode = 0; inode < numnodes; ++inode){
   int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);

   if(xLoc >= 0)
     u[xLoc] = xyz[glNums[inode]][0];

   int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);

   if(yLoc >= 0)
     u[yLoc] = xyz[glNums[inode]][1];

   int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);

   if(zLoc >= 0)
     u[zLoc] = xyz[glNums[inode]][2];

   int xRot  = c_dsa->locate(inode, DofSet::Xrot);

   if(xRot >= 0)
     u[xRot] = xyz[glNums[inode]][3];

   int yRot  = c_dsa->locate(inode, DofSet::Yrot);

   if(yRot >= 0)
     u[yRot] = xyz[glNums[inode]][4];

   int zRot  = c_dsa->locate(inode, DofSet::Zrot);

   if(zRot >= 0)
     u[zRot] = xyz[glNums[inode]][5];

   int xTemp  = c_dsa->locate(inode, DofSet::Temp);

   if(xTemp >= 0)
     u[xTemp] = xyz[glNums[inode]][6];

   int xHelm  = c_dsa->locate(inode, DofSet::Helm);

   if(xHelm >= 0)
     u[xHelm] = xyz[glNums[inode]][7];
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeAllVeloc(Scalar (*xyz)[11], Scalar *v, Scalar (*xyz_loc)[11])
{
 // xyz should be initialized to zero before being passed into this function
 int inode, nodeI;
 for(inode = 0; inode < nodes.size(); ++inode){
   if(nodes[inode] == NULL) continue;
   nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]]-1 : glNums[inode];

   int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);
   int xLoc1 =   dsa->locate(inode, DofSet::Xdisp);

   if(xLoc >= 0)
     xyz[nodeI][0] = v[xLoc];           // free
   else if(xLoc1 >= 0)
     xyz[nodeI][0] = vcx[xLoc1];        // constrained

   int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);
   int yLoc1 =   dsa->locate(inode, DofSet::Ydisp);

   if(yLoc >= 0)
     xyz[nodeI][1] = v[yLoc];
   else if(yLoc1 >= 0)
     xyz[nodeI][1] = vcx[yLoc1];

   int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);
   int zLoc1 =   dsa->locate(inode, DofSet::Zdisp);

   if(zLoc >= 0)
     xyz[nodeI][2] = v[zLoc];
   else if(zLoc1 >= 0)
     xyz[nodeI][2] = vcx[zLoc1];

   int xRot  = c_dsa->locate(inode, DofSet::Xrot);
   int xRot1 =   dsa->locate(inode, DofSet::Xrot);

   if(xRot >= 0)
     xyz[nodeI][3] = v[xRot];
   else if(xRot1 >= 0)
     xyz[nodeI][3] = vcx[xRot1];

   int yRot  = c_dsa->locate(inode, DofSet::Yrot);
   int yRot1 =   dsa->locate(inode, DofSet::Yrot);

   if(yRot >= 0)
     xyz[nodeI][4] = v[yRot];
   else if(yRot1 >= 0)
     xyz[nodeI][4] = vcx[yRot1];

   int zRot  = c_dsa->locate(inode, DofSet::Zrot);
   int zRot1 =   dsa->locate(inode, DofSet::Zrot);

   if(zRot >= 0)
     xyz[nodeI][5] = v[zRot];
   else if(zRot1 >= 0)
     xyz[nodeI][5] = vcx[zRot1];

   int xTemp  = c_dsa->locate(inode, DofSet::Temp);
   int xTemp1 =   dsa->locate(inode, DofSet::Temp);

   if(xTemp >= 0)
     xyz[nodeI][6] = v[xTemp];
   else if(xTemp1 >= 0)
     xyz[nodeI][6] = vcx[xTemp1];

   int xHelm  = c_dsa->locate(inode, DofSet::Helm);
   int xHelm1 =   dsa->locate(inode, DofSet::Helm);

   if(xHelm >= 0)
     xyz[nodeI][7] = v[xHelm];
   else if(xHelm1 >= 0)
     xyz[nodeI][7] = vcx[xHelm1];

   // transform velocities and angular velocities (if present) from DOF_FRM to basic coordinates
   if(!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
      && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
     if(xyz_loc) for(int j=0; j<11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
     bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
     transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeAllAccel(Scalar (*xyz)[11], Scalar *a, Scalar (*xyz_loc)[11])
{
 // xyz should be initialized to zero before being passed into this function
 int inode, nodeI;
 for(inode = 0; inode < nodes.size(); ++inode){
   if(nodes[inode] == NULL) continue;
   nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]]-1 : glNums[inode];

   int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);
   int xLoc1 =   dsa->locate(inode, DofSet::Xdisp);

   if(xLoc >= 0)
     xyz[nodeI][0] = a[xLoc];           // free
   else if(xLoc1 >= 0)
     xyz[nodeI][0] = acx[xLoc1];        // constrained

   int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);
   int yLoc1 =   dsa->locate(inode, DofSet::Ydisp);

   if(yLoc >= 0)
     xyz[nodeI][1] = a[yLoc];
   else if(yLoc1 >= 0)
     xyz[nodeI][1] = acx[yLoc1];

   int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);
   int zLoc1 =   dsa->locate(inode, DofSet::Zdisp);

   if(zLoc >= 0)
     xyz[nodeI][2] = a[zLoc];
   else if(zLoc1 >= 0)
     xyz[nodeI][2] = acx[zLoc1];

   int xRot  = c_dsa->locate(inode, DofSet::Xrot);
   int xRot1 =   dsa->locate(inode, DofSet::Xrot);

   if(xRot >= 0)
     xyz[nodeI][3] = a[xRot];
   else if(xRot1 >= 0)
     xyz[nodeI][3] = acx[xRot1];

   int yRot  = c_dsa->locate(inode, DofSet::Yrot);
   int yRot1 =   dsa->locate(inode, DofSet::Yrot);

   if(yRot >= 0)
     xyz[nodeI][4] = a[yRot];
   else if(yRot1 >= 0)
     xyz[nodeI][4] = acx[yRot1];

   int zRot  = c_dsa->locate(inode, DofSet::Zrot);
   int zRot1 =   dsa->locate(inode, DofSet::Zrot);

   if(zRot >= 0)
     xyz[nodeI][5] = a[zRot];
   else if(zRot1 >= 0)
     xyz[nodeI][5] = acx[zRot1];

   int xTemp  = c_dsa->locate(inode, DofSet::Temp);
   int xTemp1 =   dsa->locate(inode, DofSet::Temp);

   if(xTemp >= 0)
     xyz[nodeI][6] = a[xTemp];
   else if(xTemp1 >= 0)
     xyz[nodeI][6] = acx[xTemp1];

   int xHelm  = c_dsa->locate(inode, DofSet::Helm);
   int xHelm1 =   dsa->locate(inode, DofSet::Helm);

   if(xHelm >= 0)
     xyz[nodeI][7] = a[xHelm];
   else if(xHelm1 >= 0)
     xyz[nodeI][7] = acx[xHelm1];

   // transform accelerations and angular accelerations (if present) from DOF_FRM to basic coordinates
   if(!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
      && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
     if(xyz_loc) for(int j=0; j<11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
     bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
     transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
   }
 }
}

template<class Scalar>
void GenSubDomain<Scalar>::addUserForce(Scalar *extForce, Scalar *usrDefForce)
{
  int i;
  for(i = 0; i < claw->numUserForce; ++i)  {
    int dof = c_dsa->locate(claw->userForce[i].nnum,1<<claw->userForce[i].dofnum);
    if (dof > -1)
      extForce[dof] += usrDefForce[locToGlUserForceMap[i]]/weight[dof];
  }
}

template<class Scalar>
void GenSubDomain<Scalar>::addCtrl(Scalar *force, Scalar *ctrfrc)
{
  int i;
  for(i = 0; i < claw->numActuator; ++i) {
    int dof = c_dsa->locate(claw->actuator[i].nnum,1 << claw->actuator[i].dofnum);
    if(dof > -1)
      force[dof] += ctrfrc[locToGlActuatorMap[i]]/weight[dof];
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeStress(Scalar *locStress,  Scalar *locWeight,
                                  Scalar *globStress, Scalar *globWeight, int glNumNodes)
{
 int inode, nodeI;
 for(inode = 0; inode < numnodes; ++inode) {
   nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]]-1 : glNums[inode];
   if(nodeI >= glNumNodes) continue;
   globWeight[nodeI] += locWeight[inode];
   globStress[nodeI] += locStress[inode];
 }
}

template<class Scalar>
void GenSubDomain<Scalar>::mergeElemStress(Scalar *locStress, Scalar *globStress, Connectivity *glElemToNode)
{
  int *glOffset = glElemToNode->ptr();
  int *locOffset = elemToNode->ptr();
  for(int iElem = 0; iElem < numele; iElem++) {
    for(int iNode = 0; iNode < packedEset[iElem]->numNodes(); iNode++) {
      int glOff = glOffset[glElems[iElem]]+iNode;
      int locOff = locOffset[iElem]+iNode;
      globStress[glOff] = locStress[locOff];
    }
  }
}

inline double square(double x) { return x*x; }

template<class Scalar>
void
GenSubDomain<Scalar>::mergePrimalError(Scalar* error, Scalar* primal)
{
  for(int inode=0; inode<numnodes; ++inode) {
    Scalar nd   = 0.0;
    Scalar totP = 0.0;
    int xLoc = c_dsa->locate(inode, DofSet::Xdisp);
    if(xLoc >= 0) { totP += square(primal[xLoc]/weight[xLoc]); nd += 1.0 ; }
    int yLoc = c_dsa->locate(inode, DofSet::Ydisp);
    if(yLoc >= 0) { totP += square(primal[yLoc]/weight[yLoc]); nd += 1.0 ; }
    int zLoc = c_dsa->locate(inode, DofSet::Zdisp);
    if(zLoc >= 0) { totP += square(primal[zLoc]/weight[zLoc]); nd += 1.0 ; }
    if(nd != 0) totP = sqrt(totP/nd);
    error[glNums[inode]] += totP;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::addDMass(int glNum, int dof, double v)
{
  Domain::addDMass(glNum, dof, v);
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendDOFList(FSCommPattern<int> *pat)
{
 Connectivity &sharedNodes = *(scomm->sharedNodes);
 for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
   FSSubRecInfo<int> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
   for(int iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
     if((solInfo().isCoupled && isWetInterfaceNode(sharedNodes[iSub][iNode])) &&
        (subNumber == scomm->subNums[iSub]) && (solInfo().solvercntl->fetiInfo.fsi_corner == 0))
       sInfo.data[iNode] = wetInterfaceDofs[wetInterfaceNodeMap[sharedNodes[iSub][iNode]]].list();
     else
       sInfo.data[iNode] = (*c_dsa)[sharedNodes[iSub][iNode]].list();
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::gatherDOFList(FSCommPattern<int> *pat)
{
  int iSub, iNode, i;
  Connectivity &sharedNodes = *(scomm->sharedNodes); // contains both dry boundary nodes and wet interface nodes
  boundaryDOFs = new DofSet * [scomm->numNeighb];
  int nbdofs = 0, nwdofs = 0;
  int nbneighb = 0, nwneighb = 0;
  bool isbneighb, iswneighb;
  bool isCoupled = solInfo().isCoupled;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    FSSubRecInfo<int> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
    boundaryDOFs[iSub] = new DofSet[sharedNodes.num(iSub)];
    isbneighb = false; iswneighb = false;
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      boundaryDOFs[iSub][iNode] = DofSet(rInfo.data[iNode]); // temporarily store neighb c_dsa in boundaryDOFs
                                                             // or w dofs if neighb is self
      if(isCoupled && (isWetInterfaceNode(sharedNodes[iSub][iNode]))) { // wet interface node
        DofSet shared_wdofs = boundaryDOFs[iSub][iNode] & wetInterfaceDofs[wetInterfaceNodeMap[sharedNodes[iSub][iNode]]];
        int wcount = shared_wdofs.count();
        if(wcount > 0) { nwdofs += wcount; iswneighb = true; }
      }
      DofSet shared_rdofs = boundaryDOFs[iSub][iNode] & (*getCCDSA())[sharedNodes[iSub][iNode]];
      int bcount = shared_rdofs.count();
      nbdofs += bcount; isbneighb = true;  // make every "sharedNode" neighbor a "std sharedDOF neighb" also (simplifies augmentation)
    }
    if(iswneighb) nwneighb++;
    if(isbneighb) nbneighb++;
  }

  int *boundDofs    = new int[nbdofs];
  int *boundDofPointer = new int[nbneighb+1]; boundDofPointer[0] = 0;
  int *boundNeighbs = new int[nbneighb];
  int *wetDofs    = new int[nwdofs];
  int *wetDofPointer = new int[nwneighb+1]; wetDofPointer[0] = 0;
  int *wetNeighbs = new int[nwneighb];
  nbdofs = 0; nwdofs = 0;
  nbneighb = 0; nwneighb = 0;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    isbneighb = false; iswneighb = false;
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      if(isCoupled && (isWetInterfaceNode(sharedNodes[iSub][iNode]))) { // wet interface node
        DofSet shared_wdofs = boundaryDOFs[iSub][iNode] & wetInterfaceDofs[wetInterfaceNodeMap[sharedNodes[iSub][iNode]]];
        int dofs[7]; //, cdofs[7];
        dsa->number(sharedNodes[iSub][iNode], shared_wdofs, dofs);
        int wcount = shared_wdofs.count();
        if(wcount > 0) {
          for(int i=0; i<wcount; ++i) wetDofs[nwdofs++] = wetInterfaceMap[dofs[i]]; // WHAT ABOUT CONSTRAINTS ???
          iswneighb = true;
        }
      }
      boundaryDOFs[iSub][iNode] &= (*getCCDSA())[sharedNodes[iSub][iNode]]; // ~> shared_rdofs
      int bcount = getCCDSA()->number(sharedNodes[iSub][iNode],
                                      boundaryDOFs[iSub][iNode], boundDofs + nbdofs);
      nbdofs += bcount; isbneighb = true; // make every "sharedNode" neighbor a "std sharedDOF neighb" also (simplifies augmentation)
    }
    if(iswneighb) { wetNeighbs[nwneighb++] = scomm->subNums[iSub]; wetDofPointer[nwneighb] = nwdofs; }
    if(isbneighb) { boundNeighbs[nbneighb++] = scomm->subNums[iSub]; boundDofPointer[nbneighb] = nbdofs; }
  }

  Connectivity *stdSharedDOFs = new Connectivity(nbneighb, boundDofPointer, boundDofs);
  if(!solInfo().solvercntl->fetiInfo.bmpc) scomm->setTypeSpecificList(SComm::std, boundNeighbs, stdSharedDOFs);
  Connectivity *wetSharedDOFs = new Connectivity(nwneighb, wetDofPointer, wetDofs);
  scomm->setTypeSpecificList(SComm::wet, wetNeighbs, wetSharedDOFs);

  int ndof = localLen();
  weight = new int[ndof];
  for(i = 0; i < ndof; ++i) weight[i] = 1;
  for(i = 0; i < stdSharedDOFs->numConnect(); ++i)
    weight[(*stdSharedDOFs)[0][i]] += 1;
  if(numWIdof) {
    wweight = new Scalar[numWIdof];
    for(i=0; i<numWIdof; ++i) wweight[i] = 0.0;
    for(i = 0; i < scomm->lenT(SComm::wet); ++i)
      wweight[scomm->boundDofT(SComm::wet,i)] += 1.0;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::gatherDOFListPlus(FSCommPattern<int> *pat)
{
  int iSub, iNode, i;
  Connectivity &sharedNodes = *(scomm->sharedNodes);
  DofSet **boundaryDOFs = new DofSet * [scomm->numNeighb]; // needs to be a local variable
  int nbdofs = 0;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    FSSubRecInfo<int> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
    boundaryDOFs[iSub] = new DofSet[sharedNodes.num(iSub)];
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      boundaryDOFs[iSub][iNode] = DofSet(rInfo.data[iNode]); // temporarily store neighb c_dsa in boundaryDOFs
      DofSet shared_bdofs = boundaryDOFs[iSub][iNode] & (*getCDSA())[sharedNodes[iSub][iNode]];
      int bcount = shared_bdofs.count();
      nbdofs += bcount;
    }
  }

  int *boundDofs    = new int[nbdofs];
  int *boundDofPointer = new int[scomm->numNeighb+1]; boundDofPointer[0] = 0;
  nbdofs = 0;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      boundaryDOFs[iSub][iNode] &= (*getCDSA())[sharedNodes[iSub][iNode]]; // ~> shared_bdofs
      int bcount = getCDSA()->number(sharedNodes[iSub][iNode],
                                     boundaryDOFs[iSub][iNode], boundDofs + nbdofs);
      nbdofs += bcount;
    }
    boundDofPointer[iSub+1] = nbdofs;
  }

  Connectivity *sharedDOFsPlus = new Connectivity(scomm->numNeighb, boundDofPointer, boundDofs);
  scomm->sharedDOFsPlus = sharedDOFsPlus;

  weightPlus = new int[c_dsa->size()];
  for(i = 0; i < c_dsa->size(); ++i) weightPlus[i] = 1;
  for(i = 0; i < sharedDOFsPlus->numConnect(); ++i)
    weightPlus[(*sharedDOFsPlus)[0][i]] += 1;

  for(int i=0; i<scomm->numNeighb; ++i) delete [] boundaryDOFs[i];
  delete [] boundaryDOFs;
}

template<class Scalar>
void
GenSubDomain<Scalar>::applySplitting()
{
  // adjust discrete masses, forces and mpcs using subdomain multiplicity
  applyDmassSplitting();
  applyForceSplitting();
  applyMpcSplitting();
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyDmassSplitting()
{
  // adjust discrete masses using subdomain multiplicity
  // num = number of subdomains touching a dof
  int cdof, num;

  // discrete masses
  DMassData *cmass = firstDiMass;
  while(cmass != 0) {
    if((cdof = c_dsa->locate(cmass->node, (1 << cmass->dof))) > -1 && (num = weightPlus[cdof]) > 1)
      cmass->diMass /= double(num);
    cmass = cmass->next;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyForceSplitting()
{
  // adjust forces using subdomain multiplicity
  // num = number of subdomains touching a dof
  int cdof, num;

  // forces
  for(int i = 0; i < numNeuman; ++i) {
    if((cdof = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum))) > -1 && (num = weightPlus[cdof]) > 1)
      nbc[i].val /= double(num);
  }
  for(int i = 0; i < numComplexNeuman; ++i) {
    if((cdof = c_dsa->locate(cnbc[i].nnum, (1 << cnbc[i].dofnum))) > -1 && (num = weightPlus[cdof]) > 1) {
      cnbc[i].reval /= double(num);
      cnbc[i].imval /= double(num);
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyMpcSplitting()
{
  // adjust discrete masses, forces and mpcs using subdomain multiplicity
  // num = number of subdomains touching a dof
  int cdof, num;

  // mpcs (NOTE: optional kscaling is done later, hhs is not split)
  if(solInfo().getFetiInfo().mpc_scaling == FetiInfo::tscaling) {
    for(int iMPC = 0; iMPC < numMPC; ++iMPC) { // dual mpcs
      if(mpc[iMPC]->type == 2) continue; // bmpc
      for(int i = 0; i < mpc[iMPC]->nterms; ++i) {
        if((cdof = mpc[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
          mpc[iMPC]->terms[i].coef /= double(num);
      }
    }
  }
  // XXXX kscaling currently not supported for primal mpcs
    for(int iMPC = 0; iMPC < numMPC_primal; ++iMPC) { // primal mpcs
      for(int i = 0; i < mpc_primal[iMPC]->nterms; ++i) {
        if((cdof = mpc_primal[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
          mpc_primal[iMPC]->terms[i].coef /= double(num);
      }
    }

}

template<class Scalar>
void
GenSubDomain<Scalar>::initScaling()
{
  int i;
  if(scaling) delete [] scaling;
  scaling = new Scalar[totalInterfSize];
  for(i = 0; i < totalInterfSize; ++i) scaling[i] = Scalar(1.0);
  if(solInfo().getFetiInfo().scaling == FetiInfo::tscaling)
    for(i=0; i < scomm->lenT(SComm::std); ++i)
      scaling[scomm->mapT(SComm::std,i)] = 1.0/weight[scomm->boundDofT(SComm::std,i)];

  if(solInfo().getFetiInfo().fsi_scaling == FetiInfo::tscaling)
    for(i=0; i < scomm->lenT(SComm::wet); ++i)
      scaling[scomm->mapT(SComm::wet,i)] = 1.0/wweight[scomm->wetDofNb(i)];

  if(solInfo().getFetiInfo().mpc_precno == FetiInfo::diagCCt)
    for(i=0; i < scomm->lenT(SComm::mpc); ++i)
      scaling[scomm->mapT(SComm::mpc,i)] = 1.0;
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendDiag(GenSparseMatrix<Scalar> *s, FSCommPattern<Scalar> *vPat)
{
  int iDof = 0;
  for(int i = 0; i < scomm->numT(SComm::all); ++i) {
    FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNumber, scomm->neighbT(SComm::all,i));
    for(int  j = 0; j < scomm->lenT(SComm::all,i); ++j) {
      switch(boundDofFlag[iDof]) {
        case 0:
          scaling[iDof] = sInfo.data[j] = (s) ? s->diag(scomm->boundDofT(SComm::all,i,j)) : 1.0;
          break;
        case 1: // wet interface
          scaling[iDof] = sInfo.data[j] = (Kww) ? Kww->diag(-1-scomm->boundDofT(SComm::all,i,j)) : 1.0;
          break;
        case 2:  // dual mpc
          scaling[iDof] = sInfo.data[j] = 1.0;
          break;
      }
      iDof++;
    }
  }

  int ndof = localLen();
  // use the kweight array also for LMPCs stiffness scaling/splitting for the primal method
  if((solInfo().getFetiInfo().augment == FetiInfo::WeightedEdges) ||
     ((solInfo().getFetiInfo().mpc_scaling == FetiInfo::kscaling) && (numMPC_primal > 0))) {
    kweight = new Scalar[ndof];  // used for WeightedEdges augmentation
    for(iDof = 0; iDof < ndof; ++iDof)
      kweight[iDof] = (s) ? s->diag(iDof) : 1.0;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::fSend(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw)
{
  // Get the trace of the subdomain interfaces
  int iDof = 0;
  for(int i = 0; i < scomm->numT(SComm::all); ++i) {
    FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNumber, scomm->neighbT(SComm::all,i));
    for(int  j = 0; j < scomm->lenT(SComm::all,i); ++j) {
      int bdof = scomm->boundDofT(SComm::all,i,j);
      switch(boundDofFlag[iDof++]) {
        case 0:
          sInfo.data[j] = locF[bdof];
          break;
        case 1: // wet interface
          sInfo.data[j] = locFw[-1-bdof];
          break;
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::collectScaling(FSCommPattern<Scalar> *vPat)
{
 int offset = 0;
 int locLen = localLen();
 //Scalar *kSum = (Scalar *) dbg_alloca(sizeof(Scalar)*locLen);
 Scalar *kSum = new Scalar[locLen];
 for(int i = 0; i < locLen; ++i) kSum[i] = 0.0;
#ifdef HB_COUPLED_PRECOND
 kSumWI = new Scalar[numWIdof]; //stores the sum of the neighbourg stiffess otherwise 0.0 if not shared
#else
 Scalar* kSumWI = (Scalar *)dbg_alloca(numWIdof*sizeof(Scalar));
#endif
 bool* wflag = (bool*)dbg_alloca(numWIdof*sizeof(bool)); //HB: to avoid setting wweight multiple times
 for(int i = 0; i < numWIdof; ++i) { kSumWI[i] = 0.0; wflag[i] = false; }

 int iSub, iDof;
 for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
    if(subNumber != scomm->neighbT(SComm::all,iSub)) {
      FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all,iSub), subNumber);
      for(iDof = 0; iDof < scomm->lenT(SComm::all,iSub); ++iDof) {
        int bdof = scomm->boundDofT(SComm::all,iSub,iDof);
        if(bdof >= 0)
          kSum[bdof] += rInfo.data[iDof];
        else if(boundDofFlag[offset+iDof] == 1)
          kSumWI[-1-bdof] += rInfo.data[iDof];
      }
    }
    offset += scomm->lenT(SComm::all,iSub);
 }

 offset = 0;
 for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
    FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all,iSub), subNumber);
    for(iDof = 0; iDof < scomm->lenT(SComm::all,iSub); ++iDof) {
      int bdof = scomm->boundDofT(SComm::all,iSub,iDof);
      if(bdof >= 0)
        scaling[offset+iDof] = rInfo.data[iDof]/(scaling[offset+iDof] + kSum[bdof]);
      else if(boundDofFlag[offset+iDof] == 1) {
        scaling[offset+iDof] = scaling[offset+iDof]/(scaling[offset+iDof] + kSumWI[-1-bdof]);
        if((solInfo().getFetiInfo().fsi_scaling == FetiInfo::kscaling) & !wflag[-1-bdof]) {
          wweight[-1-bdof]= (1./scaling[offset+iDof])/wweight[-1-bdof];
          wflag[-1-bdof] = true;
        }
      }
    }
    offset += scomm->lenT(SComm::all,iSub);
 }

 // LMPCs coeff stiffness scaling/splitting for the primal method
/* XXXX kscaling currently not supported for primal mpcs
 if(solInfo().getFetiInfo().mpc_scaling == FetiInfo::kscaling) {
   for(int iMPC = 0; iMPC < numMPC_primal; iMPC++){
     for(int k = 0; k < mpc_primal[iMPC]->nterms; k++){
       int ccdof = (mpc_primal[iMPC]->terms)[k].ccdof;
       if(ccdof>=locLen) fprintf(stderr, "Strange: cdof (%d) >= locLen (%d)\n", ccdof, locLen);
       else if(ccdof >= 0) // rdof
         (mpc_primal[iMPC]->terms)[k].coef *= kweight[ccdof]/(kweight[ccdof]+kSum[ccdof]);
       else if(wetInterfaceMap) {
         int dof = (mpc_primal[iMPC]->terms)[k].dof;
         if((dof > 0) && (wetInterfaceMap[dof] > 0)) // wdof
           mpc_primal[iMPC]->terms[k].coef /= ScalarTypes::Real(wweight[wetInterfaceMap[dof]]);
       }
       // note corner dofs are always split with tscaling - see gatherDOFList
     }
   }
 }
*/
 if(solInfo().getFetiInfo().augment == FetiInfo::WeightedEdges) {
   for(iDof = 0; iDof < locLen; ++iDof)
     kweight[iDof] += kSum[iDof];
 }
 delete [] kSum;
}

template<class Scalar>
void
GenSubDomain<Scalar>::fScale(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw)
{
 bool *isShared =  (bool *) dbg_alloca(sizeof(bool)*totalInterfSize);

 int iDof, iSub;
 int offset = 0;
 for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
   FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all,iSub), subNumber);
   for(iDof = 0; iDof < scomm->lenT(SComm::all,iSub); ++iDof) {
     int bdof = scomm->boundDofT(SComm::all,iSub,iDof);
     if(bdof >= 0)
       locF[bdof] += rInfo.data[iDof];
     else if(boundDofFlag[offset+iDof] == 1) {
       if(scomm->neighbT(SComm::all,iSub) != subNumber) { // wet interface
         locFw[-1-bdof] += rInfo.data[iDof];
         isShared[offset+iDof] = true;
       }
       else isShared[offset+iDof] = false;
     }
   }
   offset += scomm->lenT(SComm::all,iSub);
 }

 Scalar *interfF = (Scalar *) dbg_alloca(sizeof(Scalar)*totalInterfSize);

 // Get the trace of the subdomain interfaces
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   if(allBoundDofs[iDof] >= 0)
     interfF[iDof] = locF[allBoundDofs[iDof]]*scaling[iDof];
   else if((boundDofFlag[iDof] == 1) && isShared[iDof]) // wet interface
     interfF[iDof] = locFw[-1-allBoundDofs[iDof]]*scaling[iDof];

 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   if(allBoundDofs[iDof] >= 0)
     locF[allBoundDofs[iDof]] -= interfF[iDof];
   else if((boundDofFlag[iDof] == 1) && isShared[iDof])  // wet interface
     locFw[-1-allBoundDofs[iDof]] -= interfF[iDof];
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendMpcScaling(FSCommPattern<Scalar> *mpcPat)
{
  int i,j;
  diagCCt = new Scalar[numMPC]; // compute weights for mpcs also
  for(i = 0; i < numMPC; ++i) {
    diagCCt[i] = 0.0;
    for(j = 0; j < mpc[i]->nterms; ++j)
      diagCCt[i] += mpc[i]->terms[j].coef * mpc[i]->terms[j].coef / mpc[i]->k[j]; //HB: for mpc kscaling;
  }

  int neighb;
  for(i = 0; i < scomm->numT(SComm::mpc); ++i) {
    if(subNumber != (neighb = scomm->neighbT(SComm::mpc, i))) {
      FSSubRecInfo<Scalar> sInfo = mpcPat->getSendBuffer(subNumber, neighb);
      for(j = 0; j < scomm->lenT(SComm::mpc,i); ++j)
        sInfo.data[j] = diagCCt[scomm->mpcNb(i,j)];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::collectMpcScaling(FSCommPattern<Scalar> *mpcPat)
{
  int i, j;
  Scalar *mpcSum = (Scalar *) dbg_alloca(sizeof(Scalar)*numMPC);
  for(i = 0; i < numMPC; ++i) mpcSum[i] = 0.0;
  int neighb;
  for(i = 0; i < scomm->numT(SComm::mpc); ++i) {
    if(subNumber != (neighb = scomm->neighbT(SComm::mpc, i))) {
      FSSubRecInfo<Scalar> rInfo = mpcPat->recData(neighb, subNumber);
      for(j = 0; j < scomm->lenT(SComm::mpc, i); ++j)
        mpcSum[scomm->mpcNb(i,j)] += rInfo.data[j];
    }
  }
  for(i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(diagCCt[locMpcNb]+mpcSum[locMpcNb] != 0.0)
      scaling[scomm->mapT(SComm::mpc,i)] /= (diagCCt[locMpcNb] + mpcSum[locMpcNb]);
  }
  delete [] diagCCt; diagCCt = 0;
}


template<class Scalar>
void
GenSubDomain<Scalar>::fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec)
{
 // Add the interface vector (interfvec) contribution to localvec
 int iDof;
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   localvec[allBoundDofs[iDof]] += interfvec[iDof];

 // solve for localvec
 if(s) s->reSolve(localvec);

 // redistribute the solution to the interface
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   interfvec[iDof] = localvec[allBoundDofs[iDof]];

}

template<class Scalar>
void
GenSubDomain<Scalar>::fetiBaseOp(Scalar *uc, GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec)
{
 // localvec += Br^T * interfvec
 multAddBrT(interfvec, localvec);

 //Scalar inorm = 0.0; for(int i=0; i<interfLen(); ++i) inorm += interfvec[i]*interfvec[i];
 //Scalar lnorm = 0.0; for(int i=0; i<localRLen(); ++i) lnorm += localvec[i]*localvec[i];
 //cerr << "sub = " << subNumber << ", before: inorm = " << inorm << ", lnorm = " << lnorm;
 //if(subNumber == 7) { cerr << "interfvec = "; for(int i=0; i<interfLen(); ++i) cerr << interfvec[i] << " "; cerr << endl; }

 // localvec = Krr^-1 * localvec
 if(s) s->reSolve(localvec);

 // interfvec = Br * localvec
 multBr(localvec, interfvec, uc);

 //inorm = 0.0; for(int i=0; i<interfLen(); ++i) inorm += interfvec[i]*interfvec[i];
 //lnorm = 0.0; for(int i=0; i<localRLen(); ++i) lnorm += localvec[i]*localvec[i];
 //cerr << ", after: inorm = " << inorm << ", lnorm = " << lnorm << endl;
 //if(subNumber == 7) { cerr << "interfvec = "; for(int i=0; i<interfLen(); ++i) cerr << interfvec[i] << " "; cerr << endl; }
}

template<class Scalar>
void
GenSubDomain<Scalar>::fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec,
                                         FSCommPattern<Scalar> *wiPat)
{
 // localvec += Br^T * interfvec
 multAddBrT(interfvec, localvec, localw);

 // solve for localvec
 if(s) s->reSolve(localvec);

 if(numWIdof) {
   int i,j;
   // compute Kww uw for this subdomain
   for(i=0; i<numWIdof; ++i) localw_copy[i] = 0.0;
   Kww->mult(localw, localw_copy);  // localw_copy = - Kww * uw

// JLchange
// Fsi elements are already added to Kww, communication via neighborKww is no longer needed
   // compute Kww uw to send to neighbors
   if (solInfo().solvercntl->fetiInfo.fsi_corner == 0)
     for(i = 0; i < scomm->numT(SComm::fsi); ++i) {
       if(subNumber != scomm->neighbT(SComm::fsi,i)) {
         FSSubRecInfo<Scalar> sInfo = wiPat->getSendBuffer(subNumber, scomm->neighbT(SComm::fsi,i));
         for(j=0; j<numNeighbWIdof[i]; ++j) sInfo.data[j] = 0.0;
         neighbKww->multAdd(localw, sInfo.data, glToLocalWImap, neighbGlToLocalWImap[i]);
       }
       else {
         neighbKww->multAdd(localw, localw_copy, glToLocalWImap);
       }
     }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::fetiBaseOpCoupled2(Scalar *uc, Scalar *localvec, Scalar *interfvec,
                                         FSCommPattern<Scalar> *wiPat, Scalar *fw)
{
 // coupled_dph
 if(numWIdof) {
   int i, j;

// JLchange
// Fsi elements are already added to Kww, communication via neighborKww is no longer needed
   if (solInfo().solvercntl->fetiInfo.fsi_corner == 0)
     for(i = 0; i < scomm->numT(SComm::fsi); ++i) {
       if(subNumber != scomm->neighbT(SComm::fsi,i)) {
         FSSubRecInfo<Scalar> rInfo = wiPat->recData(scomm->neighbT(SComm::fsi,i), subNumber);
         for(j=0; j<numWIdof; ++j) localw_copy[j] += rInfo.data[j]/wweight[j];
       }
     }

   if(Krw) Krw->transposeMultSubNew(localvec, localw_copy); // localw_copy -= Krw^T * localvec
   int numCDofs = nCoarseDofs();
   Scalar *ucLocal = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);
   for(i=0; i<numCDofs; ++i) {
     if(cornerEqNums[i] > -1) ucLocal[i] = uc[cornerEqNums[i]];
     else ucLocal[i] = 0.0;
   }
   if(Kcw) Kcw->trMult(ucLocal, localw_copy, -1.0, 1.0);  // localw_copy -= Kcw^T Bc uc
   if(Kcw_mpc) Kcw_mpc->transposeMultSubtractWI(ucLocal, localw_copy);
   if(fw) for(i=0; i<numWIdof; ++i) localw_copy[i] += fw[i];  // localw_copy += fw
 }

 // interfvec = Br * localvec
 multBr(localvec, interfvec, uc, localw_copy);
}

template<class Scalar>
void
GenSubDomain<Scalar>::fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec,
                                 Scalar *beta)

{
 // multiply Q*beta to get the MPC force contribution
 if(numMPC > 0) {
   int i,iMPC;
   for(iMPC = 0; iMPC < numMPC; ++iMPC)
     for(i = 0; i < mpc[iMPC]->nterms; ++i) {
       int cdof = mpc[iMPC]->terms[i].cdof;
       if(cdof < 0) continue;
       localvec[cdof] += mpc[iMPC]->terms[i].coef * beta[localToGlobalMPC[iMPC]];
   }
 }

 // Add the interface vector (interfvec) contribution to localvec
 int iDof;
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   localvec[allBoundDofs[iDof]] += interfvec[iDof];

 // solve for localvec
 if(s) s->reSolve(localvec);

 // redistribute the solution to the interface
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   interfvec[iDof] = localvec[allBoundDofs[iDof]];
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat)
{
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[i]);
    for(int j = 0; j < scomm->sharedNodes->num(i); ++j) {
      int n = (*(scomm->sharedNodes))[i][j];
      for(int k = 0; k < 11; ++k) sInfo.data[11*j+k] = subvec[n][k];
    }
  }
}

#ifndef MAXABS
#define MAXABS(X,Y) ((ScalarTypes::norm(X) > ScalarTypes::norm(Y)) ? X : Y)
#endif

template<class Scalar>
void
GenSubDomain<Scalar>::collectNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat)
{
  // use the one with the largest absolute value
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[i], subNumber);
    for(int j = 0; j < scomm->sharedNodes->num(i); ++j) {
      int n = (*(scomm->sharedNodes))[i][j];
      for(int k = 0; k < 11; ++k) subvec[n][k] = MAXABS(subvec[n][k],rInfo.data[11*j+k]);
    }
  }
}


template<class Scalar>
void
GenSubDomain<Scalar>::sendInterf(Scalar *interfvec, FSCommPattern<Scalar> *vPat)
{
  int iDof = 0;
  for(int i = 0; i < scomm->numT(SComm::all); ++i) {
    FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNumber, scomm->neighbT(SComm::all,i));
    for(int j = 0; j < scomm->lenT(SComm::all,i); ++j) {
      sInfo.data[j] = interfvec[iDof++];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractAndSendInterf(Scalar *subvec, FSCommPattern<Scalar> *pat)
{
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[i]);
    for(int j = 0; j < scomm->sharedDOFsPlus->num(i); ++j) {
      sInfo.data[j] = subvec[(*scomm->sharedDOFsPlus)[i][j]];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleInterf(Scalar *subvec, FSCommPattern<Scalar> *pat)
{
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[i], subNumber);
    for(int j = 0; j < scomm->sharedDOFsPlus->num(i); ++j) {
      subvec[(*scomm->sharedDOFsPlus)[i][j]] += rInfo.data[j];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::splitInterf(Scalar *subvec)
{
 Scalar *interfF = (Scalar *)alloca(scomm->sharedDOFsPlus->numConnect()*sizeof(Scalar));
 for(int iDof = 0; iDof < scomm->sharedDOFsPlus->numConnect(); ++iDof)
   interfF[iDof] = subvec[(*scomm->sharedDOFsPlus)[0][iDof]]/double(weightPlus[(*scomm->sharedDOFsPlus)[0][iDof]]);

 for(int iDof = 0; iDof < scomm->sharedDOFsPlus->numConnect(); ++iDof)
   subvec[(*scomm->sharedDOFsPlus)[0][iDof]] -= interfF[iDof];
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleInterfInvert(Scalar *subvec, FSCommPattern<Scalar> *pat)
{
  if(numMPC) cerr << "ERROR: GenSubDomain<Scalar>::assembleInterfInvert(...) not implemented for MPCs \n";
  for(int i = 0; i < numUncon(); ++i) subvec[i] = 1.0/subvec[i];
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[i], subNumber);
    for(int j = 0; j < scomm->sharedDOFsPlus->num(i); ++j) {
      subvec[(*scomm->sharedDOFsPlus)[i][j]] += 1.0/rInfo.data[j];
    }
  }
  for(int i = 0; i < numUncon(); ++i) subvec[i] = 1.0/subvec[i];
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat)
{
  int iDof = 0;
  for(int i = 0; i < scomm->numT(SComm::all); ++i) {
    FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNumber, scomm->neighbT(SComm::all,i));
    for(int  j = 0; j < scomm->lenT(SComm::all,i); ++j) {
      int bdof = scomm->boundDofT(SComm::all,i,j);
      switch(boundDofFlag[iDof]) {
        case 0: {
          if(deltaF) sInfo.data[j] = deltaF[bdof];
          else sInfo.data[j] = 0.0;
        } break;
        case 1: {  // wet interface
          int windex = -1 - bdof;
          sInfo.data[j] = deltaFwi[windex];
        } break;
        case 2: {  // dual mpc or contact
          int locMpcNb = -1 - bdof;
          sInfo.data[j] = (masterFlag[iDof]) ? deltaFmpc[locMpcNb] : -deltaFmpc[locMpcNb];
        } break;
      }
      iDof++;
    }
  }
}

template<class Scalar>
double
GenSubDomain<Scalar>::collectAndDotDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat)
{
  // if there are more than 2 subdomains sharing a mpc define the norm
  // as (f1 - f2 - f3)^2 = f1^2 + f2^2 + f3^2 - 2f1f2 - 2f1f3 + 2f2f3 --> currently implemented
  Scalar dot = 0;
  int i, iSub, jDof;

  if(deltaF)
    {
      for(i=0; i<localLen(); ++i )
	{
	  double dPrScal = 1.0/this->densProjCoeff(i);
	  dot += dPrScal*dPrScal*deltaF[i]*ScalarTypes::conj(deltaF[i]);
	}
    }

  for(i = 0; i < numMPC; ++i)
    dot += deltaFmpc[i] * ScalarTypes::conj(deltaFmpc[i]);

  for(i=0; i<numWIdof; ++i) //HB ... to be checked ...
    dot += deltaFwi[i] * ScalarTypes::conj(deltaFwi[i])/wweight[i];

  int nbdofs = 0;
  for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
    FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all,iSub), subNumber);
    for(jDof = 0; jDof < scomm->lenT(SComm::all,iSub); ++jDof) {
      int bdof = scomm->boundDofT(SComm::all,iSub,jDof);
      switch(boundDofFlag[nbdofs]) {
        case 0:
          if(deltaF)
	    {
	      double dPrScal = 1.0/this->densProjCoeff(bdof);
	      dot += dPrScal*dPrScal*deltaF[bdof] * ScalarTypes::conj(rInfo.data[jDof]);
	    }
          break;
        // do nothing for case 1 (wet interface)
        case 2: { // dual mpc
          if(subNumber != scomm->neighbT(SComm::all,iSub)) {
            int locMpcNb =  -1 - bdof;
            dot += deltaFmpc[locMpcNb] * ScalarTypes::conj(rInfo.data[jDof]);
          }
        } break;
      }
      nbdofs++;
    }
  }
  return ScalarTypes::Real(dot);
}

template<class Scalar>
void
GenSubDomain<Scalar>::interfaceJump(Scalar *interfvec, FSCommPattern<Scalar> *vPat)
{
  int i, iSub, iDof;
  int offset = 0;

  Scalar *mpcJump = (numMPC) ? new Scalar[numMPC] : 0;
  bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(i = 0; i < numMPC; ++i) mpcFlag[i] = true;

  Scalar *wiJump = (Scalar *) dbg_alloca(sizeof(Scalar)*numWIdof);
  bool *wiFlag =  (bool *) dbg_alloca(sizeof(bool)*numWIdof);
  for(i = 0; i < numWIdof; ++i) wiFlag[i] = true;

  for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
    FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all,iSub), subNumber);
    for(iDof = 0; iDof < scomm->lenT(SComm::all,iSub); ++iDof) {
      int bdof = scomm->boundDofT(SComm::all,iSub,iDof);
      switch(boundDofFlag[offset+iDof]) {
        case 0:
          interfvec[offset+iDof] -= rInfo.data[iDof];
          break;
        case 1: { // wet interface
          int windex = -1 - bdof;
          if(wiFlag[windex]) {
            wiJump[windex] = interfvec[offset+iDof];
            wiFlag[windex] = false;
          }
          if(subNumber != scomm->neighbT(SComm::all,iSub)) {
            wiJump[windex] += rInfo.data[iDof];
          }
        } break;
        case 2: {  // dual mpc
          int locMpcNb = -1 - bdof;
          if(mpcFlag[locMpcNb]) {  // do this 1st time only
            mpcJump[locMpcNb] = interfvec[offset+iDof];
            mpcFlag[locMpcNb] = false;
          }
          if(subNumber != scomm->neighbT(SComm::all,iSub))
            mpcJump[locMpcNb] += rInfo.data[iDof];
        } break;
      }
    }
    offset += scomm->lenT(SComm::all,iSub);
  }

  // add mpcJump to interfvec
  for(i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    interfvec[scomm->mapT(SComm::mpc,i)] = mpcJump[scomm->mpcNb(i)];
  }
  // add wiJump to interfvec
  for(i = 0; i < scomm->lenT(SComm::wet); ++i) {
    interfvec[scomm->mapT(SComm::wet,i)] = wiJump[scomm->wetDofNb(i)];
  }
  if(mpcJump) delete [] mpcJump;
}

template<class Scalar>
void GenSubDomain<Scalar>::extractControlData(Scalar *disp, Scalar *vel,
                                              Scalar *acc, Scalar *ctrdsp,
                                              Scalar *ctrvel, Scalar *ctracc)
{
  for(int i = 0; i < claw->numSensor; ++i) {
    int dof = c_dsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
    int gi = locToGlSensorMap[i]; // local index
    if(dof >= 0) { // free
      ctrdsp[gi] = disp[dof];
      ctrvel[gi] = vel[dof];
      ctracc[gi] = acc[dof];
    }
    else { // either constrained or non-existant
      int dof2 = dsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
      if(dof2 >= 0) { // constrained
        ctrdsp[gi] = bcx[dof2];
        ctrvel[gi] = vcx[dof2];
        ctracc[gi] = acx[dof2];
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::constructKcc()
{
  if(Kcc) delete Kcc;
  Kcc = new GenAssembledFullM<Scalar>(numCRNdof, cornerMap);
  memK += numCRNdof*numCRNdof;
}

template<class Scalar>
void
GenSubDomain<Scalar>::constructKrc()
{
 Src = new GenSparseSet<Scalar>();
 if(numCRNdof) {
   Krc = new GenCuCSparse<Scalar>(nodeToNode, dsa, cornerMap, cc_dsa->getUnconstrNum());
   Src->addSparseMatrix(Krc);
 }
 setMpcSparseMatrix();
}

template<class Scalar>
void
GenSubDomain<Scalar>::constructKrw()
{
 if(numWIdof)
   Krw = new GenCuCSparse<Scalar>(nodeToNode, dsa, wetInterfaceMap, cc_dsa->getUnconstrNum());
}

template<class Scalar>
void
GenSubDomain<Scalar>::constructKww()
{
  if(numWIdof){
    localw      = new Scalar[numWIdof];
    localw_copy = new Scalar[numWIdof];
    Kww = new GenDBSparseMatrix<Scalar>(nodeToNode, dsa, wetInterfaceMap); Kww->zeroAll();
    memK += (Kww) ? Kww->size() : 0;
    memK += (Krw) ? Krw->size() : 0;
  }

  if(!neighbKww && numWIdof)
    neighbKww = new GenFsiSparse<Scalar>(domain->getFSI(), domain->getNumFSI(), glToLocalWImap);

#ifdef HB_COUPLED_PRECOND
  if(isMixedSub & neighbKww!=0) {
    double time0, time1  =0;
    time0 = -getTime();
    // extract the "local" fsi from neighbKww (i.e. the fsi between subdomain nodes/dofs)
    fprintf(stderr," ... Create localFsiToNodes connectivity in sub %2d\n",subNumber);
    Connectivity* localFsiToNodes = neighbKww->makeLocalFsiToNodes(glToLocalWImap, glToLocalNode, &globalNMax); //need to use globalNMax ???
    double dt0 = (time0 + getTime())/1000.;
    // create a precNodeToNode for preconditioner
    // -> create an "extended" elemToNode by merging "std" elemToNode with FsiToNode
    time1 = -getTime();
    if(!elemToNode) elemToNode = new Connectivity(&packedEset);
    Connectivity* precElToNodes = elemToNode->merge(localFsiToNodes);
    Connectivity* precNodeToEls = precElToNodes->reverse();
    precNodeToNode = precNodeToEls->transcon(precElToNodes);
    double dt1 = (time1 + getTime())/1000.;
    double dtt = (time0 + getTime())/1000.;
    //fprintf(stderr," -> CPU time to make localFsiToNodes: %f s (%f %%)\n",dt0,100*dt0/dtt);
    //fprintf(stderr," -> CPUT time to make precNodeToNode: %f s (%f %%)\n",dt1,100*dt1/dtt);
    //fprintf(stderr," -> total CPU time                  : %f s\n",dtt);
    delete precElToNodes;
    delete precNodeToEls;
    delete elemToNode; elemToNode = 0;
    if(localFsiToNodes) { delete localFsiToNodes; }
  } else {
    precNodeToNode = nodeToNode;
  }
#endif
}

template<class Scalar>
void
GenSubDomain<Scalar>::scaleAndSplitKww()
{
 if (neighbKww!=0)
   neighbKww->scale(cscale_factor);

 if (neighbKww!=0)
   neighbKww->split(glToLocalWImap, wweight);

#ifdef HB_COUPLED_PRECOND
 if(solInfo().isCoupled & isMixedSub & neighbKww!=0 & KiiSparse!=0) {
   fprintf(stderr," ... Assemble localFsi into Kii in sub %2d\n",subNumber);
   if(solInfo().getFetiInfo().splitLocalFsi) {
     neighbKww->splitLocalFsi(glToLocalWImap, wweight);
     neighbKww->addLocalFsiToMatrix(KiiSparse, dsa, glToLocalNode);
   } else {
     fprintf(stderr," ... No local Fsi spliting in sub %2d\n",subNumber);
     neighbKww->addLocalFsiToMatrix(KiiSparse, dsa, glToLocalNode, kSumWI);
     //neighbKww->addLocalFsiToMatrix(KiiSparse, dsa, glToLocalNode); // for test
   }
 }
#endif
 prev_cscale_factor = cscale_factor;
}

template<class Scalar>
void
GenSubDomain<Scalar>::reScaleAndReSplitKww()
{
  double rescale_factor =  cscale_factor/prev_cscale_factor;

  if (neighbKww!=0)
    neighbKww->scale(rescale_factor);

  if (neighbKww!=0)
    if(solInfo().getFetiInfo().fsi_scaling == FetiInfo::kscaling)
      neighbKww->split(glToLocalWImap, wweight);

#ifdef HB_COUPLED_PRECOND
 if(solInfo().isCoupled & isMixedSub & neighbKww!=0) {
   fprintf(stderr," ... Assemble localFsi into Kii in sub %2d\n",subNumber);
   if(KiiSparse) neighbKww->addLocalFsiToMatrix(KiiSparse, dsa, glToLocalNode);
 }
#endif
 prev_cscale_factor = cscale_factor;
}

template<class Scalar>
void
GenSubDomain<Scalar>::constructKcw()
{
  if(Kcw) { delete Kcw; Kcw = 0; }
  if(numWIdof && nCoarseDofs()){
    Kcw = new GenCuCSparse<Scalar>(nodeToNode, dsa, wetInterfaceMap, cornerMap);
    Kcw->zeroAll();
    memK += (Kcw) ? Kcw->size() : 0;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::addSingleFsi(LMPCons *localFsi)
{
#ifndef SALINAS
  int nEle = packedEset.last();
  packedEset.fsielemadd(nEle, localFsi);
  numele = packedEset.last();
#endif
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractInterfRBMs(int numRBM, Scalar *locRBMs, Scalar *locInterfRBMs)
{
  int iDof, iRBM;
  int locLen = localLen();

  // locInterfRBMs are ordered by rbm
  for(iRBM = 0; iRBM < numRBM; ++iRBM) {
    int off = iRBM*scomm->numT(SComm::std);
    int locOff = iRBM*locLen;
    for(iDof = 0; iDof < scomm->numT(SComm::std); ++iDof)
      locInterfRBMs[off+iDof] = locRBMs[locOff+scomm->boundDofT(SComm::std,iDof)];
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendInterfRBMs(int numRBM, Scalar *locInterfRBMs, FSCommPattern<Scalar> *rbmPat)
{
  // locInterfRBMs are ordered by-RBM, need to convert to by-neighbor ordering
 for(int iSub = 0; iSub < scomm->numT(SComm::std); ++iSub) {
   FSSubRecInfo<Scalar> sInfo = rbmPat->getSendBuffer(subNumber, scomm->neighbT(SComm::std,iSub));
   int off = 0;
   for(int iRBM = 0; iRBM < nGrbm; ++iRBM) {
      int locOff = iRBM*totalInterfSize + scomm->offsetT(SComm::std,iSub);
      for(int iDof = 0; iDof < scomm->lenT(SComm::std,iSub); ++iDof) {
         sInfo.data[off + iDof] = interfaceRBMs[locOff+iDof];
      }
      off += scomm->lenT(SComm::std,iSub);
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::recvInterfRBMs(int iNeighb, int numNeighbRBM, Scalar *neighbInterfRBMs,
                                     FSCommPattern<Scalar> *rbmPat)
{
  // this is for just one neighbor
  FSSubRecInfo<Scalar> rInfo = rbmPat->recData(scomm->neighbT(SComm::std,iNeighb), subNumber);
  int len = scomm->lenT(SComm::std,iNeighb)*numNeighbRBM;
  for(int j=0; j<len; ++j) neighbInterfRBMs[j] = rInfo.data[j];
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendInterfaceGrbm(FSCommPattern<Scalar> *rbmPat)
{
 // Sub-Domain based augmented preconditioner for DP
 Grc = new GenCuCSparse<Scalar>(scomm->lenT(SComm::std), scomm->boundDofsT(SComm::std), nGrbm, interfaceRBMs);
 Src->addSparseMatrix(Grc);

 sendInterfRBMs(nGrbm, interfaceRBMs, rbmPat);
}

template<class Scalar>
void
GenSubDomain<Scalar>::receiveInterfaceGrbm(FSCommPattern<Scalar> *rbmPat)
{
 // Sub-Domain based augmented preconditioner for DP
 for(int iNeighb = 0; iNeighb < scomm->numT(SComm::std); ++iNeighb) {
   int leadingDimGs = getSComm()->lenT(SComm::std,iNeighb);
   Scalar *neighbGs = new Scalar[leadingDimGs*neighbNumGRBMs[iNeighb]];
   recvInterfRBMs(iNeighb, neighbNumGRBMs[iNeighb], neighbGs, rbmPat);
   GenCuCSparse<Scalar> *Grc = new GenCuCSparse<Scalar>(leadingDimGs, scomm->boundDofsT(SComm::std,iNeighb),
                                                        neighbNumGRBMs[iNeighb], neighbGs, leadingDimGs);
   Grc->negate();
   Src->addSparseMatrix(Grc);
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::expandRBM(Scalar *localR, VectorSet &globalR)
{
 int globalNumRBM = globalR.numVec();
 for(int iRBM=0; iRBM<globalNumRBM; ++iRBM)
   for(int inode=0; inode<numnodes; ++inode) {
     int xdof = c_dsa->locate(inode, DofSet::Xdisp);
     int ydof = c_dsa->locate(inode, DofSet::Ydisp);
     int zdof = c_dsa->locate(inode, DofSet::Zdisp);
     Scalar w1 = weight[xdof];
     Scalar w2 = weight[ydof];
     Scalar w3 = weight[zdof];
     globalR[iRBM][6*glNums[inode]+0] = localR[xdof+iRBM*localLen()]/w1;
     globalR[iRBM][6*glNums[inode]+1] = localR[ydof+iRBM*localLen()]/w2;
     globalR[iRBM][6*glNums[inode]+2] = localR[zdof+iRBM*localLen()]/w3;
   }
}

template<>
void
GenSubDomain<DComplex>::getSRMult(DComplex *lvec, DComplex *interfvec, int nRBM,
                                  double *locRBMs, DComplex *alpha);

template<>
void
GenSubDomain<double>::getSRMult(double *lvec, double *interfvec, int nRBM,
                                double *locRBMs, double *alpha);

template<class Scalar>
void
GenSubDomain<Scalar>::makeKbb(DofSetArray *dof_set_array)
{
  if(glBoundMap) delete [] glBoundMap;
  glBoundMap = makeBMaps(dof_set_array);
  if(glInternalMap) delete [] glInternalMap;
  glInternalMap = makeIMaps(dof_set_array);

  Kbb = new GenDBSparseMatrix<Scalar>(nodeToNode, dsa, glBoundMap);
  if(Kbb) memPrec += Kbb->size();
  // KHP: 3-26-98 Modified this code to always construct Kib

  if(internalLen > 0) {
#ifdef HB_COUPLED_PRECOND
    if(solInfo().isCoupled & isMixedSub & neighbKww!=0)
       Kib = new GenCuCSparse<Scalar>(precNodeToNode, dsa, glBoundMap, glInternalMap);
    else
#endif
    Kib = new GenCuCSparse<Scalar>(nodeToNode, dsa, glBoundMap, glInternalMap);
    memPrec += Kib->size();
    Kib->zeroAll();
  }
  else Kib = 0;

  if((solInfo().getFetiInfo().precno == FetiInfo::dirichlet) && (internalLen > 0)) {
#ifdef HB_COUPLED_PRECOND
    if(solInfo().isCoupled & isMixedSub & neighbKww!=0)
      KiiSolver = GenSolverFactory<Scalar>::getFactory()->createSolver(precNodeToNode, dsa, glInternalMap, *sinfo.solvercntl->fetiInfo.kii_cntl, KiiSparse);
    else
#endif
      KiiSolver = GenSolverFactory<Scalar>::getFactory()->createSolver(nodeToNode, dsa, glInternalMap, *sinfo.solvercntl->fetiInfo.kii_cntl, KiiSparse);

    KiiSolver->setPrintNullity(false);
    // Find approximate preconditioner size
    memPrec += KiiSolver->size();
  }
  else {
    KiiSparse = 0;
    KiiSolver = 0;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::rebuildKbb()
{
  if(Kbb) delete Kbb;
  if(KiiSparse) delete KiiSparse;
  if(Kib) delete Kib;
  if(numMPC) makeKbbMpc(); else makeKbb(getCCDSA());
/*
  GenSparseMatrix<Scalar> *Kas = 0;
  GenSolver<Scalar> *smat = 0;
  assemble(Kas,smat,true);
  factorKii();
*/
}


template<class Scalar>
void
GenSubDomain<Scalar>::makeKbbMpc()
{
  // make the mpc dofs be boundary dofs
  int *weight_mpc = new int[c_dsa->size()];
  for(int i = 0; i < c_dsa->size(); ++i) weight_mpc[i] = 0;
  for(int i = 0; i < cc_dsa->size(); ++i) weight_mpc[ccToC[i]] = weight[i];

  for(int i = 0; i < numMPC; i++) {
    SubLMPCons<Scalar> *m = mpc[i];
    //XXXX if(mpc[i]->active) continue;
    for(int k = 0; k < m->nterms; k++) {
      int cdof = (m->terms)[k].cdof;
      if(cdof >= 0) weight_mpc[cdof] = 2; // > 1 hence will be in boundary (see makeBmaps and makeImaps)
    }
  }

  // construct the mapping and the matrices
  int *weight_copy = weight;
  weight = weight_mpc;
  makeKbb(c_dsa);

  // recover
  delete [] weight_mpc;
  weight = weight_copy;
}

template<class Scalar>
void
GenSubDomain<Scalar>::multKbb(Scalar *u, Scalar *Pu, Scalar *deltaU, Scalar *deltaF, bool errorFlag)
{
 // KHP and DJR: 3-26-98
 // multKbb has been modified to compute subdomain primal residual using
 // deltaF (in addition to deltaU which is the displacement correction
 // and of course the lumped or dirichlet preconditioning)

 // If we are computing lumped preconditioner, we compute deltaFi using Kib
 // but deltaUi is equal to zero. For the dirichlet preconditioner, deltaUi
 // is computed but deltaFi is set equal to zero.

 // deltaFi = internal primal residual
 // deltaUi = internal displacement correction

 // boundMap    = from boundary number to subdomain number
 // internalMap = from internal number to subdomain number
 // invBoundMap  = from all subdomain dof to a unique boundary number
 // invInternalMap  =  from all subdomain dof to a unique internal number
 // allBoundDofs = indices of B, from lambda numbering to numbering of entire
 //                subdomain

 // dualToBoundary = from lambda numbering directly to boundary numbering

 // Kii = works only on the numbering of the internal dofs
 // Kbb = works only on the numbering of the boundary dofs
 // Kib = operates on the boundary numbering and returns with internal numbering

 Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);
 Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);

 int iDof;
 for(iDof = 0; iDof < boundLen; ++iDof) {
   v[iDof] = res[iDof] = 0.0;
 }
 // Karim added the following lines. I Don't think they are necessary XML
 if(deltaF && deltaU)
   for(iDof = 0; iDof < localLen(); ++iDof)
     deltaU[iDof]=deltaF[iDof]=0;

 for(iDof = 0; iDof < totalInterfSize; ++iDof) {
   v[dualToBoundary[iDof]] += u[iDof]*scaling[iDof];
   if(deltaU)
     deltaU[allBoundDofs[iDof]] = -v[dualToBoundary[iDof]];
 }

 Kbb->mult(v, res);

 Scalar *iDisp = 0; // only allocate if necessary and don't use alloca (too big)

 if((solInfo().getFetiInfo().precno == FetiInfo::dirichlet) || deltaF) {
   iDisp = new Scalar[internalLen];
   for(iDof = 0; iDof < internalLen; ++iDof) iDisp[iDof] = 0.0;
   if(Kib) Kib->transposeMultAdd(v, iDisp);
 }

 if(solInfo().getFetiInfo().precno == FetiInfo::dirichlet) {
   if(KiiSolver) KiiSolver->reSolve(iDisp);
   if(Kib) Kib->multSub(iDisp, res);
   if(deltaU)
     for(iDof = 0; iDof < internalLen; ++iDof)
       deltaU[internalMap[iDof]] = iDisp[iDof];
 }
 else {
   if(deltaF)
     for(iDof = 0; iDof<internalLen; ++iDof) {
       deltaF[internalMap[iDof]] = iDisp[iDof];
     }
 }
 if(iDisp) delete [] iDisp;

 if(deltaF) {
   for(iDof = 0; iDof < boundLen; ++iDof)
     deltaF[boundMap[iDof]] = res[iDof];
 }

 // Return preconditioned u
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
}

template<class Scalar>
void
GenSubDomain<Scalar>::multKbbCoupled(Scalar *u, Scalar *Pu, Scalar *deltaF, bool errorFlag)
{
 // modified version of multKbb for coupled_dph with primal wet interface dofs
 // included in Kii. Currently preconditioner is un-coupled, ie fluid-structure interaction
 // is ignored in the Kii solve.

 Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);
 Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);
 if(!deltaFwi) deltaFwi = new Scalar[numWIdof];

 int i, iDof;
 for(iDof = 0; iDof < boundLen; ++iDof) { v[iDof] = res[iDof] = 0.0; }
 for(iDof = 0; iDof < localLen(); ++iDof) deltaF[iDof] = 0;
 for(i=0; i<numWIdof; ++i) localw[i] = 0;

 for(iDof = 0; iDof < totalInterfSize; ++iDof) {
   if(allBoundDofs[iDof] >= 0)
     v[dualToBoundary[iDof]] += u[iDof]*scaling[iDof];
   else if(boundDofFlag[iDof] == 1) { // wet interface
     int windex = -1-allBoundDofs[iDof];
#ifdef HB_COUPLED_PRECOND
     localw[windex] = (solInfo().getFetiInfo().splitLocalFsi) ? -u[iDof]*scaling[iDof] : -u[iDof];
#else
     localw[windex] = -u[iDof]*scaling[iDof];
#endif
     if(solInfo().getFetiInfo().precno == FetiInfo::noPrec)
       localw[windex] = u[iDof];
     deltaFwi[windex] = -u[iDof];
   }
 }

 Kbb->mult(v, res);

 Scalar *iDisp = (Scalar *) dbg_alloca(sizeof(Scalar)*internalLen);
 for(iDof = 0; iDof < internalLen; ++iDof) iDisp[iDof] = 0.0;
 for(i=0; i<numWIdof; ++i) iDisp[wiInternalMap[i]] = localw[i]; // coupled_dph
 if(Kib) Kib->transposeMultAdd(v, iDisp);

 if(solInfo().getFetiInfo().precno == FetiInfo::dirichlet) {
   if(KiiSolver) KiiSolver->reSolve(iDisp);
   for(i=0; i<numWIdof; ++i) localw[i] = -iDisp[wiInternalMap[i]];
   if(Kib) Kib->multSub(iDisp, res);
 }
 else {
   for(iDof = 0; iDof<internalLen; ++iDof)
     if(internalMap[iDof] > 0) deltaF[internalMap[iDof]] = iDisp[iDof]; // not including wet interface in deltaF
 }

 for(iDof = 0; iDof < boundLen; ++iDof) { deltaF[boundMap[iDof]] = res[iDof]; }

 // Return preconditioned u
 for(iDof = 0; iDof < totalInterfSize; ++iDof) {
   if(allBoundDofs[iDof] >= 0)
     Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
   else if(boundDofFlag[iDof] == 1)
     Pu[iDof] = localw[-1-allBoundDofs[iDof]] * scaling[iDof];
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multDiagKbb(Scalar *u, Scalar *Pu)
{
 // KHP: 02-02-99
 // boundMap    = from boundary number to subdomain number
 // internalMap = from internal number to subdomain number
 // invBoundMap =  from all subdomain dof to a unique boundary number
 // allBoundDofs = indices of B, from lambda numbering to numbering of entire
 //                subdomain
 // invBoundMap[allBoundDofs[iDof]] = from lambda numbering directly to
 //                                   boundary numbering
 // Kbb = works only on the numbering of the boundary dofs

 Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);
 Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);

 // XML Karim noted that his does not work for contact...

 int iDof;
 for(iDof = 0; iDof < boundLen; ++iDof)
   v[iDof] = res[iDof] = 0.0;

 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   v[dualToBoundary[iDof]] += u[iDof]*scaling[iDof];

 // Perform diagonal multiplication
 Kbb->multDiag(v, res);

 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
}

template<class Scalar>
void
GenSubDomain<Scalar>::multMFi(GenSolver<Scalar> *s, Scalar *u, Scalar *Fiu, int nRHS)
{
 // multMFi is never called in DP. Otherwise it will crash in contact
 int iRHS, iDof;
 int ndof = localLen();

 Scalar **iDisp = new Scalar *[nRHS];
 Scalar *firstpointer = new Scalar [nRHS*ndof];

 for(iRHS=0; iRHS < nRHS; iRHS++) {
   iDisp[iRHS] = firstpointer + iRHS*ndof;
   for(iDof = 0; iDof < ndof; ++iDof)
     iDisp[iRHS][iDof] = 0.0;
 }

 // Add the interfvec contribution to localvec
 for(iRHS=0; iRHS < nRHS; iRHS++) {
   for(iDof = 0; iDof < totalInterfSize; ++iDof)
     iDisp[iRHS][allBoundDofs[iDof]] += u[iDof+iRHS*totalInterfSize];
 }

 // solve for localvec
 if(s) s->reSolve(nRHS, iDisp);

 // redistribute the solution to the interface
 for(iRHS=0; iRHS < nRHS; iRHS++) {
   for(iDof = 0; iDof < totalInterfSize; ++iDof)
     Fiu[iDof+iRHS*totalInterfSize] = iDisp[iRHS][allBoundDofs[iDof]];
 }

 delete [] firstpointer;
 delete [] iDisp;
}

template<class Scalar>
void
GenSubDomain<Scalar>::multFi(GenSolver<Scalar> *s, Scalar *u, Scalar *Fiu)
{
 // multFi is never called in DP. Otherwise it will crash in contact

 int iDof;
 int ndof = localLen();
 Scalar *iDisp = (Scalar *) dbg_alloca(sizeof(Scalar)*ndof);
 for(iDof = 0; iDof < ndof; ++iDof)
   iDisp[iDof] = 0.0;

 // Add the interface vector contribution to local vector
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   iDisp[allBoundDofs[iDof]] += u[iDof];

 // solve for local vector
 if(s) s->reSolve(iDisp);

 // redistribute solution to the interface
 for(iDof = 0; iDof < totalInterfSize; ++iDof)
   Fiu[iDof] = iDisp[allBoundDofs[iDof]];

}

template<class Scalar>
void
GenSubDomain<Scalar>::factorKii()
{
  if(KiiSolver) KiiSolver->factor();
}

template<class Scalar>
void
GenSubDomain<Scalar>::factorKrr()
{
  if(Krr) Krr->factor();
}

template<class Scalar>
void
GenSubDomain<Scalar>::getHalfInterf(Scalar *s, Scalar *t)
{
  int iTg = 0;
  int i;
  for(i=0; i<totalInterfSize; ++i)
    if(masterFlag[i]) t[iTg++] = s[i];
}

template<class Scalar>
void
GenSubDomain<Scalar>::getHalfInterf(Scalar *s, Scalar *t, Scalar *ss, Scalar *tt)
{
  int iTg = 0;
  int i;
  for(i=0; i<totalInterfSize; ++i)
    if(masterFlag[i]) {
      t[iTg] = s[i];
      tt[iTg++] = ss[i];
    }
}

template<class Scalar>
void
GenSubDomain<Scalar>::scatterHalfInterf(Scalar *s, Scalar *buffer)
{
  // note: need to send s[iTg] of mpc master to all neighbors
  int iTg = 0;
  int i;

  Scalar *mpcbuff = (Scalar *) dbg_alloca(numMPC*sizeof(Scalar));
  Scalar *wibuff = (Scalar *) dbg_alloca(numWIdof*sizeof(Scalar));

  // note: if numMPC changes (salinas) need to delete [] mpcMaster
  for(i=0; i<totalInterfSize; ++i) {
    if(masterFlag[i]) {
      switch(boundDofFlag[i]) {
        case 0:
          iTg++;
          break;
        case 1: { // wet interface
          int windex = -1 - allBoundDofs[i];
          wibuff[windex] = s[iTg++];
        } break;
        case 2: { // dual mpc
          int locMpcNb = -1 - allBoundDofs[i];
          mpcbuff[locMpcNb] = s[iTg++];
        } break;
      }
    }
  }

  iTg = 0;
  for(i=0; i<totalInterfSize; ++i) {
    if(masterFlag[i]) buffer[i] = s[iTg++];
    else {
      switch(boundDofFlag[i]) {
        case 1: { // wet interface
          int windex = -1 - allBoundDofs[i];
          buffer[i] = (wiMaster[windex]) ? wibuff[windex] : 0.0;
        } break;
        case 2: { // dual mpc
          int locMpcNb = -1 - allBoundDofs[i];
          buffer[i] = (mpcMaster[locMpcNb]) ? mpcbuff[locMpcNb] : 0.0;
        } break;
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::rebuildInterf(Scalar *v, FSCommPattern<Scalar> *vPat)
{
  int iSub, i;
  int iOff = 0;
  Scalar *mpcv = (Scalar *) dbg_alloca(sizeof(Scalar)*numMPC);
  for(i = 0; i < numMPC; ++i) mpcv[i] = 0.0;
  Scalar *wiv = (Scalar *) dbg_alloca(numWIdof*sizeof(Scalar));
  for(i=0; i<numWIdof; ++i) wiv[i] = 0.0;

  for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
    FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all,iSub), subNumber);
    for(i=0; i<scomm->lenT(SComm::all,iSub); ++i) {
      int bdof = scomm->boundDofT(SComm::all,iSub,i);
      switch(boundDofFlag[iOff+i]) {
        case 0:
          if(!masterFlag[iOff+i]) v[iOff+i] = -rInfo.data[i];
          break;
        case 1: { // wet interface
          int windex = -1 - bdof;
          if(!masterFlag[iOff+i]) {
            if(rInfo.data[i] != 0.0) wiv[windex] = rInfo.data[i];
          }
          else wiv[windex] = v[i+iOff];
        } break;
        case 2: { // dual mpc
          int locMpcNb = -1 - bdof;
          if(!masterFlag[iOff+i]) {
            if(rInfo.data[i] != 0.0) mpcv[locMpcNb] = rInfo.data[i];
          }
          else mpcv[locMpcNb] = v[i+iOff];
        } break;
      }
    }
    iOff += scomm->lenT(SComm::all,iSub);
  }

  for(i=0; i<scomm->lenT(SComm::mpc); ++i)
    v[scomm->mapT(SComm::mpc,i)] = mpcv[scomm->mpcNb(i)];

  for(i=0; i<scomm->lenT(SComm::wet); ++i)
    v[scomm->mapT(SComm::wet,i)] = wiv[scomm->wetDofNb(i)];
}


template<class Scalar>
void
GenSubDomain<Scalar>::renumberElements()
{
 for(int i=0; i < numele; ++i)
   packedEset[i]->renum(glToLocalNode);
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberElementsGlobal()
{
 for(int i=0; i < numele; ++i)
   packedEset[i]->renum(glNums);
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberSharedNodes()
{
  int *allC = scomm->sharedNodes->tgt();
  for(int i = 0; i < scomm->sharedNodes->numConnect(); ++i) {
    int gi = allC[i];
    //if(globalToLocal(gi) < 0) cerr << "error in GenSubDomain::renumberSharedNodes() \n";
    allC[i] = globalToLocal(gi);
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberDirichlet()
{
 for(int i = 0; i < numDirichlet; ++i)
   dbc[i].nnum = glToLocalNode[dbc[i].nnum];
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberBCsEtc()
{
 int i;

 for(i = 0; i < numDirichlet; ++i)
   dbc[i].nnum = glToLocalNode[dbc[i].nnum];

 for(i = 0; i < numNeuman; ++i)
   nbc[i].nnum = glToLocalNode[nbc[i].nnum];

 for(i = 0; i < numIDis; ++i)
   iDis[i].nnum = glToLocalNode[iDis[i].nnum];

 for(i = 0; i < numIDis6; ++i)
   iDis6[i].nnum = glToLocalNode[iDis6[i].nnum];

 for(i = 0; i < numIVel; ++i)
   iVel[i].nnum = glToLocalNode[iVel[i].nnum];

 for(i = 0; i < numComplexDirichlet; ++i)
   cdbc[i].nnum = glToLocalNode[cdbc[i].nnum];

 for(i = 0; i < numComplexNeuman; ++i)
   cnbc[i].nnum = glToLocalNode[cnbc[i].nnum];

 for(i=0; i < numScatter; ++i)
   scatter[i]->renum(glToLocalNode);

 for(i=0; i < numWet; ++i)
   wet[i]->renum(glToLocalNode);

 if(sinfo.ATDARBFlag==-2.0)
   for(i=0; i < numSommer; ++i)
     sommer[i]->renum(glToLocalNode);

 for(i=0; i < numNeum; ++i)
   neum[i]->renum(glToLocalNode);

 for(i=0;i<numSBoundNodes;i++) {
   int tmp = glToLocalNode[sBoundNodes[i]];
   sBoundNodes[i] = tmp;
 }

 DMassData *cmass = firstDiMass;
 while(cmass != 0) {
   cmass->node = glToLocalNode[cmass->node];
   cmass = cmass->next;
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberControlLaw()
{
  int i;
  if(claw) {
   for(i = 0; i < claw->numSensor; i++)
     claw->sensor[i].nnum = glToLocalNode[claw->sensor[i].nnum];

   for(i = 0; i < claw->numActuator; i++)
     claw->actuator[i].nnum = glToLocalNode[claw->actuator[i].nnum];

   for(i = 0; i < claw->numUserDisp; i++)
     claw->userDisp[i].nnum = glToLocalNode[claw->userDisp[i].nnum];

   for(i = 0; i < claw->numUserForce; i++)
     claw->userForce[i].nnum = glToLocalNode[claw->userForce[i].nnum];
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberMPCs()
{
  for(int iMPC = 0; iMPC < numMPC; ++iMPC) {
    for(int i = 0; i < mpc[iMPC]->nterms; ++i)
      mpc[iMPC]->terms[i].nnum = globalToLocal(mpc[iMPC]->terms[i].nnum);
  }
  for(int iMPC = 0; iMPC < numMPC_primal; ++iMPC) {
    for(int i = 0; i < mpc_primal[iMPC]->nterms; ++i)
      mpc_primal[iMPC]->terms[i].nnum = globalToLocal(mpc_primal[iMPC]->terms[i].nnum);
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeElementForce(int fileNumber, Scalar *u, int forceIndex, Scalar *elemForce)
{
  if(elemToNode==0) elemToNode = new Connectivity(&packedEset);
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs, 0.0);
  double *nodalTemperatures = getNodalTemperatures();
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes, 0.0);
  Vector elForce(maxNumNodes, 0.0);

  int iele,k;
  for(iele=0; iele<numele; ++iele) {

    packedEset[iele]->nodes(nodeNumbers);
    int NodesPerElement = packedEset[iele]->numNodes();

    for(k=0; k<allDOFs->num(iele); ++k) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
      if(cn >= 0)
        (*elDisp)[k] = ScalarTypes::Real(u[cn]);
      else
        (*elDisp)[k] = 0.0;
    }

    if(packedEset[iele]->getProperty()) {
      if(domain->solInfo().thermalLoadFlag || (domain->solInfo().thermoeFlag >= 0)) {
        for(int iNode=0; iNode<NodesPerElement; ++iNode) {
          if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
            elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
          else
            elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
        }
      }
    }

    // transform displacements from DOF_FRM to basic coordinates
    transformVectorInv(*elDisp, iele);

    if(geoSource->getOutputInfo()[fileNumber].oframe == OutputInfo::Local) {
      Vector fx(NodesPerElement,0.0), fy(NodesPerElement,0.0), fz(NodesPerElement,0.0);
      if(forceIndex == INX || forceIndex == INY || forceIndex == INZ) {
        packedEset[iele]->getIntrnForce(fx,nodes,elDisp->data(),INX,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fy,nodes,elDisp->data(),INY,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fz,nodes,elDisp->data(),INZ,elemNodeTemps.data());
      }
      else {
        packedEset[iele]->getIntrnForce(fx,nodes,elDisp->data(),AXM,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fy,nodes,elDisp->data(),AYM,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fz,nodes,elDisp->data(),AZM,elemNodeTemps.data());
      }
      Vector f(3*NodesPerElement);
      for(int iNode=0; iNode<NodesPerElement; iNode++) {
        double data[3] = { fx[iNode], fy[iNode], fz[iNode] };
        transformVector(data, nodeNumbers[iNode], false);
        switch(forceIndex) {
          case INX: case AXM: elstress[iNode] = data[0]; break;
          case INY: case AYM: elstress[iNode] = data[1]; break;
          case INZ: case AZM: elstress[iNode] = data[2]; break;
        }
      }
    }
    else {
      packedEset[iele]->getIntrnForce(elForce, nodes, elDisp->data(), forceIndex, elemNodeTemps.data());
    }

    for(k = 0; k < packedEset[iele]->numNodes(); ++k)
      elemForce[elemToNode->offset(iele) + k] = elForce[k];
  }
  delete [] nodeNumbers;
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeStressStrain(int fileNumber,
                                          Scalar *u, int stressIndex,
                                          Scalar *glStress, Scalar *glWeight)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  if(elemToNode==0) elemToNode = new Connectivity(&packedEset);

  // avgnum = 2 --> do not include stress/strain of bar/beam element in averaging
  int avgnum = oinfo[fileNumber].averageFlg;

  // ylayer and zlayer are needed when calculating the axial stress/strain in a beam element
  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;
  OutputInfo::FrameType oframe = oinfo[fileNumber].oframe;

  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();
  double *nodalTemperatures = getNodalTemperatures();

  int k;
  int surface = oinfo[fileNumber].surface;

  int iele;
  GenVector<Scalar> *elstress = new GenVector<Scalar>(maxNumNodes);
  GenVector<double> *elweight = new GenVector<double>(maxNumNodes);
  GenVector<Scalar> *elDisp = new GenVector<Scalar>(maxNumDOFs);
  GenFullM<Scalar> *p_elstress = 0;
  if(oframe == OutputInfo::Local) p_elstress = new GenFullM<Scalar>(maxNumNodes,9);

  for(iele=0; iele<numele; ++iele) {

    // Don't do anything if element is a phantom
    if (packedEset[iele]->isPhantomElement()) continue;

    // Don't include beams or bars in the averaging if nodalpartial (avgnum = 2) is requested
    if ((avgnum == 2 && packedEset[iele]->getElementType() == 6) ||
        (avgnum == 2 && packedEset[iele]->getElementType() == 7) ||
        (avgnum == 2 && packedEset[iele]->getElementType() == 1)) continue;

    elstress->zero();
    elweight->zero();
    elDisp->zero();

    int NodesPerElement = packedEset[iele]->numNodes();
    packedEset[iele]->nodes(nodeNumbers);

    if(!isFluid(iele)) {

      for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] = u[cn];
        else
          (*elDisp)[k] = Bcx((*allDOFs)[iele][k]);
      }

      int iNode;
      if(packedEset[iele]->getProperty()) {
        if(domain->solInfo().thermalLoadFlag || (domain->solInfo().thermoeFlag >= 0)) {
          for(iNode=0; iNode<NodesPerElement; ++iNode) {
            if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
              elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
            else
              elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
          }
        }

        // transform displacements from DOF_FRM to basic coordinates
        transformVectorInv(*elDisp, iele);

        // transform non-invariant stresses/strains from basic frame to DOF_FRM
        if(oframe == OutputInfo::Local && ((stressIndex >=0 && stressIndex <=5) || (stressIndex >= 7 && stressIndex <= 12))) {

          // FIRST, CALCULATE STRESS/STRAIN TENSOR FOR EACH NODE OF THE ELEMENT
          p_elstress->zero();
          int strInd = (stressIndex >=0 && stressIndex <=5) ? 0 : 1;
          packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                         *elDisp, strInd, surface,
                                         elemNodeTemps.data());

          // second, transform stress/strain tensor to nodal frame coordinates
          transformStressStrain(*p_elstress, iele);

          // third, extract the requested stress/strain value from the stress/strain tensor
          for (iNode = 0; iNode < NodesPerElement; ++iNode) {
            if(strInd == 0)
              (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex];
            else
              (*elstress)[iNode] = (*p_elstress)[iNode][stressIndex-7];
          }

        }
        else {

          packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                        *elDisp, stressIndex, surface,
                                        elemNodeTemps.data(),ylayer,zlayer,avgnum);
        }
      }
    }

    if(glWeight)  {
      for (k = 0; k<NodesPerElement; ++k) {
#ifdef DISTRIBUTED
        glStress[nodeNumbers[k]] += (*elstress)[k];
        glWeight[nodeNumbers[k]] += (*elweight)[k];
#else
        glStress[(*elemToNode)[iele][k]] += (*elstress)[k];
        glWeight[(*elemToNode)[iele][k]] += (*elweight)[k];
#endif
      }
    }
    else  // non-averaged stresses
      for(k = 0; k < packedEset[iele]->numNodes(); ++k)
        glStress[ elemToNode->offset(iele) + k ] = (*elstress)[k];

  }

  delete [] nodeNumbers;
  delete elstress;
  delete elweight;
  delete elDisp;
  if(p_elstress) delete p_elstress;
}

// -----------------------------
//
// Nonlinear Subdomain functions
//
// -----------------------------

template<class Scalar>
void
GenSubDomain<Scalar>::computeStressStrain(GeomState *gs, Corotator **allCorot,
                                          int fileNumber, int stressIndex, Scalar *glStress, Scalar *glWeight,
                                          GeomState *refState)
{

  OutputInfo *oinfo = geoSource->getOutputInfo();

  if(elemToNode==0) elemToNode = new Connectivity(&packedEset);

  int k;
  int surface = oinfo[fileNumber].surface;
  int iele;
  int *nodeNumbers = new int[maxNumNodes];
  GenVector<Scalar> *elstress = new GenVector<Scalar>(maxNumNodes);
  GenVector<double> *elweight = new GenVector<double>(maxNumNodes);
  GenVector<Scalar> *elDisp = new GenVector<Scalar>(maxNumDOFs);

  // avgnum = 2 --> do not include stress/strain of bar/beam element in averaging
  int avgnum = oinfo[fileNumber].averageFlg;

  // ylayer and zlayer are needed when calculating the axial stress/strain in a beam element
  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;

  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();
  double *nodalTemperatures = getNodalTemperatures();

  int flag;

  for(iele=0; iele<numele; ++iele) {

    // Don't do anything if element is a phantom
    if (packedEset[iele]->isPhantomElement()) continue;

    // Don't include beams or bars in the averaging if nodalpartial (avgnum = 2) is requested
    if ((avgnum == 2 && packedEset[iele]->getElementType() == 6) ||
        (avgnum == 2 && packedEset[iele]->getElementType() == 7) ||
        (avgnum == 2 && packedEset[iele]->getElementType() == 1)) continue;

    elstress->zero();
    elweight->zero();
    elDisp->zero();

    allCorot[iele]->extractDeformations(*gs, nodes, elDisp->data(), flag);

    int NodesPerElement = elemToNode->num(iele);
    packedEset[iele]->nodes(nodeNumbers);

    int iNode;
    if(packedEset[iele]->getProperty()) {
      if(domain->solInfo().thermalLoadFlag || (domain->solInfo().thermoeFlag >= 0)) {
        for(iNode=0; iNode<NodesPerElement; ++iNode) {
          if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
            elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
          else
            elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
        }
      }
    }

    if(flag == 1) {
      // USE LINEAR STRESS ROUTINE
      packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                    *elDisp, stressIndex, surface,
                                     elemNodeTemps.data(), ylayer, zlayer, avgnum);
    }
    else if(flag == 2) {
      // USE NON-LINEAR STRESS ROUTINE
      allCorot[iele]->getNLVonMises(*elstress, *elweight, *gs, refState,
                                    nodes, stressIndex, surface,
                                    elemNodeTemps.data(), ylayer, zlayer, avgnum);
    }
    else {
      // NO STRESS RECOVERY
    }


    if(glWeight)  {
      for(k = 0; k < NodesPerElement; ++k) {
#ifdef DISTRIBUTED
        glStress[nodeNumbers[k]] += (*elstress)[k];
        glWeight[nodeNumbers[k]] += (*elweight)[k];
#else
        glStress[(*elemToNode)[iele][k]] += (*elstress)[k];
        glWeight[(*elemToNode)[iele][k]] += (*elweight)[k];
#endif
      }
    }
    else  // non-averaged stresses
      for(k = 0; k < packedEset[iele]->numNodes(); ++k)
        glStress[ elemToNode->offset(iele) + k ] = (*elstress)[k];

  }
  delete elstress;
  delete elweight;
  delete elDisp;
  delete [] nodeNumbers;
}

template<class Scalar>
void
GenSubDomain<Scalar>::reBuildKbb(FullSquareMatrix *kel)
{
 // Zero sparse matrices
 if(Kcc) Kcc->zero();
 if(Krc) Krc->zeroAll();
 if(Kbb) Kbb->zeroAll();
 if(Kib) Kib->zeroAll();
 if(KiiSparse) KiiSparse->zeroAll();

 // Assemble new subdomain sparse matrices
 int iele;
 for(iele=0; iele < numele; ++iele) {
   if(Kcc) Kcc->add(kel[iele],(*allDOFs)[iele]);
   if(Krc) Krc->add(kel[iele],(*allDOFs)[iele]);
   if(Kbb) Kbb->add(kel[iele],(*allDOFs)[iele]);
   if(Kib) Kib->add(kel[iele],(*allDOFs)[iele]);
   if(KiiSparse) KiiSparse->add(kel[iele],(*allDOFs)[iele]);
 }

 // Factor Kii if necessary (used for dirichlet preconditioner)
 if(KiiSolver) KiiSolver->factor();
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDisp(Scalar (*xyz)[11], GeomState* u)//DOfSet::max_known_nonL_dof
{
 // Loop over all nodes, Subtract the initial configuration from the
 // Current configuration to get the displacement
 int i, nodeI;
 for(i = 0; i < numnodes; ++i) {
   if(!nodes[i]) continue;
   nodeI = (domain->outFlag) ? domain->nodeTable[glNums[i]]-1 : glNums[i];

   if(dynamic_cast<TemperatureState*>(u)) { xyz[nodeI][0] = (*u)[i].x; continue; }
   xyz[nodeI][0] = ((*u)[i].x - nodes[i]->x);
   xyz[nodeI][1] = ((*u)[i].y - nodes[i]->y);
   xyz[nodeI][2] = ((*u)[i].z - nodes[i]->z);
   double rot[3];
   mat_to_vec((*u)[i].R,rot);
   xyz[nodeI][3] = rot[0];
   xyz[nodeI][4] = rot[1];
   xyz[nodeI][5] = rot[2];
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDistributedNLDisp(Scalar (*xyz)[11], GeomState* u)//DofSet::max_known_nonL_dof
{
 // Loop over all nodes, Subtract the initial configuration from the
 // Current configuration to get the displacement
 int i;
 for(i = 0; i < numnodes; ++i) {
   if(!nodes[i]) continue;
   if(dynamic_cast<TemperatureState*>(u)) { xyz[i][0] = (*u)[i].x; continue; }
   xyz[i][0] = ((*u)[i].x - nodes[i]->x);
   xyz[i][1] = ((*u)[i].y - nodes[i]->y);
   xyz[i][2] = ((*u)[i].z - nodes[i]->z);
   double rot[3];
   mat_to_vec((*u)[i].R,rot);
   xyz[i][3] = rot[0];
   xyz[i][4] = rot[1];
   xyz[i][5] = rot[2];
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeForces(Scalar (*mergedF)[6], Scalar *subF)
{
 for(int inode = 0; inode < numnodes; ++inode) {
   int nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]]-1 : glNums[inode];
   for(int jdof = 0; jdof < 6; ++jdof) {
     int cdof  = c_dsa->locate(inode, 1 << jdof);
     if(cdof >= 0)
       mergedF[nodeI][jdof] += subF[cdof]; // free: assemble into global array
     else
       mergedF[nodeI][jdof] = 0.0;         // constrained or doesn't exist
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeReactions(Scalar (*mergedF)[11], Scalar *subF)
{
 // TODO: just loop over dbc, cdbc
 for(int inode = 0; inode < numnodes; ++inode) {
   int nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]]-1 : glNums[inode];
   for(int jdof = 0; jdof < 11; ++jdof) {
     int dof = dsa->locate(inode, 1 << jdof);
     if(dof >= 0)  {
       int cdof = c_dsa->invRCN(dof);
       if(cdof >= 0) mergedF[nodeI][jdof] += subF[cdof]; // constrained
     }
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDistributedForces(Scalar (*mergedF)[6], Scalar *subF)
{
 for(int inode = 0; inode < numnodes; ++inode) {
   for(int jdof = 0; jdof < 6; ++jdof) {
     int cdof  = c_dsa->locate(inode, 1 << jdof);
     if(cdof >= 0)
       mergedF[inode][jdof] = subF[cdof];  // free
     else
       mergedF[inode][jdof] = 0.0;         // constrained or doesn't exist
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDistributedReactions(Scalar (*mergedF)[11], Scalar *subF)
{
 // TODO: just loop over dbc, cdbc
 for(int inode = 0; inode < numnodes; ++inode) {
   for(int jdof = 0; jdof < 11; ++jdof) {
     int dof = dsa->locate(inode, 1 << jdof);
     if(dof >= 0)  {
       int cdof = c_dsa->invRCN(dof);
       if(cdof >= 0) mergedF[inode][jdof] += subF[cdof]; // constrained
     }
   }
 }
}


template<class Scalar>
void
GenSubDomain<Scalar>::updatePrescribedDisp(GeomState *geomState, Scalar deltaLambda)
{
  if(numDirichlet > 0)
    geomState->updatePrescribedDisplacement(dbc, numDirichlet, deltaLambda);
}

template<class Scalar>
void
GenSubDomain<Scalar>::updatePrescribedDisp(GeomState *geomState)
{
  if(domain->solInfo().initialTime == 0.0) {
    // note 1: "if both IDISP and IDISP6 are present in the input file, FEM selects IDISP6 to initialize the displacement field"
    if((domain->numInitDisp() > 0) && (domain->numInitDisp6() == 0)) {
      if(this->numInitDisp() > 0) {
        geomState->updatePrescribedDisplacement(this->getInitDisp(), this->numInitDisp(), nodes);
      }
    }

    if(domain->numInitDisp6() > 0) {
      if(this->numInitDisp6() > 0)
        geomState->updatePrescribedDisplacement(this->getInitDisp6(), this->numInitDisp6(), nodes);
    }

    if(numDirichlet > 0)
      geomState->updatePrescribedDisplacement(dbc, numDirichlet, nodes);
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::setMpcSparseMatrix()
{
 // should always be the second term added to Src
 if(numMPC_primal > 0) {
   MPCsparse = new GenMpcSparse<Scalar>(numMPC_primal, mpc_primal, cc_dsa);
   Src->addSparseMatrix(MPCsparse);

   // to take into account wet interface / mpc interaction
   if(numWIdof) {
     int mpcOffset = numCRNdof;
     Kcw_mpc = new GenMpcSparse<Scalar>(numMPC_primal, mpc_primal, cc_dsa, dsa, wetInterfaceMap, mpcOffset);
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleMpcIntoKcc()
{
  // compute mpcOffset for my Subdomain
  int mpcOffset = numCRNdof;

  int i,iMPC;
  for(iMPC = 0; iMPC < numMPC_primal; ++iMPC) {
    for(i = 0; i < mpc_primal[iMPC]->nterms; ++i) {
      int   d = mpc_primal[iMPC]->terms[i].dof;
      int dof = mpc_primal[iMPC]->terms[i].ccdof;
      if((dof < 0) && (d >= 0) && !isWetInterfaceDof(d)) {
        int column = cornerMap[d];
        int row    = mpcOffset + iMPC;
        if(row > Kcc->dim())
          cout << " *** ERROR: Dimension Error Row = " << row    << " > " << Kcc->dim() << endl;
        if(column > Kcc->dim())
          cout << " *** ERROR: Dimension Error Col = " << column << " > " << Kcc->dim() << endl;
        // i.e. an mpc touches a corner node that also has DBCs
        if(column >= 0) {
          (*Kcc)[row][column] += mpc_primal[iMPC]->terms[i].coef;
          (*Kcc)[column][row] += mpc_primal[iMPC]->terms[i].coef;
        }
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multKcc()
{
 // resize Kcc if necessary to add space for "primal" mpcs and augmentation
 if(Src->numCol() != Kcc->dim()) {
   GenAssembledFullM<Scalar> *Kcc_copy = Kcc;
   Kcc = new GenAssembledFullM<Scalar>(Src->numCol(), cornerMap);
   for(int i=0; i<numCRNdof; ++i) for(int j=0; j<numCRNdof; ++j) (*Kcc)[i][j] = (*Kcc_copy)[i][j];
   delete Kcc_copy;
 }

 // add in MPC coefficient contributions for Kcc^(s).
 if(numMPC_primal > 0)
   assembleMpcIntoKcc();

 if(Krr == 0) return;

 // Kcc* -> Kcc - Krc^T Krr^-1 Krc

 // first, perform Krc I = iDisp
 // which extracts the correct rhs vectors for the forward/backwards

 int nRHS = Src->numCol();
 Scalar **iDisp = new Scalar *[nRHS];
 Scalar *firstpointer = new Scalar [nRHS*nRHS];

 int numEquations = Krr->neqs();
 if(BKrrKrc) {
   if(BKrrKrc[0]) { delete [] BKrrKrc[0]; BKrrKrc[0] = 0; }
   delete [] BKrrKrc; BKrrKrc = 0;
 }
 BKrrKrc = (nRHS > 0) ? new Scalar *[nRHS] : 0;
 Scalar *secondpointer = new Scalar [nRHS*totalInterfSize];
 Scalar *thirdpointer = new Scalar [nRHS*numEquations];
 Scalar **KrrKrc = (Scalar **) dbg_alloca(nRHS*sizeof(Scalar *));
 //if(nRHS*numEquations == 0)
 //  fprintf(stderr, "We have a zero size %d %d %d\n",numEquations,totalInterfSize,nRHS);

 int iRHS, iDof;
 for(iRHS=0; iRHS < nRHS; ++iRHS) {
   iDisp[iRHS]  = firstpointer  + iRHS*nRHS;
   BKrrKrc[iRHS] = secondpointer+iRHS*totalInterfSize;
   KrrKrc[iRHS] = thirdpointer + iRHS*numEquations;
   for(iDof=0; iDof<numEquations; ++iDof)
     KrrKrc[iRHS][iDof] = 0.0;
 }
 if(Src) Src->multIdentity(KrrKrc);

 // KrrKrc <- Krr^-1 Krc
 if(Krr) Krr->reSolve(nRHS, KrrKrc); // this can be expensive when nRHS is large eg for coupled 

 // -Krc^T KrrKrc
 for(iRHS=0; iRHS < nRHS; ++iRHS)
   for(iDof=0; iDof<nRHS; ++iDof)
     iDisp[iRHS][iDof] = 0.0;

 // Multiple RHS version of multSub: iDisp <- -Krc^T KrrKrc
 if(Src) Src->multSub(nRHS, KrrKrc, iDisp);

 if(Kcc) Kcc->add(iDisp);

 delete [] iDisp; delete [] firstpointer;

 int k;
 for(iRHS = 0; iRHS < nRHS; ++iRHS) {
   bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
   for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
   bool *wiFlag = (bool *) dbg_alloca(sizeof(bool)*numWIdof);
   for(int i = 0; i < numWIdof; ++i) wiFlag[i] = true;

   if(Krw) Krw->transposeMultNew(KrrKrc[iRHS], localw);

   for(iDof = 0; iDof < totalInterfSize; iDof++) {
     switch(boundDofFlag[iDof]) {
       case 0:
         BKrrKrc[iRHS][iDof] = KrrKrc[iRHS][allBoundDofs[iDof]];
         break;
       case 1: { // wet interface
         int windex = -1-allBoundDofs[iDof];
         if(wiFlag[windex]) {
           BKrrKrc[iRHS][iDof] = -localw[-1-allBoundDofs[iDof]];
           wiFlag[windex] = false;
         }
         else BKrrKrc[iRHS][iDof] = 0.0;
       } break;
       case 2: { // dual mpc
         int locMpcNb = -1-allBoundDofs[iDof];
         SubLMPCons<Scalar> *m = mpc[locMpcNb];
         BKrrKrc[iRHS][iDof] = 0.0;
         if(mpcFlag[locMpcNb]) {
           for(k = 0; k < m->nterms; k++) {
             int cc_dof = (m->terms)[k].ccdof;
             if(cc_dof >= 0) BKrrKrc[iRHS][iDof] += KrrKrc[iRHS][cc_dof] * (m->terms)[k].coef;
/* experimental code for mpc / wet interface interaction
             else {
               int c_dof = c_dsa->locate((m->terms)[k].nnum, (1 << (m->terms)[k].dofnum));
               if(c_dof > -1) {
                 int dof = dsa->locate((m->terms)[k].nnum, (1 << (m->terms)[k].dofnum));
                 if(wetInterfaceMap && (dof > -1) && (wetInterfaceMap[dof] > -1))
                   BKrrKrc[iRHS][iDof] -= localw[wetInterfaceMap[dof]];
               }
             }
*/
           }
           mpcFlag[locMpcNb] = false;
         }
       } break;
     }
   }
 }
 delete [] thirdpointer;
}


template<class Scalar>
void
GenSubDomain<Scalar>::reMultKcc()
{
 if(Krr == 0) return;

 // Kcc* -> Kcc - Krc^T Krr^-1 Krc

 // first, perform Krc I = iDisp
 // which extracts the correct rhs vectors for the forward/backwards

 int nRHS = Src->numCol();

 int numEquations = Krr->neqs();
 if(BKrrKrc) {
   if(BKrrKrc[0]) { delete [] BKrrKrc[0]; BKrrKrc[0] = 0; }
   delete [] BKrrKrc; BKrrKrc = 0;
 }
 BKrrKrc = (nRHS > 0) ? new Scalar *[nRHS] : 0;
 Scalar *secondpointer = new Scalar [nRHS*totalInterfSize];
 Scalar *thirdpointer = new Scalar [nRHS*numEquations];
 Scalar **KrrKrc = (Scalar **) dbg_alloca(nRHS*sizeof(Scalar *));
 if(nRHS*numEquations == 0)
   fprintf(stderr, "We have a zero size %d %d %d\n",numEquations,totalInterfSize,nRHS);

 int iRHS, iDof;
 for(iRHS=0; iRHS < nRHS; ++iRHS) {
   BKrrKrc[iRHS] = secondpointer+iRHS*totalInterfSize;
   KrrKrc[iRHS] = thirdpointer + iRHS*numEquations;
   for(iDof=0; iDof<numEquations; ++iDof)
     KrrKrc[iRHS][iDof] = 0.0;
 }
 if(Src) Src->multIdentity(KrrKrc);

 // KrrKrc <- Krr^-1 Krc
 if(Krr) Krr->reSolve(nRHS, KrrKrc); // this can be expensive when nRHS is large eg for coupled freq sweep

 int k;
 for(iRHS = 0; iRHS < nRHS; ++iRHS) {
   bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
   for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
   bool *wiFlag = (bool *) dbg_alloca(sizeof(bool)*numWIdof);
   for(int i = 0; i < numWIdof; ++i) wiFlag[i] = true;

   if(Krw) Krw->transposeMultNew(KrrKrc[iRHS], localw);

   for(iDof=0; iDof<totalInterfSize; iDof++) {
     switch(boundDofFlag[iDof]) {
       case 0:
         BKrrKrc[iRHS][iDof] = KrrKrc[iRHS][allBoundDofs[iDof]];
         break;
       case 1: { // wet interface
         int windex = -1-allBoundDofs[iDof];
         if(wiFlag[windex]) {
           BKrrKrc[iRHS][iDof] = -localw[-1-allBoundDofs[iDof]];
           wiFlag[windex] = false;
         }
         else BKrrKrc[iRHS][iDof] = 0.0;
       } break;
       case 2: { // dual mpc
         int locMpcNb = -1-allBoundDofs[iDof];
         SubLMPCons<Scalar> *m = mpc[locMpcNb];
         BKrrKrc[iRHS][iDof] = 0.0;
         if(mpcFlag[locMpcNb]) {
           for(k = 0; k < m->nterms; k++) {
             int cc_dof = (m->terms)[k].ccdof;
             if(cc_dof >= 0) BKrrKrc[iRHS][iDof] += KrrKrc[iRHS][cc_dof] * (m->terms)[k].coef;
/* experimental code for mpc / wet interface interaction
             else {
               int c_dof = c_dsa->locate((m->terms)[k].nnum, (1 << (m->terms)[k].dofnum));
               if(c_dof > -1) {
                 int dof = dsa->locate((m->terms)[k].nnum, (1 << (m->terms)[k].dofnum));
                 if(wetInterfaceMap && (dof > -1) && (wetInterfaceMap[dof] > -1))
                   BKrrKrc[iRHS][iDof] -= localw[wetInterfaceMap[dof]];
               }
             }
*/
           }
           mpcFlag[locMpcNb] = false;
         }
       } break;
     }
   }
 }
 delete [] thirdpointer;
}

template<class Scalar>
void
GenSubDomain<Scalar>::multKrc(Scalar *fr, Scalar *uc)
{
 int numCDofs = Src->numCol();
 Scalar *ucLocal = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);

 for(int i=0; i<numCDofs; ++i) {
   if(cornerEqNums[i] > -1) ucLocal[i] = uc[cornerEqNums[i]];
   else ucLocal[i] = 0.0;
 }

 // Perform fr = fr - Krc uc
 if(Src) Src->transposeMultSubtract(ucLocal, fr);
}

template<class Scalar>
void
GenSubDomain<Scalar>::multfc(Scalar *fr, /*Scalar *fc,*/ Scalar *lambda)
{
 int i, iDof, k;

 Scalar *force = new Scalar[localLen()];
 for(iDof = 0; iDof < localLen(); ++iDof) force[iDof] = -fr[iDof];

 //add the lambda contribution to fr, ie: -fr + Br^(s)^T lambda
 bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
 for(i = 0; i < numMPC; ++i) mpcFlag[i] = true;
 for(iDof = 0; iDof < totalInterfSize; ++iDof) {
   switch(boundDofFlag[iDof]) {
     case 0:
       force[allBoundDofs[iDof]] -= lambda[iDof];
       break;
     case 1:  // wet interface
       localw[-1-allBoundDofs[iDof]] = lambda[iDof];
       break;
     case 2: { // dual mpc or contact
       int locMpcNb = -1-allBoundDofs[iDof];
       if(mpcFlag[locMpcNb]) {
         SubLMPCons<Scalar> *m = mpc[locMpcNb];
         for(k = 0; k < m->nterms; k++) {
           int ccdof = (m->terms)[k].ccdof;
           if(ccdof >= 0) force[ccdof] -= lambda[iDof] * (m->terms)[k].coef;
         }
         mpcFlag[locMpcNb] = false;
       }
     } break;
   }
 }

 if(numWIdof) Krw->multAddNew(localw, force);  // coupled_dph: force += Krw * uw

 if(Krr) Krr->reSolve(force);

 // re-initialization required for mpc/contact
 if(fcstar) delete [] fcstar;  fcstar = new Scalar[Src->numCol()];
 for(i=0; i<Src->numCol(); ++i) fcstar[i] = 0.0;

 // fcstar = - (Krr^-1 Krc)^T fr
 //        = - Krc^T (Krr^-1 fr)
 //        = Src force
 if(Src) Src->multAdd(force, fcstar);
 delete [] force;

 // for coupled_dph add fcstar -= Kcw Bw uw
 if(numWIdof) {
   if(Kcw) Kcw->mult(localw, fcstar, -1.0, 1.0);
   if(Kcw_mpc) Kcw_mpc->multSubWI(localw, fcstar);
 }

 // add Bc^(s)^T lambda
 for(i = 0; i < numMPC; ++i) mpcFlag[i] = true;
 for(i = 0; i < scomm->lenT(SComm::mpc); ++i) {
   int locMpcNb = scomm->mpcNb(i);
   if(mpcFlag[locMpcNb]) {
     SubLMPCons<Scalar> *m = mpc[locMpcNb];
     for(k = 0; k < m->nterms; k++) {
       int dof = (m->terms)[k].dof;
       if((dof >= 0) && (cornerMap[dof] >= 0))
         fcstar[cornerMap[dof]] += lambda[scomm->mapT(SComm::mpc,i)] * (m->terms)[k].coef;
     }
     mpcFlag[locMpcNb] = false;
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multFcB(Scalar *p)
{
  int i, k;
  if(Src->numCol() == 0) return;
  if((totalInterfSize == 0) || (localLen() == 0)) {
    for(i=0; i<Src->numCol(); ++i) fcstar[i]=0.0;
    return;
  }

  // fcstar = - (Krr^-1 Krc)^T p
  //        = - Krc^T Krr^-1 p
  //        = - Acr p
  GenStackFullM<Scalar> Acr(Src->numCol(), totalInterfSize, BKrrKrc[0]);
  if(!fcstar) fcstar = new Scalar[Src->numCol()];
  Acr.mult(p, fcstar, -1.0, 0.0);

  // for coupled_dph add fcstar += Kcw Bw uw
  if(numWIdof) {
    for(i = 0; i < scomm->lenT(SComm::wet); ++i)
      localw[scomm->wetDofNb(i)] = p[scomm->mapT(SComm::wet,i)];
    if(Kcw) Kcw->mult(localw, fcstar, -1.0, 1.0);
    if(Kcw_mpc) Kcw_mpc->multSubWI(localw, fcstar);
  }

  // fcstar += Bc^(s)^T p
  bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(i = 0; i < numMPC; ++i) mpcFlag[i] = true;
  for(i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpcFlag[locMpcNb]) {
      SubLMPCons<Scalar> *m = mpc[locMpcNb];
      for(k = 0; k < m->nterms; k++) {
        int dof = (m->terms)[k].dof;
        if((dof >= 0) && (cornerMap[dof] >= 0))
          fcstar[cornerMap[dof]] += p[scomm->mapT(SComm::mpc,i)] * (m->terms)[k].coef;
      }
      mpcFlag[locMpcNb] = false;
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::getFc(Scalar *f, Scalar *Fc)
{
  int i, j;
  int dNum[DofSet::max_known_dof];
  int iOff = 0;
  for(i=0; i<numCRN; ++i) {
    int nd = c_dsa->number(cornerNodes[i], cornerDofs[i], dNum);
    for(j = 0; j < nd; ++j) Fc[iOff+j] = f[dNum[j]];
    iOff += nd;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::getFw(Scalar *f, Scalar *fw)
{
 if(numWIdof) {
  int i,j;
  int dofs[DofSet::max_known_dof];
  int cdofs[DofSet::max_known_dof];
  for(i = 0; i < numWInodes; ++i) {
    DofSet thisDofSet = wetInterfaceDofs[i];
    int nd = thisDofSet.count();
    dsa->number(wetInterfaceNodes[i], thisDofSet, dofs);
    c_dsa->number(wetInterfaceNodes[i], thisDofSet, cdofs);
    for(j = 0; j < nd; ++j)
      fw[wetInterfaceMap[dofs[j]]] = f[cdofs[j]];
  }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::getFr(Scalar *f, Scalar *fr)
{
  int i,iNode;
  int rDofs[DofSet::max_known_dof];
  int oDofs[DofSet::max_known_dof];
  for(iNode = 0; iNode < numnodes; ++iNode) {
    DofSet thisDofSet = (*cc_dsa)[iNode];
    int nd = thisDofSet.count();
    c_dsa->number(iNode, thisDofSet, oDofs);
    cc_dsa->number(iNode, thisDofSet, rDofs);
    for(i = 0; i < nd; ++i) {
      fr[rDofs[i]] = f[oDofs[i]];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeUr(Scalar *ur, Scalar *uc, Scalar *u, Scalar *lambda)
{
 int i, iNode;
 int rDofs[DofSet::max_known_dof];
 int oDofs[DofSet::max_known_dof];
 for(iNode = 0; iNode < numnodes; ++iNode) {
    DofSet thisDofSet = (*cc_dsa)[iNode];
    int nd = thisDofSet.count();
    c_dsa->number(iNode, thisDofSet, oDofs);
    cc_dsa->number(iNode, thisDofSet, rDofs);
    for(i = 0; i < nd; ++i) {
      u[oDofs[i]] = ur[rDofs[i]];
    }
 }
 int j;
 int iOff = 0;
 for(i=0; i<numCRN; ++i) {
    int nd = c_dsa->number(cornerNodes[i], cornerDofs[i], oDofs);
    for(j = 0; j < nd; ++j) {
       if(cornerEqNums[iOff+j] > -1)
         u[oDofs[j]] = uc[cornerEqNums[iOff+j]];
       else
         u[oDofs[j]] = 0.0;
    }
    iOff += nd;
 }

 // extract uw
 Scalar *uw = (Scalar *) dbg_alloca(numWIdof*sizeof(Scalar));
 for(i=0; i < scomm->lenT(SComm::wet); ++i)
   uw[scomm->wetDofNb(i)] = lambda[scomm->mapT(SComm::wet,i)];

 for(i=0; i<numWInodes; ++i) {
   DofSet thisDofSet = wetInterfaceDofs[i]; // (*c_dsa)[wetInterfaceNodes[i]];
   int nd = thisDofSet.count();
   dsa->number(wetInterfaceNodes[i], thisDofSet, rDofs);
   c_dsa->number(wetInterfaceNodes[i], thisDofSet, oDofs);
   for(j = 0; j < nd; ++j) {
     u[oDofs[j]] = uw[wetInterfaceMap[rDofs[j]]];
   }
 }

 // keep a local copy of the lagrange multipliers
 setLocalLambda(lambda);
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeQ()
{
 switch(solInfo().getFetiInfo().augment) {
   default:
     break;
   case FetiInfo::Gs:
   {
     nGrbm = myMin(rigidBodyModes->numRBM(), solInfo().getFetiInfo().nGs);
     int iLen = scomm->lenT(SComm::std);
     interfaceRBMs = new Scalar[nGrbm*iLen];
     rbms = new Scalar[nGrbm*iLen];
     int *boundToCDSA = (int *) dbg_alloca(sizeof(int)*iLen);
     int i;
     for(i = 0; i < iLen; ++i)
       if(scomm->boundDofT(SComm::std,i) >= 0)
         boundToCDSA[i] = ccToC[scomm->boundDofT(SComm::std,i)];
       else
         boundToCDSA[i] = -1;
     int offset = 0;
     if((rigidBodyModes->numRBM() == 6) && (solInfo().getFetiInfo().rbmType == FetiInfo::rotation))
       offset = 3;
     rigidBodyModes->getRBMs(interfaceRBMs, iLen, boundToCDSA, nGrbm, offset);
   }
   break;
   case FetiInfo::WeightedEdges:
   case FetiInfo::Edges:
   {
     switch(solInfo().getFetiInfo().rbmType )
     {
       default:
       case FetiInfo::translation:
       case FetiInfo::rotation:
       case FetiInfo::all:
       case FetiInfo::None:
         if(isMixedSub) {
           makeEdgeVectorsPlus(true); // build augmentation for fluid dofs
           makeEdgeVectorsPlus(false);  // build augmentation for structure dofs
         }
         else {
           makeEdgeVectorsPlus(isFluid(0));
         }
         break;
       case FetiInfo::averageTran:
       case FetiInfo::averageRot:
       case FetiInfo::averageAll:
         makeAverageEdgeVectors();
         break;
     }
   }
   break;
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::weightEdgeGs()
{
  // for WeightedEdge augmentation
  if(solInfo().getFetiInfo().scaling == FetiInfo::tscaling)
    Grc->doWeighting(weight);
  else if(solInfo().getFetiInfo().scaling == FetiInfo::kscaling)
    Grc->doWeighting(kweight);
}

// I've changed this routine to compute the following Q vectors:
//
// averageTran = [Qx + Qy]  i.e. [ 1 1 0 ]^T at a node
// averageRot  = [Qy + Qz]  i.e. [ 0 1 1 ]^T at a node
// averageAll  = computes both of these Q vectors
//
// These 2 vectors are better than [Qx Qy Qz] together, which seems too
// "weak" of a constraint to help the convergence very much. Maybe there
// are other "stronger" vectors that will help more.
//

template<class Scalar>
void
GenSubDomain<Scalar>::makeAverageEdgeVectors()
{
     int i;
     int numR=1;
     if(solInfo().getFetiInfo().rbmType == FetiInfo::averageAll)
       numR=2;
     int totalLengthGrc=0;

     Connectivity &sharedNodes = *(scomm->sharedNodes);
     edgeDofSize = new int[scomm->numNeighb];

     // 1. first count number of edge dofs
     int iSub, iNode, nE=0;
     for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
       edgeDofSize[iSub]=0;
       for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) { edgeDofSize[iSub]=numR;  }
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) { edgeDofSize[iSub]=numR;  }
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) { edgeDofSize[iSub]=numR;  }
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Xrot))  { edgeDofSize[iSub]=numR;  }
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Yrot))  { edgeDofSize[iSub]=numR;  }
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Zrot))  { edgeDofSize[iSub]=numR;  }
       }
       if(edgeDofSize[iSub] > 0) nE++;
     }

     int *xyzCount = new int[nE*numR];
     for(i=0; i<numR*nE; ++i)
       xyzCount[i] = 0;

     nE=0;
     int nDof = 0;
     if(numR > 1) nDof=1;

   for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
     if(edgeDofSize[iSub]==0) continue;
     int xCount=0, yCount=0, zCount=0, xrCount=0, yrCount=0, zrCount=0;
     for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) xCount++;
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) yCount++;
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) zCount++;
        if(solInfo().getFetiInfo().rbmType == FetiInfo::averageRot || (numR > 1)) {
          if(boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) xrCount++;
          if(boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) yrCount++;
          if(boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) zrCount++;
        }
     }
     int edgeLength = 0;
     if(solInfo().getFetiInfo().rbmType == FetiInfo::averageTran || (numR > 1)) {
       xyzCount[numR*nE+0] = xCount + yCount + zCount;
       edgeLength         += xCount + yCount + zCount;
     }
     if(solInfo().getFetiInfo().rbmType == FetiInfo::averageRot || (numR > 1)) {
       // KHP
       //xyzCount[numR*nE+nDof]  = xrCount + yrCount + zrCount;
       //edgeLength             += xrCount + yrCount + zrCount;
         xyzCount[numR*nE+nDof]  = xCount + yCount + zCount;
         edgeLength             += xCount + yCount + zCount;
     }
     totalLengthGrc += edgeLength;

     nE++;
  }

// --------------------------------------------------------
  int *xyzList = new int[totalLengthGrc];
  Scalar *xyzCoefs = new Scalar[totalLengthGrc];
  int xOffset=0;
  int xrOffset=0;
  nE=0;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    if(edgeDofSize[iSub]==0) continue;
    if(numR > 1) {
      xrOffset  =  xOffset + xyzCount[numR*nE+0];
    }
    Scalar sign = (scomm->subNums[iSub] < subNumber) ? 1.0 : -1.0;
    int used=0;
    int middleNode = sharedNodes.num(iSub)/2;
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
          int xDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Xdisp);
          xyzList[xOffset]    = xDof;
          xyzCoefs[xOffset++] = (middleNode == iNode) ? sign*1.0 : 0.0;
          used=1;
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
          int yDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Ydisp);
          xyzList[xOffset]    = yDof;
          xyzCoefs[xOffset++] = (middleNode == iNode) ? sign*0.0 : 0.0;
          used = 1;
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
          int zDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Zdisp);
          xyzList[xOffset]    = zDof;
          xyzCoefs[xOffset++] = (middleNode == iNode) ? sign*1.0 : 0.0;
          used=1;
        }

       if(solInfo().getFetiInfo().rbmType == FetiInfo::averageRot || (numR > 1)) {
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
          int xDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Xdisp);
          xyzList[xrOffset]    = xDof;
          xyzCoefs[xrOffset++] = (middleNode == iNode) ? sign*1.0 : 0.0;
          used=1;
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
          int yDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Ydisp);
          xyzList[xrOffset]    = yDof;
          xyzCoefs[xrOffset++] = (middleNode == iNode) ? sign*1.0 : 0.0;
          used = 1;
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
          int zDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Zdisp);
          xyzList[xrOffset]    = zDof;
          xyzCoefs[xrOffset++] = (middleNode == iNode) ? sign*0.0 : 0.0;
          used=1;
        }
/*
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) {
           xyzList[xrOffset] = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Xrot);
           xyzCoefs[xrOffset++] = sign; used=1;
         }
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) {
           xyzList[xrOffset] = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Yrot);
           xyzCoefs[xrOffset++] = sign; used=1;
         }
         if(boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) {
           xyzList[xrOffset] = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Zrot);
           xyzCoefs[xrOffset++] = sign; used=1;
         }
*/
       }

    }
   int off = 0;
   if(numR > 1)
     off += xyzCount[numR*nE+1];

   xOffset += off;

   if(used) nE++;
   }
   Grc = new GenCuCSparse<Scalar>(numR*nE, cc_dsa->size(),
                                  xyzCount, xyzList, xyzCoefs);
   // Src->setSparseMatrices(1, Grc);
   Src->addSparseMatrix(Grc);
   delete [] xyzCount;
}

template<>
void
GenSubDomain<DComplex>::precondGrbm();

template<>
void
GenSubDomain<double>::precondGrbm();

template<class Scalar>
void
GenSubDomain<Scalar>::orthoWithOrigR(Scalar *origV, Scalar *V, int numV, int length)
{
  int i,j;
  for(j=0; j<numV; ++j) {
    GenStackVector<Scalar> v1(origV+j*length,length);
    for(i=0; i<numV; ++i) {
      GenStackVector<Scalar> v2(V+i*length,length);
      Scalar s = v1*v2;
      v2 -= s*v1; // Vector -= operation
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::ortho(Scalar *vectors, int numVectors, int length)
{
  int i,j;
  for(j=0; j<numVectors; ++j) {
    GenStackVector<Scalar> v1(vectors+j*length,length);
    for(i=0; i<j; ++i) {
      GenStackVector<Scalar> v2(vectors+i*length,length);
      Scalar s = v1*v2;
      v1 -= s*v2; // Vector -= operation
    }
    v1 *= 1.0/v1.norm();
  }
}

template<class Scalar>
void GenSubDomain<Scalar>::initMpcScaling()
{
  // sets scaling = 1.0 for all mpc virtual dofs, for use with generalized preconditioner
  for(int i=0; i<scomm->lenT(SComm::mpc); ++i)
    scaling[scomm->mapT(SComm::mpc,i)] = 1.0;
}

template<class Scalar>
void GenSubDomain<Scalar>::setUserDefBC(double *usrDefDisp, double *usrDefVel, double *usrDefAcc, bool nlflag)
{
  // modify boundary condition values for output purposes
  int i;
  for(i = 0; i < claw->numUserDisp; ++i) {
    int dof = dsa->locate(claw->userDisp[i].nnum,1 << claw->userDisp[i].dofnum);
    if(dof >= 0) {
      bcx[dof] = usrDefDisp[locToGlUserDispMap[i]];
      if(bcx_scalar) bcx_scalar[dof] = bcx[dof];
      vcx[dof] = usrDefVel[locToGlUserDispMap[i]];
      acx[dof] = usrDefAcc[locToGlUserDispMap[i]];
    }
  }
  updateUsddInDbc(usrDefDisp, locToGlUserDispMap); // CHECK

  if(nlflag) return; // don't need to adjust rhs for nonlinear dynamics

  for(int i=0; i < numMPC; ++i) {
    mpc[i]->rhs = mpc[i]->original_rhs;
    for(int j=0; j < mpc[i]->nterms; ++j) {
      int mpc_node = mpc[i]->terms[j].nnum;
      int mpc_dof = mpc[i]->terms[j].dofnum;
      for(int k=0; k<numDirichlet; ++k) {
        int dbc_node = dbc[k].nnum;
        int dbc_dof = dbc[k].dofnum;
        if((dbc_node == mpc_node) && (dbc_dof == mpc_dof)) {
          mpc[i]->rhs -= mpc[i]->terms[j].coef * dbc[k].val;
        }
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeKccDofsExp(ConstrainedDSA *cornerEqs, int augOffset,
                                     Connectivity *subToEdge, int mpcOffset, GlobalToLocalMap& nodeMap)
{
  int numC = numCoarseDofs();
  if(cornerEqNums) delete [] cornerEqNums;
  cornerEqNums = new int[numC];

  // numbers the corner equations 
  int offset = 0;
  for(int i=0; i<numCRN; ++i) {
    offset += cornerEqs->number(nodeMap[glCornerNodes[i]], cornerDofs[i].list(), cornerEqNums+offset);          
  }
  //std::cerr << "here in GenSubDomain<Scalar>::makeKccDofsExp, cornerEqNums = ";
  //for(int i=0; i<numC; ++i) std::cerr << cornerEqNums[i] << " "; std::cerr << std::endl;
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeKccDofsExp2(int nsub, GenSubDomain<Scalar> **sd)
{
  nCDofs = 0;
  for(int i=0; i<nsub; ++i) nCDofs += sd[i]->getCDSA()->size();
  if(cornerEqNums) delete [] cornerEqNums;
  cornerEqNums = new int[nCDofs];

  // numbers the corner equations 
  int offset = 0;
  for(int i=0; i<numCRN; ++i) {
    int offset2 = 0;
    for(int j=0; j<nsub; ++j) {
      GlobalToLocalMap& nodeMap = sd[j]->getGlobalToLocalNodeMap();
      ConstrainedDSA *cornerEqs = sd[j]->getCDSA();
      if(nodeMap[glCornerNodes[i]] > -1) {
         int count = cornerEqs->number(nodeMap[glCornerNodes[i]], cornerDofs[i].list(), cornerEqNums+offset);
         for(int k=0; k<count; ++k) cornerEqNums[offset+k] += offset2;
         offset += count;
      }
      offset2 += cornerEqs->size();
    }
  }
  //std::cerr << "here in GenSubDomain<Scalar>::makeKccDofsExp, cornerEqNums = ";
  //for(int i=0; i<nCDofs; ++i) std::cerr << cornerEqNums[i] << " "; std::cerr << std::endl;
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeKccDofs(DofSetArray *cornerEqs, int augOffset,
                                  Connectivity *subToEdge, int mpcOffset)
{
  int numC = numCoarseDofs();
  if(cornerEqNums) delete [] cornerEqNums;
  cornerEqNums = new int[numC];

  // numbers the corner equations
  int offset = 0;
  for(int i=0; i<numCRN; ++i) {
    offset += cornerEqs->number(glCornerNodes[i], cornerDofs[i].list(), cornerEqNums+offset);
  }

  // number the mpc equations
  for(int i = 0; i < numMPC_primal; ++i) {
    int fDof = cornerEqs->firstdof(mpcOffset+localToGlobalMPC_primal[i]);
    cornerEqNums[offset++] = fDof;
  }

  // number the augmentation equations
  if(solInfo().getFetiInfo().augment == FetiInfo::Gs) {
    int fDof = cornerEqs->firstdof(augOffset+subNumber);
    for(int i=0; i<nGrbm; ++i)
      cornerEqNums[offset++] = fDof+i;
    for(int iNeighb=0; iNeighb<scomm->numNeighb; ++iNeighb) {
      int fDof = cornerEqs->firstdof(augOffset+scomm->subNums[iNeighb]);
      for(int i=0; i<neighbNumGRBMs[iNeighb]; ++i)
        cornerEqNums[offset++] = fDof+i;
    }
  }
  else if(solInfo().getFetiInfo().isEdgeAugmentationOn()) {
    int iEdgeN = 0;
    for(int iNeighb=0; iNeighb<scomm->numNeighb; ++iNeighb) {
      if(scomm->isEdgeNeighb[iNeighb]) {
        int fDof = cornerEqs->firstdof(augOffset+(*subToEdge)[subNumber][iEdgeN]);
        if(isMixedSub) {
          for(int i=0; i<edgeDofSize[iNeighb]-edgeDofSizeTmp[iNeighb]; ++i)  // fluid
            cornerEqNums[offset++] = fDof+i;
        }
        else {
          for(int i=0; i<edgeDofSize[iNeighb]; ++i)
            cornerEqNums[offset++] = fDof+i;
        }
        iEdgeN++;
      }
    }
    iEdgeN = 0;
    if(isMixedSub) {
      for(int iNeighb=0; iNeighb<scomm->numNeighb; ++iNeighb) {
        if(scomm->isEdgeNeighb[iNeighb]) {
          int fDof = cornerEqs->firstdof(augOffset+(*subToEdge)[subNumber][iEdgeN]) +
            edgeDofSize[iNeighb]-edgeDofSizeTmp[iNeighb];
          for(int i=0; i<edgeDofSizeTmp[iNeighb]; ++i)  // structure
            cornerEqNums[offset++] = fDof+i;
          iEdgeN++;
        }
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleKccStar(GenSparseMatrix<Scalar> *KccStar)
{
  KccStar->add(*Kcc, cornerEqNums);
}

template<class Scalar>
void
GenSubDomain<Scalar>::deleteKcc()
{
  delete Kcc; Kcc = 0;
}

template<class Scalar>
void
GenSubDomain<Scalar>::multKbbMpc(Scalar *u, Scalar *Pu, Scalar *deltaU, Scalar *deltaF, bool errorFlag)
{
  // KHP and DJR: 3-26-98
  // multKbb has been modified to compute subdomain primal residual using
  // deltaF (in addition to deltaU which is the displacement correction
  // and of course the lumped or dirichlet preconditioning)

  // If we are computing lumped preconditioner, we compute deltaFi using Kib
  // but deltaUi is equal to zero. For the dirichlet preconditioner, deltaUi
  // is computed but deltaFi is set equal to zero.

  // deltaFi = internal primal residual
  // deltaUi = internal displacement correction
  // boundMap = from boundary number to subdomain number
  // internalMap = from internal number to subdomain number
  // invBoundMap = from all subdomain dof to a unique boundary number
  // invInternalMap =  from all subdomain dof to a unique internal number
  // allBoundDofs = indices of B, from lambda numbering to numbering of entire subdomain
  // dualToBoundary = from lambda numbering directly to boundary numbering

  // Kii works only on the numbering of the internal dofs
  // Kbb works only on the numbering of the boundary dofs
  // Kib operates on the boundary numbering and returns with internal numbering

  Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);
  Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar)*boundLen);
  if(!deltaFwi) deltaFwi = new Scalar[numWIdof]; // coupled_dph

  int i, iDof;
  for(iDof = 0; iDof < boundLen; ++iDof) v[iDof] = res[iDof] = 0.0;
  for(iDof = 0; iDof < localLen(); ++iDof) {
    if(deltaU) deltaU[iDof] = 0.0;
    if(deltaF) deltaF[iDof] = 0.0;
  }
  for(i=0; i<numWIdof; ++i) localw[i] = 0;

  applyBtransposeAndScaling(u, v, deltaU, localw);

  //Scalar norm = 0; for(i=0; i<boundLen; ++i) norm += v[i]*v[i]; cerr << "1. norm = " << sqrt(norm) << endl;
  if((solInfo().getFetiInfo().precno == FetiInfo::lumped) ||
     (solInfo().getFetiInfo().precno == FetiInfo::dirichlet) || errorFlag) Kbb->mult(v, res);  // res = Kbb * v
  //norm = 0; for(i=0; i<boundLen; ++i) norm += res[i]*res[i]; cerr << "2. norm = " << sqrt(norm) << endl;

  if((solInfo().getFetiInfo().precno == FetiInfo::dirichlet) || errorFlag) {
    Scalar *iDisp = new Scalar[internalLen];
    for(iDof = 0; iDof < internalLen; ++iDof) iDisp[iDof] = 0.0;
    for(i=0; i<numWIdof; ++i) iDisp[wiInternalMap[i]] = localw[i]; // coupled_dph
    if(Kib) Kib->transposeMultAdd(v, iDisp); // iDisp += Kib^T * v

    if(solInfo().getFetiInfo().precno == FetiInfo::dirichlet) {
      //norm = 0; for(i=0; i<internalLen; ++i) norm += iDisp[i]*iDisp[i]; cerr << "3. norm = " << sqrt(norm) << endl;
      if(KiiSolver) KiiSolver->reSolve(iDisp);
      //norm = 0; for(i=0; i<internalLen; ++i) norm += iDisp[i]*iDisp[i]; cerr << "4. norm = " << sqrt(norm) << endl;
      for(i=0; i<numWIdof; ++i) localw[i] = iDisp[wiInternalMap[i]]; // coupled_dph
      if(Kib) Kib->multSub(iDisp, res); // res -= Kib*iDisp
    }
    else if(deltaF) { // improves estimate of error
      for(iDof = 0; iDof<internalLen; ++iDof) {
        if(internalMap[iDof] > -1) {
          int ccDof = cToCC[internalMap[iDof]];
          if(ccDof > -1) deltaF[ccDof] = iDisp[iDof];
        }
      }
    }
    delete [] iDisp;

    if(deltaF) {
      for(iDof = 0; iDof < totalInterfSize; ++iDof)
        if(allBoundDofs[iDof] >= 0) { // deltaF of ctc and mpc nodes computed in applyScalingAndB below
          deltaF[allBoundDofs[iDof]] = res[dualToBoundary[iDof]];
        }
    }
  }

  if(solInfo().getFetiInfo().precno == FetiInfo::identity) {
    for(iDof = 0; iDof < boundLen; ++iDof) res[iDof] = v[iDof];
  }

  // Return preconditioned u
  applyScalingAndB(res, Pu, localw);
}

template<class Scalar>
void
GenSubDomain<Scalar>::getQtKQ(GenSolver<Scalar> *s)
{
 if(numMPC == 0) return;

 int numDOFs = localLen();
 locKpQ = new Scalar[numMPC*numDOFs];
 int i;
 for(i = 0; i < numMPC*numDOFs; ++i)
   locKpQ[i] = 0.0;

 // loop over mpc structure and fill coefficients
 int iMPC;
 for(iMPC = 0; iMPC < numMPC; ++iMPC) {
   for(i = 0; i < mpc[iMPC]->nterms; ++i) {
     int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
                            (1 << mpc[iMPC]->terms[i].dofnum));
     if(dof >= 0) {
       locKpQ[dof + iMPC*numDOFs] = mpc[iMPC]->terms[i].coef;
     }
   }
 }

 for(iMPC = 0; iMPC < numMPC; ++iMPC)
   if(s) s->reSolve(locKpQ + iMPC*numDOFs);

 QtKpBt = new Scalar [numMPC*totalInterfSize];

 for(i = 0; i < numMPC*totalInterfSize; ++i)
   QtKpBt[i] = 0.0;

 int j;
 for(i = 0; i < numMPC; ++i)
   for(j=0; j<totalInterfSize; ++j)
     QtKpBt[j+i*totalInterfSize] = locKpQ[allBoundDofs[j]+i*numDOFs];

 qtkq = new GenFullM<Scalar>(numMPC,numMPC);
 qtkq->zero();

 int jMPC;
 for(iMPC = 0; iMPC < numMPC; ++iMPC)
   for(jMPC = 0; jMPC < numMPC; ++jMPC)
     for(i=0; i<mpc[iMPC]->nterms; ++i) {
       int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
                                (1 << mpc[iMPC]->terms[i].dofnum));
       if(dof < 0) continue;

       (*qtkq)[iMPC][jMPC] += mpc[iMPC]->terms[i].coef * locKpQ[dof+jMPC*numDOFs];
     }

 delete [] locKpQ; locKpQ = 0;
}

template<class Scalar>
void
GenSubDomain<Scalar>::getQtKQ(int glMPCnum, Scalar *QtKQ)
{
 int thisMPC=  globalToLocalMPC[glMPCnum];
 int iMPC;
 for(iMPC = 0; iMPC < numMPC; ++iMPC)
   QtKQ[iMPC] = (*qtkq)[thisMPC][iMPC];
}

template<class Scalar>
void
GenSubDomain<Scalar>::multQtKBt(int glMPCnum, Scalar *G, Scalar *QtKBtG,
                                Scalar alpha, Scalar beta)
{
 int iMPC = globalToLocalMPC[glMPCnum];
 // QtKBtG = QtKBt*G
 Tgemv('T', totalInterfSize, 1, alpha, QtKpBt+iMPC*totalInterfSize,
       totalInterfSize, G, 1, beta, QtKBtG, 1);
}

template<class Scalar>
void
GenSubDomain<Scalar>::multQt(int glMPCnum, Scalar *V, int numV, Scalar *QtV)
{
 int numDofs = localLen();
 int iMPC = globalToLocalMPC[glMPCnum];
 int i,n;
 for(n = 0; n < numV; ++n) {
   for(i = 0; i < mpc[iMPC]->nterms; ++i) {
     int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
                              (1 << mpc[iMPC]->terms[i].dofnum));
     if(dof < 0) continue;
     Scalar *beta = V+n*numDofs;
     QtV[n] += mpc[iMPC]->terms[i].coef * beta[dof];
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multQt(int glMPCnum, Scalar *beta, Scalar *result)
{
 int iMPC = globalToLocalMPC[glMPCnum];
 int i;
 for(i = 0; i < mpc[iMPC]->nterms; ++i) {
   int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
                            (1 << mpc[iMPC]->terms[i].dofnum));
   if(dof < 0) continue;
   result[0] += mpc[iMPC]->terms[i].coef * beta[dof];
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::clean_up()
{
 int i;
 for(i=0; i < numele; ++i) packedEset[i]->renum(glNums);
 if(weight) { delete [] weight; weight = 0; }
 if(scaling) { delete [] scaling; scaling = 0; }
 if(bcx) { delete [] bcx; bcx = 0; }
 if(vcx) { delete [] vcx; vcx = 0; }
 if(acx) { delete [] acx; acx = 0; }
 if(boundMap) { delete [] boundMap; boundMap = 0; }
 if(dualToBoundary) { delete [] dualToBoundary; dualToBoundary = 0; }
 if(internalMap) { delete [] internalMap; internalMap = 0; }

 long m1 = 0;

 if(Kbb) {
   m1 = -memoryUsed();
   Kbb->clean_up();
   m1 += memoryUsed();
 }

 if(Kib) {
   m1 = -memoryUsed();
   Kib->clean_up();
   m1 += memoryUsed();
 }

 if(KiiSparse) {
   m1 = -memoryUsed();
   KiiSparse->clean_up();
   m1 += memoryUsed();
 }

 if(Krr) {
   m1 = -memoryUsed();
   Krr->clean_up();
   m1 += memoryUsed();
 }

 if(Krc) {
   m1 = -memoryUsed();
   Krc->clean_up();
   m1 += memoryUsed();
 }

 if(Kcc) {
   m1 = -memoryUsed();
   Kcc->clean_up();
   m1 += memoryUsed();
 }

 if(cornerMap) { delete [] cornerMap; cornerMap = 0; }
 if(cornerEqNums) { delete [] cornerEqNums; cornerEqNums = 0; }
 if(cornerNodes) { delete [] cornerNodes; cornerNodes = 0; }
 if(glCornerNodes) { delete [] glCornerNodes; glCornerNodes = 0; }
 if(ccToC) { delete [] ccToC; ccToC = 0; }
 if(fcstar) { delete [] fcstar; fcstar = 0; }
 if(dsa) dsa->clean_up();
 if(c_dsa) c_dsa->clean_up();
 if(cc_dsa) cc_dsa->clean_up();
 if(neighbNumGRBMs) { delete [] neighbNumGRBMs; neighbNumGRBMs = 0; }

 // pade
  if(a) { for(int i=0; i<ia; ++i) delete a[i]; delete [] a; a = 0; }
  if(b) { for(int i=0; i<ib; ++i) delete b[i]; delete [] b; b = 0;}
  if(P) { delete P; P = 0; }
  if(Q) { delete Q; Q = 0; }
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums)
{
  // extracts the mpc residual components from subv and inserts them in global vector (mpcv)
  for(int i=0; i<scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    mpcv[mpcEqNums->firstdof(localToGlobalMPC[locMpcNb])] = subv[scomm->mapT(SComm::mpc,i)];
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::insertMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums)
{
  // extracts the mpc residual components from mpcv and inserts them in the interface vector subv
  for(int i=0; i<scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    subv[scomm->mapT(SComm::mpc,i)] = mpcv[mpcEqNums->firstdof(localToGlobalMPC[locMpcNb])];
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::setMpcRhs(Scalar *interfvec, double _t)
{
  // set the rhs of inequality mpcs to the geometric gap and reset the rhs of the equality mpcs to the original rhs
  // (used in nonlinear analysis)
  // idea: initalize to dual-active if gap is open (+ve) // XXXX
  for(int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpc[locMpcNb]->getSource() != mpc::ContactSurfaces) {
      mpc[locMpcNb]->rhs = mpc[locMpcNb]->original_rhs - interfvec[scomm->mapT(SComm::mpc,i)];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::updateMpcRhs(Scalar *interfvec)
{
  bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
  for(int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpcFlag[locMpcNb] == true) {
      mpc[locMpcNb]->rhs = mpc[locMpcNb]->original_rhs + interfvec[scomm->mapT(SComm::mpc,i)];
      mpcFlag[locMpcNb] = false;
    }
  }
}

template<class Scalar>
double
GenSubDomain<Scalar>::getMpcError()
{
  double ret = 0;
  for(int i = 0; i < numMPC; ++i) {
    if(mpcMaster[i]) {
      if(mpc[i]->type == 0) {
        ret += ScalarTypes::sqNorm(mpc[i]->rhs);
      }
      else if(mpc[i]->type == 1) ret += ScalarTypes::sqNorm(MIN(0.0, mpc[i]->rhs));
    }
  }
  return ret;
} 

template<class Scalar>
void
GenSubDomain<Scalar>::solveLocalCCt(Scalar *subv)
{
  // used for block diagonal CCt preconditioner for MPCs
  if(numMPC > 0) {
    int i;

    // Step 1: extract the mpc residual components from subv and inserts them in a local vector (mpcv)
    GenVector<Scalar> mpcv(numMPC, 0.0);
    for(i=0; i<scomm->lenT(SComm::mpc); ++i)
      mpcv[scomm->mpcNb(i)] = subv[scomm->mapT(SComm::mpc,i)];

    // Step 2: solve CCt^-1 * mpcv
    localCCtsolver->reSolve(mpcv);

    // Step 3: redistribute mpcv to the interface vector subv
    for(i=0; i<scomm->lenT(SComm::mpc); ++i)
      subv[scomm->mapT(SComm::mpc,i)] = mpcv[scomm->mpcNb(i)];
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractBlockMpcResidual(int block, Scalar *subv, GenVector<Scalar> *mpcv,
                                              SimpleNumberer *blockMpcEqNums)
{
  // extracts the mpc residual components from subv and inserts them in global vector (mpcv)
  // modified to loop only over the subdomain lmpcs belonging to block
  // Make sure that the input block (global) Id is in the range of blocks Id that have contribution
  // from this subdomain
  int i,j;
  for(j=0; j<blockToLocalMpc->num(block); ++j) {
    i = (*blockToLocalMpc)[block][j];
    int iDof = (*mpcToBoundDof)[i][0];
    int bij = (*blockToBlockMpc)[block][j];
    (*mpcv)[blockMpcEqNums->firstdof(bij)] = subv[iDof];
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::insertBlockMpcResidual(Scalar *subv, GenVector<Scalar> **mpcv, Connectivity *mpcToBlock,
                                             SimpleNumberer **blockMpcEqNums)
{
  // extracts the mpc residual components from mpcv and inserts them in the interface vector subv
  int i, j, k;
  for(i = 0; i < numMPC; ++i) {
    int iDof = (*mpcToBoundDof)[i][0];
    int gi = localToGlobalMPC[i];
    double w = double(mpcToBlock->num(gi));
    subv[iDof] = 0.0;
    for(j=0; j<localMpcToBlock->num(i); ++j) {
      int jBlock = (*localMpcToBlock)[i][j];
      int bij = (*localMpcToBlockMpc)[i][j];
      subv[iDof] += (*mpcv[jBlock])[blockMpcEqNums[jBlock]->firstdof(bij)] / w;
    }
    for(k=1; k<mpcToBoundDof->num(i); ++k) {
      int kDof = (*mpcToBoundDof)[i][k];
      subv[kDof] = subv[(*mpcToBoundDof)[i][0]];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::setMpcCommSize(FSCommPattern<Scalar> *mpcPat)
{
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc,i);
    int len = (subNumber != neighb) ? scomm->lenT(SComm::mpc,i) : 0;
    mpcPat->setLen(subNumber, neighb, len);
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec)
{
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc,i);
    if(subNumber != neighb) {
      FSSubRecInfo<Scalar> sInfo = mpcPat->getSendBuffer(subNumber, neighb);
      for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j)
        sInfo.data[j] = interfvec[scomm->mapT(SComm::mpc,i,j)];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::combineMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec)
{
 // this function combines the interfvec values corresponding to mpcs
 // that are shared between subdomains. Used with approximated blockDiag and diag CCt preconditioners
 int i,j;
 Scalar *mpcCombo = (Scalar *) dbg_alloca(sizeof(Scalar)*numMPC);
 int *mpcCount =  (int *) dbg_alloca(sizeof(int)*numMPC);
 for(i = 0; i < numMPC; ++i) mpcCount[i] = 1;
 for(i = 0; i < scomm->numT(SComm::mpc); ++i) {
   int neighb = scomm->neighbT(SComm::mpc, i);
   if(subNumber != neighb) {
     FSSubRecInfo<Scalar> rInfo = mpcPat->recData(neighb, subNumber);
     for(j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
       int locMpcNb = scomm->mpcNb(i,j);
       if(mpcCount[locMpcNb] == 1)
         mpcCombo[locMpcNb] = interfvec[scomm->mapT(SComm::mpc,i,j)];
       mpcCount[locMpcNb]++;
       mpcCombo[locMpcNb] += rInfo.data[j];
     }
   }
 }
 // now add mpcCombo to interfvec
 for(i = 0; i < scomm->numT(SComm::mpc); ++i) {
   int neighb = scomm->neighbT(SComm::mpc, i);
   if(subNumber != neighb) {
     for(j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
       int locMpcNb = scomm->mpcNb(i,j);
       interfvec[scomm->mapT(SComm::mpc,i,j)] = mpcCombo[locMpcNb]/double(mpcCount[locMpcNb]);
     }
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::setMpcCommSize(FSCommPattern<int> *mpcPat)
{
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc,i);
    int len = (subNumber != neighb) ? scomm->lenT(SComm::mpc,i) : 0;
    mpcPat->setLen(subNumber, neighb, len);
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendMpcStatus(FSCommPattern<int> *mpcPat, int flag)
{
  // if flag = -1 only send status of mpcMaster otherwise send -1
  // note: could use SComm::ieq list
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc,i);
    if(subNumber != neighb) {
      FSSubRecInfo<int> sInfo = mpcPat->getSendBuffer(subNumber, neighb);
      for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
        int locMpcNb = scomm->mpcNb(i,j);
        if(flag == -1) sInfo.data[j] = (mpcMaster[locMpcNb]) ? int(!mpc[locMpcNb]->active) : -1; else
        sInfo.data[j] = int(!mpc[locMpcNb]->active);
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::recvMpcStatus(FSCommPattern<int> *mpcPat, int flag, bool &statusChange)
{
 // this function is to make sure that the status of an mpc is the same in all subdomains which share it
 // needed to due to possible roundoff error
 // if flag == 1 then make dual constraint not active in all subdomains if not active in at least one (use in proportioning step)
 // if flag == 0 then make dual constraint active in all subdomains if it is active in at least one (use in expansion step)
 // if flag == -1 use mpc master status in all subdomains
 // note: could use SComm::ieq list
 int i,j;
 bool *tmpStatus = (bool *) alloca(sizeof(bool)*numMPC);
 for(int i=0; i<numMPC; ++i) tmpStatus[i] = !mpc[i]->active;
 for(i = 0; i < scomm->numT(SComm::mpc); ++i) {
   int neighb = scomm->neighbT(SComm::mpc, i);
   if(subNumber != neighb) {
     FSSubRecInfo<int> rInfo = mpcPat->recData(neighb, subNumber);
     for(j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
       int locMpcNb = scomm->mpcNb(i,j);
       if(flag == -1) tmpStatus[locMpcNb] = (rInfo.data[j] > -1) ? bool(rInfo.data[j]) : tmpStatus[locMpcNb]; else // XXXX
       tmpStatus[locMpcNb] = (flag == 1) ? (tmpStatus[locMpcNb] || bool(rInfo.data[j])) : (tmpStatus[locMpcNb] && bool(rInfo.data[j]));
     }
   }
 }

 bool print_debug = false;
 statusChange = false;
 for(i = 0; i < numMPC; ++i) {
   if(solInfo().getFetiInfo().contactPrintFlag && mpcMaster[i]) {
     if(!mpc[i]->active && !tmpStatus[i]) { cerr << "-"; if(print_debug) cerr << " recvMpcStatus: sub = " << subNumber << ", mpc = " << localToGlobalMPC[i] << endl; }
     else if(mpc[i]->active && tmpStatus[i]) { cerr << "+"; if(print_debug) cerr << " recvMpcStatus: sub = " << subNumber << ", mpc = " << localToGlobalMPC[i] << endl; }
   }
   mpc[i]->active = !tmpStatus[i];
   if(mpcStatus2[i] == mpc[i]->active) statusChange = true;
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::printMpcStatus()
{
 for(int i = 0; i < numMPC; ++i) {
   if(mpc[i]->type == 1) {
     cerr<< (mpc[i]->active ? 'o' : 'x');
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::initMpcStatus()
{
 for(int i = 0; i < numMPC; ++i) {
   mpc[i]->active = false;
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::saveMpcStatus()
{
 // this saves the status before first update iteration so it can be reset if nonmonotic
 if(!mpcStatus) mpcStatus = new int[numMPC];
 for(int i = 0; i < numMPC; ++i) {
   mpcStatus[i] = int(!mpc[i]->active);
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::restoreMpcStatus()
{
 for(int i = 0; i < numMPC; ++i) {
   if(solInfo().getFetiInfo().contactPrintFlag && mpcMaster[i]) {
     if(!mpc[i]->active && !mpcStatus[i]) cerr << "-";
     else if(mpc[i]->active && mpcStatus[i]) cerr << "+";
   }
   mpc[i]->active = bool(!mpcStatus[i]);
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::saveMpcStatus1()
{
  if(!mpcStatus1) mpcStatus1 = new bool[numMPC];
  for(int i = 0; i < numMPC; ++i) mpcStatus1[i] = !mpc[i]->active;
}

template<class Scalar>
void
GenSubDomain<Scalar>::saveMpcStatus2()
{
  if(!mpcStatus2) mpcStatus2 = new bool[numMPC];
  for(int i = 0; i < numMPC; ++i) mpcStatus2[i] = !mpc[i]->active;
}

template<class Scalar>
void
GenSubDomain<Scalar>::cleanMpcData()
{
  if(mpcStatus) { delete [] mpcStatus; mpcStatus = 0; }
  if(mpcStatus1) { delete [] mpcStatus1; mpcStatus1 = 0; }
  if(mpcStatus2) { delete [] mpcStatus2; mpcStatus2 = 0; }
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyBtransposeAndScaling(Scalar *u, Scalar *v, Scalar *deltaU, Scalar *localw)
{
  int i, iDof, k;
  bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(i = 0; i < numMPC; ++i) mpcFlag[i] = true;

  for(iDof = 0; iDof < totalInterfSize; ++iDof) {
    switch(boundDofFlag[iDof]) {
      case 0:
        v[dualToBoundary[iDof]] += u[iDof] * scaling[iDof];
        if(deltaU) deltaU[allBoundDofs[iDof]] = -v[dualToBoundary[iDof]];
        break;
      case 1: { // wet interface
        int windex = -1-allBoundDofs[iDof];
        localw[windex] = u[iDof]*scaling[iDof];
        deltaFwi[windex] = u[iDof];
      } break;
      case 2: { // dual mpc
        int locMpcNb = -1-allBoundDofs[iDof];
        if(mpcFlag[locMpcNb]) {
          SubLMPCons<Scalar> *m = mpc[locMpcNb];
          if(!mpc[locMpcNb]->active) {
            for(k = 0; k < m->nterms; k++) {
              int cdof = (m->terms)[k].cdof;
              if(cdof >= 0) { // mpc dof that exists
                Scalar coef = (m->terms)[k].coef / m->k[k]; // 1/m->k[k] = A, see generalized preconditioner
                if(invBoundMap[cdof] < 0) cerr << "error here in GenSubDomain<Scalar>::applyBtransposeAndScaling\n";
                v[invBoundMap[cdof]] += u[iDof] * coef * scaling[iDof];
              }
            }
          }
          mpcFlag[locMpcNb] = false;
        }
      } break;
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyScalingAndB(Scalar *res, Scalar *Pu, Scalar *localw)
{
  int i, iDof, k;

  if(numMPC && !deltaFmpc)
    deltaFmpc = new Scalar[numMPC]; // only need to allocate 1st time (unless numMPC changes)
  bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(i = 0; i < numMPC; ++i) mpcFlag[i] = true;

  // Return preconditioned u
  for(iDof = 0; iDof < totalInterfSize; ++iDof) {
    switch(boundDofFlag[iDof]) {
      case 0:
        Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
        break;
      case 1: // wet interface
        Pu[iDof] = localw[-1-allBoundDofs[iDof]] * scaling[iDof];
        break;
      case 2: { // dual mpc or contact
        int locMpcNb = -1-allBoundDofs[iDof];
        SubLMPCons<Scalar> *m = mpc[locMpcNb];
        Pu[iDof] = 0.0;
        if(mpcFlag[locMpcNb]) deltaFmpc[locMpcNb] = 0.0;
        if(!mpc[locMpcNb]->active) {
          for(k = 0; k < m->nterms; k++) {
            int cdof = (m->terms)[k].cdof;
            if(cdof > -1) { // mpc dof that exists
              Scalar coef = (m->terms)[k].coef;
              Pu[iDof] += res[invBoundMap[cdof]] * coef * scaling[iDof] / m->k[k]; // 1/m->k[k] = A, see generalized preconditioner
              if(mpcFlag[locMpcNb]) deltaFmpc[locMpcNb] += res[invBoundMap[cdof]] * coef;
            }
          }
        }
        mpcFlag[locMpcNb] = false;
      } break;
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::subtractMpcRhs(Scalar *interfvec)
{
  for(int i=0; i < scomm->lenT(SComm::mpc); ++i) {
    interfvec[scomm->mapT(SComm::mpc,i)] -= mpc[scomm->mpcNb(i)]->rhs;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::setLocalLambda(Scalar *_localLambda)
{
  if(localLambda) delete [] localLambda;
  localLambda = new double[totalInterfSize];
  for(int i=0; i<totalInterfSize; ++i) localLambda[i] = ScalarTypes::Real(_localLambda[i]);
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeContactPressure(Scalar *globStress, Scalar *globWeight)
{
  for(int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpc[locMpcNb]->type == 1) { // inequality constraint
      for(int j = 0; j < mpc[locMpcNb]->nterms; ++j) {
        int node = mpc[locMpcNb]->terms[j].nnum;
        int glNode = (domain->outFlag) ? domain->nodeTable[glNums[node]]-1 : glNums[node];
        globStress[glNode] += abs(localLambda[scomm->mapT(SComm::mpc,i)]);
        globWeight[glNode] += 1.0;
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeLocalContactPressure(Scalar *stress, Scalar *weight)
{
  for(int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpc[locMpcNb]->type == 1) { // inequality constraint
      for(int j = 0; j < mpc[locMpcNb]->nterms; ++j) {
        int node = mpc[locMpcNb]->terms[j].nnum;
        stress[node] += abs(localLambda[scomm->mapT(SComm::mpc,i)]);
        weight[node] += 1.0;
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::getLocalMpcForces(double *mpcLambda, DofSetArray *cornerEqs,
                                        int mpcOffset, GenVector<Scalar> &uc)
{
// XXXX needs some work to map both dual and primal into single mpcLambda array
  if(numMPC > 0 && numMPC_primal > 0) cerr << "unsupported feature in GenSubDomain::getLocalMpcForces \n";
  for(int i = 0; i < scomm->lenT(SComm::mpc); ++i) { // dual mpcs
    int locMpcNb = scomm->mpcNb(i);
    if(localLambda) mpcLambda[locMpcNb] = localLambda[scomm->mapT(SComm::mpc,i)];
    else mpcLambda[locMpcNb] = 0;
  }
  for(int i = 0; i < numMPC_primal; ++i) {
    int glMpcNb = localToGlobalMPC_primal[i];
    int dof = cornerEqs->firstdof(mpcOffset+glMpcNb);
    mpcLambda[i] = -ScalarTypes::Real(uc[dof]);
  }
  if(salinasFlag) for(int i = 0; i < numMPC+numMPC_primal; ++i) mpcLambda[i] = -mpcLambda[i];  // different sign convention
}

template<class Scalar>
void
GenSubDomain<Scalar>::getConstraintMultipliers(std::map<std::pair<int,int>,double> &mu, std::vector<double> &lambda)
{
  bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

  for(int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
    int l = scomm->mpcNb(i);
    if(mpcFlag[l]) {
      if(mpc[l]->getSource() == mpc::ContactSurfaces) {
        mu[mpc[l]->id] = (localLambda) ? localLambda[scomm->mapT(SComm::mpc,i)] : 0;
      } else {
        lambda.push_back( (localLambda) ? localLambda[scomm->mapT(SComm::mpc,i)] : 0 );
      }
      mpcFlag[l] = false;
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::initialize()
{
  kweight = 0; deltaFmpc = 0; scaling = 0; 
  rigidBodyModes = 0; rigidBodyModesG = 0;
  Src = 0; Krr = 0; KrrSparse = 0; BKrrKrc = 0; Kcc = 0; Krc = 0;
  Grc = 0; rbms = 0; interfaceRBMs = 0; qtkq = 0; KiiSparse = 0;
  KiiSolver = 0; Kib = 0; MPCsparse = 0; Kbb = 0; corotators = 0;
  fcstar = 0; QtKpBt = 0; locKpQ = 0; glInternalMap = 0; glBoundMap = 0;
  mpcForces = 0; mpc = 0; mpc_primal = 0; localCCtsolver = 0; localCCtsparse = 0; diagCCt = 0;
  Kww = 0; Kcw = 0, Krw = 0; neighbKww = 0; localw = 0; localw_copy = 0; Kcw_mpc = 0;
  deltaFwi = 0; M = 0; Muc = 0;
  wweight = 0;
  lengthCCtData = 0; CCtrow = 0; CCtcol = 0; CCtval = 0;
  bcx_scalar = 0;
  a = 0; b = 0; P = 0; Q = 0; rebuildPade = true;  // pade
  numC_deriv = 0; C_deriv = 0; Cuc_deriv = 0;
#ifdef HB_COUPLED_PRECOND
  kSumWI = 0;
  precNodeToNode = 0;
#endif
  weightPlus = 0;
  mpcStatus = 0; mpcStatus1 = 0; mpcStatus2 = 0;
  G = 0; neighbG = 0;
  sharedRstar_g = 0; tmpRstar_g = 0;
  l2g = 0;
}

template<class Scalar>
GenSubDomain<Scalar>::~GenSubDomain()
{
  if(KiiSparse) { delete KiiSparse; KiiSparse = 0; KiiSolver = 0; }
  if(Krr) { delete Krr; Krr = 0; KrrSparse = 0; }
  if(Kbb) { delete Kbb; Kbb = 0; }
  if(scaling) { delete [] scaling; scaling = 0; }
  if(kweight) { delete [] kweight; kweight = 0; }
  if(Kcc) { delete Kcc; Kcc = 0; }
  if(fcstar) { delete [] fcstar; fcstar = 0; }
  if(Kib) { delete Kib; Kib = 0; }
  if(deltaFmpc) { delete [] deltaFmpc; deltaFmpc = 0; }
  if(BKrrKrc) {
    if(BKrrKrc[0]) { delete [] BKrrKrc[0]; BKrrKrc[0] = 0; }
    delete [] BKrrKrc; BKrrKrc = 0;
  }
  if(interfaceRBMs) { delete [] interfaceRBMs; interfaceRBMs = 0; }
  if(Grc) { delete Grc; Grc = 0; }
  if(locKpQ) { delete [] locKpQ; locKpQ = 0; }
  if(QtKpBt) { delete [] QtKpBt; QtKpBt = 0; }
  if(qtkq) { delete qtkq; qtkq = 0; }
  if(Krc) { delete Krc; Krc = 0; }
  if(rigidBodyModes) { delete rigidBodyModes; rigidBodyModes = 0; }
  if(rigidBodyModesG) { delete rigidBodyModesG; rigidBodyModesG = 0; }
  if(Src) { delete Src; Src = 0; }
  if(MPCsparse) { delete MPCsparse; MPCsparse = 0; }
  if(glInternalMap) { delete [] glInternalMap; glInternalMap = 0; }
  if(glBoundMap) { delete [] glBoundMap; glBoundMap = 0; }
  if(mpcForces) { delete [] mpcForces; mpcForces = 0; }
  if(mpc) {
    for(int i=0; i<numMPC; ++i)
      if(mpc[i]) { delete mpc[i]; mpc[i] = 0; }
    delete [] mpc; mpc = 0;
  }
  if(mpc_primal) {
    for(int i=0; i<numMPC_primal; ++i)
      if(mpc_primal[i]) { delete mpc_primal[i]; mpc_primal[i] = 0; }
    delete [] mpc_primal; mpc_primal = 0;
  }
  if(diagCCt) { delete [] diagCCt; diagCCt = 0; }
  if(localw) { delete [] localw; localw = 0; }
  if(localw_copy) { delete [] localw_copy; localw_copy = 0; }
  if(deltaFwi) { delete [] deltaFwi; deltaFwi = 0; }
  if(Kww) { delete Kww; Kww = 0; }
  if(Kcw) { delete Kcw; Kcw = 0; }
  if(Krw) { delete Krw; Krw = 0; }
  if(neighbKww) { delete neighbKww; neighbKww = 0; }
  if(wweight) { delete [] wweight; wweight = 0; }
  if(M) { delete M; M = 0; }
  if(Muc) { delete Muc; Muc = 0; }
  if(CCtrow) { delete [] CCtrow; CCtrow = 0; }
  if(CCtcol) { delete [] CCtcol; CCtcol = 0; }
  if(CCtval) { delete [] CCtval; CCtval = 0; }
  lengthCCtData = 0;
  // don't delete boundaryDOFs, rbms
  if(bcx_scalar) { delete [] bcx_scalar; bcx_scalar = 0; }
  if(Kcw_mpc) { delete Kcw_mpc; }
  // pade
  if(a) { for(int i=0; i<ia; ++i) delete a[i]; delete [] a; }
  if(b) { for(int i=0; i<ib; ++i) delete b[i]; delete [] b; }
  if(P) delete P; if(Q) delete Q;
  if(C_deriv) { for(int i=0; i<numC_deriv; ++i) delete C_deriv[i]; delete [] C_deriv; }
  if(Cuc_deriv) { for(int i=0; i<numC_deriv; ++i) delete Cuc_deriv[i]; delete [] Cuc_deriv; }
#ifdef HB_COUPLED_PRECOND
  if(kSumWI) { delete kSumWI; kSumWI = 0; }
  if(isMixedSub && precNodeToNode) { delete precNodeToNode; precNodeToNode = 0; }
#endif
  if(weightPlus) { delete [] weightPlus; weightPlus = 0; }
  if(mpcStatus) { delete [] mpcStatus; mpcStatus = 0; }
  if(mpcStatus1) { delete [] mpcStatus1; mpcStatus1 = 0; }
  if(mpcStatus2) { delete [] mpcStatus2; mpcStatus2 = 0; }

  if(sharedRstar_g) { delete sharedRstar_g; sharedRstar_g = 0; }
  if(tmpRstar_g) { delete tmpRstar_g; tmpRstar_g = 0; }
  int numPCN = scomm->numT(SComm::mpc);
  if(G) {
    for(int i=0; i<numPCN; ++i)
      if(G[i]) { delete G[i]; G[i] = 0; }
     delete [] G; G = 0;
  }
  if(neighbG) {
    for(int i=0; i<numPCN; ++i)
      if(neighbG[i]) { delete neighbG[i]; neighbG[i] = 0; }
     delete [] neighbG; neighbG = 0;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::deleteMPCs()
{
  if(mpc) {
    for(int i = 0; i < numMPC; ++i)
      if(mpc[i]) { delete mpc[i]; mpc[i] = 0; }
    delete [] mpc; mpc = 0;
  }
  if(localToGlobalMPC) { delete [] localToGlobalMPC; localToGlobalMPC = 0; }
  deleteG();
  numMPC = 0;
  scomm->deleteTypeSpecificList(SComm::mpc);

  scomm->deleteTypeSpecificList(SComm::all);
  if(mpcStatus) { delete [] mpcStatus; mpcStatus = 0; }
  if(mpcStatus1) { delete [] mpcStatus1; mpcStatus1 = 0; }
  if(mpcStatus2) { delete [] mpcStatus2; mpcStatus2 = 0; }
  if(deltaFmpc) { delete [] deltaFmpc; deltaFmpc = 0; }
  //if(mpcForces) { delete [] mpcForces; mpcForces = 0; }
  if(diagCCt) { delete [] diagCCt; diagCCt = 0; }
  if(mpcMaster) { delete [] mpcMaster; mpcMaster = 0; }
  if(CCtrow) { delete [] CCtrow; CCtrow = 0; }
  if(CCtcol) { delete [] CCtcol; CCtcol = 0; }
  if(CCtval) { delete [] CCtval; CCtval = 0; }

  if(mpcToDof) { delete mpcToDof; mpcToDof = 0; }
  if(localMpcToMpc) { delete localMpcToMpc; localMpcToMpc = 0; }

}

template<class Scalar>
void
GenSubDomain<Scalar>::makeEdgeVectorsPlus(bool isFluidSub)
{
  int i;
  //bool isCoupled = solInfo().isCoupled;
  int spaceDim = solInfo().solvercntl->fetiInfo.spaceDimension;
  // Number of directions in the coarse problem, choose from 0, 1, ..., 13
  int numDirec = solInfo().getFetiInfo().numdir;
  int numWaves = spaceDim;      // Number of long and trans waves
  if(isFluidSub)
    numWaves = 1;
  int numCS = 2;  // Choose cos or/and sin mode

  // If no direction added, put numCS = 0, and numWaves = 0.
  if(numDirec == 0) {
    numWaves = 0;
    numCS    = 0;
  }

  // To store the Directions of the long and transverse waves
  double *d_x  = (double *) dbg_alloca(sizeof(double)*numDirec);
  double *d_y  = (double *) dbg_alloca(sizeof(double)*numDirec);
  double *d_z  = (double *) dbg_alloca(sizeof(double)*numDirec);
  double *t_x  = (double *) dbg_alloca(sizeof(double)*numDirec);
  double *t_y  = (double *) dbg_alloca(sizeof(double)*numDirec);
  double *t_z  = (double *) dbg_alloca(sizeof(double)*numDirec);

  double *wDir_x = 0;
  double *wDir_y = 0;
  double *wDir_z = 0;

  double pi = 3.141592653589793;
  if(spaceDim == 2) {
    for(int i = 0; i < numDirec; i++) {
      d_x[i] = cos( (pi*i)/(numDirec) );
      d_y[i] = sin( (pi*i)/(numDirec) );
      d_z[i] = 0.0;
      t_x[i] = -d_y[i];
      t_y[i] = d_x[i];
      t_z[i] = 0.0;
    }
  }
  else {
    getDirections(numDirec, numWaves, wDir_x, wDir_y, wDir_z);
    for(i=0; i<numDirec; i++) {
      d_x[i] = wDir_x[numWaves*i+0];
      d_y[i] = wDir_y[numWaves*i+0];
      d_z[i] = wDir_z[numWaves*i+0];
    }
  }

  int numR = solInfo().getFetiInfo().nGs;
  int totalLengthGrc=0;
  DofSet desired;
  if(isFluidSub) {
    if(numR > 0) numR = 1;
    desired = DofSet::Helm;
  }
  else {
    if(solInfo().getFetiInfo().rbmType == FetiInfo::rotation || (numR > 3)) {
      numR = 6;
      desired = DofSet::XYZdisp | DofSet::XYZrot;
    }
    else desired = DofSet::XYZdisp;
  }

  Connectivity &sharedNodes = *(scomm->sharedNodes);
  // edgeDofSize: number of augmentation degree of freedom per edge
  if(!edgeDofSize) {
    edgeDofSize = new int[scomm->numNeighb];
    for(i=0; i<scomm->numNeighb; ++i) edgeDofSize[i] = 0;
  }
  if(!edgeDofSizeTmp) edgeDofSizeTmp = new int[scomm->numNeighb];
  // total: total number of augmentaions per sub
  int total = 0;

  int numdofperNode = 0;
  if(isFluidSub)
    numdofperNode = 1;
  else if(solInfo().getFetiInfo().waveType == FetiInfo::solid)
    numdofperNode = 3;
  else numdofperNode = 6;

  int numInterfNodes = sharedNodes.numConnect();
  int lQ = (numR + numDirec*numWaves*numCS)*numdofperNode*numInterfNodes;

  int nQPerNeighb = numR + numDirec*numWaves*numCS;
  double *Q = new double[lQ];
  bool *isUsed = (bool *) dbg_alloca(sizeof(bool)*(nQPerNeighb*scomm->numNeighb));
  for(i = 0; i < nQPerNeighb*scomm->numNeighb; ++i)
    isUsed[i] = false;
  for(i = 0; i < lQ; ++i)
    Q[i] = 0;

  // 1. first count number of edge dofs
  int iSub, iNode;
  DofSet *found = new DofSet[scomm->numNeighb];
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    edgeDofSizeTmp[iSub]=0;
    if(scomm->subNums[iSub] == subNumber) continue;
    int nhelm = 0;
    int nx = 0, ny = 0, nz = 0;
    int nxr = 0, nyr = 0, nzr = 0;
    int label = 0;
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      if((isFluidSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) ||
         (!isFluidSub && boundaryDOFs[iSub][iNode].contains(DofSet::Helm))) continue;

      if(boundaryDOFs[iSub][iNode].contains(DofSet::Helm))
        nhelm++;
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp))
        nx++;
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp))
        ny++;
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp))
        nz++;
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Xrot))
        nxr++;
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Yrot))
        nyr++;
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Zrot))
        nzr++;

      if(boundaryDOFs[iSub][iNode].containsAllDisp(spaceDim)) label++;

      found[iSub] |= boundaryDOFs[iSub][iNode] & desired;
    }

    // Check if we should add rotation for problems that do not
    // define rotational degrees of freedom
    if(desired.contains(DofSet::XYZrot)) {
      // if ny +nz is bigger than 2 (at least 3) there is no danger in puttin g a Xrot
      if(ny  + nz > 2) found[iSub] |= DofSet::Xrot;
      // if nx+nz > 2 and there are some y AND some z, we add the Yrot
      if(nx + nz > 2 && (ny*nx) != 0) found[iSub] |= DofSet::Yrot;
      if(nx + ny > 2 && (nx*ny*nz) != 0) found[iSub] |= DofSet::Zrot;
    }

    if(numR > 0) {
      if(nhelm > 0)
        isUsed[iSub*nQPerNeighb] = true;
      if(nx > 0)
        isUsed[0+iSub*nQPerNeighb] = true;
      if(ny > 0)
        isUsed[1+iSub*nQPerNeighb] = true;
      if(nz > 0)
        isUsed[2+iSub*nQPerNeighb] = true;
      if(desired.contains(DofSet::XYZrot)) {
        if((nxr > 0) || ( ny + nz > 2))
          isUsed[3+iSub*nQPerNeighb] = true;
        if((nyr > 0) || ( nx + nz > 2))
          isUsed[4+iSub*nQPerNeighb] = true;
        if((nzr > 0) || ( nx + ny > 2))
          isUsed[5+iSub*nQPerNeighb] = true;
      }
    }

    if((label > 0) || (nhelm > 0))
      for(i=0; i<numDirec*numWaves*numCS; i++)
        isUsed[numR+i+iSub*nQPerNeighb] = true;
    if((label == 0) && (nhelm == 0))
      edgeDofSizeTmp[iSub] = found[iSub].count();
    else {
      if(numR > 0)
        edgeDofSizeTmp[iSub] = found[iSub].count() + numDirec * numWaves * numCS;
      else
        edgeDofSizeTmp[iSub] = numDirec * numWaves * numCS;
    }
    total += edgeDofSizeTmp[iSub];
  }

  // Form Q matrix
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    if(scomm->subNums[iSub] == subNumber) continue;
    int sOffset = sharedNodes.offset(iSub);

    // compute the wave numbers for EdgeWs augmentation
    double k_pSolid = 0.0, k_sSolid = 0.0, k_pShell = 0.0, k_sShell = 0.0, k_pFluid = 0.0;
    //if(!isFluidSub && (numDirec > 0)) {
    if(numDirec > 0) { // support for multiple fluids
      if(solInfo().solvercntl->fetiInfo.waveMethod == FetiInfo::uniform) {
        k_pSolid = k_pShell = k_p;
        k_sSolid = k_sShell = k_s;
        k_pFluid = k_f;
      }
      if(solInfo().solvercntl->fetiInfo.waveMethod == FetiInfo::averageK) {
        k_pSolid = k_pShell = neighbK_p[iSub];
        k_sSolid = k_sShell = neighbK_s[iSub];
        k_pFluid = neighbK_f[iSub];
      }
      else if(solInfo().solvercntl->fetiInfo.waveMethod == FetiInfo::averageMat) {
        double omega2 = geoSource->shiftVal();
        double lambda = (neighbPrat[iSub]*neighbYmod[iSub])/
                        (1.0+neighbPrat[iSub])/(1.0-2.0*neighbPrat[iSub]);
        double mu = neighbYmod[iSub]/2.0/(1.0+neighbPrat[iSub]);
        k_pSolid = sqrt(omega2 * neighbDens[iSub])/sqrt(lambda + 2*mu);
        k_sSolid = sqrt(omega2 * neighbDens[iSub])/sqrt(mu);
        k_pFluid = sqrt(omega2)/neighbSspe[iSub]; // should use bulk modulus and density?
        if(neighbThih[iSub] > 0.0) {
          double di = neighbYmod[iSub]*neighbThih[iSub]*neighbThih[iSub]*neighbThih[iSub]/
                      (12.0*(1.0-neighbPrat[iSub]*neighbPrat[iSub]));
          double beta4 = omega2*neighbDens[iSub]*neighbThih[iSub]/di;
          k_pShell = k_sShell = sqrt(sqrt(beta4));
        }
        else {
         k_pShell = k_sShell = 0.0;
        }
      }
    }
    if(salinasFlag && isFluidSub) { cerr << "salinas check\n"; k_pFluid = neighbK_p[iSub]; }

    // Find the center of the edge/face for EdgeGs augmentation
    double xc = 0, yc = 0, zc = 0;
    if(numR > 0) {
      for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
        xc += nodes[sharedNodes[iSub][iNode]]->x;
        yc += nodes[sharedNodes[iSub][iNode]]->y;
        zc += nodes[sharedNodes[iSub][iNode]]->z;
      }
      xc /= sharedNodes.num(iSub);
      yc /= sharedNodes.num(iSub);
      zc /= sharedNodes.num(iSub);
    }

    // loop over all the nodes on the interface and fill Q
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      if((isFluidSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) ||
         (!isFluidSub && boundaryDOFs[iSub][iNode].contains(DofSet::Helm))) continue;

      double x = nodes[sharedNodes[iSub][iNode]]->x;
      double y = nodes[sharedNodes[iSub][iNode]]->y;
      double z = nodes[sharedNodes[iSub][iNode]]->z;

      if(numR > 0) {
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Helm))
          Q[numdofperNode*(iNode+sOffset) + 0] = 1.0;
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
          Q[numdofperNode*(iNode+sOffset) + 0] = 1;
          if(found[iSub].contains(DofSet::Yrot) ) {
            Q[numdofperNode*(iNode+sOffset+4*numInterfNodes) + 0] = (z-zc);
          }
          if(found[iSub].contains(DofSet::Zrot) ) {
            Q[numdofperNode*(iNode+sOffset+5*numInterfNodes) + 0] = -(y-yc);
          }
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
          Q[numdofperNode*(iNode+sOffset+numInterfNodes) + 1] = 1;
          if(found[iSub].contains(DofSet::Xrot) ) {
            Q[numdofperNode*(iNode+sOffset+3*numInterfNodes) + 1] = -(z-zc);
          }
          if(found[iSub].contains(DofSet::Zrot) ) {
            Q[numdofperNode*(iNode+sOffset+5*numInterfNodes) + 1] = (x-xc);
          }
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
          Q[numdofperNode*(iNode+sOffset+2*numInterfNodes) + 2] = 1;
          if(found[iSub].contains(DofSet::Xrot) ) {
            Q[numdofperNode*(iNode+sOffset+3*numInterfNodes) + 2] = (y-yc);
          }
          if(found[iSub].contains(DofSet::Yrot) ) {
            Q[numdofperNode*(iNode+sOffset+4*numInterfNodes) + 2] = -(x-xc);
          }
        }
        if(numdofperNode==6) {
          if((boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) && found[iSub].contains(DofSet::Xrot))
            Q[numdofperNode*(iNode+sOffset+3*numInterfNodes) + 3] = 1;
          if((boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) && found[iSub].contains(DofSet::Yrot))
            Q[numdofperNode*(iNode+sOffset+4*numInterfNodes) + 4] = 1;
          if((boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) && found[iSub].contains(DofSet::Zrot))
            Q[numdofperNode*(iNode+sOffset+5*numInterfNodes) + 5] = 1;
        }
      }

      if(boundaryDOFs[iSub][iNode].containsAllDisp(spaceDim) ||
         boundaryDOFs[iSub][iNode].contains(DofSet::Helm))
        for (int iDir=0; iDir<numDirec; iDir++) {
          double k_p, k_s;
          if(boundaryDOFs[iSub][iNode].containsAnyRot()) {
            k_p = k_pShell;
            k_s = k_sShell;
          }
          else {
            k_p = k_pSolid;
            k_s = k_sSolid;
          }
          double cosValP, cosValS, sinValP, sinValS, cosValH, sinValH;
          double ddotx = d_x[iDir]*x + d_y[iDir]*y + d_z[iDir]*z;
          cosValP = cos(k_p * ddotx);
          cosValS = cos(k_s * ddotx);
          sinValP = sin(k_p * ddotx);
          sinValS = sin(k_s * ddotx);
          sinValH = sin(k_pFluid * ddotx);
          cosValH = cos(k_pFluid * ddotx);
          double *cosVal = new double[numWaves];
          double *sinVal = new double[numWaves];
          for(int iW=0; iW<numWaves; iW++) {
            if(iW == 0) {
              cosVal[iW] = cosValP;
              sinVal[iW] = sinValP;
            } else {
              cosVal[iW] = cosValS;
              sinVal[iW] = sinValS;
            }
            for(int iCS=0; iCS<numCS; iCS++) {
              int waveOffset = (numR + iDir*numWaves*numCS+numCS*iW+iCS)*numInterfNodes+sOffset;
              if(iCS == 0) {
                if(numdofperNode == 1) Q[numdofperNode*(waveOffset+iNode)+0] = cosValH;
                else {
                  if(spaceDim == 3) {
                    Q[numdofperNode*(waveOffset+iNode)+0] = cosVal[iW]*wDir_x[iDir*numWaves+iW];
                    Q[numdofperNode*(waveOffset+iNode)+1] = cosVal[iW]*wDir_y[iDir*numWaves+iW];
                    Q[numdofperNode*(waveOffset+iNode)+2] = cosVal[iW]*wDir_z[iDir*numWaves+iW];
                  }
                  else {
                    if(iW == 0) {
                      Q[numdofperNode*(waveOffset+iNode)+0] = cosVal[iW]*d_x[iDir];
                      Q[numdofperNode*(waveOffset+iNode)+1] = cosVal[iW]*d_y[iDir];
                    }
                    else {
                      Q[numdofperNode*(waveOffset+iNode)+0] = cosVal[iW]*t_x[iDir];
                      Q[numdofperNode*(waveOffset+iNode)+1] = cosVal[iW]*t_y[iDir];
                    }
                  }
                }
              }
              if(iCS == 1 ) {
                if(numdofperNode == 1) Q[numdofperNode*(waveOffset+iNode)+0] = sinValH;
                else {
                  if(spaceDim == 3) {
                    Q[numdofperNode*(waveOffset+iNode)+0] = sinVal[iW]*wDir_x[iDir*numWaves+iW];
                    Q[numdofperNode*(waveOffset+iNode)+1] = sinVal[iW]*wDir_y[iDir*numWaves+iW];
                    Q[numdofperNode*(waveOffset+iNode)+2] = sinVal[iW]*wDir_z[iDir*numWaves+iW];
                  }
                  else {
                    if(iW ==0) {
                      Q[numdofperNode*(waveOffset+iNode)+0] = sinVal[iW]*d_x[iDir];
                      Q[numdofperNode*(waveOffset+iNode)+1] = sinVal[iW]*d_y[iDir];
                    }
                    else {
                      Q[numdofperNode*(waveOffset+iNode)+0] = sinVal[iW]*t_x[iDir];
                      Q[numdofperNode*(waveOffset+iNode)+1] = sinVal[iW]*t_y[iDir];
                    }
                  }
                }
              }
            }
          }
          delete [] cosVal;
          delete [] sinVal;
        }
    }
  }

  GramSchmidt(Q, isUsed, numdofperNode, nQPerNeighb); // see Driver.d/BaseSub.C

  int *nQAdd = (int *) dbg_alloca(sizeof(int)*(scomm->numNeighb));
  int *nQAddWaves = (int *) dbg_alloca(sizeof(int)*(scomm->numNeighb));

  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    nQAdd[iSub] = 0;
    nQAddWaves[iSub] = 0;
  }

  // Count the number of Qs that each neighbor will have, and reset total, edgeDofSizeTmp[]
  total = 0;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    if(scomm->subNums[iSub] == subNumber) continue;
    for(int iQ = 0; iQ < nQPerNeighb; ++iQ)
      if(isUsed[iQ+iSub*nQPerNeighb]) {
        nQAdd[iSub]++;
        if ( iQ >= numR )
          nQAddWaves[iSub]++;
      }
    edgeDofSizeTmp[iSub] = nQAdd[iSub];
    total += edgeDofSizeTmp[iSub];
  }

  int *HelmCount = new int[total];
  int *xyzCount = new int[total];

  int oldTot = total;
  for(i=0; i<total; ++i) {
    HelmCount[i] = 0;
    xyzCount[i] = 0;
  }

  total = 0;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    if(edgeDofSizeTmp[iSub]==0) continue;
    int hCount = 0, xCount=0, yCount=0, zCount=0, xrCount=0, yrCount=0, zrCount=0;
    int waveCount=0;
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      if((isFluidSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) ||
         (!isFluidSub && boundaryDOFs[iSub][iNode].contains(DofSet::Helm))) continue;

      if(boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) hCount++;
      if(numR > 0) {
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) xCount++;
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) yCount++;
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) zCount++;
        if(numdofperNode==6) {
          if(boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) xrCount++;
          if(boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) yrCount++;
          if(boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) zrCount++;
        }
      }
      if(boundaryDOFs[iSub][iNode].containsAllDisp(spaceDim)) waveCount++;
    }

    int edgeLength = 0;
    if(numR > 0) {
      if(desired.contains(DofSet::Helm))
        if(found[iSub].contains(DofSet::Helm) && isUsed[0+iSub*nQPerNeighb] )
          edgeLength += HelmCount[total++] = hCount;
      if(desired.contains(DofSet::XYZdisp)) {
        if(found[iSub].contains(DofSet::Xdisp) && isUsed[0+iSub*nQPerNeighb] )
          edgeLength += xyzCount[total++] = xCount;
        if(found[iSub].contains(DofSet::Ydisp) && isUsed[1+iSub*nQPerNeighb] )
          edgeLength += xyzCount[total++] = yCount;
        if(found[iSub].contains(DofSet::Zdisp) && isUsed[2+iSub*nQPerNeighb] )
          edgeLength += xyzCount[total++] = zCount;
      }
      if(desired.contains(DofSet::XYZrot)) {
        if(found[iSub].contains(DofSet::Xrot) && isUsed[3+iSub*nQPerNeighb] )
          edgeLength += xyzCount[total++] = xrCount + yCount + zCount;
        if(found[iSub].contains(DofSet::Yrot) && isUsed[4+iSub*nQPerNeighb] )
          edgeLength += xyzCount[total++] = yrCount + xCount + zCount;
        if(found[iSub].contains(DofSet::Zrot) && isUsed[5+iSub*nQPerNeighb] )
          edgeLength += xyzCount[total++] = zrCount + xCount + yCount;
      }
    }
    if(desired.contains(DofSet::XYZdisp)) {
      if(waveCount > 0) {
        if(found[iSub].containsAllDisp(spaceDim)) {
          for(i=0; i< numDirec * numWaves * numCS; i++)
            if(isUsed[numR +i+iSub*nQPerNeighb])
              edgeLength += xyzCount[total++] = numWaves*waveCount;
        }
      }
    }
    if(desired.contains(DofSet::Helm))
      if((hCount > 0) && (found[iSub].contains(DofSet::Helm)))
        for(i=0; i< numDirec * numWaves * numCS; i++)
          if(isUsed[numR +i+iSub*nQPerNeighb])
            edgeLength += HelmCount[total++] = 1*hCount;

    totalLengthGrc += edgeLength;
  }
  if(total != oldTot) fprintf(stderr, "Non match %d %d\n", total,oldTot);

  int *HelmList = new int[totalLengthGrc];
  Scalar *HelmCoefs = new Scalar[totalLengthGrc];
  int *xyzList = new int[totalLengthGrc];
  Scalar *xyzCoefs = new Scalar[totalLengthGrc];

  for(i=0; i<totalLengthGrc; ++i) {
    HelmList[i] = 0;
    HelmCoefs[i] = 0.0;
    xyzList[i] = 0;
    xyzCoefs[i] = 0.0;
  }

  int hOffset=0, xOffset=0, yOffset=0, zOffset=0, waveOffset=0, xrOffset=0, yrOffset=0, zrOffset=0;
  total = 0;
  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    if(edgeDofSizeTmp[iSub]==0) continue;

    int index = 0;
    int off;
    if((desired.contains(DofSet::XYZdisp)) || (desired.contains(DofSet::Helm))) {
      if(numR > 0) {
        yOffset = xOffset +
                  ((found[iSub].contains(DofSet::Xdisp) && isUsed[0+iSub*nQPerNeighb]) ? xyzCount[total++] : 0);
        zOffset = yOffset +
                  ((found[iSub].contains(DofSet::Ydisp) && isUsed[1+iSub*nQPerNeighb]) ? xyzCount[total++] : 0);
        xrOffset= zOffset +
                  ((found[iSub].contains(DofSet::Zdisp) && isUsed[2+iSub*nQPerNeighb]) ? xyzCount[total++] : 0);
        yrOffset = xrOffset +
                   ((found[iSub].contains(DofSet::Xrot)  && isUsed[3+iSub*nQPerNeighb]) ? xyzCount[total++] : 0);
        zrOffset = yrOffset +
                   ((found[iSub].contains(DofSet::Yrot)  && isUsed[4+iSub*nQPerNeighb]) ? xyzCount[total++] : 0);
        if(numdofperNode == 1) {
          waveOffset = hOffset +
                       ((found[iSub].contains(DofSet::Helm) && isUsed[0+iSub*nQPerNeighb]) ? HelmCount[total++] : 0);
        }
        else {
          waveOffset = zrOffset +
                       ((found[iSub].contains(DofSet::Zrot)  && isUsed[5+iSub*nQPerNeighb]) ? xyzCount[total++] : 0);
        }
      }
      if(nQAddWaves[iSub] > 0) {
        if(numdofperNode == 1)
          index = HelmCount[(total+=nQAddWaves[iSub])-1];
        else
          index = xyzCount[(total+=nQAddWaves[iSub])-1];
      }
      else index = 0;
      off = waveOffset + ((nQAddWaves[iSub] > 0) ? nQAddWaves[iSub]*index : 0);
    }

    double sign = (scomm->subNums[iSub] < subNumber) ? 1.0 : -1.0;
    int sOffset = sharedNodes.offset(iSub);
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      if((isFluidSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) ||
         (!isFluidSub && boundaryDOFs[iSub][iNode].contains(DofSet::Helm))) continue;

      int qOff = sOffset + iNode;
      int hDof = -1, xDof = -1, yDof = -1, zDof = -1;
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) {
        hDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Helm);
        if(numR > 0 ) {
          if(found[iSub].contains(DofSet::Helm) && isUsed[0+iSub*nQPerNeighb]) {
            HelmList[hOffset] = hDof;
            HelmCoefs[hOffset++] = sign * Q[numdofperNode*qOff+0];
          }
        }
      }
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
        xDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Xdisp);
        if(numR > 0) {
          if(found[iSub].contains(DofSet::Xdisp) && isUsed[0+iSub*nQPerNeighb]) {
            xyzList[xOffset] = xDof;
            xyzCoefs[xOffset++] = sign*Q[numdofperNode*qOff];
          }
          if(found[iSub].contains(DofSet::Yrot)  && isUsed[4+iSub*nQPerNeighb]) {
            xyzList[yrOffset] = xDof;
            xyzCoefs[yrOffset++] = sign* Q[numdofperNode*(qOff+4*numInterfNodes)];
          }
          if(found[iSub].contains(DofSet::Zrot)  && isUsed[5+iSub*nQPerNeighb]) {
            xyzList[zrOffset] = xDof;
            xyzCoefs[zrOffset++] = sign* Q[numdofperNode*(qOff+5*numInterfNodes)];
          }
        }
      }
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
        yDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Ydisp);
        if(numR > 0) {
          if(found[iSub].contains(DofSet::Ydisp) && isUsed[1+iSub*nQPerNeighb]) {
            xyzList[yOffset] = yDof;
            xyzCoefs[yOffset++] = sign*Q[numdofperNode*(qOff+numInterfNodes)+1];
          }
          if(found[iSub].contains(DofSet::Xrot)  && isUsed[3+iSub*nQPerNeighb]) {
            xyzList[xrOffset] = yDof;
            xyzCoefs[xrOffset++] = sign* Q[numdofperNode*(qOff+3*numInterfNodes)+1];
          }
          if(found[iSub].contains(DofSet::Zrot)  && isUsed[5+iSub*nQPerNeighb]) {
            xyzList[zrOffset] = yDof;
            xyzCoefs[zrOffset++] = sign* Q[numdofperNode*(qOff+5*numInterfNodes)+1];
          }
        }
      }
      if(boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
        zDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Zdisp);
        if(numR > 0) {
          if(found[iSub].contains(DofSet::Zdisp) && isUsed[2+iSub*nQPerNeighb]) {
            xyzList[zOffset] = zDof;
            xyzCoefs[zOffset++] = sign*Q[numdofperNode*(qOff+2*numInterfNodes)+2];
          }
          if(found[iSub].contains(DofSet::Xrot)  && isUsed[3+iSub*nQPerNeighb]) {
            xyzList[xrOffset] = zDof;
            xyzCoefs[xrOffset++] = sign* Q[numdofperNode*(qOff+3*numInterfNodes)+2];
          }
          if(found[iSub].contains(DofSet::Yrot)  && isUsed[4+iSub*nQPerNeighb]) {
            xyzList[yrOffset] = zDof;
            xyzCoefs[yrOffset++] = sign* Q[numdofperNode*(qOff+4*numInterfNodes)+2];
          }
        }
      }
      if((numR > 3) && (numdofperNode==6)) {
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) {
          int xrDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Xrot);
          if(found[iSub].contains(DofSet::Xrot) && isUsed[3+iSub*nQPerNeighb]) {
            xyzList[xrOffset] = xrDof;
            xyzCoefs[xrOffset++] = sign*Q[numdofperNode*(qOff+3*numInterfNodes)+3];
          }
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) {
          int yrDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Yrot);
          if(found[iSub].contains(DofSet::Yrot) && isUsed[4+iSub*nQPerNeighb]) {
            xyzList[yrOffset] = yrDof;
            xyzCoefs[yrOffset++] = sign*Q[numdofperNode*(qOff+4*numInterfNodes)+4];
          }
        }
        if(boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) {
          int zrDof = cc_dsa->locate(sharedNodes[iSub][iNode],DofSet::Zrot);
          if(found[iSub].contains(DofSet::Zrot) && isUsed[5+iSub*nQPerNeighb]) {
            xyzList[zrOffset] = zrDof;
            xyzCoefs[zrOffset++] = sign*Q[numdofperNode*(qOff+5*numInterfNodes)+5];
          }
        }
      }
      if(nQAddWaves[iSub] > 0) {
      if(boundaryDOFs[iSub][iNode].containsAllDisp(spaceDim) ||
         boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) {
        int UsedWaves = 0;
        for(int iDir=0; iDir<numDirec; iDir++)
          for(int iW=0; iW<numWaves; iW++)
            for(int iCS=0; iCS<numCS; iCS++) {
              if(isUsed[(numR+iDir*numWaves*numCS+iW*numCS+iCS)+iSub*nQPerNeighb] == false)
                continue;
              qOff = (numR+iDir*numWaves*numCS+iW*numCS+iCS)*numInterfNodes+sOffset+iNode;
              if(UsedWaves == 0) {
                if(numdofperNode == 1) {
                  if(hDof > -1) {
                    HelmList[waveOffset]= hDof;
                    HelmCoefs[waveOffset++]= sign*Q[numdofperNode*qOff+0];
                  }
                }
                else {
                  if(xDof > -1) {
                    xyzList[waveOffset]= xDof;
                    xyzCoefs[waveOffset++]= sign*Q[numdofperNode*qOff+0];
                  }
                  if(yDof > -1) {
                    xyzList[waveOffset]= yDof;
                    xyzCoefs[waveOffset++]= sign*Q[numdofperNode*qOff+1];
                  }
                  if(zDof > -1) {
                    xyzList[waveOffset]= zDof;
                    xyzCoefs[waveOffset++]= sign*Q[numdofperNode*qOff+2];
                  }
                }
              }
              else {
                if(numdofperNode == 1) {
                  if(hDof > -1) {
                    HelmList[UsedWaves*index+waveOffset-1]  = hDof;
                    HelmCoefs[UsedWaves*index+waveOffset-1] = sign*Q[numdofperNode*qOff+0];
                  }
                }
                else {
                  if(xDof > -1) {
                    xyzList[UsedWaves*index+waveOffset-3]  = xDof;
                    xyzCoefs[UsedWaves*index+waveOffset-3] = sign*Q[numdofperNode*qOff+0];
                  }
                  if(yDof > -1) {
                    xyzList[UsedWaves*index+waveOffset-2]  = yDof;
                    xyzCoefs[UsedWaves*index+waveOffset-2] = sign*Q[numdofperNode*qOff+1];
                  }
                  if(zDof > -1) {
                    xyzList[UsedWaves*index+waveOffset-1]  = zDof;
                    xyzCoefs[UsedWaves*index+waveOffset-1] = sign*Q[numdofperNode*qOff+2];
                  }
                }
              }
              UsedWaves++;
            }
        if(UsedWaves != nQAddWaves[iSub])
          cerr << " Something is wrong for the number of used waves " << endl;
      }
      }
    }
    if(numR > 0)
      if(numdofperNode == 1)
        hOffset = off;
      else
        xOffset = off;
    else
      waveOffset = off;
  }
  if(oldTot != total)
    fprintf(stderr, " *** ERROR: total is incorrect %d %d\n", oldTot, total);

  if(numdofperNode == 1) {
    Grc = new GenCuCSparse<Scalar>(total, cc_dsa->size(), HelmCount, HelmList, HelmCoefs);
    delete [] xyzList;
    delete [] xyzCoefs;
  }
  else {
    Grc = new GenCuCSparse<Scalar>(total, cc_dsa->size(), xyzCount, xyzList, xyzCoefs);
    delete [] HelmList;
    delete [] HelmCoefs;
  }

  int ii = (isFluidSub) ? 1 : 0;
  if(edgeQindex[ii] == -1)
    edgeQindex[ii] = Src->addSparseMatrix(Grc);  // store index for possible rebuild (multiple LHS freq sweep)
  else
     Src->setSparseMatrix(edgeQindex[ii], Grc);

  for(i=0; i<scomm->numNeighb; ++i) edgeDofSize[i] += edgeDofSizeTmp[i];

  if(wDir_x) delete [] wDir_x;
  if(wDir_y) delete [] wDir_y;
  if(wDir_z) delete [] wDir_z;
  delete [] HelmCount;
  delete [] xyzCount;
  delete [] found;
  if(Q) delete [] Q;
}

template<class Scalar>
void
GenSubDomain<Scalar>::setMpcDiagCommSize(FSCommPattern<Scalar> *mpcDiagPat)
{
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc,i);
    int len = 0;
    if(subNumber != neighb) {
      for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j)
        len += mpc[scomm->mpcNb(i,j)]->gsize;
    }
    mpcDiagPat->setLen(subNumber, neighb, len);
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendMpcDiag(FSCommPattern<Scalar> *mpcDiagPat)
{
  for(int i = 0; i < numMPC; ++i) mpc[i]->initKsum();

  // Get the trace of the subdomain interfaces
  int iNeighb, iDof, j;
  for(iNeighb = 0; iNeighb < scomm->numT(SComm::mpc); ++iNeighb) {
    int neighb = scomm->neighbT(SComm::mpc, iNeighb);
    FSSubRecInfo<Scalar> sInfo = mpcDiagPat->getSendBuffer(subNumber, neighb);
    int nOff = 0;
    for(iDof = 0; iDof < scomm->lenT(SComm::mpc,iNeighb); ++iDof) {
      int locMpcNb = scomm->mpcNb(iNeighb,iDof);
      if(subNumber != neighb)
        for(j = 0; j < mpc[locMpcNb]->gsize; ++j) sInfo.data[nOff+j] = 0.0;
      for(j = 0; j < mpc[locMpcNb]->nterms; ++j) {
        int c_dof = mpc[locMpcNb]->terms[j].cdof;
        if(c_dof > -1) {
          int b_dof = invBoundMap[c_dof];
          mpc[locMpcNb]->k[j] = (Kbb) ? Kbb->diag(b_dof) : 1.0;
          if(ScalarTypes::norm(mpc[locMpcNb]->k[j]) < 1.0e-12) cerr << " *** WARNING: Kbb diagonal < 1.0e-12 \n";
          if(subNumber != neighb)
            sInfo.data[nOff+mpc[locMpcNb]->gi[j]] = mpc[locMpcNb]->k[j];
          mpc[locMpcNb]->ksum[j] = mpc[locMpcNb]->k[j];
        }
        else {
          if(subNumber != neighb)
            sInfo.data[nOff+mpc[locMpcNb]->gi[j]] = 0.0;
          mpc[locMpcNb]->k[j] = 0.0;
          mpc[locMpcNb]->ksum[j] = 0.0;
        }
      }
      nOff += mpc[locMpcNb]->gsize;
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::collectMpcDiag(FSCommPattern<Scalar> *mpcDiagPat)
{
 int iNeighb, iDof, j;
 // 1) gets & sum the stiffness contribution from neighbourg subds
 for(iNeighb = 0; iNeighb < scomm->numT(SComm::mpc); ++iNeighb) {
   // ksum already contains its own contribution (i.e. initialized in sendMpcDiag)
   // -> loop only on neighboring subds
   int neighb = scomm->neighbT(SComm::mpc,iNeighb);
   if(subNumber != neighb) {
     FSSubRecInfo<Scalar> rInfo = mpcDiagPat->recData(neighb, subNumber);
     int nOff = 0;
     for(iDof = 0; iDof < scomm->lenT(SComm::mpc,iNeighb); ++iDof) {
       int locMpcNb = scomm->mpcNb(iNeighb,iDof);
       for(j = 0; j < mpc[locMpcNb]->nterms; ++j)  {
         mpc[locMpcNb]->ksum[j] += rInfo.data[nOff+mpc[locMpcNb]->gi[j]];
       }
       nOff += mpc[locMpcNb]->gsize;
     }
   }
 }

 // 2) adjust dual mpcs using subdomain multiplicity
 //    c(i) -> c(i).k(i,i)/sum[k(j,j)]
 if(solInfo().getFetiInfo().mpc_scaling == FetiInfo::kscaling) {
   int iMPC, i;
#ifdef DEBUG_MPC
   cerr << "before k scaling: \n";
   for(iMPC = 0; iMPC < numMPC; ++iMPC) mpc[iMPC]->print();
#endif
   for(iMPC = 0; iMPC < numMPC; ++iMPC) {
     if(mpc[iMPC]->type == 2) continue; // bmpc
     for(i = 0; i < mpc[iMPC]->nterms; ++i) {
       if(ScalarTypes::norm(mpc[iMPC]->ksum[i]) < 1.0e-12) {
         //cerr << " *** WARNING: ksum = " << mpc[iMPC]->ksum[i] << ", cdof = " << mpc[iMPC]->terms[i].cdof << ", coef = " << mpc[iMPC]->terms[i].coef << endl;
         mpc[iMPC]->ksum[i] = 1.0;
         mpc[iMPC]->k[i] = 1.0;
       }
       mpc[iMPC]->terms[i].coef *= (mpc[iMPC]->k[i]/mpc[iMPC]->ksum[i]);
     }
   }
#ifdef DEBUG_MPC
   cerr << "after k scaling: \n";
   for(iMPC = 0; iMPC < numMPC; ++iMPC) mpc[iMPC]->print();
#endif
 }
}

template<class Scalar>
Scalar
GenSubDomain<Scalar>::getMpcRhs(int glMPCnum)
{
  return mpc[globalToLocalMPC[glMPCnum]]->rhs;
}

template<class Scalar>
Scalar
GenSubDomain<Scalar>::getMpcRhs_primal(int glMPCnum)
{
  return mpc_primal[globalToLocalMPC_primal[glMPCnum]]->rhs;
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractMPCs(int glNumMPC, ResizeArray<LMPCons *> &lmpc)
{
  // PASS 1 count number of local mpcs
  numMPC = 0;
  for(int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
    if(lmpc[iMPC]->isPrimalMPC()) continue;
    for(int i = 0; i < lmpc[iMPC]->nterms; ++i) {
      if(globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
        numMPC++;
        break;
      }
    }
  }
  if(numMPC == 0) return;

  mpc = new SubLMPCons<Scalar> * [numMPC]; // use specific class for subdomain lmpcs
  for(int iMPC = 0; iMPC < numMPC; ++iMPC) mpc[iMPC] = 0;
  localToGlobalMPC = new int[numMPC];

  // PASS 2: Get the mpc values
  numMPC = 0;
  int count = 0;
  for(int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
    if(lmpc[iMPC]->isPrimalMPC()) continue;
    int used = 0;
    for(int i = 0; i < lmpc[iMPC]->nterms; ++i) {
      if(globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
        //if(c_dsa->locate((lmpc[iMPC]->terms)[i].nnum, (1 << (lmpc[iMPC]->terms)[i].dofnum)) < 0) continue; // XXXX
        if(mpc[numMPC] == 0) {
          Scalar rhs = lmpc[iMPC]->template getRhs<Scalar>();
          GenLMPCTerm<Scalar> term0 = lmpc[iMPC]->template getTerm<Scalar>(i);
          if(lmpc[iMPC]->isBoundaryMPC()) {
            if(lmpc[iMPC]->psub == subNumber) term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? 0.707106781 :*/ 1.0;
            else if(lmpc[iMPC]->nsub == subNumber) term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? -0.707106781 :*/ -1.0;
          }
          mpc[numMPC] = new SubLMPCons<Scalar>(lmpc[iMPC]->lmpcnum, rhs, term0, lmpc[iMPC]->nterms, i);
          mpc[numMPC]->type = lmpc[iMPC]->type; // this is to be phased out
          mpc[numMPC]->setType(lmpc[iMPC]->getType());
          mpc[numMPC]->setSource(lmpc[iMPC]->getSource());
          mpc[numMPC]->id = lmpc[iMPC]->id;
          used = 1;
        }
        else {
          GenLMPCTerm<Scalar> term = lmpc[iMPC]->template getTerm<Scalar>(i);
          mpc[numMPC]->addterm(term, i);
        }
      }
    }
    if(used==1) {
      localToGlobalMPC[numMPC] = count;
      numMPC++;
    }
    count++;
  }
  globalToLocalMPC.initialize(numMPC, localToGlobalMPC);
  //globalToLocalMPC.print();

#ifdef DEBUG_MPC
  cerr << "DUAL MPCs: \n";
  for(int iMPC = 0; iMPC < numMPC; ++iMPC) mpc[iMPC]->print();
#endif
}

template<class Scalar>
void
GenSubDomain<Scalar>::printLMPC()
{
  cerr << "sub = " << subNumber << ", numMPC = " << numMPC << endl;
  for(int iMPC=0; iMPC<numMPC; ++iMPC) {
    cerr << "local mpc # " << iMPC << ", global mpc # " << localToGlobalMPC[iMPC] << endl;
    mpc[iMPC]->print();
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractMPCs_primal(int glNumMPC, ResizeArray<LMPCons *> &lmpc)
{
  // PASS 1 count number of local mpcs
  numMPC_primal = 0;
  for(int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
    if(!lmpc[iMPC]->isPrimalMPC()) continue;
    for(int i = 0; i < lmpc[iMPC]->nterms; ++i) {
      if(globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
        numMPC_primal++;
        break;
      }
    }
  }
  if(numMPC_primal == 0) return;

  mpc_primal = new SubLMPCons<Scalar> * [numMPC_primal]; // use specific class for subdomain lmpcs
  for(int iMPC = 0; iMPC < numMPC_primal; ++iMPC) mpc_primal[iMPC] = 0;
  localToGlobalMPC_primal = new int[numMPC_primal];

  // PASS 2: Get the mpc values
  numMPC_primal = 0;
  int count = 0;
  for(int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
    if(!lmpc[iMPC]->isPrimalMPC()) continue;
    int used = 0;
    for(int i = 0; i < lmpc[iMPC]->nterms; ++i) {
      if(globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
        if(mpc_primal[numMPC_primal] == 0) {
          Scalar rhs = lmpc[iMPC]->template getRhs<Scalar>();
          GenLMPCTerm<Scalar> term0 = lmpc[iMPC]->template getTerm<Scalar>(i);
          if(lmpc[iMPC]->isBoundaryMPC()) {
            if(lmpc[iMPC]->psub == subNumber) term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? 0.707106781 :*/ 1.0;
            else if(lmpc[iMPC]->nsub == subNumber) term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? -0.707106781 :*/ -1.0;
          }
          mpc_primal[numMPC_primal] = new SubLMPCons<Scalar>(numMPC_primal, rhs, term0, lmpc[iMPC]->nterms, i);
          mpc_primal[numMPC_primal]->type = lmpc[iMPC]->type;
          mpc_primal[numMPC_primal]->setType(lmpc[iMPC]->getType());
          mpc_primal[numMPC_primal]->setSource(lmpc[iMPC]->getSource());
          mpc_primal[numMPC_primal]->id = lmpc[iMPC]->id;
          used = 1;
        }
        else {
          GenLMPCTerm<Scalar> term = lmpc[iMPC]->template getTerm<Scalar>(i);
          mpc_primal[numMPC_primal]->addterm(term, i);
        }
      }
    }
    if(used==1) {
      localToGlobalMPC_primal[numMPC_primal] = count;
      numMPC_primal++;
    }
    count++;
  }
  globalToLocalMPC_primal.initialize(numMPC_primal, localToGlobalMPC_primal);
  //globalToLocalMPC_primal.print();

#ifdef DEBUG_MPC
  cerr << "PRIMAL MPCs: \n";
  for(int iMPC = 0; iMPC < numMPC_primal; ++iMPC) mpc_primal[iMPC]->print();
#endif
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeLocalMpcToDof()
{
  if(mpcToDof) delete mpcToDof; mpcToDof = 0;

  // step 1: make mpcToDof
  int size = numMPC;
  // step 1.1: find size of target: total number of coefficients involving a different dof
  int numtarget = 0;
  int i, j, jj;
  for(i = 0; i < size; i++) {
    for(j = 0; j < mpc[i]->nterms; j++) {
      int dofj = mpc[i]->terms[j].dof;
      for(jj=0; jj<j; jj++) {
        int dofjj = mpc[i]->terms[jj].dof;
        if(dofj == dofjj) break;
      }
      if((jj==j) && (dofj >= 0)) numtarget++;
    }
  }
  // step 1.2: fill target with coefficient dofs
  int *pointer = new int[size+1];
  int *target  = new int[numtarget];
  int count = 0;
  for(i = 0; i < size; i++) {
    pointer[i] = count;
    for(j = 0; j < mpc[i]->nterms; j++) {
      int dofj = mpc[i]->terms[j].dof;
      for(jj=0; jj<j; jj++) {
        int dofjj = mpc[i]->terms[jj].dof;
        if(dofj == dofjj) break;
      }
      if((jj==j) && (dofj >= 0)) {
        target[count] = dofj;
        count++;
      }
    }
  }
  pointer[i] = numtarget;
  // step 1.3: construct mpcToDof connectivity
  mpcToDof = new Connectivity(size, pointer, target);
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeLocalMpcToMpc()
{
  // step 1: make mpcToDof
  makeLocalMpcToDof();

  // step 2: make localMpcToMpc connectivity
  Connectivity *dofToMpc = mpcToDof->reverse();
  localMpcToMpc = mpcToDof->transcon(dofToMpc);
  delete dofToMpc;
}

template<class Scalar>
void
GenSubDomain<Scalar>::locateMpcDofs()
{
  for(int i = 0; i < numMPC; ++i)
    for(int k = 0; k < mpc[i]->nterms; ++k) {
      (mpc[i]->terms)[k].dof = dsa->locate((mpc[i]->terms)[k].nnum, (1 << (mpc[i]->terms)[k].dofnum));
      (mpc[i]->terms)[k].cdof = c_dsa->locate((mpc[i]->terms)[k].nnum, (1 << (mpc[i]->terms)[k].dofnum));
      (mpc[i]->terms)[k].ccdof = cc_dsa->locate((mpc[i]->terms)[k].nnum, (1 << (mpc[i]->terms)[k].dofnum));
    }
  for(int i = 0; i < numMPC_primal; ++i)
    for(int k = 0; k < mpc_primal[i]->nterms; ++k) {
      (mpc_primal[i]->terms)[k].dof = dsa->locate((mpc_primal[i]->terms)[k].nnum, (1 << (mpc_primal[i]->terms)[k].dofnum));
      (mpc_primal[i]->terms)[k].cdof = c_dsa->locate((mpc_primal[i]->terms)[k].nnum, (1 << (mpc_primal[i]->terms)[k].dofnum));
      (mpc_primal[i]->terms)[k].ccdof = cc_dsa->locate((mpc_primal[i]->terms)[k].nnum, (1 << (mpc_primal[i]->terms)[k].dofnum));
    }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver, SimpleNumberer *mpcEqNums)
{
  int i, j, k, l;
  for(i=0; i<numMPC; ++i) {
    Scalar dotii = 0.0;
    int gi = localToGlobalMPC[i];
    int renum_gi = mpcEqNums->firstdof(gi);
    if(mpc[i]->active) { CCtsolver->addone(1.0, renum_gi, renum_gi); continue; } // trick to prevent singularities in CCt when rebuit for contact
    for(k = 0; k < mpc[i]->nterms; ++k) {
      int dof = (mpc[i]->terms)[k].cdof;
      if(dof >= 0)
        dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling
    }
    CCtsolver->addone(dotii, renum_gi, renum_gi);
    for(j=0; j < localMpcToMpc->num(i); ++j) {
      Scalar dotij = 0.0;
      int lj = (*localMpcToMpc)[i][j];
      if(mpc[lj]->active) continue;
      int gj = localToGlobalMPC[lj];
      int renum_gj = mpcEqNums->firstdof(gj);
      if(renum_gj > renum_gi) {  // work with upper symmetric half
        // now find matching dof/s
        for(k = 0; k < mpc[i]->nterms; ++k) {
          int dofk = (mpc[i]->terms)[k].cdof;
          if(dofk >= 0) {
            for(l = 0; l < mpc[lj]->nterms; ++l) {
              int dofl = (mpc[lj]->terms)[l].cdof;
              if(dofk == dofl) {
                dotij += mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef / mpc[lj]->k[l]; // for mpc kscaling
              }
            }
          }
        }
        CCtsolver->addone(dotij, renum_gi, renum_gj);
      }
    }
  }
}

//HB: attempt to improve parallel efficiency of assembling the global CCt matrix
//    This method ONLY compute the subdomain contributions to the global CCt matrix,
//    and store them (values and their respective row & col position in the global CCt matrix)
//    LOCALLY. Thus this method can be executed CONCURRENTLY.
//    The assembling in the global CCt matrix is done SEQUENTIALLY (see method below this one).
template<class Scalar>
void
GenSubDomain<Scalar>::computeSubContributionToGlobalCCt(SimpleNumberer *mpcEqNums)
{
  int i, j, k, l;
  // Step 1. Determine the size of the array & allocate array
  lengthCCtData=0;
  // this is an upper estimate of the required array size
  // -> nearly 2x the required size

  // this is the exact required array size
  for(i = 0; i < numMPC; ++i){
    if(mpc[i]->active) continue;
    int gi = localToGlobalMPC[i];
    int renum_gi = mpcEqNums->firstdof(gi);
    lengthCCtData++; //diagonal term
    for(j=0; j < localMpcToMpc->num(i); ++j) {
      int lj = (*localMpcToMpc)[i][j];
      if(mpc[lj]->active) continue;
      int gj = localToGlobalMPC[lj];
      int renum_gj = mpcEqNums->firstdof(gj);
      if(renum_gj > renum_gi)   // work with upper symmetric half
         lengthCCtData++;
    }
  }
  CCtrow = new int[lengthCCtData];
  CCtcol = new int[lengthCCtData];
  CCtval = new Scalar[lengthCCtData];

  // Step 2. Fill the array
  lengthCCtData = 0; // use it as counter (at the end, it should be the exact number of contributions)
  for(i = 0; i < numMPC; ++i) {
    if(mpc[i]->active) continue;
    Scalar dotii = 0.0;
    int gi = localToGlobalMPC[i];
    int renum_gi = mpcEqNums->firstdof(gi);
    for(k = 0; k < mpc[i]->nterms; ++k) {
      int dof = (mpc[i]->terms)[k].cdof;
      if(dof >= 0)
        dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling
    }
    CCtrow[lengthCCtData] = renum_gi;
    CCtcol[lengthCCtData] = renum_gi;
    CCtval[lengthCCtData] = dotii;
    lengthCCtData++;
    for(j=0; j < localMpcToMpc->num(i); ++j) {
      Scalar dotij = 0.0;
      int lj = (*localMpcToMpc)[i][j];
      if(mpc[lj]->active) continue;
      int gj = localToGlobalMPC[lj];
      int renum_gj = mpcEqNums->firstdof(gj);
      if(renum_gj > renum_gi) {  // work with upper symmetric half
        // now find matching dof/s
        for(k = 0; k < mpc[i]->nterms; ++k) {
          int dofk = (mpc[i]->terms)[k].cdof;
          if(dofk >= 0) {
            for(l = 0; l < mpc[lj]->nterms; ++l) {
              int dofl = (mpc[lj]->terms)[l].cdof;
              if(dofk == dofl) {
                dotij += mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef / mpc[lj]->k[l]; // for mpc kscaling
              }
            }
          }
        }
        CCtrow[lengthCCtData] = renum_gi;
        CCtcol[lengthCCtData] = renum_gj;
        CCtval[lengthCCtData] = dotij;
        lengthCCtData++;
      }
    }
  }
}
// HB: this method add/assemble the locally stored contributions into the global CCt matrix.
//     MUST be called SEQUENTIALLY to avoid writting concurrently at the same memory location.
template<class Scalar>
void
GenSubDomain<Scalar>::assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver)
{
  for(int i=0; i<lengthCCtData; i++)
    CCtsolver->addone(CCtval[i], CCtrow[i], CCtcol[i]);

  // delete local contributions
  //(may need to delete this in parallel as it has been allocated in parallel ...)
  if(CCtrow) { delete [] CCtrow; CCtrow = 0; }
  if(CCtcol) { delete [] CCtcol; CCtcol = 0; }
  if(CCtval) { delete [] CCtval; CCtval = 0; }
  lengthCCtData = 0;
}

template<class Scalar>
void
GenSubDomain<Scalar>::constructLocalCCtsolver()
{
  // Step 1. initialize solver object
  SimpleNumberer *mpcEqNums = new SimpleNumberer(numMPC);
  for(int i = 0; i < numMPC; ++i) mpcEqNums->setWeight(i, 1);
  mpcEqNums->makeOffset();
  localCCtsolver = GenSolverFactory<Scalar>::getFactory()->createSolver(localMpcToGlobalMpc, mpcEqNums, *sinfo.solvercntl->fetiInfo.cct_cntl, localCCtsparse);
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleLocalCCtsolver()
{
  // Step 2. add local mpc CC^t terms to solver
  int i, j, k, l;
  for(i = 0; i<numMPC; ++i) {
    if(mpc[i]->active) continue;
    Scalar dotii = 0.0;
    for(k = 0; k < mpc[i]->nterms; ++k) {
      int dof = (mpc[i]->terms)[k].cdof;
      if(dof >= 0)
        dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling;
    }
    localCCtsolver->addone(dotii, i, i);
    for(j=0; j < localMpcToMpc->num(i); ++j) {
      Scalar dotij = 0.0;
      int lj = (*localMpcToMpc)[i][j];
      if(mpc[lj]->active) continue;
      if(lj > i) { // work with upper symmetric half only
        // now find matching dof/s
        for(k = 0; k < mpc[i]->nterms; ++k) {
          int dofk = (mpc[i]->terms)[k].cdof;
          if(dofk >= 0) {
            for(l = 0; l < mpc[lj]->nterms; ++l) {
              int dofl = (mpc[lj]->terms)[l].cdof;
              if(dofk == dofl) {
                dotij += mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef / mpc[lj]->k[l]; // for mpc kscaling
                //HB: can probably be optimized by assuming that no dupplicate term (i.e. dof) exist in a lmpc
                //    so that we can stop the l loop (break) when dofk == dofl ?? The non dupplicate assumption
                //    is or can be enforced in a preprocess step
              }
            }
          }
        }
        localCCtsolver->addone(dotij, i, lj);
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::setCCtCommSize(FSCommPattern<Scalar> *cctPat)
{
  // note: this is an upper bound, the required length will be less
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i)
    cctPat->setLen(subNumber, scomm->neighbT(SComm::mpc, i), localMpcToGlobalMpc->numConnect());
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, Connectivity *mpcToSub)
{
  int i,j,k;
  for(i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc,i);
    if(subNumber != neighb) {
      int count = 0;
      FSSubRecInfo<Scalar> sInfo = cctPat->getSendBuffer(subNumber, neighb);
      for(j = 0; j < numMPC; ++j) {
        int gj = localToGlobalMPC[j];
        if(mpcToSub->offset(gj, neighb) > -1) {
          sInfo.data[count++] = localCCtsolver->getone(j, j);
          for(k = j+1; k<numMPC; ++k) {
            int gk = localToGlobalMPC[k];
            if((mpcToSub->offset(gk, neighb) > -1) && (localMpcToGlobalMpc->offset(j, k) > -1)) {
              sInfo.data[count] = localCCtsolver->getone(j, k);
              count++;
            }
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::recNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, Connectivity *mpcToSub)
{
  int i, j, k;
  for(i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc, i);
    if(subNumber != neighb) {
      int count = 0;
      FSSubRecInfo<Scalar> rInfo = cctPat->recData(neighb, subNumber);
      for(j = 0; j < numMPC; ++j) {
        int gj = localToGlobalMPC[j];
        if(mpcToSub->offset(gj, neighb) > -1) {
          localCCtsolver->addone(rInfo.data[count++], j, j);
          for(k = j+1; k<numMPC; ++k) {
            int gk = localToGlobalMPC[k];
            if((mpcToSub->offset(gk, neighb) > -1) && (localMpcToGlobalMpc->offset(j, k) > -1)) {
              localCCtsolver->addone(rInfo.data[count], j, k);
              count++;
            }
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::factorLocalCCtsolver()
{
  localCCtsolver->factor();
  int numCCtSing = localCCtsolver->numRBM();
  if(numCCtSing > 0)
    cerr << "sub = " << subNumber << ", Number of singularities in CCt = "
         << numCCtSing << endl;
}

template<class Scalar>
void
GenSubDomain<Scalar>::zeroLocalCCtsolver()
{
  localCCtsolver->zeroAll();
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleBlockCCtsolver(int iBlock, GenSolver<Scalar> *CCtsolver, SimpleNumberer *blockMpcEqNums)
{
  // optimize by looping ONLY over the lmpc eqs contributed to this iBlock
  // Make sure that the input block (global) Id iBlock is in the range of blocks this subdomain
  // contributes to
  int i, j, k, l,p;
  for(p=0; p<blockToLocalMpc->num(iBlock); ++p) {
    i = (*blockToLocalMpc)[iBlock][p];
    if(mpc[i]->active) continue;
    int bi = (*blockToBlockMpc)[iBlock][p];
    int renum_bi = blockMpcEqNums->firstdof(bi);
    Scalar dotii = 0.0;
    for(k = 0; k < mpc[i]->nterms; ++k) {
      int dof = (mpc[i]->terms)[k].cdof;
      if(dof >= 0)
        dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling;
    }
    CCtsolver->addone(dotii, renum_bi, renum_bi);
    for(j=0; j < localMpcToMpc->num(i); ++j) {
      int lj = (*localMpcToMpc)[i][j];
       if(mpc[lj]->active) continue;
      int jb = localMpcToBlock->cOffset(lj, iBlock);
      if(jb > -1) {
        int bj = (*localMpcToBlockMpc)[lj][jb];
        int renum_bj = blockMpcEqNums->firstdof(bj);
        if(renum_bj > renum_bi) { // work with upper symmetric part only
          Scalar dotij = 0.0;
          // now find matching dof/s
          for(k = 0; k < mpc[i]->nterms; ++k) {
            int dofk = (mpc[i]->terms)[k].cdof;
            if(dofk >= 0) {
              for(l = 0; l < mpc[lj]->nterms; ++l) {
                int dofl = (mpc[lj]->terms)[l].cdof;
                if(dofk == dofl)
                  dotij += mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef / mpc[lj]->k[l]; // for mpc kscaling;
              }
            }
          }
          CCtsolver->addone(dotij, renum_bi, renum_bj);
        }
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::constraintProduct(int num_vect, const double* R[], Scalar** V, int trans)
{
  int i, n, iMPC;

  if(trans) {
    for(iMPC = 0; iMPC < numMPC; ++iMPC) {
      for(n = 0; n < num_vect; ++n) {
        for(i = 0; i < mpc[iMPC]->nterms; ++i) {
          int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
                                  (1 << mpc[iMPC]->terms[i].dofnum));
          if(dof < 0) continue;
          V[n][dof] += mpc[iMPC]->terms[i].coef * R[n][dof];
        }
      }
    }
  }
  else {
    for(iMPC = 0; iMPC < numMPC; ++iMPC) {
      for(i = 0; i < mpc[iMPC]->nterms; ++i) {
        int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
                                (1 << mpc[iMPC]->terms[i].dofnum));
        if(dof < 0) continue;
        for(n = 0; n < num_vect; ++n)
          V[iMPC][n] += mpc[iMPC]->terms[i].coef * R[iMPC][dof];
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::addConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &f)
{
  bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

  vector<double>::iterator it2 = lambda.begin();
  for(int l = 0; l < scomm->lenT(SComm::mpc); ++l) {
    int i = scomm->mpcNb(l);
    if(!mpcFlag[i]) continue;
    if(mpc[i]->getSource() == mpc::ContactSurfaces) { // contact
      std::map<std::pair<int,int>, double>::iterator it1 = mu.find(mpc[i]->id);
      if(it1 != mu.end()) {
        for(int j = 0; j < mpc[i]->nterms; ++j) {
          int dof = c_dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
          if(dof < 0) continue;
          f[dof] += mpc[i]->terms[j].coef*it1->second;
        }
      }
    }
    else {
      if(it2 != lambda.end()) {
        for(int j = 0; j < mpc[i]->nterms; ++j) {
          int dof = c_dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
          if(dof < 0) continue;
          f[dof] += mpc[i]->terms[j].coef*(*it2);
        }
        it2++;
      }
    }
    mpcFlag[i] = false;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::addCConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &fc, double s)
{
  bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

  vector<double>::iterator it2 = lambda.begin();
  for(int l = 0; l < scomm->lenT(SComm::mpc); ++l) {
    int i = scomm->mpcNb(l);
    if(!mpcFlag[i]) continue;
    if(mpc[i]->getSource() == mpc::ContactSurfaces) { // contact
      std::map<std::pair<int,int>, double>::iterator it1 = mu.find(mpc[i]->id);
      if(it1 != mu.end()) {
        for(int j = 0; j < mpc[i]->nterms; ++j) {
          int dof = dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
          if(dof < 0) continue;
          int cdof = c_dsa->invRCN(dof);
          if(cdof >= 0)
            fc[cdof] += s*mpc[i]->terms[j].coef*it1->second;
        }
      }
    }
    else {
      if(it2 != lambda.end()) {
        for(int j = 0; j < mpc[i]->nterms; ++j) {
          int dof = dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
          if(dof < 0) continue;
          int cdof = c_dsa->invRCN(dof);
          if(cdof >= 0)
            fc[cdof] += s*mpc[i]->terms[j].coef*(*it2);
        }
        it2++;
      }
    }
    mpcFlag[i] = false;
  }
}

template<class Scalar> template<class Scalar1>
void
GenSubDomain<Scalar>::dispatchNodalData(FSCommPattern<Scalar> *pat, DistVec<Scalar1> *vec)
{
  Scalar1 *v = vec->subData(localSubNumber);
  for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
    for(int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
      sInfo.data[iNode] = v[(*scomm->sharedNodes)[iSub][iNode]];
  }
}

template<class Scalar> template<class Scalar1>
void
GenSubDomain<Scalar>::addNodalData(FSCommPattern<Scalar> *pat, DistVec<Scalar1> *vec)
{
  Scalar1 *v = vec->subData(localSubNumber);
  for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
    for(int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
      v[(*scomm->sharedNodes)[iSub][iNode]] += rInfo.data[iNode];
  }
}

template<> template<>
inline void
GenSubDomain<DComplex>::dispatchNodalData(FSCommPattern<DComplex> *pat, DistVec<double> *vec)
{
  double *v = vec->subData(localSubNumber);
  for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    FSSubRecInfo<DComplex> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
    for(int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
      sInfo.data[iNode] = v[(*scomm->sharedNodes)[iSub][iNode]];
  }
}

template<> template<>
inline void
GenSubDomain<DComplex>::addNodalData(FSCommPattern<DComplex> *pat, DistVec<double> *vec)
{
  double *v = vec->subData(localSubNumber);
  for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    FSSubRecInfo<DComplex> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
    for(int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
      v[(*scomm->sharedNodes)[iSub][iNode]] += ScalarTypes::Real(rInfo.data[iNode]);
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::setWICommSize(FSCommPattern<Scalar> *pat)
{
  for(int i = 0; i < scomm->numT(SComm::fsi); ++i)
    pat->setLen(subNumber,  scomm->neighbT(SComm::fsi,i), numNeighbWIdof[i]);
}

template<> double GenSubDomain<double>::Bcx(int i);
template<> DComplex GenSubDomain<DComplex>::Bcx(int i);

template<class Scalar>
void
GenSubDomain<Scalar>::multM(Scalar *localrhs, GenStackVector<Scalar> **u, int k)
{
// RT: 11/17/08
//  double omega2 = (isFluid(0) && !sinfo.isCoupled) ? packedEset[0]->helmCoef() : geoSource->shiftVal();
  double omega2 = geoSource->shiftVal();
  double omega = sqrt(omega2);

  // assemble localvec for MatVec product
  Scalar *localvec = (Scalar *) dbg_alloca(sizeof(Scalar)*c_dsa->size());
  if (u==0) {
    for(int i=0; i<c_dsa->size(); ++i)
      localrhs[i] = 0.0;
    makeFreqSweepLoad(localrhs, k, omega);
    return;
  }
  for(int i=0; i<c_dsa->size(); ++i) {
    localvec[i] = double(k)*(double(k-1)*(*u[k-1])[i] + 2.0*omega*(*u[k])[i]);
  }
  M->mult(localvec, localrhs);


  if(C_deriv) {
    for(int j=0; j<=k-1; ++j) {
      if(C_deriv[k-j-1]) {
        double ckj = DCombination(k,j);
        for(int i=0; i<c_dsa->size(); ++i) localvec[i] = -ckj*(*u[j+1])[i];
        C_deriv[k-j-1]->multAdd(localvec, localrhs);
      }
    }
  }

  makeFreqSweepLoad(localrhs, k, omega);  // this adds the residual effects of any prescribed displacements
                                          // and other loads that have non-zero kth derivative wrt omega
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeFreqSweepLoad(Scalar *d, int iRHS, double omega)
{
 int numUncon = c_dsa->size();
 GenStackVector<Scalar> force(numUncon, d);

 // Compute Right Hand Side Force = Fext + Fgravity + Fnh + Fpressure
 buildFreqSweepRHSForce<Scalar>(force, Muc, Cuc_deriv, iRHS, omega);

 for(int i = 0; i < numMPC; ++i) mpc[i]->rhs = 0.0; // set mpc rhs to zero for higher-order derivative solves
}


template<class Scalar>
void
GenSubDomain<Scalar>::multMCoupled1(Scalar *localrhs, GenStackVector<Scalar> **u, int k,
                                    FSCommPattern<Scalar> *wiPat)
{
  double omega2 = geoSource->shiftVal();
  double omega = sqrt(omega2);

  // assemble localvec for MatVec product
  Scalar *localvec = (Scalar *) dbg_alloca(sizeof(Scalar)*c_dsa->size());
  for(int i=0; i<c_dsa->size(); ++i) {
    localvec[i] = double(k)*(double(k-1)*(*u[k-1])[i] + 2.0*omega*(*u[k])[i]);
  }

  if(numWIdof) {
    int i, j;
    for(i=0; i<numWIdof; ++i) { localw[i] = localw_copy[i] = 0; }

    for(i = 0; i < numWInodes; ++i) {
      //DofSet thisDofSet = wetInterfaceDofs[i]^DofSet(DofSet::Helm);
      DofSet thisDofSet = wetInterfaceDofs[i]^(wetInterfaceDofs[i] & DofSet(DofSet::Helm)); // unmark Helm
      int nd = thisDofSet.count();
      int cdofs[6], dofs[6];
      c_dsa->number(wetInterfaceNodes[i], thisDofSet, cdofs);
      dsa->number(wetInterfaceNodes[i], thisDofSet, dofs);
      for(j = 0; j < nd; ++j) {
        localw[wetInterfaceMap[dofs[j]]] = localvec[cdofs[j]];
      }
    }

    // compute C uw to send to neighbors
    for(i = 0; i < scomm->numT(SComm::fsi); ++i) {
      if(subNumber != scomm->neighbT(SComm::fsi,i)) {
        FSSubRecInfo<Scalar> sInfo = wiPat->getSendBuffer(subNumber, scomm->neighbT(SComm::fsi,i));
        for(j=0; j<numNeighbWIdof[i]; ++j) sInfo.data[j] = 0.0;
        neighbKww->multAdd(localw, sInfo.data, glToLocalWImap, neighbGlToLocalWImap[i], true);
      }
      else {
        neighbKww->multAdd(localw, localw_copy, glToLocalWImap, true);
      }
    }
  }

  M->mult(localvec, localrhs);

  if(C_deriv) {
    for(int j=0; j<=k-1; ++j) {
      if(C_deriv[k-j-1]) {
        double ckj = DCombination(k,j);
        for(int i=0; i<c_dsa->size(); ++i) localvec[i] = -ckj*(*u[j+1])[i];
        C_deriv[k-j-1]->multAdd(localvec, localrhs);
      }
    }
  }

  makeFreqSweepLoad(localrhs, k, omega);
}

template<class Scalar>
void
GenSubDomain<Scalar>::multMCoupled2(Scalar *localrhs, FSCommPattern<Scalar> *wiPat)
{
  int i, j;
  if(numWIdof) {
    for(i = 0; i < scomm->numT(SComm::fsi); ++i) {
      if(subNumber != scomm->neighbT(SComm::fsi,i)) {
        FSSubRecInfo<Scalar> rInfo = wiPat->recData(scomm->neighbT(SComm::fsi,i), subNumber);
        for(j=0; j<numWIdof; ++j) localw_copy[j] += rInfo.data[j]/wweight[j];
      }
    }

    // put localw_copy into localrhs
    double omega2 = geoSource->shiftVal();
    for(i = 0; i < numWInodes; ++i) {
      //DofSet thisDofSet = wetInterfaceDofs[i]&DofSet(DofSet::Helm);
      //if(thisDofSet.count()) {
      if(wetInterfaceDofs[i].contains(DofSet::Helm)) {
        int thisNode = wetInterfaceNodes[i];
        int cdof= c_dsa->locate(thisNode, DofSet::Helm);
        int dof = dsa->locate(thisNode, DofSet::Helm);
        localrhs[cdof] -= localw_copy[wetInterfaceMap[dof]]/omega2;
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeBcx_scalar()
{
  int numdofs = dsa->size();
  bcx_scalar = new Scalar[numdofs];
  for(int i=0; i<numdofs; ++i) bcx_scalar[i] = Bcx(i);
}

template<class Scalar>
void
GenSubDomain<Scalar>::pade(GenStackVector<Scalar> *sol,  GenStackVector<Scalar> **u, double *h, double x)
{
  int len = sol->size();
  int i, j, k, n, r;
  // general N-point pade extrapolation
  if(rebuildPade) {
    // first time, allocate storage
    int l = domain->solInfo().getSweepParams()->padeL;
    int m = domain->solInfo().getSweepParams()->padeM;
    int padeN = domain->solInfo().getSweepParams()->padeN;
    int nRHS = domain->solInfo().getSweepParams()->nFreqSweepRHS;
    if(subNumber == 0) fprintf(stderr, " ... Computing %d-point Pade coefficients (l = %d, m = %d) ... \n", padeN, l, m);
    ia = l+1;
    ib = m+1;
    // allocate storage first time only
    ia = l+1;
    if(a == 0) {
      a = new GenVector<Scalar> * [ia];
      for(i = 0; i < ia; ++i) a[i] = new GenVector<Scalar>(len);
      b = new GenVector<Scalar> * [ib];
      for(i = 0; i < ib; ++i) b[i] = new GenVector<Scalar>(len);
      P = new GenVector<Scalar>(len);
      Q = new GenVector<Scalar>(len);
    }
    // compute P, Q coefficients by solving system Ax = b
    int dim = l+m+1;
    GenFullM<Scalar> A(dim); // LHS
    Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar)*dim);
    for(i=0; i<sol->size(); ++i) {
      A.zero();
       // assemble
      for(j=0; j<padeN; ++j) {
        int offset = j*(nRHS+1)+1;
        int I;
        for(n=0; n<nRHS; ++n) {
          if((I = j*nRHS+n) >= dim) break;
          // fill left block (a coefficients) of RHS
          for(k=n; k<=l; ++k)
            A[I][k] = -double(DFactorial(k)/DFactorial(k-n)*pow(h[j],k-n));
          // fill right block (b coefficients) of RHS
          for(r=0; r<=n; ++r)
            for(k=r; k<=m; ++k)
              if(k>0) A[I][l+k] += double(DCombination(n,r)*DFactorial(k)/DFactorial(k-r)*pow(h[j],k-r))*(*u[offset+n-r])[i];
          // fill LHS
          v[I] = -(*u[offset+n])[i];
        }
      }
      // factorize
      A.factor();
      // solve
      A.reSolve(v);
      // extract coefficients
      for(j=0; j<=l; ++j) (*a[j])[i] = v[j];
      (*b[0])[i] = 1.0; // b0 = 1
      for(j=1; j<=m; ++j) (*b[j])[i] = v[l+j];
    }
    rebuildPade = false;
  }
  // assemble P, Q and solution P/Q
  P->zero();
  Q->zero();
  for(i = 0; i < ia; ++i) P->linAdd(pow(x,i),*a[i]);
  for(i = 0; i < ib; ++i) Q->linAdd(pow(x,i),*b[i]);
  for(j=0; j<sol->size(); ++j) (*sol)[j] = (*P)[j] / (*Q)[j];
}

template<class Scalar>
void
GenSubDomain<Scalar>::updateActiveSet(Scalar *v, double tol, int flag, bool &statusChange)
{
  // flag = 0 : dual planing
  // flag = 1 : primal planing
  int *chgstatus = (int *) alloca(numMPC*sizeof(int));
  for(int i = 0; i < numMPC; ++i) chgstatus[i] = -1;  // set to 0 to remove, 1 to add

  for(int i = 0; i<scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpc[locMpcNb]->type == 1) { // inequality constraint requiring planing

      if(flag == 0) { // dual planing
        if(mpcStatus1[locMpcNb] == 1) { // active set expansion only: if constraint was initially active then it will not change status
          // if active and lambda < 0 then remove from active set
          if(mpc[locMpcNb]->active && ScalarTypes::lessThan(v[scomm->mapT(SComm::mpc, i)], tol)) chgstatus[locMpcNb] = 1; 
          // if not active and lambda >= 0 then add to the active set
          if(!mpc[locMpcNb]->active && ScalarTypes::greaterThanEq(v[scomm->mapT(SComm::mpc, i)], tol)) chgstatus[locMpcNb] = 0;
        }
      }
      else { // primal planing
        if(mpcStatus1[locMpcNb] == 0) { // active set contraction only: if constraint was initially inactive then it will not change status
          // if not active and w <= 0 then add to active set
          if(!mpc[locMpcNb]->active && ScalarTypes::lessThanEq(v[scomm->mapT(SComm::mpc, i)], tol)) chgstatus[locMpcNb] = 0;
          // if active and w > 0 then remove from the active set
          if(mpc[locMpcNb]->active && ScalarTypes::greaterThan(v[scomm->mapT(SComm::mpc, i)], tol)) chgstatus[locMpcNb] = 1;
        }
      }

    }
  }

  statusChange = false;
  for(int i = 0; i < numMPC; ++i) {
    if(chgstatus[i] > -1) {
      statusChange = true;
      mpc[i]->active = !bool(chgstatus[i]);
      if(solInfo().getFetiInfo().contactPrintFlag && mpcMaster[i]) { if(chgstatus[i] == 0) cerr << "-"; else cerr << "+"; }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::projectActiveIneq(Scalar *v)
{
  for(int i = 0; i<scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpc[locMpcNb]->type == 1 && mpc[locMpcNb]->active)
      v[scomm->mapT(SComm::mpc, i)] = 0.0;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::split(Scalar *v, Scalar *v_f, Scalar *v_c)
{
  // split v into free (v_f) and chopped (v_c) components
  for(int i = 0; i<totalInterfSize; ++i) { v_f[i] = v[i]; v_c[i] = 0.0; }
  for(int i = 0; i<scomm->lenT(SComm::mpc); ++i) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpc[locMpcNb]->type == 1 && mpc[locMpcNb]->active) {
      int iDof = scomm->mapT(SComm::mpc, i);
      if(ScalarTypes::greaterThan(v[iDof], 0.0)) v_c[iDof] = v[iDof];
      v_f[iDof] = 0.0;
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::bmpcQualify(vector<LMPCons *> *bmpcs, int *pstatus, int *nstatus)
{
  for(int i=0; i<bmpcs->size(); ++i) {
    LMPCons *bmpc = (*bmpcs)[i];
    if(bmpc->psub == subNumber) {
      int ccdof = cc_dsa->locate(globalToLocal((bmpc->terms)[0].nnum), (1 << (bmpc->terms)[0].dofnum));
      pstatus[i] = (ccdof > -1) ? 1 : 0;
    }
    if(bmpc->nsub == subNumber) {
      int ccdof = cc_dsa->locate(globalToLocal((bmpc->terms)[0].nnum), (1 << (bmpc->terms)[0].dofnum));
      nstatus[i] = (ccdof > -1) ? 1 : 0;
    }
  }
}

// ref: Dostal, Horak and Stefanica (IMACS 2005) "A scalable FETI-DP algorithm for coercive variational inequalities"
// every row of C matrix should have a unit norm to improve the condition number of the Feti operator
template<class Scalar>
void
GenSubDomain<Scalar>::normalizeCstep1(Scalar *cnorm)
{
  for(int i = 0; i < numMPC; ++i)
    for(int j = 0; j < mpc[i]->nterms; ++j)
      cnorm[localToGlobalMPC[i]] += mpc[i]->terms[j].coef*mpc[i]->terms[j].coef;
}

template<class Scalar>
void
GenSubDomain<Scalar>::normalizeCstep2(Scalar *cnorm)
{
  for(int i = 0; i < numMPC; ++i)
    for(int j = 0; j < mpc[i]->nterms; ++j)
      mpc[i]->terms[j].coef /= cnorm[localToGlobalMPC[i]];
}

template<class Scalar>
void GenSubDomain<Scalar>::mergeElemProps(double* props,
					  double* weights,
					  int propType)
{
#ifdef DISTRIBUTED
  int *nodeNumbers = (int *) dbg_alloca(sizeof(int)*maxNumNodes);
#endif
  for(int iele=0; iele<numele; ++iele)
    {
      const StructProp* sProp = packedEset[iele]->getProperty();
      double eleattr = 0.0;
      switch(propType)
	{
	case YOUNG:
	  eleattr = sProp->E;
	  break;
	case MDENS:
	  eleattr = sProp->rho;
	  break;
	case THICK:
	  eleattr = sProp->eh;
	  break;
	default:
	  assert(0);
	}

      int NodesPerElement = elemToNode->num(iele);
#ifdef DISTRIBUTED
      packedEset[iele]->nodes(nodeNumbers);
#endif
      for(int k=0; k<NodesPerElement; ++k)
	{
#ifdef DISTRIBUTED
	  int glbNodeNum = nodeNumbers[k];
#else
	  int locNodeNum = (*elemToNode)[iele][k];
	  int glbNodeNum = this->localToGlobal(locNodeNum);
#endif
	  // not thread safe!!!
	  props  [glbNodeNum] += eleattr;
	  weights[glbNodeNum] += 1.0;
	}
    }
  return;
}

template<class Scalar>
void
GenSubDomain<Scalar>::addSommer(SommerElement *ele)
{
 ele->dom = this;
 sommer[numSommer++] = ele;
 //if(sinfo.ATDARBFlag != -2.0) packedEset.elemadd(numele++,ele);  // XDEBUG
 if(sinfo.ATDARBFlag != -2.0)
  { ele->renum(glToLocalNode); packedEset.elemadd(numele++,ele); }
}

