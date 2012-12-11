#include <iostream>
using std::ifstream;
#include <sstream>
using std::istringstream;
#include <ctime>
#include <Utils.d/dbg_alloca.h>

#include <map>
#include <list>
using std::map;
using std::list;

#include <unistd.h>
#include <Timers.d/GetTime.h>
#include <Threads.d/PHelper.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>

#include <Utils.d/dofset.h>
#include <Driver.d/Domain.h>
#include <Element.d/Element.h>
#include <Utils.d/ModeData.h>
#include <Math.d/mathUtility.h>
#include <Driver.d/GeoSource.h>
#include <Feti.d/DistrVector.h>
#include <Corotational.d/DistrGeomState.h>

#include <Element.d/Rigid.d/RigidBeam.h>
#include <Element.d/Rigid.d/RigidThreeNodeShell.h>
#include <Element.d/Rigid.d/RigidFourNodeShell.h>
#include <Element.d/Rigid.d/RigidSolid6Dof.h>

#include<Sfem.d/Sfem.h>

extern int verboseFlag;
extern Sfem *sfem;
extern GeoSource *geoSource;
extern int totalNewtonIter;

// Global variable for mode data
ModeData modeData;

//----------------------------------------------------------------------------------

Domain::Domain(Domain &d, int nele, int *eles, int nnodes, int *nnums)
  : nodes(*new CoordSet(nnodes)), lmpc(0), fsi(0), ymtt(0), ctett(0),
    SurfEntities(0), MortarConds(0)
{
 initialize();

 int iele;
 numele = nele;        // number of elements
 for(int i=0; i < numele; ++i)
   packedEset.elemadd(i, d.packedEset[eles[i]]);

 numnodes = nnodes; // number of nodes
 for(int i=0; i < numnodes; ++i)
   if(d.nodes[nnums[i]] != NULL)
     nodes.nodeadd(i, *d.nodes[nnums[i]]);

 if(d.gravityFlag() ) {
   gravityAcceleration = new double [3];
   gravityAcceleration[0] = d.gravityAcceleration[0];
   gravityAcceleration[1] = d.gravityAcceleration[1];
   gravityAcceleration[2] = d.gravityAcceleration[2];
 }

 mftval = d.mftval;
 mptval = d.mptval;
 hftval = d.hftval;

 if(verboseFlag == 0) setSilent();
 else setVerbose();
}

Domain::Domain(Domain &d, Elemset *_elems, CoordSet *_nodes)
  : nodes(*_nodes), lmpc(0), fsi(0), ymtt(0), ctett(0), SurfEntities(0), MortarConds(0)
{
 initialize();

 int iele;
 numele = _elems->last();
 for(iele=0; iele < numele; ++iele)
   packedEset.elemadd(iele, (*_elems)[iele]);

 numnodes = _nodes->size();

 if(d.gravityFlag() ) {
   gravityAcceleration = new double [3];
   gravityAcceleration[0] = d.gravityAcceleration[0];
   gravityAcceleration[1] = d.gravityAcceleration[1];
   gravityAcceleration[2] = d.gravityAcceleration[2];
 }

 mftval = d.mftval;
 mptval = d.mptval;
 hftval = d.hftval;

 if(verboseFlag == 0) setSilent();
 else setVerbose();
}

Domain::Domain(int iniSize) : nodes(*(new CoordSet(iniSize*16))), packedEset(iniSize*16), lmpc(0,iniSize),
   fsi(0,iniSize), ymtt(0,iniSize), ctett(0,iniSize), SurfEntities(0,iniSize), MortarConds(0,iniSize)
{
 initialize();

 if(verboseFlag == 0) setSilent();
 else setVerbose();

 matrixTimers = new MatrixTimers;
}

void
Domain::makeAllDOFs() // build the dof connectivity
{
 // Test if allDOFs has been made already
 //if(allDOFs) return;

 int numele = packedEset.last(); // PJSA 5-2-05: include phantoms here

 int iele;
 int *pointers = new int[numele+1];
 pointers[0] = 0;
 maxNumDOFs  = 0;
 for(iele=0; iele < numele; ++iele) {
   int numDOFs = packedEset[iele]->numDofs();
   if(numDOFs > maxNumDOFs) maxNumDOFs = numDOFs;
   pointers[iele+1] = pointers[iele] + numDOFs;
 }
 int *targets = new int[ pointers[numele] ];
 for(iele=0; iele < numele; ++iele)
   packedEset[iele]->dofs(*dsa, targets + pointers[iele]);
 if(allDOFs) delete allDOFs;
 allDOFs = new Connectivity(numele, pointers, targets);

 maxNumNodes = 0;
 // compute maximum number of nodes per element
 for(iele=0; iele < numele; ++iele) {
   int numNodesPerElement = packedEset[iele]->numNodes();
   maxNumNodes = myMax(maxNumNodes, numNodesPerElement);
 }
}

// build the dof connectivity for fluid
// ADDED FOR HEV PROBLEM, EC, 20070820
void
Domain::makeAllDOFsFluid()
{
 // Test if allDOFsFluid has been made already
 if(allDOFsFluid) return;
 int numele = (*(geoSource->getPackedEsetFluid())).last(); // include phantoms here

 int iele;
 int *pointers = new int[numele+1];
 pointers[0] = 0;
 maxNumDOFsFluid  = 0;
 for(iele=0; iele < numele; ++iele) {
   int numDOFs = (*(geoSource->getPackedEsetFluid()))[iele]->numDofs();
   if(numDOFs > maxNumDOFsFluid) maxNumDOFsFluid = numDOFs;
   pointers[iele+1] = pointers[iele] + numDOFs;
 }
 int *targets = new int[ pointers[numele] ];
 for(iele=0; iele < numele; ++iele)
   (*(geoSource->getPackedEsetFluid()))[iele]->dofs(*dsaFluid, targets + pointers[iele]);
 allDOFsFluid = new Connectivity(numele, pointers, targets);

 maxNumNodesFluid = 0;
 // compute maximum number of nodes per element
 for(iele=0; iele < numele; ++iele) {
   int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
   maxNumNodesFluid = myMax(maxNumNodesFluid, numNodesPerElement);
 }
}

// This routine creates the array of boundary condition
// in the global DOF vector
void
Domain::make_bc(int *bc, double *bcx)
{
 // Initialize all boundary conditions to free
 // and all boundary condition values to zero
 int numdof = dsa->size();
 int i;
 for(i=0; i<numdof; ++i) {
   bc[i]  = BCFREE;
   bcx[i] = 0.0;
 }

 // Set the Neuman boundary conditions
 for(i=0; i<numNeuman; ++i) {
   int dof  = dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found FORCE on non-existant dof: node %d dof %d\n",
     //          nbc[i].nnum+1,nbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCLOAD) {
     //filePrint(stderr," *** WARNING: check input, multiple FORCEs defined at node %d"
     //          ", dof %d\n",nbc[i].nnum+1,nbc[i].dofnum+1);
     bcx[dof] += nbc[i].val;
   }
   else {
     bc[dof] = BCLOAD;
     bcx[dof] = nbc[i].val;
   }
 }

 // Set the real part of the Complex Neumann boundary conditions
 for(i=0; i<numComplexNeuman; ++i) {
   int dof  = dsa->locate(cnbc[i].nnum, 1 << cnbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found FORCE on non-existant dof: node %d dof %d\n",
     //          cnbc[i].nnum+1,cnbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCLOAD) {
     //filePrint(stderr," *** WARNING: check input, multiple FORCEs defined at node %d"
     //          ", dof %d\n",cnbc[i].nnum+1,cnbc[i].dofnum+1);
     bcx[dof] += cnbc[i].reval;
   }
   else {
     bc[dof] = BCLOAD;
     bcx[dof] = cnbc[i].reval;
   }
 }

 // Set the dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
   int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found DISP on non-existant dof: node %d dof %d\n",
     //                dbc[i].nnum+1,dbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCFIXED) {
     //filePrint(stderr," *** WARNING: check input, found repeated DISP"
     //               " (node %d, dof %d)\n",dbc[i].nnum+1,dbc[i].dofnum+1);
   }

   bc[dof] = BCFIXED;
   bcx[dof] = dbc[i].val;
 }

 // Set the real part of the Complex Dirichlet boundary conditions
 for(i=0; i<numComplexDirichlet; ++i) {
   int dof  = dsa->locate(cdbc[i].nnum, 1 << cdbc[i].dofnum);
   if(dof < 0) {
     //filePrint(stderr," *** WARNING: Found DISP on non-existant dof: node %d dof %d\n",
     //                cdbc[i].nnum+1,cdbc[i].dofnum+1);
     continue;
   }
   if(bc[dof] == BCFIXED) {
     //filePrint(stderr," *** WARNING: check input, found repeated Complex DISP"
     //               " (node %d, dof %d)\n",cdbc[i].nnum+1,cdbc[i].dofnum+1);
   }
   bc[dof] = BCFIXED;
   bcx[dof] = cdbc[i].reval;
 }
}

void
Domain::make_bc(int *bc, DComplex *bcx)
{
 int numdof = dsa->size();
 int i;
 for(i=0; i<numdof; ++i) {
   bc[i]  = BCFREE;
   bcx[i] = 0.0;
 }

// Set the real Neumann boundary conditions
 for(i=0; i<numNeuman; ++i) {
   int dof  = dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCLOAD) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " FORCE (node %d, dof %d)\n",nbc[i].nnum,nbc[i].dofnum);
   }
   bc[dof] = BCLOAD;
   bcx[dof] = DComplex(nbc[i].val,0.0);
 }

// Set the Complex Neuman boundary conditions
 for(i=0; i<numComplexNeuman; ++i) {
   int dof  = dsa->locate(cnbc[i].nnum, 1 << cnbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCLOAD) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " HFORCE (node %d, dof %d)\n",
        //                 cnbc[i].nnum,cnbc[i].dofnum);
   }
   bc[dof] = BCLOAD;
   bcx[dof] = DComplex(cnbc[i].reval,cnbc[i].imval);
 }

// Set the real Dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
   int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCFIXED) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " DISP (node %d, dof %d)\n",dbc[i].nnum,dbc[i].dofnum);
   }
   bc[dof] = BCFIXED;
   bcx[dof] = DComplex(dbc[i].val,0.0);
 }

// Set the Complex Dirichlet boundary condtions
 for(i=0; i<numComplexDirichlet; ++i) {
   int dof  = dsa->locate(cdbc[i].nnum, 1 << cdbc[i].dofnum);
   if(dof < 0) continue;
   if(bc[dof] == BCFIXED) {
        //fprintf(stderr," *** WARNING: check input, found repeated"
        //               " HDISP (node %d, dof %d)\n",
        //               cdbc[i].nnum,cdbc[i].dofnum);
   }
   bc[dof] = BCFIXED;
   bcx[dof] = DComplex(cdbc[i].reval, cdbc[i].imval);
 }
}


int
Domain::addDMass(int n, int d, double mass, int jdof)
{
 DMassData *dmass = new DMassData;
 dmass->next = firstDiMass;
 firstDiMass = dmass;

 firstDiMass->node   = n;
 firstDiMass->dof    = d;
 firstDiMass->diMass = mass;
 firstDiMass->jdof   = jdof;

 nDimass++;

 return 0;
}

// set Linear Multipoint Constrains
int
Domain::addLMPC(LMPCons *_MPC, bool checkflag)
{
 // as currently implemented all LMPCs (both real & complex) must use unique id numbers
 // eg if LMPC (real) id numbers 1 to 10, CLMPC (complex) should be 11 to 20
 // however you can use CLMPC to add a complex term to a previously real LMPC with the same id number
 // also inequality constraint id numbers can be the same as equality constraint id numbers
 // NEW IDEA: if checkflag is set to false, then always create new LMPC, used for various internal calls but not Parser
 if(!checkflag) { lmpc[numLMPC++] = _MPC; return numLMPC-1; } // PJSA 7-19-06

 // Verify if lmpc was already defined
 int i=0;
 while((i < numLMPC) && ((lmpc[i]->lmpcnum != _MPC->lmpcnum) || lmpc[i]->type != _MPC->type)) i++;

 // if LMPC not previously defined create new
 if(i==numLMPC) {
   lmpc[numLMPC++] = _MPC;
   if(_MPC->isComplex) numComplexLMPC++;
   if(_MPC->type == 1) numCTC++;
 }
 // if LMPC already defined overwrite rhs and add terms
 else {
  //filePrint(stderr," *** WARNING: rhs of LMPC number %d overwritten\n", lmpc[i]->lmpcnum);
  //filePrint(stderr," ***          previous terms are added\n");
  int j;
  for(j=0;j<lmpc[i]->nterms;j++)
    _MPC->addterm(&lmpc[i]->terms[j]);
  delete lmpc[i];
  lmpc[i] = _MPC;
 }
 return i; // PJSA: return global MPC id
}

void Domain::printLMPC()
{
 filePrint(stderr," ... LMPC list :                    ...\n");
 int i;
 for(i=0; i<numLMPC; i++) {
   filePrint(stderr," ...    lmpc number %d (i = %d) : ", lmpc[i]->lmpcnum, i);
   if(lmpc[i]->isComplex) {
     filePrint(stderr,"rhs = (%f,%f)\n", lmpc[i]->rhs.c_value.real(), lmpc[i]->rhs.c_value.imag());
     filePrint(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       filePrint(stderr,"        %d        %d        (%f,%f)\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.c_value.real(), lmpc[i]->terms[j].coef.c_value.imag());
   }
   else {
     filePrint(stderr,"rhs = %f\n", lmpc[i]->rhs.r_value);
     filePrint(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       filePrint(stderr,"        %d        %d        %20.16f\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.r_value);
   }
 }
}

void Domain::printLMPC2()
{
 // use fprintf instead of filePrint
 fprintf(stderr," ... LMPC list :                    ...\n");
 int i;
 for(i=0; i<numLMPC; i++) {
   fprintf(stderr," ...    lmpc number %d : ", lmpc[i]->lmpcnum);
   if(lmpc[i]->isComplex) {
     fprintf(stderr,"rhs = (%f,%f)\n", lmpc[i]->rhs.c_value.real(), lmpc[i]->rhs.c_value.imag());
     fprintf(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       fprintf(stderr,"        %d        %d        (%f,%f)\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.c_value.real(), lmpc[i]->terms[j].coef.c_value.imag());
   }
   else {
     fprintf(stderr,"rhs = %f\n", lmpc[i]->rhs.r_value);
     fprintf(stderr," ...    node      dof       coef\n");
     int j;
     for(j=0; j<lmpc[i]->nterms; j++)
       fprintf(stderr,"        %d        %d        %20.16f\n",
                 lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum,
                 lmpc[i]->terms[j].coef.r_value);
   }
 }
}


void Domain::normalizeLMPC()
{
 // PJSA 5-24-06
 for(int i=0; i<numLMPC; i++) {
   if(lmpc[i]->isComplex) {
     cerr << " *** WARNING: normalizeLMPC not implemented for complex coefficients \n";
   }
   else {
     double cnorm = 0.0;
     for(int j=0; j<lmpc[i]->nterms; j++)
       cnorm += lmpc[i]->terms[j].coef.r_value * lmpc[i]->terms[j].coef.r_value;
     cnorm = sqrt(cnorm);
     for(int j=0; j<lmpc[i]->nterms; j++)
       lmpc[i]->terms[j].coef.r_value /= cnorm;
     lmpc[i]->rhs.r_value /= cnorm;
   }
 }
}

void Domain::setPrimalLMPCs(int& numDual, int &numPrimal)
{
 numDual = numLMPC; numPrimal = 0;
 if(solInfo().fetiInfo.mpcflag == 2) { // convert all dual mpcs (type 0) to primal
   for(int i=0; i<numLMPC; i++)
     if(lmpc[i]->type == 0) {
       lmpc[i]->type = 3;
       numPrimal++;
       numDual--;
     }
 }
/*
 if(solInfo().fetiInfo.cmpc) { // convert corner bmpcs (type 2) to primal
   if(cornerWeight) {
     for(int i=0; i<numLMPC; i++)
       if(lmpc[i]->type == 2) {
         if(cornerWeight[9*lmpc[i]->terms[0].nnum+lmpc[i]->terms[0].dofnum] > 0) {
           lmpc[i]->type = 5;
           numPrimal++;
           numDual--;
         }
       }
   }
 }
*/
 //cerr << "numDual = " << numDual << ", numPrimal = " << numPrimal << endl;
 // note other strategies may be implemented here
}

Connectivity *
Domain::makeLmpcToNode()
{
 int size = 0;
 // find size of target: total number of coefficients
 // involving a different node
 // Note: done by double loop because assume number of terms is small
 int numtarget = 0;
 int i,j,jj;
 for(i=0; i<numLMPC; i++){
   if(lmpc[i]->isPrimalMPC()) continue;  // skip primal mpcs
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j) numtarget++;
   }
   size++;
 }
 // fill target with coefficient nodes
 int *pointer = new int[size+1];
 int *target  = new int[numtarget];
 size = 0; numtarget = 0;
 for(i=0; i<numLMPC; i++){
   if(lmpc[i]->isPrimalMPC()) continue;  // skip primal mpcs
   pointer[size]=numtarget;
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j){ target[numtarget]=lmpc[i]->terms[j].nnum; numtarget++;}
   }
   size++;
 }
 pointer[size]=numtarget;
 return new Connectivity(size,pointer,target);
}

Connectivity *
Domain::makeLmpcToNode_primal()
{
 int size = 0;
 // find size of target: total number of coefficients
 // involving a different node
 // Note: done by double loop because assume number of terms is small
 int numtarget = 0;
 int i,j,jj;
 for(i=0; i<numLMPC; i++){
   if(!lmpc[i]->isPrimalMPC()) continue;  // skip dual mpcs
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j) numtarget++;
   }
   size++;
 }
 // fill target with coefficient nodes
 int *pointer = new int[size+1];
 int *target  = new int[numtarget];
 size = 0; numtarget = 0;
 for(i=0; i<numLMPC; i++){
   if(!lmpc[i]->isPrimalMPC()) continue;  // skip dual mpcs
   pointer[size]=numtarget;
   for(j=0; j<lmpc[i]->nterms; j++){
     for(jj=0; jj<j; jj++)
       if(lmpc[i]->terms[j].nnum==lmpc[i]->terms[jj].nnum) break;
     if (jj==j){ target[numtarget]=lmpc[i]->terms[j].nnum; numtarget++;}
   }
   size++;
 }
 pointer[size]=numtarget;
 return new Connectivity(size,pointer,target);
}

void
Domain::makeFsiToNode()
{
 int size = numFSI;
 // find size of target: total number of coefficients
 // involving a different node
 // Note: done by double loop because assume number of terms is small
 int numtarget = 0;
 int i,j,jj;
 for(i=0; i<size; i++){
   numtarget++;  // for fluid node
   for(j=0; j<fsi[i]->nterms; j++){
     if(fsi[i]->terms[j].nnum==fsi[i]->lmpcnum) continue; //HB: case of same fluid & structure node
     for(jj=0; jj<j; jj++) // look if not already found among structure nodes
       if(fsi[i]->terms[j].nnum==fsi[i]->terms[jj].nnum) break;
     if(jj==j) numtarget++;
   }
 }
 // fill target with coefficient nodes
 int *pointer = new int[size+1];
 int *target  = new int[numtarget];
 int count = 0;
 for(i=0; i<numFSI; i++){
   pointer[i]=count;
   target[count++] = fsi[i]->lmpcnum;  // fluid node
   for(j=0; j<fsi[i]->nterms; j++){
     if(fsi[i]->terms[j].nnum==fsi[i]->lmpcnum) continue; //HB: case of same fluid & structure node
     for(jj=0; jj<j; jj++) // look if not already found among structure nodes
       if(fsi[i]->terms[j].nnum==fsi[i]->terms[jj].nnum) break;
     if(jj==j) { target[count]=fsi[i]->terms[j].nnum; count++; }
   }
 }
 pointer[i]=numtarget;
 if (count!=pointer[i]) filePrint(stderr,"*** ERROR in Domain::makeFsiToNode");  // Check
 fsiToNode = new Connectivity(size,pointer,target);
 nodeToFsi = fsiToNode->reverse();
}

void Domain::getInterestingDofs(DofSet &ret, int glNode)
{
  if(glNode > nodeToFsi->csize()) return;
  for(int i=0; i<nodeToFsi->num(glNode); ++i) {
    int iFSI = (*nodeToFsi)[glNode][i];
    if(fsi[iFSI]->lmpcnum == glNode) ret.mark(DofSet::Helm) ;// check fluid dof
    for(int j=0; j<fsi[iFSI]->nterms; j++) // check structure dofs
      if(fsi[iFSI]->terms[j].nnum == glNode) ret.mark(1 << fsi[iFSI]->terms[j].dofnum);
  }
}

double ** Domain::getCMatrix()
//returns the coupling matrix for hydroelastic vibration problems
//ADDED FOR HEV PROBLEM, EC, 20070820
{
  if(C_condensed) return C_condensed;
  int np = c_dsaFluid->size();
  int nu = c_dsa->size();
  nuNonZero = 0;
  double ** C = new double * [nu];
  int* umap_inv = new int[nu]; // map from u wet interface numbering to c_dsa numbering
  int* umap_add_temp = new int[nu];
  int** pmap_inv = new int* [nu];
  int* npNonZero_full = new int[nu];

  for(int i=0; i<nu; ++i)   {
    C[i] = 0;
    umap_inv[i] = -1;
    umap_add_temp[i] = -1;
    pmap_inv[i] = new int[np];
    npNonZero_full[i] = 0;
    for (int j=0; j<np; ++j)  {
      pmap_inv[i][j] = -1;
    }
  }

  for(int iFSI=0; iFSI<numFSI; ++iFSI) {

    int pindex = c_dsaFluid->locate(fsi[iFSI]->fluid_node, DofSet::Potential);

    if (pindex < 0) {
      continue;
    }

    for(int j=0; j<fsi[iFSI]->nterms; j++) {
      int uindex = c_dsa->locate(fsi[iFSI]->terms[j].nnum, 1 << fsi[iFSI]->terms[j].dofnum);

      if (uindex < 0)  {
        continue;
      }

      double C_up = fsi[iFSI]->terms[j].coef.r_value;
      if (C_up != 0.) {
        if(C[uindex] == 0) {
          C[uindex] = new double[np];
          for (int jp=0; jp<np; ++jp)
            C[uindex][jp]=0;
          umap_inv[uindex] = nuNonZero;
          umap_add_temp[nuNonZero++] = dsa->locate(fsi[iFSI]->terms[j].nnum, 1 << fsi[iFSI]->terms[j].dofnum);
        }
        C[uindex][pindex] = C_up;
        pmap_inv[uindex][pindex] = (npNonZero_full[uindex])++;
      }
    }
  }

  umap = new int[nuNonZero];
  umap_add = new int[nuNonZero];
  pmap = new int*[nuNonZero];
  npNonZero = new int[nuNonZero];
  C_condensed = new double * [nuNonZero];

  for(int i=0; i<nuNonZero; ++i)  {
    C_condensed[i] = 0;
    umap[i] = -1;
    npNonZero[i] = 0;
  }

  for(int i=0; i<nu; ++i)  {
    int iu = umap_inv[i];
    if (iu >= 0 )  {
      umap[iu] = i;
      pmap[iu] = new int[npNonZero_full[i]];
      for (int j=0; j < npNonZero_full[i]; ++j)  {
        pmap[iu][j] = -1;
      }
      for (int jp=0; jp < np; ++jp)  {
        int K = pmap_inv[i][jp];
        if (K>=0)  {
          pmap[iu][K] = jp;
        }
      }
    }
  }

  for(int i=0; i<nuNonZero; ++i)  {
    umap_add[i] = umap_add_temp[i];
    C_condensed[i] = C[umap[i]];
    npNonZero[i] = npNonZero_full[umap[i]];
  }
  delete [] C;
  return C_condensed;
}

void
Domain::trMultC(const Vector& x, Vector& y)
{
/*
  // computes y = C^T*x
  double** C_NZrows = getCMatrix();

  int nnp = domain->numUnconFluid();
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  int *unconstrNum = c_dsa->getUnconstrNum();

  y.zero();
  for(int i=0; i < domain->nuNonZero; ++i) {
    for(int jp=0; jp < nnp; ++jp) {
      y[jp] += rhoFluid*C_NZrows[i][jp]*x[unconstrNum[umap_add[i]]];
    }
  }
*/
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  y.zero();
  for(int i=0; i<numFSI; ++i) {
    int pindex = c_dsaFluid->locate(fsi[i]->fluid_node, DofSet::Potential);
    if (pindex < 0) continue;
    for(int j=0; j<fsi[i]->nterms; j++) {
      int uindex = c_dsa->locate(fsi[i]->terms[j].nnum, 1 << fsi[i]->terms[j].dofnum);
      if (uindex < 0) continue;
      double coef = fsi[i]->terms[j].coef.r_value;
      y[pindex] += rhoFluid*coef*x[uindex];
    }
  }
}

void
Domain::multC(const Vector& x, Vector& y)
{
/*
  // computes y = C*x
  double** C_NZrows = getCMatrix();

  int nnp = domain->numUnconFluid();
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  int *unconstrNum = c_dsa->getUnconstrNum();

  y.zero();
  for(int i=0; i<domain->nuNonZero; ++i) {
    for(int k=0; k<domain->npNonZero[i]; ++k) {
      int K = domain->pmap[i][k];
      y[unconstrNum[umap_add[i]]] += rhoFluid*C_NZrows[i][K]*x[K];
    }
  }
*/
  double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;
  y.zero();
  for(int i=0; i<numFSI; ++i) {
    int pindex = c_dsaFluid->locate(fsi[i]->fluid_node, DofSet::Potential);
    if (pindex < 0) continue;
    for(int j=0; j<fsi[i]->nterms; j++) {
      int uindex = c_dsa->locate(fsi[i]->terms[j].nnum, 1 << fsi[i]->terms[j].dofnum);
      if (uindex < 0) continue;
      double coef = fsi[i]->terms[j].coef.r_value;
      y[uindex] += rhoFluid*coef*x[pindex];
    }
  }
}

// set dirichlet boundary conditions
int Domain::setDirichlet(int _numDirichlet, BCond *_dbc)
{
  // dbc = dirichlet boundary conditions
  if (dbc) {
    // Allocate memory for correct number of dbc
    BCond *nd = new BCond[numDirichlet+_numDirichlet];

    // copy old dbc
    int i;
    for (i = 0; i < numDirichlet; ++i)
      nd[i] = dbc[i];

    // copy new dbc
    for (i = 0; i<_numDirichlet; ++i)
      nd[i+numDirichlet] = _dbc[i];

    // set correct number of dbc
    numDirichlet += _numDirichlet;

    // delete old array of dbc
    delete [] dbc;

    // set new pointer to correct number of dbc
    dbc = nd;
  }
  else {
    numDirichlet = _numDirichlet;
    dbc          = _dbc;
  }

  return 0;
}

// set dirichlet boundary conditions in Fluid
int Domain::setDirichletFluid(int _numDirichletFluid, BCond *_dbcFluid)
{
  if (dbcFluid) {
    // Allocate memory for correct number of dbcFluid
    BCond *ndFluid = new BCond[numDirichletFluid+_numDirichletFluid];

    // copy old dbcFluid
    int i;
    for (i = 0; i < numDirichletFluid; ++i)
      ndFluid[i] = dbcFluid[i];

    // copy new dbcFluid
    for (i = 0; i<_numDirichletFluid; ++i)
      ndFluid[i+numDirichletFluid] = _dbcFluid[i];

    // set correct number of dbcFluid
    numDirichletFluid += _numDirichletFluid;

    // delete old array of dbcFluid
    delete [] dbcFluid;

    // set new pointer to correct number of dbcFluid
    dbcFluid = ndFluid;
  }
  else {
    numDirichletFluid = _numDirichletFluid;
    dbcFluid          = _dbcFluid;
  }

  return 0;
}

int Domain::setNeuman(int _numNeuman, BCond *_nbc)
{
 if(nbc) {
   // Allocate memory for correct number of dbc
   BCond *nd = new BCond[numNeuman+_numNeuman];

   // copy old nbc
   int i;
   for(i=0; i < numNeuman; ++i)
      nd[i] = nbc[i];

   // copy new nbc
   for(i=0; i<_numNeuman; ++i)
     nd[i+numNeuman] = _nbc[i];

   // set correct number of dbc
   numNeuman += _numNeuman;

   // delete old array of dbc
   delete [] nbc;

   // set new pointer to correct number of dbc
   nbc = nd;
 }
 else {
   numNeuman = _numNeuman;
   nbc          = _nbc;
 }
 return 0;
}

int
Domain::setMFTT(MFTTData *_mftval)
{
 mftval = _mftval;
 return 0;
}

int
Domain::setMPTT(MFTTData *_mptval)
{
 mptval = _mptval;
 return 0;
}

int
Domain::setHFTT(MFTTData *_hftval)
{
 hftval = _hftval;
 return 0;
}

int
Domain::addYMTT(MFTTData *_ymtt)
{
 //--- Verify if ymtt was already defined
 int i = 0;
 while(i < numYMTT && ymtt[i]->getID() != _ymtt->getID()) i++;

 // if YMTT not previously defined create new
 if(i == numYMTT) ymtt[numYMTT++] = _ymtt;

 // if YMTT already defined print warning message
 else
   filePrint(stderr," *** WARNING: YMTT %d has already been defined \n", _ymtt->getID());

 return 0;
}

void Domain::printYMTT()
{
  if(numYMTT > 0) filePrint(stderr," ... YMTT list :                    ...\n");
  int i;
  for(i = 0; i < numYMTT; i++) {
    filePrint(stderr," ...    ymtt %d : \n", ymtt[i]->getID());
    filePrint(stderr," ...    temperature   Young's Modulus \n");
    int j;
    for(j = 0; j < ymtt[i]->getNumPoints(); j++)
      filePrint(stderr,"         %f     %f\n",
              ymtt[i]->getT(j), ymtt[i]->getV(j));
  }
}

int
Domain::addCTETT(MFTTData *_ctett)
{
 //--- Verify if ctett was already defined
 int i = 0;
 while(i < numCTETT && ctett[i]->getID() != _ctett->getID()) i++;

 // if CTETT not previously defined create new
 if(i == numCTETT) ctett[numCTETT++] = _ctett;

 // if CTETT already defined print warning message
 else
   filePrint(stderr," *** WARNING: CTETT %d has already been defined \n", _ctett->getID());

 return 0;
}

void Domain::printCTETT()
{
  if(numCTETT > 0) filePrint(stderr," ... TETT list :                    ...\n");
  int i;
  for(i = 0; i < numCTETT; i++) {
    filePrint(stderr," ...    tett %d : \n", ctett[i]->getID());
    filePrint(stderr," ...    temperature   Coeff of Thermal Expansion \n");
    int j;
    for(j = 0; j < ctett[i]->getNumPoints(); j++)
      filePrint(stderr,"         %f     %f\n",
              ctett[i]->getT(j), ctett[i]->getV(j));
  }
}

int
Domain::setIDis6(int _numIDis6, BCond *_iDis6)
{
 numIDis6 = _numIDis6;
 iDis6    = _iDis6;
 return 0;
}

int
Domain::setIDis(int _numIDis, BCond *_iDis)
{
 numIDis = _numIDis;
 iDis    = _iDis;
 return 0;
}

int
Domain::setIDisModal(int _numIDisModal, BCond *_iDisModal)
{
 numIDisModal = _numIDisModal;
 iDisModal    = _iDisModal;
 return 0;
}

int
Domain::setIVel(int _numIVel, BCond *_iVel)
{
 numIVel = _numIVel;
 iVel    = _iVel;
 return 0;
}

int
Domain::setIVelModal(int _numIVelModal, BCond *_iVelModal)
{
 numIVelModal = _numIVelModal;
 iVelModal    = _iVelModal;
 return 0;
}

int Domain::setIAcc(int , BCond *)
{
// INITIAL ACCELERATION IS NOT NEEDED CURRENTLY
// SO IT IS READ FROM THE INPUT FILE BUT NOT STORED
// numIAcc = _numIAcc;
// iAcc    = _iAcc;
 return 0;
}

void
Domain::setGravity(double ax, double ay, double az)
{
 gravityAcceleration = new double [3];
 gravityAcceleration[0] = ax;
 gravityAcceleration[1] = ay;
 gravityAcceleration[2] = az;
}

void
Domain::setGravitySloshing(double gg)
{
 gravitySloshing = gg;
}

void
Domain::setUpData()
{
  startTimerMemory(matrixTimers->setUpDataTime, matrixTimers->memorySetUp);

  Elemset eset_tmp;
  if(sinfo.type == 0) {
    if(numLMPC) geoSource->getNonMpcElems(eset_tmp);
  }

  geoSource->setUpData();

  if(sinfo.type == 0) {
    if(numLMPC) {
      Connectivity *elemToNode_tmp = new Connectivity(&eset_tmp);
      Connectivity *nodeToElem_tmp = elemToNode_tmp->reverse();
      Connectivity *nodeToNode_tmp = nodeToElem_tmp->transcon(elemToNode_tmp);
      renumb_nompc = nodeToNode_tmp->renumByComponent(0);
      int numnodes_tmp = nodeToNode_tmp->csize();
      int *order = new int[numnodes_tmp];
      for(int i=0; i<numnodes_tmp; ++i) order[i] = -1;
      for(int i=0; i<numnodes_tmp; ++i)
        if(renumb_nompc.renum[i] >= 0)
          order[renumb_nompc.renum[i]] = i;
      renumb_nompc.order = order;
      delete elemToNode_tmp; delete nodeToElem_tmp; delete nodeToNode_tmp;
    }
  }


  if(!haveNodes) { numnodes = geoSource->getNodes(nodes); haveNodes = true; }
  else numnodes = geoSource->totalNumNodes();
  numele = geoSource->getElems(packedEset);
  numele = packedEset.last(); // XXXX

  // set boundary conditions
  int numBC;
  int numTextBC;
  BCond *bc;
  BCond *textBC;

  // set dirichlet
  numBC = geoSource->getDirichletBC(bc);
  numTextBC = geoSource->getTextDirichletBC(textBC);
  if (numTextBC)
    geoSource->augmentBC(numTextBC, textBC, bc, numBC);
  if (solInfo().HEV) {
    int numBCFluid;
    BCond *bcFluid;
    numBCFluid = geoSource->getDirichletBCFluid(bcFluid);
    setDirichletFluid(numBCFluid, bcFluid);
  }
  setDirichlet(numBC, bc);

  // set neuman
  numBC = geoSource->getNeumanBC(bc);
  setNeuman(numBC, bc);

  // set initial displacements
  numBC = geoSource->getIDis(bc);
  setIDis(numBC, bc);
  numBC = geoSource->getIDisModal(bc);
  if(solInfo().modal) { // for modal dynamics, keep the modal idisp and non-modal idisp separate
    setIDisModal(numBC, bc);
  }
  else { // for non-modal dynamics convert the modal idisp into non-modal idisp
    if(numBC) {
      filePrint(stderr, " ... Compute initial displacement from given modal basis ...\n");
      int numIDisModal = modeData.numNodes*6;
      BCond *iDis_new = new BCond[numIDis+numIDisModal];
      for(int i=0; i<numIDis; ++i) iDis_new[i] = iDis[i]; 
      modeData.addMultY(numBC, bc, iDis_new+numIDis, 6);
      numIDis += numIDisModal;
      delete [] iDis;
      iDis = iDis_new;
    }
  }

  // set init disp6
  numBC = geoSource->getIDis6(bc);
  setIDis6(numBC, bc);

  // set initial velocities
  numBC = geoSource->getIVel(bc);
  setIVel(numBC, bc);
  numBC = geoSource->getIVelModal(bc);
  if(solInfo().modal) { // for modal dynamics, keep the modal ivel and non-modal ivel separate
    setIVelModal(numBC, bc);
  }
  else {
    if(numBC) {
      filePrint(stderr, " ... Compute initial velocity from given modal basis ...\n");
      int numIVelModal = modeData.numNodes*6;
      BCond *iVel_new = new BCond[numIVel+numIVelModal];
      for(int i=0; i<numIVel; ++i) iVel_new[i] = iVel[i];
      modeData.addMultY(numBC, bc, iVel_new+numIVel, 6);
      numIVel += numIVelModal;
      delete [] iVel;
      iVel = iVel_new;
    }
  }

  // set Control Law
  claw = geoSource->getControlLaw();

  // PJSA: compute temperature dependent material properties
  computeTDProps();

  initNodalTemperatures();

/*
  -- Include here checks on LMPC --
     check for redundant LMPC by applying a stabilized
     Gram-Schmit on LMPC global constrain matrix (domain->lmpc)
        --> yields the left nullspace L in Lu=g
            => if rhs, i.e. g, not orth to left nullspace: ERROR
            => otherwise conserve only a set of non-redundant MPC
        (D.Rixen 04-28-99)
*/

  stopTimerMemory(matrixTimers->setUpDataTime, matrixTimers->memorySetUp);

}

#ifndef OUTPUTMESSAGE
#define OUTPUTMESSAGE
const char* OutputMessage[] = {
" ... Outputing Displacement         ... \n",
" ... Outputing Temperature          ... \n",
" ... Outputing StressXX             ... \n",
" ... Outputing StressYY             ... \n",
" ... Outputing StressZZ             ... \n",
" ... Outputing StressXY             ... \n",
" ... Outputing StressYZ             ... \n",
" ... Outputing StressXZ             ... \n",
" ... Outputing StrainXX             ... \n",
" ... Outputing StrainYY             ... \n",
" ... Outputing StrainZZ             ... \n",
" ... Outputing StrainYZ             ... \n",
" ... Outputing StrainXY             ... \n",
" ... Outputing StrainXZ             ... \n",
" ... Outputing HeatFlXX             ... \n",
" ... Outputing HeatFlXY             ... \n",
" ... Outputing HeatFlXZ             ... \n",
" ... Outputing GrdTempX             ... \n",
" ... Outputing GrdTempY             ... \n",
" ... Outputing GrdTempZ             ... \n",
" ... Outputing Von Mises Stress     ... \n",
" ... Outputing Principal Stress #1  ... \n",
" ... Outputing Principal Stress #2  ... \n",
" ... Outputing Principal Stress #3  ... \n",
" ... Outputing Principal Strain #1  ... \n",
" ... Outputing Principal Strain #2  ... \n",
" ... Outputing Principal Strain #3  ... \n",
" ... Outputing InXForce             ... \n",
" ... Outputing InYForce             ... \n",
" ... Outputing InZForce             ... \n",
" ... Outputing AXMoment             ... \n",
" ... Outputing AYMoment             ... \n",
" ... Outputing AZMoment             ... \n",
" ... Outputing Energies             ... \n",
" ... Outputing AeroForce            ... \n",
" ... Outputing EigenPair            ... \n",
" ... Outputing Von Mises Strain     ... \n",
" ... Outputing Pressure             ... \n",
" ... Outputing All 6 dofs           ... \n",
" ... Outputing Aero X Force         ... \n",
" ... Outputing Aero Y Force         ... \n",
" ... Outputing Aero Z Force         ... \n",
" ... Outputing Aero X Moment        ... \n",
" ... Outputing Aero Y Moment        ... \n",
" ... Outputing Aero Z Moment        ... \n",
" ... Outputing Velocity             ... \n",
" ... Outputing Acceleration         ... \n",
" ... Outputing Youngs Modulus       ... \n",
" ... Outputing Material Density     ... \n",
" ... Outputing Thickness            ... \n",
" ... Outputing Shape Attributes     ... \n",
" ... Outputing Shape Attributes     ... \n",
" ... Outputing Composite            ... \n",
" ... Outputing Dx                   ... \n",
" ... Outputing Dy                   ... \n",
" ... Outputing Dz                   ... \n",
" ... Outputing Rx                   ... \n",
" ... Outputing Ry                   ... \n",
" ... Outputing Rz                   ... \n",
" ... Outputing DispMod              ... \n",
" ... Outputing RotMod               ... \n",
" ... Outputing TotMod               ... \n",
" ... Outputing Rigid Body Modes     ... \n",
" ... Outputing Elem To Node Table   ... \n",
" ... Outputing Node To Elem Table   ... \n",
" ... Outputing Node To Node Table   ... \n",
" ... Outputing Aero Flux            ... \n",
" ... Outputing Heat Flux            ... \n",
" ... Outputing Grd Temp             ... \n",
" ... Outputing Velocity (6 dofs)    ... \n",
" ... Outputing Acceleration (6 dofs)... \n",
" ... Outputing Alphas Modes         ... \n",
" ... Outputing Error for Modes Decomp. ... \n",
" ... Outputing Sloshing Eigen Modes ... \n",
" ... Outputing Sloshing x-disp.     ... \n"
" ... Outputing Sloshing y-disp.     ... \n"
" ... Outputing Sloshing z-disp.     ... \n"
};
#endif

void
Domain::resProcessing(Vector &totRes, int index, double t)
{
  int iNode;
  int i;
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (i = 0; i < numOutInfo; ++i)  {

    if (oinfo[i].interval == 1) {
      fprintf(stderr,"%s",OutputMessage[oinfo[i].type]);
      switch(oinfo[i].type)  {

        default:
          break;
        case OutputInfo::Reactions:
          double (*xyz)[3] = new double[numnodes][3];
          for(iNode=0; iNode<numnodes; ++iNode) {

            if(nodes[iNode] == 0) continue;
            DofSet dof[] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp };
            int k;
            for(k = 0; k < 3; ++k) {
              int loc =   dsa->locate( iNode, dof[k].list());

              double v;
              if(loc >= 0)       v = totRes[loc];
              else               v = 0.0;
              xyz[iNode][k] = v;
            }
          }
          geoSource->outputNodeVectors(i, xyz, numnodes, t);
          delete [] xyz;
          break;
      }
    }
  }
}

Connectivity *
Domain::makeSommerToNode() 
{
 int size = numSommer;
 // Find out the number of targets we will have
 int *pointer = new int[size+1] ;
 int pp = 0;
 for(int i=0; i < size; ++i) {
   pointer[i] = pp;
   pp += sommer[i] ? sommer[i]->numNodes() : 0;
 }
 pointer[size] = pp;
 //int numtarget = pp;
 // Create the target array
 int *target = new int[pp];
 // Fill it in
 for(int i=0; i < size; ++i) {
   if(sommer[i]) sommer[i]->nodes(target+pointer[i]);
 }

 Connectivity *ret = new Connectivity(size, pointer, target);
 return ret;
}

Connectivity *
Domain::prepDirectMPC()
{
  // get all of the lmpcs from the packedEset
  // TODO don't extract elements which use lagrange multipliers or penalty
  for(int i=0; i<numLMPC; ++i) if(lmpc[i]) delete lmpc[i];
  numLMPC = 0;
  lmpc.deleteArray();
  int nEle = packedEset.last();
  int nRigid = 0;
  Elemset *rigidSet = new Elemset;
  for(int i = 0; i < nEle; ++i) {
    Element *ele = packedEset[i];
    if(ele != 0) {
      if(false /*dynamic_cast<RigidBeam*>(ele) || dynamic_cast<RigidThreeNodeShell*>(ele)
         || dynamic_cast<RigidFourNodeShell*>(ele) || dynamic_cast<RigidSolid6Dof*>(ele)*/) // TODO
        rigidSet->elemadd(nRigid++, ele);
      else {
        int n = ele->getNumMPCs();
        if(n > 0) {
          LMPCons **l = ele->getMPCs();
          for(int j = 0; j < n; ++j)
            lmpc[numLMPC++] = l[j];
          delete [] l;
        }
      }
    }
  }
  //printLMPC();
  if(nRigid > 0) {
    //std::cerr << "found " << nRigid << " rigid beam/shell/solid6 elements in the element set\n";
    std::set<int> blockedNodes;
    for(int i = 0; i < numDirichlet; ++i) blockedNodes.insert(dbc[i].nnum);
    rigidSet->collapseRigid6(blockedNodes);
    nRigid = rigidSet->last();
    StructProp p; // dummy property
    for(int i = 0; i < nRigid; ++i) {
      int n = (*rigidSet)[i]->getNumMPCs();
      if(n > 0) {
        if((*rigidSet)[i]->getProperty() == 0) { // this element was instantiated in Elemset::collapseRigid6
          (*rigidSet)[i]->buildFrame(nodes);     // need to call buildFrame and setProp to prep it.
          (*rigidSet)[i]->setProp(&p);
        }
        LMPCons **l = (*rigidSet)[i]->getMPCs();
        for(int j = 0; j < n; ++j)
          lmpc[numLMPC++] = l[j];
        delete [] l;
      }
    }
  }
  //cerr << "extracted " << numLMPC << " lmpcs from the element set\n";

  geoSource->makeDirectMPCs(numLMPC, lmpc);
  // MPC Connectivity treatment in direct way.
  std::multimap<int, int> mpcConnect;
  std::multimap<int, int>::iterator it;
  std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> ret;
  int ndMax = nodeToElem->csize();
  for(int i = 0; i < numLMPC; ++i) {
    for(int j = 1; j < lmpc[i]->nterms; ++j)
      if(lmpc[i]->terms[0].nnum != lmpc[i]->terms[j].nnum) {

        ret = mpcConnect.equal_range(lmpc[i]->terms[0].nnum);
        bool found = false;
        for(it = ret.first; it != ret.second; ++it)
          if((*it).second == lmpc[i]->terms[j].nnum) { found = true; break; }

        if(!found)
          mpcConnect.insert(std::pair<int,int>(lmpc[i]->terms[0].nnum, lmpc[i]->terms[j].nnum));
      }
  }

  int *ptr = new int[ndMax+1];
  for(int i = 0; i < ndMax; ++i)
    ptr[i] = 1;
  ptr[ndMax] = 0;
  for(it = mpcConnect.begin(); it != mpcConnect.end(); ++it)
    ptr[it->first]++;
  for(int i = 0; i < ndMax; ++i)
    ptr[i+1] += ptr[i];
  int *tg = new int[ptr[ndMax]];

  for(it = mpcConnect.begin(); it != mpcConnect.end(); ++it)
    tg[--ptr[it->first]] = it->second;
  for(int i = 0; i < ndMax; ++i)
    tg[--ptr[i]] = i;

  Connectivity renumToNode(ndMax, ptr, tg);
  Connectivity *elemToRenum = elemToNode->transcon(&renumToNode);
  Connectivity *renumToElem = elemToRenum->reverse();
  Connectivity *nodeToNodeDirect = renumToElem->transcon(elemToRenum);
  delete renumToElem;
  delete elemToRenum;
  delete rigidSet;
  return nodeToNodeDirect;
}

Renumber
Domain::getRenumbering()
{
 // create node to element connectivity from element to node connectivity
 if(nodeToElem) delete nodeToElem;
 nodeToElem = elemToNode->reverse();

 // create node to node connectivity
 if(nodeToNode) delete nodeToNode;
 nodeToNode = nodeToElem->transcon(elemToNode);

 //ADDED FOR HEV PROBLEM, EC, 20070820
 if(solInfo().HEV == 1 && solInfo().addedMass == 1) {
   int numWetNodes = 0;
   int *wetIFNodes = getAllWetInterfaceNodes(numWetNodes);
   Connectivity* nodeToNodeTemp = nodeToNode->combineAll(numWetNodes,wetIFNodes);
   delete nodeToNode;
   nodeToNode = nodeToNodeTemp->copy();
   delete nodeToNodeTemp;
 }

 // get number of nodes
 numnodes = nodeToNode->csize();

 if(solInfo().type == 0 || solInfo().type == 1) makeNodeToNode_sommer(); // single domain solvers

 // delete any previously allocated memory
 if(renumb.order) { delete [] renumb.order; renumb.order=0; }
 if(renumb.renum) { delete [] renumb.renum; renumb.renum=0; }
 if(renumb.xcomp) { delete [] renumb.xcomp; renumb.xcomp=0; }

 // renumber the nodes
 renumb = nodeToNode->renumByComponent(sinfo.renum);

 int *order = new int[numnodes];

 int i;
 for(i=0; i<numnodes; ++i)
   order[i] = -1;

 for(i=0; i<numnodes; ++i)
   if(renumb.renum[i] >= 0)
     order[renumb.renum[i]] = i;

 //if(sinfo.renum > 0 && verboseFlag && renumb.numComp > 1)
 //  filePrint(stdout," ... Number of components =%2d	...\n",renumb.numComp);

 Renumber ret;
 ret.order    = order;
 ret.renumb   = renumb.renum;
 renumb.order = order;

 return ret;
}

Renumber*
Domain::getRenumberingFluid()
{
 // create node to element connectivity from element to node connectivity
 if(nodeToElemFluid) delete nodeToElemFluid;
 nodeToElemFluid = elemToNodeFluid->reverse();

 // create node to node connectivity
 if(nodeToNodeFluid) delete nodeToNodeFluid;
 nodeToNodeFluid = nodeToElemFluid->transcon(elemToNodeFluid);

 // get number of nodes
 numnodesFluid = nodeToNodeFluid->csize();

 // renumber the nodes
#ifdef TFLOP
 alloca(0);
#endif
 renumbFluid = nodeToNodeFluid->renumByComponent(sinfo.renum);

 int *order = new int[numnodesFluid];

 int i;
 for(i=0; i<numnodesFluid; ++i)
   order[i] = -1;

 for(i=0; i<numnodesFluid; ++i)
   if(renumbFluid.renum[i] >= 0)
     order[renumbFluid.renum[i]] = i;

 Renumber* ret = new Renumber;
 ret->order    = order;
 ret->renumb   = renumbFluid.renum;
 renumbFluid.order = order;

 // Just for Gnu alloca function
#ifdef TFLOP
 alloca(0);
#endif

 return ret;

}
void
Domain::makeNodeToNode_sommer()
{
 if(solInfo().doFreqSweep && (numSommer > 0)) {
   Connectivity *sommerToNode = makeSommerToNode();
   Connectivity *nodeToSommer = sommerToNode->reverse();
   nodeToNode_sommer = nodeToSommer->transcon(sommerToNode);
   delete nodeToSommer; delete sommerToNode;
 }

 if( (solInfo().newmarkBeta==0.0)&&(solInfo().isAcoustic()) ) {
   Connectivity *sommerToNode = makeSommerToNode();
   Connectivity *nodeToSommer = sommerToNode->reverse();
   //Watch nodeToNode_sommer may not be of the right size !
   int numnodes = nodeToNode->csize();
   Connectivity *temp0=new Connectivity(numnodes);
   Connectivity *temp1 = temp0->modify();
   Connectivity *temp2 = nodeToSommer->transcon(sommerToNode);
   nodeToNode_sommer = temp1->transcon(temp2);
   delete nodeToSommer; delete sommerToNode;
   delete temp0;
   delete temp1;
   delete temp2;
 }
}

void
Domain::readInModes(char* modeFileName)
{
 filePrint(stderr," ... Read in Modes from file: %s ...\n",modeFileName);

 // Open file containing mode shapes and frequencies.
 FILE *f;
 if((f=fopen(modeFileName,"r"))==(FILE *) 0 ){
   filePrint(stderr," *** ERROR: Cannot open %s ... exiting\n",modeFileName);
   exit(0);
 }
 fflush(f);

 // Read in number of modes and number of nodes
 int count = fscanf(f, "%d%d", &modeData.numModes, &modeData.numNodes);

 // Allocation of memory for frequencies and mode shapes
 modeData.frequencies  = new double[modeData.numModes];

 typedef double (*T6)[6];
 modeData.modes        = new T6[modeData.numModes];
 modeData.nodes        = new int[modeData.numNodes];

 // Read frequencies and mode shapes
 int iMode, iNode, numsubstring;
 char input[500], *substring;
 double tmpinput[7];

 //int jj;
 for(iMode=0; iMode<modeData.numModes; ++iMode) {
   modeData.modes[iMode] = new double[modeData.numNodes][6];

   count = fscanf(f,"%lf\n",&modeData.frequencies[iMode]);

   //int nodeNum;
   for(iNode=0; iNode<modeData.numNodes; ++iNode) {

     char *c = fgets(input, 500, f);
     substring = strtok(input, " ");
     numsubstring = 0;
     do {
       tmpinput[numsubstring] = strtod(substring, NULL);
       if (strncmp(substring, "\n", 1) == 0)
         break;
       substring = strtok(NULL, " ");
       ++numsubstring;

     } while (substring != NULL);

     switch(numsubstring){
       case 3:
         modeData.nodes[iNode] = iNode;
         modeData.modes[iMode][iNode][0] = tmpinput[0];
         modeData.modes[iMode][iNode][1] = tmpinput[1];
         modeData.modes[iMode][iNode][2] = tmpinput[2];
         modeData.modes[iMode][iNode][3] = 0.0;
         modeData.modes[iMode][iNode][4] = 0.0;
         modeData.modes[iMode][iNode][5] = 0.0;
         break;
       case 4:
         modeData.nodes[iNode] = (int) tmpinput[0];
         modeData.modes[iMode][iNode][0] = tmpinput[1];
         modeData.modes[iMode][iNode][1] = tmpinput[2];
         modeData.modes[iMode][iNode][2] = tmpinput[3];
         modeData.modes[iMode][iNode][3] = 0.0;
         modeData.modes[iMode][iNode][4] = 0.0;
         modeData.modes[iMode][iNode][5] = 0.0;
         break;
       case 6:
         modeData.nodes[iNode] = iNode;
         modeData.modes[iMode][iNode][0] = tmpinput[0];
         modeData.modes[iMode][iNode][1] = tmpinput[1];
         modeData.modes[iMode][iNode][2] = tmpinput[2];
         modeData.modes[iMode][iNode][3] = tmpinput[3];
         modeData.modes[iMode][iNode][4] = tmpinput[4];
         modeData.modes[iMode][iNode][5] = tmpinput[5];
         break;
       case 7:
         modeData.nodes[iNode] = (int) tmpinput[0];
         modeData.modes[iMode][iNode][0] = tmpinput[1];
         modeData.modes[iMode][iNode][1] = tmpinput[2];
         modeData.modes[iMode][iNode][2] = tmpinput[3];
         modeData.modes[iMode][iNode][3] = tmpinput[4];
         modeData.modes[iMode][iNode][4] = tmpinput[5];
         modeData.modes[iMode][iNode][5] = tmpinput[6];
         break;
       default:
         filePrint(stderr, " *** ERROR: Check input for mode data (numstrings = %d) in file %s ... exiting\n", numsubstring, modeFileName);
         exit(0);
     }
   }
 }
}

void
Domain::getElementAttr(int fileNumber,int iAttr,double time) {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  double eleattr;

  // ... OUTPUT precision
  int p = oinfo[fileNumber].precision;

  // ... Allocate Vectors for Nodal Attributes
  Vector nodattr(numnodes,0.0);
  Vector ndWeight(numnodes,0.0);

  int iele,k;
  for(iele=0; iele<numele; ++iele) {

     int NodesPerElement = elemToNode->num(iele);

     switch (iAttr) {
     case YOUNG:
       eleattr=packedEset[iele]->getProperty()->E;
       break;
     case MDENS:
       eleattr=packedEset[iele]->getProperty()->rho;
       break;
     case THICK:
       eleattr=packedEset[iele]->getProperty()->eh;
       break;
     default:
       assert(0);
     }

     for(k=0; k<NodesPerElement; ++k) {
       nodattr[(*elemToNode)[iele][k]]  += eleattr;
       ndWeight[(*elemToNode)[iele][k]] += 1.0;
     }

     if(avgnum == 0 ) {
       for(k=0; k<NodesPerElement; ++k) {
        filePrint(oinfo[fileNumber].filptr,"% .*E",p,eleattr);
       }
       filePrint(oinfo[fileNumber].filptr,"\n");
       fflush(oinfo[fileNumber].filptr);
     }
  }

  if(avgnum == 1 ) {

    for (k=0; k<numnodes; ++k)  {
      if (ndWeight[k] == 0.0)
	{ nodattr[k] = 0.0; }
      else
	{ nodattr[k] /= ndWeight[k]; }
    }

    geoSource->outputNodeScalars(fileNumber, nodattr.data(), numnodes, time);
  }

}

void
Domain::getCompositeData(int iInfo,double time) {

  OutputInfo *oinfo = geoSource->getOutputInfo();

  int w   = oinfo[iInfo].width;
  int p   = oinfo[iInfo].precision;
  int lay = oinfo[iInfo].nodeNumber;

  // ... For time = 0.0: Output of pseudo-structure

  if ( time == 0.0 ) {

    int iele,icp;
    int numCPoints = numele*4;

    // ... Create Arrays

    if ( ! MidPoint ) {

      MidPoint = new double*[numele];

      int maxLayer=0;
      for(iele=0; iele<numele; ++iele) {
        MidPoint[iele] = packedEset[iele]->getMidPoint(nodes);
        maxLayer = myMax(maxLayer,packedEset[iele]->getCompositeLayer());
      }

      CPoint = new double**[maxLayer];

      int icl;
      for (icl=0;icl<maxLayer;icl++) {
        CPoint[icl] = new double*[numCPoints];
        for (icp=0;icp<numCPoints;icp++)
          CPoint[icl][icp] = new double[3];
      }
    }

    // ... Open top-file for pseudo-structure ".cross"

    char * basename = getBasename(oinfo[iInfo].filename);
    int fnamesize   = strlen(basename)+7;
    char * crossfile  = (char*) malloc(sizeof(char*)*(fnamesize));
    sprintf(crossfile,"%s%s%d",basename,".cross",lay+1);
    FILE * crossout = fopen(crossfile,"w");

    // ... Global Scaling factor for length of cross bars

    double minDist=10000000000.0;
    double dx,dy,dz,dd;
    int jele;
    for (iele=0; iele<numele; ++iele) {
      for (jele=iele+1; jele<numele; ++jele) {
         dx = MidPoint[iele][0]-MidPoint[jele][0];
         dy = MidPoint[iele][1]-MidPoint[jele][1];
         dz = MidPoint[iele][2]-MidPoint[jele][2];
	 dd = sqrt(dx*dx+dy*dy+dz*dz);
	 minDist = myMin(minDist,dd);
      }
    }

    crossScale = minDist * 1.25;

    // ... Determine end-points of cross bars

    for(iele=0; iele<numele; ++iele) {

       // ... get E1, E2, Phi of layer & cframe

       double * layData = packedEset[iele]->getCompositeData(lay);
       double * cFrame  = packedEset[iele]->getCompositeFrame();

       if ( ! layData || ! cFrame ) continue;

       double E1  = layData[0];
       double E2  = layData[1];
       double Phi = layData[8];

       double length1 = crossScale * E1 / myMax(E1,E2) / 2.0;
       double length2 = crossScale * E2 / myMax(E1,E2) / 2.0;

       // .... Adjust Phi

       if ( Phi < 0.0 || Phi > 360.0 ) {
         int irot = int(Phi / 360.0);
         Phi = Phi - double (irot) * 360.0;
         if ( Phi < 0.0 ) Phi = 360.0 - Phi;
       }
       Phi = (3.141593*Phi)/180.0;

       // .... determine cross points in cframe

       double * cP1 = new double[3]; double * cP2 = new double[3];
       double * cP3 = new double[3]; double * cP4 = new double[3];

       cP1[0] =  length1 * cos(Phi);
       cP1[1] =  length1 * sin(Phi);
       cP1[2] =  0.0;

       cP2[0] = -length1 * cos(Phi);
       cP2[1] = -length1 * sin(Phi);
       cP2[2] =  0.0;

       cP3[0] =  length2 * cos(Phi+1.5708);
       cP3[1] =  length2 * sin(Phi+1.5708);
       cP3[2] =  0.0;

       cP4[0] = -length2 * cos(Phi+1.5708);
       cP4[1] = -length2 * sin(Phi+1.5708);
       cP4[2] =  0.0;

       // ... transform into global frame

       int n1 = iele*4; int n2 = n1+1; int n3 = n1+2; int n4 = n1+3;

       CPoint[lay][n1][0] = MidPoint[iele][0] + cFrame[0]*cP1[0] +  cFrame[3]*cP1[1] + cFrame[6]*cP1[2];
       CPoint[lay][n1][1] = MidPoint[iele][1] + cFrame[1]*cP1[0] +  cFrame[4]*cP1[1] + cFrame[7]*cP1[2];
       CPoint[lay][n1][2] = MidPoint[iele][2] + cFrame[2]*cP1[0] +  cFrame[5]*cP1[1] + cFrame[8]*cP1[2];

       CPoint[lay][n2][0] = MidPoint[iele][0] + cFrame[0]*cP2[0] +  cFrame[3]*cP2[1] + cFrame[6]*cP2[2];
       CPoint[lay][n2][1] = MidPoint[iele][1] + cFrame[1]*cP2[0] +  cFrame[4]*cP2[1] + cFrame[7]*cP2[2];
       CPoint[lay][n2][2] = MidPoint[iele][2] + cFrame[2]*cP2[0] +  cFrame[5]*cP2[1] + cFrame[8]*cP2[2];

       CPoint[lay][n3][0] = MidPoint[iele][0] + cFrame[0]*cP3[0] +  cFrame[3]*cP3[1] + cFrame[6]*cP3[2];
       CPoint[lay][n3][1] = MidPoint[iele][1] + cFrame[1]*cP3[0] +  cFrame[4]*cP3[1] + cFrame[7]*cP3[2];
       CPoint[lay][n3][2] = MidPoint[iele][2] + cFrame[2]*cP3[0] +  cFrame[5]*cP3[1] + cFrame[8]*cP3[2];

       CPoint[lay][n4][0] = MidPoint[iele][0] + cFrame[0]*cP4[0] +  cFrame[3]*cP4[1] + cFrame[6]*cP4[2];
       CPoint[lay][n4][1] = MidPoint[iele][1] + cFrame[1]*cP4[0] +  cFrame[4]*cP4[1] + cFrame[7]*cP4[2];
       CPoint[lay][n4][2] = MidPoint[iele][2] + cFrame[2]*cP4[0] +  cFrame[5]*cP4[1] + cFrame[8]*cP4[2];

       delete [] cP1;
       delete [] cP2;
       delete [] cP3;
       delete [] cP4;
    }

    // ... print pseudo structure

    filePrint(crossout,"Nodes cnodes%d\n",lay+1);
    for (icp=0;icp<numCPoints;icp++)
       filePrint(crossout,"%d\t % 14.6f\t% 14.6f\t % 14.6f\n",
                              icp+1,CPoint[lay][icp][0],CPoint[lay][icp][1],CPoint[lay][icp][2]);


    filePrint(crossout,"Elements cele%d using cnodes%d\n",lay+1,lay+1);
    for (iele=0;iele<numele;iele++) {
       filePrint(crossout,"%d 1 %d %d \n",iele*2+1,iele*4+1,iele*4+2);
       filePrint(crossout,"%d 1 %d %d \n",iele*2+2,iele*4+3,iele*4+4);
    }

    // ... close top-file for pseudo structure

    fclose(crossout);

    // ... print header

    filePrint(oinfo[iInfo].filptr,
            "Vector COMPOSIT_LAYER%d under %s for %s%d\n%d\n",
            lay+1,"Attributes","cnodes",lay+1,4*numele);

    // .... print zero deformation

    filePrint(oinfo[iInfo].filptr,"  %f\n",0.0);
    for (icp=0;icp<numCPoints;icp++)
       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                                    ,w,p,0.0,w,p,0.0,w,p,0.0);

    fflush(oinfo[iInfo].filptr);
  }
  else {

    filePrint(oinfo[iInfo].filptr,"  %f\n",time);

    int iele;
    for(iele=0; iele<numele; ++iele) {

       // ... get E1, E2, Phi of layer & cframe

       double * layData = packedEset[iele]->getCompositeData(lay);
       double * cFrame  = packedEset[iele]->getCompositeFrame();

       if ( ! layData || ! cFrame ) continue;

       double E1  = layData[0];
       double E2  = layData[1];
       double Phi = layData[8];

       double length1 = crossScale * E1 / myMax(E1,E2) / 2.0;
       double length2 = crossScale * E2 / myMax(E1,E2) / 2.0;

       // .... Adjust Phi

       if ( Phi < 0.0 || Phi > 360.0 ) {
         int irot = int(Phi / 360.0);
         Phi = Phi - double (irot) * 360.0;
         if ( Phi < 0.0 ) Phi = 360.0 - Phi;
       }
       Phi = (3.141593*Phi)/180.0;

       // .... determine cross points in cframe

       double * cP1 = new double[3]; double * cP2 = new double[3];
       double * cP3 = new double[3]; double * cP4 = new double[3];

       cP1[0] =  length1 * cos(Phi);
       cP1[1] =  length1 * sin(Phi);
       cP1[2] =  0.0;

       cP2[0] = -length1 * cos(Phi);
       cP2[1] = -length1 * sin(Phi);
       cP2[2] =  0.0;

       cP3[0] =  length2 * cos(Phi+1.5708);
       cP3[1] =  length2 * sin(Phi+1.5708);
       cP3[2] =  0.0;

       cP4[0] = -length2 * cos(Phi+1.5708);
       cP4[1] = -length2 * sin(Phi+1.5708);
       cP4[2] =  0.0;

       // ... transform into global frame

       int n1 = iele*4; int n2 = n1+1; int n3 = n1+2; int n4 = n1+3;

       double x,y,z;


       x = MidPoint[iele][0] + cFrame[0]*cP1[0] +  cFrame[3]*cP1[1] + cFrame[6]*cP1[2] - CPoint[lay][n1][0];
       y = MidPoint[iele][1] + cFrame[1]*cP1[0] +  cFrame[4]*cP1[1] + cFrame[7]*cP1[2] - CPoint[lay][n1][1];
       z = MidPoint[iele][2] + cFrame[2]*cP1[0] +  cFrame[5]*cP1[1] + cFrame[8]*cP1[2] - CPoint[lay][n1][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       x = MidPoint[iele][0] + cFrame[0]*cP2[0] +  cFrame[3]*cP2[1] + cFrame[6]*cP2[2] - CPoint[lay][n2][0];
       y = MidPoint[iele][1] + cFrame[1]*cP2[0] +  cFrame[4]*cP2[1] + cFrame[7]*cP2[2] - CPoint[lay][n2][1];
       z = MidPoint[iele][2] + cFrame[2]*cP2[0] +  cFrame[5]*cP2[1] + cFrame[8]*cP2[2] - CPoint[lay][n2][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       x = MidPoint[iele][0] + cFrame[0]*cP3[0] +  cFrame[3]*cP3[1] + cFrame[6]*cP3[2] - CPoint[lay][n3][0];
       y = MidPoint[iele][1] + cFrame[1]*cP3[0] +  cFrame[4]*cP3[1] + cFrame[7]*cP3[2] - CPoint[lay][n3][1];
       z = MidPoint[iele][2] + cFrame[2]*cP3[0] +  cFrame[5]*cP3[1] + cFrame[8]*cP3[2] - CPoint[lay][n3][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       x = MidPoint[iele][0] + cFrame[0]*cP4[0] +  cFrame[3]*cP4[1] + cFrame[6]*cP4[2] - CPoint[lay][n4][0];
       y = MidPoint[iele][1] + cFrame[1]*cP4[0] +  cFrame[4]*cP4[1] + cFrame[7]*cP4[2] - CPoint[lay][n4][1];
       z = MidPoint[iele][2] + cFrame[2]*cP4[0] +  cFrame[5]*cP4[1] + cFrame[8]*cP4[2] - CPoint[lay][n4][2];

       filePrint(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",w,p,x,w,p,y,w,p,z);

       fflush(oinfo[iInfo].filptr);

       delete [] cP1;
       delete [] cP2;
       delete [] cP3;
       delete [] cP4;
    }
  }
}

//------------------------------------------------------------------------------

char * Domain::getBasename(char * fname)
{
   char * bname;

   int    ssize  = strlen(fname);
   char * fpoint = strstr(fname,".");
   int    psize  = fpoint - fname;

   if ( psize < ssize )
   {
     bname  = new char[psize+1];
     int i;
     for (i=0;i<psize;i++) bname[i]=fname[i];
     bname[psize] = 0;
   }
   else
   {
     bname = fname;
   }

   return bname;
}

// New functions added for Helmholtz-Axisymmetric

Connectivity *
Domain::getNodeToNode() {
 if(nodeToNode) return nodeToNode;
 if(elemToNode == 0)
   elemToNode = new Connectivity(&packedEset) ;
 if(nodeToElem == 0)
   nodeToElem = elemToNode->reverse();
 nodeToNode = nodeToElem->transcon(elemToNode);
 return nodeToNode;
}


ConstrainedDSA *
Domain::makeCDSA(int nbc, BCond *bcs) {
  //MODIFIED FOR HEV PROBLEM, EC, 20070820
  //if (c_dsa) return c_dsa;

  if(!c_dsa) c_dsa = new ConstrainedDSA(*dsa, nbc, bcs);
  if(solInfo().HEV) {
    fprintf(stderr," *** c_dsaFluid NOT built in this version of makeCDSA ***\n");
    //if(!c_dsaFluid) c_dsaFluid = new ConstrainedDSA(*dsaFluid, nbcFluid, bcs);
  }

  return c_dsa;
}

void Domain::getNormal2D(int node1, int node2, double &nx, double &ny) {

 int nconele = nodeToElem->num(node1);

 int icon;

 double r=0;
 double z=0;

 for(icon=0; icon<nconele; ++icon) {
   int iele = (*nodeToElem)[node1][icon];
   int NodesPerElement = elemToNode->num(iele);
   int k;
   int test=0;

   r=0;
   z=0;

   for (k=0; k<NodesPerElement; k++) {
     if ((*elemToNode)[iele][k]==node2)
        test = 1;

     Node nd = nodes.getNode((*elemToNode)[iele][k]);

     r += nd.x/NodesPerElement;
     z += nd.y/NodesPerElement;
   }

   if (test==1)
      break;
 }
 Node nd = nodes.getNode(node1);

 double x1 = nd.x;
 double y1 = nd.y;

 nd = nodes.getNode(node2);

 double x2 = nd.x;
 double y2 = nd.y;

 nx = -(y2-y1);
 ny = x2-x1;

 double dot= nx*(r-0.5*(x2+x1))+ny*(z-0.5*(y1+y2));

 if (dot>=0) {
    nx = -nx;
    ny = -ny;
 }

}

void Domain::computeTDProps()
{
  if((numYMTT > 0) || (numCTETT > 0)) {
    // PJDS: used to calculate temperature dependent material properties
    //double defaultTemp = -10000000.0;

    // compute maximum number of nodes per element
    // note this is also done in Domain::makeAllDOFs()
    maxNumNodes = 0;
    int iele;
    for(iele=0; iele < numele; ++iele) {
      int numNodesPerElement = packedEset[iele]->numNodes();
      maxNumNodes = myMax(maxNumNodes, numNodesPerElement);
    }

    int i;
    // create id mapping for ymtt
    int maxid = 0;
    for(i = 0; i < numYMTT; ++i)
      if(ymtt[i]->getID() > maxid) maxid = ymtt[i]->getID();
    int *ymttmap = new int[maxid + 1];
    for(i = 0; i < numYMTT; ++i) ymttmap[ymtt[i]->getID()] = i;

    // create id mapping for ctett
    maxid = 0;
    for(i = 0; i < numCTETT; ++i)
      if(ctett[i]->getID() > maxid) maxid = ctett[i]->getID();
    int *ctettmap = new int[maxid + 1];
    for(i = 0; i < numCTETT; ++i) ctettmap[ymtt[i]->getID()] = i;

    // compute average temperatures and material properties
    int *nodeNumbers = new int[maxNumNodes];
    Vector elemNodeTemps(maxNumNodes);
    elemNodeTemps.zero();
    double *nodalTemperatures = getNodalTemperatures();
    for(iele = 0; iele < numele; ++iele) {
      // note: packedEset[iele]->numNodes() > 2 is temp fix to avoid springs
      // this means that until fixed, beams can't have temp-dependent material props
      if((packedEset[iele]->numNodes() > 2) && !packedEset[iele]->isPhantomElement()) {
        if((packedEset[iele]->getProperty()->E < 0) ||
           (packedEset[iele]->getProperty()->W < 0)) { // iele has temp-dependent E or W
          int NodesPerElement = packedEset[iele]->numNodes();
          packedEset[iele]->nodes(nodeNumbers);

          // compute average temperature in element
          double avTemp = 0.0;
          int iNode;
          for(iNode = 0; iNode < NodesPerElement; ++iNode) {
            if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
              elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
            else
              elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
            avTemp += elemNodeTemps[iNode];
          }
          avTemp /= NodesPerElement;
          // if(avTemp > 0) cerr << "element = " << iele << ", avTemp = " << avTemp << endl;

          StructProp *newProp = new StructProp(*packedEset[iele]->getProperty());
          // compute E using interp table
          if(packedEset[iele]->getProperty()->E < 0) {
            int id = (int) -packedEset[iele]->getProperty()->E;
            newProp->E = ymtt[ymttmap[id]]->getValAlt(avTemp);
          }

          // compute coeff of thermal expansion using interp table
          if(packedEset[iele]->getProperty()->W < 0) {
            int id  = (int) -packedEset[iele]->getProperty()->W;
            newProp->W = ctett[ctettmap[id]]->getValAlt(avTemp);
          }
          packedEset[iele]->setProp(newProp);
        }
      }
    }
    delete [] nodeNumbers;
    delete [] nodalTemperatures;
  }
}

// ************************************************************************
// HB: Temporary methods for testing tied Mortar condition by generating
//     the LMPC from the Mortar tied conditions
// ************************************************************************
//                           EXPERIMENTAL
// ************************************************************************
// HB: Add Surface Entity to the SurfEntities array
int
Domain::AddSurfaceEntity(SurfaceEntity* _SurfEntity)
{ //--- Verify if _SurfEntity was already defined
 int i=0;
 while (i<nSurfEntity &&
        SurfEntities[i]->ID()!=_SurfEntity->ID())
        i++;
 // if _SurfEntity not previously defined create new
 if (i==nSurfEntity)
  SurfEntities[nSurfEntity++] = _SurfEntity;
 else 
  if(_SurfEntity->ID() != 0) filePrint(stderr," *** WARNING: Surface Entity of Id %d has already been defined !!!\n", SurfEntities[i]->ID());

 return 0;
}

//HB: Add Surface Entity at the position isurf in the SurfEntities array
int
Domain::AddSurfaceEntity(SurfaceEntity* _SurfEntity, int isurf)
{
 //--- Verify if _SurfEntity was already defined
 if(SurfEntities[isurf]) {
   if(SurfEntities[isurf]->ID()!=_SurfEntity->ID()) {
     filePrint(stderr," *** ERROR in Domain::AddSurfaceEntity(SurfaceEntity*, int):\n");
     filePrint(stderr," -> try to overwrite the surface stored at index %d in the SurfEntities array by a different surface\n",isurf);
     exit(-1);
    } else {
      filePrint(stderr," *** WARNING in Domain::AddSurfaceEntity(SurfaceEntity*, int):\n");
      filePrint(stderr," -> surface %d already exist in the SurfEntities array at position %d\n",isurf);
    }
 } else {
   SurfEntities[isurf] = _SurfEntity;
 }
 return 0;
}

// HB: Print informations about the Surface Entity
//     in the SurfEntities array
void Domain::PrintSurfaceEntities()
{
 if(nSurfEntity>0) {
   filePrint(stderr," ... Model has %2d Surface Entities  ...\n",nSurfEntity);
#ifdef MORTAR_DEBUG
   for(int i=0; i<nSurfEntity; i++)
      SurfEntities[i]->Print();
#endif
 } else {
   filePrint(stderr," ... No Surface Entity in the pb    ...\n");
 }
}

// Add a wet surface Id.
int 
Domain::AddAeroEmbedSurfaceId(int Id)
{
  aeroEmbeddedSurfaceId.insert(Id);
  return 0;
}
// HB: Add Mortar Conditions to the MortarConds array
// Warning: CURRENTLY NO CHECK FOR MULTIPLE DEFINITION OF THE SAME MORTAR
//          CONDITION
int
Domain::AddMortarCond(MortarHandler* _MortarCond)
{
  // By default, set the  MortarCond Id = the order in which
  //we read them (see no check !!)
  _MortarCond->SetId(nMortarCond);

  // PJSA dirty fix for self contact
  if(_MortarCond->GetMasterEntityId() == _MortarCond->GetSlaveEntityId()) {
    SurfaceEntity *dummy = new SurfaceEntity(0); 
    AddSurfaceEntity(dummy);
    _MortarCond->SetMasterEntityId(0);
    _MortarCond->SetSelfContact(true);
  }

  // Add the MortarCond to the MortarConds array
  MortarConds[nMortarCond] = _MortarCond;

  // Count number of MortarCond
  nMortarCond++;

  if(_MortarCond->GetInteractionType() == MortarHandler::CTC) nContactSurfacePairs++;

  return 0;
}

void
Domain::DeleteMortarConds()
{
  if(nMortarCond) {
    for(int i=0; i<nMortarCond; ++i)
      if(MortarConds[i]) delete MortarConds[i];
  }
  nMortarCond = 0;
}

// HB: Print informations about Mortar Conditions
//     in the MortarConds array
void
Domain::PrintMortarConds()
{
  if(nMortarCond>0){
   filePrint(stderr," ... Model has %2d Mortar conditions ...\n",nMortarCond);
#ifdef MORTAR_DEBUG
  for(int i=0; i<nMortarCond; i++)
    MortarConds[i]->Print();
#endif
  } else
    filePrint(stderr," ... No Mortar condition in the pb  ...\n");

}

// HB: Set Ptr to Surface entities in the MortarCond objects
void
Domain::SetMortarPairing()
{
  if(nMortarCond>0){
#ifdef MORTAR_DEBUG
  filePrint(stderr," ... Set Ptr to Surface entities in the MortarCond objects ...\n");
#endif
  std::map<int, SurfaceEntity*> SurfIdToPtrSurfMap;
  for(int iSurf=0; iSurf<nSurfEntity; iSurf++){
    SurfIdToPtrSurfMap[SurfEntities[iSurf]->ID()] = SurfEntities[iSurf];
    //SurfIdToPtrSurfMap.insert(map<int, SurfaceEntity*>::value_type(SurfId, SurfEntities[iSurf]));
  }

  for(int iMortar=0; iMortar<nMortarCond; iMortar++){
    MortarConds[iMortar]->SetPtrSurfaceEntity(SurfIdToPtrSurfMap[MortarConds[iMortar]->GetMasterEntityId()],
                                              SurfIdToPtrSurfMap[MortarConds[iMortar]->GetSlaveEntityId()]);

#ifdef MORTAR_DEBUG
    // debug: check pairing
    filePrint(stderr," -> Pairing master surf. Id %d with slave surf. Id %d in Mortar cond Id %d\n",
              MortarConds[iMortar]->GetMasterEntityId(),MortarConds[iMortar]->GetSlaveEntityId(),MortarConds[iMortar]->ID());
    filePrint(stderr," -> check throught ptr: master surf Id = %d, slave surf Id = %d\n",
              MortarConds[iMortar]->GetPtrMasterEntity()->ID(),
              MortarConds[iMortar]->GetPtrSlaveEntity()->ID());
#endif
  }
  }
}

// HB: to setup internal data & renumber surfaces
void Domain::SetUpSurfaces(CoordSet* cs)
{
#if defined(MORTAR_LOCALNUMBERING) && defined(MORTAR_DEBUG)
  // if we want to work with local numbering of the surface entities (Salinas)
  if(nSurfEntity) filePrint(stderr," ... Use local numbering (and local nodeset) in the surface entities\n");
#endif
  for(int iSurf=0; iSurf<nSurfEntity; iSurf++){

    // For acme shells it is necessary to include both sides of the element in the face block
    // currently we only use acme shells when we are using the tdenforcement module of acme
    // to compute the contact forces. In other cases we don't use acme shells since acme doesn't
    // support the extraction of interactions for shells
    if(SurfEntities[iSurf]->GetIsShellFace() && tdenforceFlag()) {
      int nFaceElems = SurfEntities[iSurf]->GetnFaceElems();
      for(int i = 0; i < nFaceElems; ++i) {
        int etype = SurfEntities[iSurf]->GetFaceElemSet()[i]->GetFaceElemType();
        if(etype != 1 && etype != 3) {
          cerr << " *** ERROR: Surface element type " << etype << " not supported with SHELL_THICKNESS option\n";
          exit(-1);
        }
        int nNodes = SurfEntities[iSurf]->GetFaceElemSet()[i]->nNodes();
        int *nodes = new int[nNodes];
        for(int j=0; j<nNodes; ++j) nodes[nNodes-1-j] = SurfEntities[iSurf]->GetFaceElemSet()[i]->GetNode(j);
        SurfEntities[iSurf]->AddFaceElement(nFaceElems+i, etype, nNodes, nodes);
        delete [] nodes;
      }
    }

#ifdef MORTAR_DEBUG
    //filePrint(stderr," ------------------------------------------------------------------------\n");
    //filePrint(stderr,"  average normal of face element of surface %2d BEFORE local renumbering\n",SurfEntities[iSurf]->ID());
    //filePrint(stderr," ------------------------------------------------------------------------\n");
    //if(cs) SurfEntities[iSurf]->PrintFaceNormal(cs);
#endif
    SurfEntities[iSurf]->SetUpData(cs);
#ifdef MORTAR_LOCALNUMBERING
    SurfEntities[iSurf]->Renumber();
#endif
#ifdef MORTAR_DEBUG
  #ifdef HB_NODALNORMAL
    SurfEntities[iSurf]->ComputeNodalNormals(SurfEntities[iSurf]->GetNodeSet());
    SurfEntities[iSurf]->PrintNodalNormals();
  #endif
    filePrint(stderr," ------------------------------------------------------------------------\n");
    if(SurfEntities[iSurf]->IsRenumbered())
      filePrint(stderr,"  average normal of face element of surface %2d AFTER local renumbering\n",SurfEntities[iSurf]->ID());
    else
      filePrint(stderr,"  average normal of face element of surface %2d AFTER\n",SurfEntities[iSurf]->ID());
    filePrint(stderr," ------------------------------------------------------------------------\n");
    SurfEntities[iSurf]->PrintFaceNormal(SurfEntities[iSurf]->GetNodeSet());
    SurfEntities[iSurf]->Print();
#endif
  }
}

// PJSA 5-30-08 *********************************************************************************************************
// These functions use ACME's dynamic 2-configuation search algorithm and contact enforcement model for explicit dynamics
void Domain::InitializeDynamicContactSearch(int numSub, SubDomain **sd)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    CurrentMortarCond->SetDistAcme(sinfo.dist_acme);
    CurrentMortarCond->build_search(true, numSub, sd);
    CurrentMortarCond->build_td_enforcement();
    CurrentMortarCond->set_search_data(1); // interaction_type = 1 (NodeFace) 
    CurrentMortarCond->SetNoSecondary(solInfo().no_secondary);
    CurrentMortarCond->set_search_options();
    if(numSub == 0 && !sd) CurrentMortarCond->set_node_constraints(numDirichlet, dbc);
    else CurrentMortarCond->set_node_constraints(numSub, sd);
  }
}

void Domain::RemoveGap(Vector &g)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    CurrentMortarCond->remove_gap(g);
  }
}

void Domain::UpdateSurfaces(GeomState *geomState, int config_type) // config_type = 1 for current, 2 for predicted
{
  for(int iSurf=0; iSurf<nSurfEntity; iSurf++) {
    SurfEntities[iSurf]->UpdateNodeData(geomState);
  }
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    CurrentMortarCond->set_node_configuration(config_type);
  }
}

void Domain::UpdateSurfaces(DistrGeomState *geomState, int config_type, SubDomain **sd) // config_type = 1 for current, 2 for predicted
{
  if(solInfo().dist_acme != 2) {
    for(int iSurf=0; iSurf<nSurfEntity; iSurf++) {
      SurfEntities[iSurf]->UpdateNodeData(geomState, sd);
    }
  }
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(solInfo().dist_acme == 2) {
      CurrentMortarCond->GetPtrMasterEntity()->UpdateNodeData(geomState, sd);
      CurrentMortarCond->GetPtrSlaveEntity()->UpdateNodeData(geomState, sd);
    }
    CurrentMortarCond->set_node_configuration(config_type, geomState->getNumSub(), sd);
  }
}

void Domain::PerformDynamicContactSearch(double dt_old, double dt)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    int search_algorithm = 4; // augmented dynamic 2-configuration
    double dt_old = dt;
    CurrentMortarCond->perform_search(search_algorithm, dt_old, dt);
  }
}

void Domain::AddContactForces(double dt_old, double dt, Vector &f)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->compute_td_contact_force(dt_old, dt, f);
    }
  }
}

void Domain::AddContactForces(double dt_old, double dt, DistrVector &f)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->compute_td_contact_force(dt_old, dt, f);
    }
  }
}

void Domain::MakeNodalMass(SparseMatrix *M, SparseMatrix *Mcc)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->make_nodal_mass(M, c_dsa, Mcc);
      CurrentMortarCond->make_kinematic_partitioning(packedEset, nodeToElem);
    }
  }
}

#include <Paral.d/SubDOp.h>
void Domain::MakeNodalMass(SubDOp *M, SubDomain **sd)
{
  for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC || CurrentMortarCond->GetInteractionType() == MortarHandler::TIED) {
      CurrentMortarCond->make_nodal_mass(M, sd);
      CurrentMortarCond->make_kinematic_partitioning(M->getNumSub(), sd);
    }
  }
}
// **********************************************************************************************************************

// These functions use ACME's static 1-configuation search algorithm and face-face interaction, and Henri's mortar LMPCs for statics and implicit dynamics
void Domain::InitializeStaticContactSearch(MortarHandler::Interaction_Type t, int numSub, SubDomain **sd)
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      CurrentMortarCond->SetDistAcme(sinfo.dist_acme);
      CurrentMortarCond->SetMortarScaling(sinfo.mortar_scaling);
      CurrentMortarCond->SetMortarIntegrationRule(sinfo.mortar_integration_rule);
      CurrentMortarCond->build_search(false, numSub, sd);
      CurrentMortarCond->set_search_data(4); // interaction_type = 4 (FaceFace) 
      CurrentMortarCond->set_node_configuration(1);
      CurrentMortarCond->SetNoSecondary(solInfo().no_secondary);
      CurrentMortarCond->set_search_options();
      if(numSub == 0 && !sd) CurrentMortarCond->set_node_constraints(numDirichlet, dbc);
      else {
        CurrentMortarCond->set_node_constraints(numSub, sd);
        CurrentMortarCond->make_share(numSub, sd);
      }
    }
  }
}

void Domain::UpdateSurfaces(MortarHandler::Interaction_Type t, GeomState *geomState)
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      CurrentMortarCond->GetPtrMasterEntity()->UpdateNodeData(geomState);
      CurrentMortarCond->GetPtrSlaveEntity()->UpdateNodeData(geomState);
      CurrentMortarCond->set_node_configuration(1); // 1 --> config_type current
    }
  }
}

void Domain::UpdateSurfaces(MortarHandler::Interaction_Type t, DistrGeomState *geomState, SubDomain **sd) 
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      CurrentMortarCond->GetPtrMasterEntity()->UpdateNodeData(geomState, sd);
      CurrentMortarCond->GetPtrSlaveEntity()->UpdateNodeData(geomState, sd);
      CurrentMortarCond->set_node_configuration(1, geomState->getNumSub(), sd); // 1 --> config_type current
    }
  }
}

void Domain::PerformStaticContactSearch(MortarHandler::Interaction_Type t)
{
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
      int search_algorithm = 1; // static 1-configuration
      CurrentMortarCond->perform_search(search_algorithm);
      int interaction_type = 4;
      CurrentMortarCond->get_interactions(interaction_type);
    }
  }
}

void Domain::ExpComputeMortarLMPC(MortarHandler::Interaction_Type t, int nDofs, int *dofs)
{
  int num_interactions = 0;
  for(int iMortar = 0; iMortar < nMortarCond; iMortar++) {
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
    if(CurrentMortarCond->GetInteractionType() == t) {
#if (MAX_MORTAR_DERIVATIVES > 0)
      if(CurrentMortarCond->GetInteractionType() == MortarHandler::CTC && !tdenforceFlag())
        CurrentMortarCond->CreateFFIPolygon<ActiveDouble>();
      else
#endif
      CurrentMortarCond->CreateFFIPolygon<double>();
      CurrentMortarCond->AddMortarLMPCs(&lmpc, numLMPC, numCTC, nDofs, dofs);
      nMortarLMPCs += CurrentMortarCond->GetnMortarLMPCs();
      num_interactions += CurrentMortarCond->GetnFFI();
      CurrentMortarCond->DeleteFFIData();
    }
  }
  if(verboseFlag) filePrint(stderr," ... Built %d Mortar Surface/Surface Interactions ...\n", nMortarLMPCs);
#ifdef HB_ACME_FFI_DEBUG
  if(sinfo.ffi_debug && num_interactions > 0) {
    char fname[16];
    sprintf(fname,"FFI.top.%d",totalNewtonIter);
    filePrint(stderr," -> Write FFI top file: %s for %d interactions\n", fname, num_interactions);
    FILE* FFITopFile = fopen(fname,"w");
    WriteFFITopFile(FFITopFile);
    fclose(FFITopFile);
  }
#endif
}

// HB: compute & add to the standard LMPC array (lmpc) the Mortar LMPCs:
void Domain::ComputeMortarLMPC(int nDofs, int *dofs)
{
  if(nMortarCond>0){
  double time00=-getTime();
  for(int iMortar=0; iMortar<nMortarCond; iMortar++){
    MortarHandler* CurrentMortarCond = MortarConds[iMortar];
#ifdef MORTAR_TIMINGS
    double time0= -getTime();
    double time = -getTime();
#endif
#ifdef MORTAR_DEBUG
    filePrint(stderr," ... Treat Mortar condition Id %d\n",CurrentMortarCond->ID());
    CurrentMortarCond->Print();
#endif
    switch(int(CurrentMortarCond->GetGeomType())){
     case MortarHandler::EQUIVALENCED:
       CurrentMortarCond->CreateACMEFFIData();
       break;
     case MortarHandler::NON_MATCHING:
     default:
       CurrentMortarCond->PerformACMEFFISearch();
       break;
    }
#ifdef MORTAR_TIMINGS
    time += getTime();
    filePrint(stderr," -> time spent in the ACME search: %e s\n",time/1000);
    time = -getTime();
#endif
#if (MAX_MORTAR_DERIVATIVES > 0)
    if(CurrentMortarCond->GetInteractionType()==MortarHandler::CTC && !tdenforceFlag())
      CurrentMortarCond->CreateFFIPolygon<ActiveDouble>();
    else
#endif
    CurrentMortarCond->CreateFFIPolygon<double>();
#ifdef MORTAR_TIMINGS
    time += getTime();
    filePrint(stderr," -> time spent in building the Mortar B matrices: %e s\n",time/1000);
#endif
    if(CurrentMortarCond->GetInteractionType()!=MortarHandler::FSI){
#ifdef MORTAR_TIMINGS
      time = -getTime();
#endif
      CurrentMortarCond->AddMortarLMPCs(&lmpc, numLMPC, numCTC, nDofs, dofs);
#ifdef MORTAR_TIMINGS
      time += getTime();
      time0 += getTime();
      filePrint(stderr," -> time spent in creating the Mortar LMPCs: %e s\n",time/1000);
      filePrint(stderr," -> total time spent in building those Mortar LMPCs: %e s\n",time0/1000);
#endif
      //printLMPC();
      // count the total number of Mortar LMPCs generated
      nMortarLMPCs += CurrentMortarCond->GetnMortarLMPCs();
    } else
      CurrentMortarCond->AddWetFSI(&fsi, numFSI);

// CONTACT DEBUG    CurrentMortarCond->DeleteFFIData();
  }

  time00 += getTime();
  if(verboseFlag) filePrint(stderr," ... Built %d Mortar Surface/Surface Interactions ...\n", nMortarLMPCs+numFSI);
#ifdef MORTAR_TIMINGS
  filePrint(stderr," ... CPU time for building mortar surface/surface interactions: %e s\n",time00/1000);
#endif
  //if(numFSI){
  //  printFSI();
  //  long memFSI = 0;
  //  for(int iFSI=0; iFSI<numFSI; iFSI++)
  //    memFSI += fsi[iFSI]->mem();
  //  filePrint(stderr," ### memnory used by fsi array %14.3f Mb \n",memFSI/(1024.*1024.));
  //}
 }
#ifdef HB_ACME_FFI_DEBUG
 filePrint(stderr," -> Write FFI top file: FFI.top\n");
 FILE* FFITopFile = fopen("FFI.top","w");
 WriteFFITopFile(FFITopFile); //valid only if cs!=0 ...
 fclose(FFITopFile);
#endif
}

#ifdef HB_ACME_FFI_DEBUG
void Domain::WriteFFITopFile(FILE* file)
{
  if(nMortarCond){
    int firstNodeId = 1;
    fprintf(file,"Nodes FFINodes\n");
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* MortarCond = MortarConds[iMortar];
      MortarCond->PrintFFIPolyVertices(file, firstNodeId);
    }
    firstNodeId = 1;
    int firstElemId = 1;
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* MortarCond = MortarConds[iMortar];
      fprintf(file,"Elements SlaveFFI_%d using FFINodes\n",MortarCond->GetSlaveEntityId());
      MortarCond->PrintFFIPolyTopo(file, firstElemId, firstNodeId);
      fprintf(file,"Elements MasterFFI_%d using FFINodes\n",MortarCond->GetMasterEntityId());
      MortarCond->PrintFFIPolyTopo(file, firstElemId, firstNodeId);
    }
  }
}
#endif

// HB: create (Tied) MortarToMPC connectivity
// WARNING: (1) WE ASSUME THAT THE MORTAR LMPCs ARE AFTER THE REGULAR LMPCs
//              IN THE STANDARD LMPC ARRAY (lmpc)
//          (2) WE ASSUME THAT IN EACH MORTAR CONDITION, THE MORTAR LMPCs
//              ARE NUMBERED SEQUENTIALY
// We currently do not separate a Tied Mortar LMPC into 2 or 3 (2D/3D)
// Tied Mortar LMPC (one for each axis: x, y, z)
// this splitting can be done at the CCt level (see component renumbering)
// DONE: changed the order in which the Mortar LMPCs are added in the
// standard LMPC array (lmpc):
// before    [(x,y,z),...,(x,y,z)]
// currently [(x,...,x),(y,...,y),(z,...,z)]
void
Domain::CreateMortarToMPC()
{
  int* pointermap = 0;
  int* target     = 0;

  if(mortarToMPC) { delete mortarToMPC; mortarToMPC = 0; }

  if(nMortarCond>0){
#ifdef MORTAR_DEBUG
    filePrint(stderr," ... Create (Tied) Mortar To LMPCs connectivity ...\n");
#endif

    // 1) Set the map pointer array & count number of
    //    targets = total number of (Tied) Mortar LMPCs
    pointermap = new int[nMortarCond+1];
    int ntargets = 0;
    pointermap[0] = 0;
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* CurrentMortarCond = MortarConds[iMortar];
      int nLMPCs = CurrentMortarCond->GetnMortarLMPCs();
      pointermap[iMortar+1] = pointermap[iMortar] + nLMPCs;
      ntargets += nLMPCs;
    }

    // 2) Fill the target (map) array
    target = new int[ntargets];
    for(int iMortar=0; iMortar<nMortarCond; iMortar++){
      MortarHandler* CurrentMortarCond = MortarConds[iMortar];
      int nLMPCs       = CurrentMortarCond->GetnMortarLMPCs();
      int IdFirstLMPC  = CurrentMortarCond->GetIdFirstMortarLMPC();
      for(int i=0; i<nLMPCs; i++)
        target[pointermap[iMortar]+i] = IdFirstLMPC+i;
    }

    mortarToMPC = new Connectivity(nMortarCond, pointermap, target);
    //mortarToMPC->print();
  }
}
// HB: return the Mortar to LMPCs connectivity and the total number of
// Mortar LMPCs
Connectivity*
Domain::GetMortarToMPC() { return mortarToMPC; }
int
Domain::GetnMortarLMPCs() { return nMortarLMPCs; }

// HB: write in a file ALL the generated Mortar LMPCs
// make the assumptions that the Mortar LMPC id is NEGATIVE (<0)
// basicaly, loop over ALL the lmpc and look for the negative lmpc ids
// Another possibility will be to use the mortarToMPC connectivity
// !!! TO DO !!!
/*
void Domain::WriteToFileMortarLMPCs(FILE *file)
{

  int i,ilmpc, nmotarlmpc=0;
  int firstmortarlmpcId = numLMPC-nMortarLMPCs;
  if(firstmortarlmpcId<0) firstmortarlmpcId = 0;
  for(ilmpc=0; ilmpc<numLMPC; ilmpc++){
    int lmpcId = lmpc[ilpmc].lmpcnum;
    if(lmpcId<0){
     nmotarlmpc++;
     double rhs = lmpc[ilpmc].rhs;
     lmpcId = firstmortarlmpcId+nmotarlmpc;
     filePrint(file," %4d  %f\n",lmpcId,rhs);
     for(i=0;i<lmpc[ilpmc].nterms;i++){
       int node  = lmpc[ilpmc].terms[i].nnum;
       int dof   = lmpc[ilpmc].terms[i].dofnum;
       double val= lmpc[ilpmc].terms[i].val;
       filePrint(file,"    %4d  %d   %f\n",node,dof,val);
     }
    }
  }
  filePrint(stderr," -> has written %d motar lmpcs to file\n",nmotarlmpc);
}
*/

// ************************************************************************

// returns the value of the pressure force flag
int
Domain::pressureFlag() { return geoSource->pressureFlag(); }

// function that returns composite layer info
LayInfo *Domain::getLayerInfo(int num) { return geoSource->getLayerInfo(num); }



void
Domain::initialize()
{
 numdofs = 0; numDispDirichlet = 0; numContactPairs = 0;
 numIDis = 0; numIDisModal = 0; numIVel = 0; numIVelModal = 0; numDirichlet = 0; numNeuman = 0; numSommer = 0;
 numComplexDirichlet = 0; numComplexNeuman = 0; 
 firstDiMass = 0; numIDis6 = 0; gravityAcceleration = 0;
 allDOFs = 0; stress = 0; weight = 0; elstress = 0; elweight = 0; claw = 0;
 numLMPC = 0; numYMTT = 0; numCTETT = 0; MidPoint = 0; temprcvd = 0;
 heatflux = 0; elheatflux = 0; elTemp = 0; dbc = 0; nbc = 0; 
 iDis = 0; iDisModal = 0; iVel = 0; iVelModal = 0; iDis6 = 0; elemToNode = 0; nodeToElem = 0;
 nodeToNode = 0; dsa = 0; c_dsa = 0; cdbc = 0; cnbc = 0;
 dsaFluid = 0; c_dsaFluid = 0; allDOFsFluid = 0; dbcFluid = 0;
 elemToNodeFluid = 0; nodeToElemFluid = 0; nodeToNodeFluid = 0;
 nSurfEntity = 0; nMortarCond = 0; nMortarLMPCs= 0; mortarToMPC = 0;
 solver = 0; csolver = 0; mftval = 0; mptval = 0; hftval = 0;
 flExchanger = 0; outFile = 0; elDisp = 0; p_stress = 0; p_elstress = 0; stressAllElems = 0;
 previousExtForce = 0; previousAeroForce = 0; previousDisp = 0; previousCq = 0;
 temprcvd = 0; optinputfile = 0;
 nSurfEntity = 0; nMortarLMPCs = 0; mortarToMPC = 0; matrixTimers = 0;
 allDOFs = 0; haveNodes = false; nWetInterface = 0; wetInterfaces = 0;
 numFSI = 0; firstOutput = true; nodeToNode_sommer = 0; sowering = false; nDimass = 0;
 fluidDispSlosh = 0; elPotSlosh =0; elFluidDispSlosh =0; //ADDED FOR SLOSHING PROBLEM, EC, 20070723
 elFluidDispSloshAll =0; //ADDED FOR SLOSHING PROBLEM, EC, 20081101
 nodeToFsi = 0; //HB
 numCTC=0;
 output_match_in_top = false;//TG
 C_condensed = 0;
 nContactSurfacePairs = 0;
 outFlag = 0;
 nodeTable = 0;
 MpcDSA = 0; nodeToNodeDirect = 0;
 p = 0;
 g_dsa = 0;
}

Domain::~Domain()
{
 // delete &nodes;
 if(nodeToNode)   { delete nodeToNode;   nodeToNode=0;   }
 if(nodeToNodeFluid)   { delete nodeToNodeFluid;   nodeToNodeFluid=0;   } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(renumb.order) { delete [] renumb.order; renumb.order=0; }
 if(renumb.renum) { delete [] renumb.renum; renumb.renum=0; }
 if(renumb.xcomp) { delete [] renumb.xcomp; renumb.xcomp=0; }
 if(renumbFluid.order) { delete [] renumbFluid.order; renumbFluid.order=0; } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(renumbFluid.renum) { delete [] renumbFluid.renum; renumbFluid.renum=0; } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(gravityAcceleration) { delete [] gravityAcceleration; gravityAcceleration = 0; }
 if(matrixTimers) { delete matrixTimers; matrixTimers = 0; }
 if(allDOFs) { delete allDOFs; allDOFs = 0; }
 if(dsa) { delete dsa; dsa = 0; }
 if(c_dsa) { delete c_dsa; c_dsa = 0; }
 if(elemToNode) { delete elemToNode; elemToNode = 0; }
 if(nodeToElem) { delete nodeToElem; nodeToElem = 0; }
 if(allDOFsFluid) { delete allDOFsFluid; allDOFsFluid = 0; } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(dsaFluid) { delete dsaFluid; dsaFluid = 0; } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(c_dsaFluid) { delete c_dsaFluid; c_dsaFluid = 0; } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(elemToNodeFluid) { delete elemToNodeFluid; elemToNodeFluid = 0; } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(nodeToElemFluid) { delete nodeToElemFluid; nodeToElemFluid = 0; } //ADDED FOR HEV PROBLEM, EC, 20070820
 if(mortarToMPC) { delete mortarToMPC; mortarToMPC = 0; }
 if(dbc) { delete [] dbc; dbc = 0; }
 if(nbc) { delete [] nbc; nbc = 0; }
 //if(cvbc) { delete [] cvbc; cvbc = 0; }
 if(iDis) { delete [] iDis; iDis = 0; }
 if(iDisModal) { delete [] iDisModal; iDisModal = 0; }
 if(iVel) { delete [] iVel; iVel = 0; }
 if(iVelModal) { delete [] iVelModal; iVelModal = 0; }
 if(iDis6) { delete [] iDis6; iDis6 = 0; }
 if(cdbc) { delete [] cdbc; cdbc = 0; }
 if(cnbc) { delete [] cnbc; cnbc = 0; }
 if(stress) { delete stress; stress = 0; }
 if(weight) { delete weight; weight = 0; }
 if(elDisp) { delete elDisp; elDisp = 0; }
 if(elTemp) { delete elTemp; elTemp = 0; }
 if(fluidDispSlosh) { delete fluidDispSlosh; fluidDispSlosh = 0; } //ADDED FOR SLOSHING PROBLEM, EC, 20070723
 if(elPotSlosh) { delete elPotSlosh; elPotSlosh = 0; }  //ADDED FOR SLOSHING PROBLEM, EC, 20070723
 if(elFluidDispSlosh) { delete elFluidDispSlosh; elFluidDispSlosh = 0; } //ADDED FOR SLOSHING PROBLEM, EC, 20070723
 if(elstress) { delete elstress; elstress = 0; }
 if(elweight) { delete elweight; elweight = 0; }
 if(p_stress) { delete p_stress; p_stress = 0; }
 if(p_elstress) { delete p_elstress; p_elstress = 0; }
 if(stressAllElems) { delete stressAllElems; stressAllElems = 0; }
 if(claw) { delete claw; claw = 0; }
 //if(mftval) { delete mftval; mftval = 0; }
 //if(mptval) { delete mptval; mptval = 0; }
 //if(hftval) { delete hftval; hftval = 0; }
 if(firstDiMass) { delete firstDiMass; firstDiMass = 0; }
 if(previousExtForce) { delete previousExtForce; previousExtForce = 0; }
 if(previousAeroForce) { delete previousAeroForce; previousAeroForce = 0; }
 if(previousDisp) { delete previousDisp; previousDisp = 0; }
 if(previousCq) { delete previousCq; previousCq = 0; }
 if(heatflux) { delete heatflux; heatflux = 0; }
 if(elheatflux) { delete elheatflux; elheatflux = 0; }
 delete &nodes;
 if(wetInterfaces) { delete wetInterfaces; wetInterfaces = 0; }
 if(nodeToNode_sommer) { delete nodeToNode_sommer; nodeToNode_sommer = 0; }
 if(nodeToFsi) { delete nodeToFsi; nodeToFsi = 0; } //HB
 //HB: try to avoid memory leaks but trouble with ~BlockAlloc() ...
 for(int i=0; i<SurfEntities.max_size(); i++)
   if(SurfEntities[i]) { SurfEntities[i]->~SurfaceEntity(); SurfEntities[i] = 0; }
 for(int i=0; i<nMortarCond; ++i) // HB
   if(MortarConds[i]) delete MortarConds[i];
 for(int i=0; i<numLMPC; ++i) // HB
   if(lmpc[i]) delete lmpc[i];
 if(nodeTable) delete [] nodeTable;
 if(MpcDSA) delete MpcDSA; if(nodeToNodeDirect) delete nodeToNodeDirect;
 if(p) delete p;
 for(int i=0; i<contactSurfElems.size(); ++i)
   packedEset.deleteElem(contactSurfElems[i]);
 if(g_dsa) delete g_dsa;
}

#include <Element.d/Helm.d/HelmElement.h>

int
Domain::isFluidElement(int i)
{
 HelmElement *he = dynamic_cast<HelmElement *>(packedEset[i]);
 if(he==0) return 0; // non-Helmholtz element found
 else return he->isFluid();
}

int
Domain::isStructureElement(int i)
{
  return (!isFluidElement(i) && !packedEset[i]->isConstraintElement());
}

bool
Domain::isHomogeneous()
{
  if((geoSource->getNumProps() == 1) && (numYMTT == 0) && (numCTETT == 0)) return true;
  else return false;
}

void
Domain::addWetInterface(int _fluidSurfaceID, int _structureSurfaceID, double normal_tol, double tangential_tol)
{
  if(nWetInterface == 0) wetInterfaces = new ResizeArray<WetInterface *>(0, 1);  // initial size is 1
  WetInterface *newWI = new WetInterface();
  newWI->fluidSurfaceID = _fluidSurfaceID;
  newWI->structureSurfaceID = _structureSurfaceID;
  (*wetInterfaces)[nWetInterface++] = newWI;

  // PJSA: temporary fix, use mortar code to build LMPCs to make interface
  MortarHandler* _MortarCond = new MortarHandler(_structureSurfaceID, _fluidSurfaceID, normal_tol, tangential_tol);
  _MortarCond->SetId(nMortarCond);
  _MortarCond->SetMortarType(MortarHandler::STD); //HB: trigger std mortar basis
  _MortarCond->SetInteractionType(MortarHandler::FSI);
  //filePrint(stderr," _fluidSurfaceID = %d, _structureSurfaceID = %d\n",_fluidSurfaceID,_structureSurfaceID); //HB
  if(_fluidSurfaceID == _structureSurfaceID)
    _MortarCond->SetGeomType(MortarHandler::EQUIVALENCED);
  MortarConds[nMortarCond] = _MortarCond;
  //MortarConds[nMortarCond]->Print();
  nMortarCond++;

  return;
}

//HB: I think it may be more safe to get them from the fsi array in case ACME misses
//some interactions (non-matching wet interfaces) -> in which case only some of the nodes of the
//wet surfaces given in the input file are currently trully "wet interface nodes".
int *
Domain::getAllWetInterfaceNodes(int &count)
{
  count = 0;
#define USE_FSI_NODES
  if(numFSI) {
#ifdef USE_FSI_NODES //HB: assume the FSI has already been computed. to be validated ...
    int *nodeMap = new int[numnodes];
    for(int i=0; i<numnodes; ++i) nodeMap[i] = -1;
    for(int i=0; i<numFSI; i++) {
      if(nodeMap[fsi[i]->lmpcnum]<0) nodeMap[fsi[i]->lmpcnum] = count++; // fluid node
      for(int j=0; j<fsi[i]->nterms; j++) // Loops over structure's nodes
        if(nodeMap[fsi[i]->terms[j].nnum]<0) nodeMap[fsi[i]->terms[j].nnum] = count++;
    }
#else
  if(nWetInterface > 0) {
    std::map<int, SurfaceEntity*> SurfIdToPtrSurfMap;
    for(int i=0; i<nSurfEntity; i++){
      int SurfId = SurfEntities[i]->ID();
      SurfIdToPtrSurfMap[SurfId] = SurfEntities[i];
    }

    int *nodeMap = new int[numnodes];
    for(int i=0; i<numnodes; ++i) nodeMap[i] = -1;
    for(int i=0; i<nWetInterface; i++) {
      int fluidId = (*wetInterfaces)[i]->fluidSurfaceID;
      int structureId = (*wetInterfaces)[i]->structureSurfaceID;
      // fluid
      SurfaceEntity* fluidEntity = SurfIdToPtrSurfMap[fluidId];
      int nFluidElem = fluidEntity->nFaceElements();
      FaceElemSet *fluidElemSet = fluidEntity->GetPtrFaceElemSet();
      int *glFluidNodes = fluidEntity->GetPtrGlNodeIds();
      int nFluidNodes = fluidEntity->GetnNodes();
      for(int j=0; j<nFluidNodes; ++j)
        if(nodeMap[glFluidNodes[j]] == -1) nodeMap[glFluidNodes[j]] = count++;
      if(structureId == fluidId) continue;
      // structure
      SurfaceEntity* structureEntity  = SurfIdToPtrSurfMap[structureId] ;
      int *glStructureNodes = structureEntity->GetPtrGlNodeIds();
      int nStructureNodes = structureEntity->GetnNodes();
      for(int j=0; j<nStructureNodes; ++j)
        if(nodeMap[glStructureNodes[j]] == -1) nodeMap[glStructureNodes[j]] = count++;
      }
    }
#endif
    glWetNodeMap = nodeMap;
    int *allWetInterfaceNodes = new int[count];
    for(int i=0; i<numnodes; ++i)
      if(nodeMap[i] != -1) allWetInterfaceNodes[nodeMap[i]] = i;

    //if(nodeMap) delete [] nodeMap;
    return allWetInterfaceNodes;
  }
  else {
    //cerr << " *** WARNING: No Wet Interfaces have been defined \n";
    return 0;
  }
}

void Domain::printFSI(FILE* file)
{
 filePrint(file," ... FSI list :                    ...\n");
 int i;
 for(i=0; i<numFSI; i++) {
   filePrint(file," %4d ... fluid node %6d : \n", i+1,fsi[i]->lmpcnum+1);
   if(fsi[i]->isComplex) {
     filePrint(file," ... structure node   dof       coef\n");
     for(int j=0; j<fsi[i]->nterms; j++)
       filePrint(file,"        %6d            %d        (%6e,%6e)\n",
                 fsi[i]->terms[j].nnum+1, fsi[i]->terms[j].dofnum,
                 fsi[i]->terms[j].coef.c_value.real(), fsi[i]->terms[j].coef.c_value.imag());
   }
   else {
     filePrint(stderr," ... structure node   dof       coef\n");
     for(int j=0; j<fsi[i]->nterms; j++)
       filePrint(file,"        %6d            %d        %6e\n",
                 fsi[i]->terms[j].nnum+1, fsi[i]->terms[j].dofnum,
                 fsi[i]->terms[j].coef.r_value);
   }
 }
}

double
Domain::getFrequencyOrWavenumber()
{
  double ret = 0.0;
  if(geoSource->isShifted()) {
    if(domain->solInfo().doFreqSweep) {
      if(isCoarseGridSolve) {
        ret = domain->coarse_frequencies->front();
        ret = getSavedFreq();
      }
      else ret = domain->frequencies->front();
    }
    else {
      ret = geoSource->omega();
    }
    if(domain->solInfo().isAcousticHelm()) ret /= domain->fluidCelerity;
    else ret /= (2.0*PI);
  }
  return ret;
}

void
Domain::computeCoupledScaleFactors()
{
  // use global average properties to compute Lame constants & coupled scaling factor
  // (equation proposed by Jan Mandel)
  double ymod = geoSource->global_average_E;
  double prat = geoSource->global_average_nu;
  //double rhofavg = geoSource->global_average_rhof;
  double lame_lambda = prat*ymod/((1.0+prat)*(1.0-2.0*prat));
  double lame_mu_times2 = 2.0*ymod/(2.0*(1.0+prat));
  double lame_max = (lame_lambda > lame_mu_times2) ? lame_lambda : lame_mu_times2;
/* old scaling
  double rhoFomega2 = (geoSource->isShifted()) ? domain->fluidDensity * geoSource->shiftVal() : 1000.0;
  coupledScaling = ((lame_max > 0.0) && (rhoFomega2 > 0.0)) ? 1.0/sqrt(rhoFomega2*lame_max) : 1.0;
  coupledScaling *= domain->solInfo().coupled_scale; // PJSA 10-21-05 can be used to adjust scaling up or down
  cscale_factor = rhoFomega2*coupledScaling;
*/
  double omega2 = (geoSource->isShifted()) ? geoSource->shiftVal() : 1000.0;
  coupledScaling = ((lame_max > 0.0) && (geoSource->global_average_rhof> 0.0)) ?
               1.0/sqrt(geoSource->global_average_rhof*omega2*lame_max) : 1.0; // RADEK
  coupledScaling *= domain->solInfo().coupled_scale; // PJSA 10-21-05 can be used to adjust scaling up or down
  cscale_factor = omega2*coupledScaling; // RADEK

  cscale_factor2 = cscale_factor*coupledScaling;
  //cerr << "coupledScaling = " << coupledScaling << ", cscale_factor = " << cscale_factor
  //     << ", cscale_factor2 = " << cscale_factor2 << endl;
}

void
Domain::initSfem()
{
#ifndef SALINAS
  if(domain->solInfo().noninpc || domain->solInfo().inpc) {
 //   geoSource->printGroups();
    sfem->computeLP();
  }
#endif
}

// add nodal contact mpcs
int
Domain::addNodalCTC(int n1, int n2, double nx, double ny, double nz,
                    double normalGap, int _mode, int lagrangeMult, double penalty)
{
 // contact note: if normalGapPresent is false perhaps the default should be to compute the geometric gap
 // using the nodal coordinates
 int mode = (_mode > -1) ? _mode : domain->solInfo().contact_mode;  // 0 -> normal tied + tangents free, 1 -> normal contact + tangents free
                                                                    // 2 -> normal+tangents tied, 3 -> normal contact + tied tangents
/*CoordSet &cs = geoSource->GetNodes();
 double normal[3] = { cs[n2]->x - cs[n1]->x, cs[n2]->y - cs[n1]->y, cs[n2]->z - cs[n1]->z };
 double gap = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
 normalize(normal);
 cerr << n1+1 << " " << n2+1 << "  " << std::setprecision(12) << -normal[0] << " " << -normal[1] << " " << -normal[2] << " GAP " << std::setprecision(2) << -gap << endl;
*/
 int lmpcnum = 0;

 // normal constraint
 LMPCons *_CTC = new LMPCons(lmpcnum, normalGap);
 double norm = sqrt(nx*nx + ny*ny + nz*nz);
 if(nx != 0.0) {
   nx /= norm;
   LMPCTerm *term1 = new LMPCTerm(n1, 0, nx);
   _CTC->addterm(term1);
   LMPCTerm *term2 = new LMPCTerm(n2, 0, -nx);
   _CTC->addterm(term2);
 }
 if(ny != 0.0) {
   ny /= norm;
   LMPCTerm *term1 = new LMPCTerm(n1, 1, ny);
   _CTC->addterm(term1);
   LMPCTerm *term2 = new LMPCTerm(n2, 1, -ny);
   _CTC->addterm(term2);
 }
 if(nz != 0.0) {
   nz /= norm;
   LMPCTerm *term1 = new LMPCTerm(n1, 2, nz);
   _CTC->addterm(term1);
   LMPCTerm *term2 = new LMPCTerm(n2, 2, -nz);
   _CTC->addterm(term2);
 }
 _CTC->type = (mode == 1 || mode == 3); // this is to be phased out
 if(mode == 1 || mode == 3) _CTC->setType(mpc::Inequality);
 _CTC->setSource(mpc::NodalContact);
 addLMPC(_CTC,false);
 if(_CTC->type == 1) numCTC++; // inequality constraint

 // PJSA 7-12-2007 tangential constraints ... note this may lead to redundant constraints (singularity in CCt)
 // for 2D you need to set spacedim 2 in the input file, also it is currently assumed that 2D model is defined in XY plane
 if(mode == 2 || mode == 3) {
   // 1. select direction for initial cross product... use z unless normal is parallel or close to parallel to z axis, in which case use y
   double n[3] = { nx, ny, nz };
   double t1[3], t2[3];
   if(!(nx < 0.1 && ny < 0.1) || solInfo().fetiInfo.spaceDimension == 2) { t2[0] = 0.0; t2[1] = 0.0; t2[2] = 1.0; }
   else if(!(nx < 0.1 && nz < 0.1)) { t2[0] = 0.0; t2[1] = 1.0; t2[2] = 0.0; }
   else { t2[0] = 1.0; t2[1] = 0.0; t2[2] = 0.0; }
   crossprod(n,t2,t1);
   normalize(t1);
   LMPCons *_TGT1 = new LMPCons(0, 0.0);
   if(t1[0] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 0, t1[0]); _TGT1->addterm(term1);
                                    LMPCTerm *term2 = new LMPCTerm(n2, 0, -t1[0]); _TGT1->addterm(term2); }
   if(t1[1] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 1, t1[1]); _TGT1->addterm(term1);
                                    LMPCTerm *term2 = new LMPCTerm(n2, 1, -t1[1]); _TGT1->addterm(term2); }
   if(t1[2] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 2, t1[2]); _TGT1->addterm(term1);
                                    LMPCTerm *term2 = new LMPCTerm(n2, 2, -t1[2]); _TGT1->addterm(term2); }
   _TGT1->setSource(mpc::NodalContact);
   addLMPC(_TGT1,false);
   if(solInfo().fetiInfo.spaceDimension == 3) {
     crossprod(n,t1,t2);
     normalize(t2);
     LMPCons *_TGT2 = new LMPCons(0, 0.0);
     if(t2[0] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 0, t2[0]); _TGT2->addterm(term1);
                                      LMPCTerm *term2 = new LMPCTerm(n2, 0, -t2[0]); _TGT2->addterm(term2); }
     if(t2[1] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 1, t2[1]); _TGT2->addterm(term1);
                                      LMPCTerm *term2 = new LMPCTerm(n2, 1, -t2[1]); _TGT2->addterm(term2); }
     if(t2[2] != 0.0) { LMPCTerm *term1 = new LMPCTerm(n1, 2, t2[2]); _TGT2->addterm(term1);
                                      LMPCTerm *term2 = new LMPCTerm(n2, 2, -t2[2]); _TGT2->addterm(term2); }
     _TGT2->setSource(mpc::NodalContact);
     addLMPC(_TGT2,false);
   }
 }

 return numCTC;
}

/*
void
Domain::addNodeToNodeLMPCs(int lmpcnum, int n1, int n2, double face_normal[3], double gap_vector[3], int itype)
{
  // type = 0: generate lmpcs for tied nodes with gap vector using global cartesian coordinate frame
  //           note: addNodalCTC can be used to tie nodes in normal and two tangential directions, however only the normal gap is specified therefore tangential gaps are always zero
  // type = 1: generate lmpcs for normal contact
  if(itype == 0) { // tie
    for(int i=0; i<solInfo().fetiInfo.spaceDimension; ++i) {
      LMPCons *_MPC = new LMPCons(lmpcnum+i, gap_vector[i]);
      LMPCTerm *term1 = new LMPCTerm(n1, i, 1.0);
      _MPC->addterm(term1);
      LMPCTerm *term2 = new LMPCTerm(n2, i, -1.0);
      _MPC->addterm(term2);
      addLMPC(_MPC,true);
    }
  }
  else if(itype == 1) { // normal contact
    double norm = sqrt(face_normal[0]*face_normal[0]+face_normal[1]*face_normal[1]+face_normal[2]*face_normal[2]);
    if(norm == 0.0) cerr << " *** ERROR in Domain::addNodeToNodeLMPCs, face_normal has length 0.0." << endl;
    double gap = (face_normal[0]*gap_vector[0] + face_normal[1]*gap_vector[1] + face_normal[2]*gap_vector[2])/norm;
    addNodalCTC(n1, n2, face_normal[0], face_normal[1], face_normal[2], gap, true, 1, true, lmpcnum);
  }
}
*/

void
Domain::addDirichletLMPCs(int _numDirichlet, BCond *_dbc)
{
  for(int i=0; i<_numDirichlet; ++i) {
    LMPCons *dmpc = new LMPCons(0, _dbc[i].val,
         new LMPCTerm(_dbc[i].nnum,_dbc[i].dofnum,1.0));
    addLMPC(dmpc,false);
  }
}

void
Domain::checkLMPCs(Connectivity *nodeToSub)
{
  if(numLMPC > 0) {
    // TODO consider the case where there is an mpc involving a node/dof that is not connected to any other elements. For example this may
    // be used to connect node with a force or a lumped mass to the structure
    if(nodeToSub) {
      for(int i=0; i < numLMPC; ++i) {
        for(int j=0; j < lmpc[i]->nterms; ++j) {
          int node = lmpc[i]->terms[j].nnum;
          if(node > -1 && nodeToSub->num(node) <= 0) // salinas mpcs can have node = -1 (indicates that node is not in any subdomains on this cpu)
            fprintf(stderr," *** WARNING: MPC %d involves bad node %d \n", lmpc[i]->lmpcnum, node+1);
        }
      }
    }
  }
}


void Domain::computeMatchingWetInterfaceLMPC() {

 if(numWet==0) return; // HB
// Create coupling LMPc's
 double tWI = 0.0;

 tWI -= getTime();

 //fprintf(stderr," ... In computeMatchingWetInterfaceLMPC(Domain.C), numWet is %d ...\n", numWet);
 //fprintf(stderr," ... wet[0] is %d ...\n", wet[0]);

 Connectivity *wetElemToNode = new Connectivity(&(wet[0]), numWet);
 Connectivity *nodeToWetElem = wetElemToNode->reverse();

 Connectivity *nodeToNode = nodeToWetElem->transcon(wetElemToNode);

 //fprintf(stderr," ... Printing WetLMPC nodeToNode ...\n");
 //nodeToNode->print();

 int maxElNodes = 0;

 int iele;
 for(iele=0;iele<numWet;iele++) {
   int nn = wet[iele]->numNodes();
   if (nn>maxElNodes) maxElNodes = nn;
 }

// double *marray = (double *) alloca(sizeof(double)*maxElNodes*maxElNodes*3);

 int iNode = 0;
 for (iNode=0; iNode<nodeToWetElem->csize(); iNode++) {
   //fprintf(stderr," ... iNode = %d, nodeToWetElem->num(iNode) is %d ...\n", iNode, nodeToWetElem->num(iNode));
   if (nodeToWetElem->num(iNode) > 0) {
     LMPCons *wetFSI = new LMPCons(iNode,0.0);
     int iEle;
     for(iEle = 0; iEle < nodeToWetElem->num(iNode); iEle++) {
       int jEle = (*nodeToWetElem)[iNode][iEle];
//       wet[jEle]->wetInterfaceMatrix(geoSource->GetNodes(),marray);
       wet[jEle]->wetInterfaceLMPC(geoSource->GetNodes(),wetFSI,iNode);
     }
     fsi[numFSI++] = wetFSI;
     //wetFSI->print();
   }
 }

 //fprintf(stderr," ... numFSI is %d ...\n", numFSI);

 tWI += getTime();
 //fprintf(stderr,"Time to compute wet interface LMPCs: %f\n",tWI/1000.00);
 if(wetElemToNode) delete wetElemToNode; // HB: added to avoid memory leaks
 if(nodeToWetElem) delete nodeToWetElem; // HB: added to avoid memory leaks
 if(nodeToNode)    delete nodeToNode;    // HB: added to avoid memory leaks
}

int Domain::glToPackElem(int e)
{
 return geoSource->glToPackElem(e);
}

void
Domain::ProcessSurfaceBCs()
{
  BCond *surface_dbc;
  int numSurfaceDirichletBC = geoSource->getSurfaceDirichletBC(surface_dbc);
  for(int i=0; i<numSurfaceDirichletBC; ++i) {
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_dbc[i].nnum) {
        int *glNodes = SurfEntities[j]->GetPtrGlNodeIds();
        int nNodes = SurfEntities[j]->GetnNodes();
        switch(surface_dbc[i].type) {
          default:
          case BCond::Displacements : {
            BCond *bc = new BCond[nNodes];
            for(int k=0; k<nNodes; ++k) { 
              bc[k].nnum = glNodes[k];
              bc[k].dofnum = surface_dbc[i].dofnum;
              bc[k].val = surface_dbc[i].val; 
              bc[k].type = surface_dbc[i].type;
            }
            int numDirichlet_copy = geoSource->getNumDirichlet();
            geoSource->setDirichlet(nNodes, bc);
            if(numDirichlet_copy != 0) delete [] bc;
          } break;
          case BCond::Lmpc : {
            for(int k=0; k<nNodes; ++k) {
              LMPCTerm term(glNodes[k], surface_dbc[i].dofnum, 1);
              LMPCons *lmpc = new LMPCons(-1,surface_dbc[i].val,&term);
              lmpc->setSource(mpc::Lmpc);
              addLMPC(lmpc,false);
            }
          } break;
        }
      }
    }
  }

  BCond *surface_nbc;
  int numSurfaceNeumanBC = geoSource->getSurfaceNeumanBC(surface_nbc);
  for(int i=0; i<numSurfaceNeumanBC; ++i) {
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_nbc[i].nnum) {
        int *glNodes = SurfEntities[j]->GetPtrGlNodeIds();
        int nNodes = SurfEntities[j]->GetnNodes();
        BCond *bc = new BCond[nNodes];
        for(int k=0; k<nNodes; ++k) { bc[k].nnum = glNodes[k]; bc[k].dofnum = surface_nbc[i].dofnum; bc[k].val = surface_nbc[i].val; bc[k].type = surface_nbc[i].type; bc[k].caseid = surface_nbc[i].caseid; }
        int numNeuman_copy = geoSource->getNumNeuman();
        geoSource->setNeuman(nNodes, bc);
        if(numNeuman_copy != 0) delete [] bc;
      }
    }
  }

  BCond *surface_pres;
  int numSurfacePressure = geoSource->getSurfacePressure(surface_pres);
  int nEle = geoSource->getElemSet()->last();
  for(int i=0; i<numSurfacePressure; ++i) {
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_pres[i].nnum) {
        FaceElemSet &faceElemSet = SurfEntities[j]->GetFaceElemSet();
        for(int iele = 0; iele < faceElemSet.last(); ++iele) {
          int nVertices = faceElemSet[iele]->nVertices();
          if(nVertices == 3 || nVertices == 4 || nVertices == 6 || nVertices == 8) {
             int *nodes = new int[nVertices];
             for(int inode = 0; inode < nVertices; ++inode)
               nodes[inode] = SurfEntities[j]->GetPtrGlVertexIds()[faceElemSet[iele]->GetVertex(inode)];
             int type;
             switch(nVertices) {
               case 3: type = 15; break;
               case 4: type = 16; break;
               case 6: type = 17; break;
               case 8: type = 18; break;
             }
             addNeumElem(-1, type, surface_pres[i].val, nVertices, nodes);
             delete [] nodes;
          }
          else {
            cerr << " *** ERROR: can't set pressure for surface " << SurfId << " element " << iele << " nVertices = " << nVertices << endl;
            continue;
          }
        }
      }
    }
  }

  BCond *surface_cfe;
  int numSurfaceConstraint = geoSource->getSurfaceConstraint(surface_cfe);
  for(int i=0; i<numSurfaceConstraint; ++i) {
    double frame_data[9];
    EFrame *ef;
    for(int j=0; j<geoSource->getNumCSframes(); ++j) {
      EFrameData &efd = geoSource->getCSframes()[j];
      if(int(surface_cfe[i].val) == efd.elnum) {
        ef = &efd.frame;
        for(int k=0; k<3; ++k) for(int l=0; l<3; ++l) frame_data[3*k+l] = efd.frame[k][l];
        break;
      }
    }
    for(int j=0; j<nSurfEntity; j++) {
      int SurfId = SurfEntities[j]->ID();
      if(SurfId-1 == surface_cfe[i].nnum) {
        int *glNodes = SurfEntities[j]->GetPtrGlNodeIds();
        int nNodes = SurfEntities[j]->GetnNodes();
        switch(surface_cfe[i].type) {
          case BCond::PointPointDistance : {
            for(int k=0; k<nNodes; ++k) {
              geoSource->getElemSet()->elemadd(nEle, 77, 1, glNodes+k);
              geoSource->setAttrib(nEle, surface_cfe[i].dofnum);
              geoSource->setFrame(nEle, frame_data);
              nEle++;
            }
          } break;
          case BCond::PointLineDistance : {
            for(int k=0; k<nNodes; ++k) {
              geoSource->getElemSet()->elemadd(nEle, 78, 1, glNodes+k);
              geoSource->setAttrib(nEle, surface_cfe[i].dofnum);
              geoSource->setFrame(nEle, frame_data);
              nEle++;
            }
          } break;
          case BCond::PointPlaneDistance : {
            for(int k=0; k<nNodes; ++k) {
              geoSource->getElemSet()->elemadd(nEle, 79, 1, glNodes+k);
              geoSource->setAttrib(nEle, surface_cfe[i].dofnum);
              geoSource->setFrame(nEle, frame_data);
              nEle++;
            }
          } break;
        }
      }
    }
  }

}


void Domain::setNewProperties(int s)
{
  int na = geoSource->getNumAttributes();
  map<int, Attrib> &attr = geoSource->getAttributes();
  if(s==0) {
//   elems_copy.setMyData(true);
   for(int j=0; j<packedEset.last(); ++j) elems_fullcopy.elemadd(j,packedEset[j]);
   for(map<int, Group >::iterator it = geoSource->group.begin(); it != geoSource->group.end(); ++it) {   // loop over all of the groups
     for(int i = 0; i < int(it->second.attributes.size()); ++i) { // loop over attributes
        int iattr = it->second.attributes[i];

        for(int j=0; j<na; ++j) {
          if(attr[j].attr != iattr) continue;
          int jele = glToPackElem(attr[j].nele);
          if(jele == -1) continue;

          StructProp *newProp = new StructProp(*packedEset[jele]->getProperty());
          switch(it->second.randomProperties[0].rprop)  { // switch on the rprop type
            case 0:
               newProp->A = it->second.randomProperties[0].mean; // assuming one group has one random property
              break;
            case 1:
               newProp->E = it->second.randomProperties[0].mean;
              break;
            case 5:
               newProp->kx = it->second.randomProperties[0].mean;
              break;
            case 6:
               newProp->ky = it->second.randomProperties[0].mean;
              break;
            case 7:
               newProp->kz = it->second.randomProperties[0].mean;
              break;
            default:
               cerr << "case " << it->second.randomProperties[0].rprop << " not handled\n";
              break;
          }

          packedEset[jele]->setProp(newProp,true);
          elems_copy.elemadd(jele, packedEset[jele]);
        }
      }
    }
  }
  else {
    packedEset.deleteElems();    // clear elemset
    int count = 0;
    for(int i = 0; i < int(geoSource->group[s-1].attributes.size()); ++i) { // loop over attributes
      int iattr = geoSource->group[s-1].attributes[i];
      for(int j=0; j<na; ++j) {
          if(attr[j].attr != iattr) continue;
          int jele = glToPackElem(attr[j].nele);
          if(jele == -1) continue;
          StructProp *newProp = new StructProp(*elems_copy[jele]->getProperty());
          switch(geoSource->group[s-1].randomProperties[0].rprop)  { // switch on the rprop type
            case 0:
               if(sfem->Gauss) newProp->A = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->A = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 1:
               if(sfem->Gauss) newProp->E = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->E = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 5:
               if(sfem->Gauss) newProp->kx = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->kx = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 6:
               if(sfem->Gauss) newProp->ky = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->ky = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            case 7:
               if(sfem->Gauss) newProp->kz = geoSource->group[s-1].randomProperties[0].std_dev;
               else newProp->kz = geoSource->group[s-1].randomProperties[0].std_dev/sqrt(2.0);
              break;
            default:
               cerr << "case " << geoSource->group[s-1].randomProperties[0].rprop << " not handled\n";
              break;
          }
          elems_copy[jele]->setProp(newProp,true);

          Element* elem = elems_copy[jele];
          packedEset.elemadd(count++,elem);
        }
      }
      setNumElements(packedEset.last());
      deleteAllDOFs();
      makeAllDOFs();
  }
}


void Domain::assignRandMat()  // Equivalent to the non-intrusive version, but used for intrusive
{
  int ndim = sfem->getndim();
  double *xitemp = sfem->getxi();
  int na = geoSource->getNumAttributes(); // Returning the total number of elements, another option domain->numElements();
  map<int, Attrib> &attr = geoSource->getAttributes();
  int jele, count, iattr;
  for(map<int, Group >::iterator it = geoSource->group.begin(); it != geoSource->group.end(); ++it) {   // loop over all of the groups
     for(int i = 0; i < int(it->second.attributes.size()); ++i) { // loop over attributes
       iattr = it->second.attributes[i];
       for(int j=0; j<na; ++j) {
         if(attr[j].attr != iattr) continue;
            jele = glToPackElem(attr[j].nele);
          if(jele == -1)  continue;
         StructProp *newProp = new StructProp(*packedEset[jele]->getProperty());
         switch(it->second.randomProperties[0].rprop)  { // switch on the rprop type
           case 0:
              newProp->A = it->second.randomProperties[0].mean; // assuming one group has one random property
             break;
           case 1:
              newProp->E = it->second.randomProperties[0].mean;
             break;
           case 5:
              newProp->kx = it->second.randomProperties[0].mean;
             break;
           case 6:
              newProp->ky = it->second.randomProperties[0].mean;
             break;
           case 7:
              newProp->kz = it->second.randomProperties[0].mean;
             break;
           default:
              cerr << "case " << it->second.randomProperties[0].rprop << " not handled\n";
             break;
         }
         packedEset[jele]->setProp(newProp,true);
       }
     }
   }

   for (int s=1;s<=ndim;s++) {
    count = 0;
    for(int i = 0; i < int(geoSource->group[s-1].attributes.size()); ++i) { // loop over attributes
      iattr = geoSource->group[s-1].attributes[i];
      for(int j=0; j<na; ++j) {
          if(attr[j].attr != iattr) continue;
          jele = glToPackElem(attr[j].nele);
          if(jele == -1) continue;
          StructProp *newProp =  new StructProp(*packedEset[jele]->getProperty());
          switch(geoSource->group[s-1].randomProperties[0].rprop)  { // switch on the rprop type
            case 0:
               if(sfem->Gauss) newProp->A = newProp->A +  geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->A = newProp->A +  geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 1:
               if(sfem->Gauss) newProp->E = newProp->E + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->E = newProp->E + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 5:
               if(sfem->Gauss) newProp->kx = newProp->kx + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->kx = newProp->kx + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 6:
               if(sfem->Gauss) newProp->ky = newProp->ky + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->ky = newProp->ky + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            case 7:
               if(sfem->Gauss) newProp->kz = newProp->kz + geoSource->group[s-1].randomProperties[0].std_dev*xitemp[s-1];
               else newProp->kz = newProp->kz + geoSource->group[s-1].randomProperties[0].std_dev*(pow(xitemp[s-1],2)-1)/sqrt(2.0);
              break;
            default:
               cerr << "case " << geoSource->group[s-1].randomProperties[0].rprop << " not handled\n";
              break;
          }
          packedEset[jele]->setProp(newProp,true);
        }
      }
  }

}


void Domain::retrieveElemset()
{
 int na = geoSource->getNumAttributes();
 map<int, Attrib> &attr = geoSource->getAttributes();
 packedEset.deleteElems();    // clear elemset
 int count = 0;
 for(int j=0; j<na; ++j) {
   int jele = glToPackElem(attr[j].nele);
   if(jele == -1) continue;
   Element* elem = elems_fullcopy[jele];
   packedEset.elemadd(count++,elem);
 }
 setNumElements(packedEset.last());
 deleteAllDOFs();
 makeAllDOFs();
}

//CBM: new stuff
void
Domain::setEigenValue(double _lbound, int _nshifts, int _maxArnItr)
{
 // for multiple shifts eigen anaylsis
  if((_lbound < 0.0)) {
    filePrint(stderr, " *** WARNING: lbound is negative \n");
    return;
  } else {
    sinfo.lbound = _lbound;
    sinfo.nshifts = _nshifts;
  }

  sinfo.maxArnItr = _maxArnItr;
  sinfo.doEigSweep = true;
}

void
Domain::setEigenValues(double _lbound , double _ubound, int _neigps, int _maxArnItr)
{
  // for multiple shifts eigen anaylsis
  if(_lbound < 0.0 || _ubound < 0.0) {
    filePrint(stderr, " *** WARNING: lbound or ubound is negative \n");
    return;
  } else if (_lbound > _ubound) {
    filePrint(stderr, " *** WARNING: lbound > ubound, resetting interval to [ubound lbound]\n");
    sinfo.lbound = _ubound;
    sinfo.ubound = _lbound;
  } else {
    sinfo.lbound = _lbound;
    sinfo.ubound = _ubound;
  }
  sinfo.neigps = _neigps;
  sinfo.maxArnItr = _maxArnItr;
  sinfo.doEigSweep = true;
}

void
Domain::deleteAllLMPCs()
{
  lmpc.deleteArray();
  lmpc.restartArray();
  numLMPC = 0;
  nMortarLMPCs = 0;
  numCTC = 0;
  if(mortarToMPC) {
    delete mortarToMPC;
    mortarToMPC = 0;
  }
}

void
Domain::deleteSomeLMPCs(mpc::ConstraintSource s)
{ 
  int j = 0;
  for(int i = 0; i < numLMPC; ++i) {
    if(lmpc[i]->getSource() == s) {
      if(lmpc[i]->getType() == mpc::Inequality) numCTC--;
      if(s == mpc::ContactSurfaces || s == mpc::TiedSurfaces) nMortarLMPCs--;
      delete lmpc[i];
      lmpc[i] = 0;
    }
    else {
      if(i != j) {
        lmpc[j] = lmpc[i];
        lmpc[i] = 0;
      }
      j++; 
    }
  }
  //cerr << "deleted " << numLMPC-j << " mpcs, ";
  numLMPC = j;
  if(mortarToMPC && (s == mpc::ContactSurfaces || s == mpc::TiedSurfaces)) {
    delete mortarToMPC;
    mortarToMPC = 0;
  }
}

void
Domain::UpdateContactSurfaceElements(GeomState *geomState)
{
  // copy the lagrange multipliers from geomState
  std::vector<double> mu;
  if(sinfo.lagrangeMult) {
    for(int i = 0; i < numLMPC; ++i) {
      if(lmpc[i]->getSource() == mpc::ContactSurfaces) {
        mu.push_back(geomState->getMultiplier(lmpc[i]->id));
      }
    }
    geomState->clearMultiplierNodes();
  }

  if(!p) {
    p = new StructProp(); 
    p->lagrangeMult = sinfo.lagrangeMult;
    p->initialPenalty = p->penalty = sinfo.penalty;
    p->type = StructProp::Constraint;
    p->constraint_hess = sinfo.constraint_hess;
    p->constraint_hess_eps = sinfo.constraint_hess_eps;
  }
  int count = 0;
  int nEle = packedEset.last();
  int count1 = 0;
  int nNode = geomState->numNodes();
  for(int i = 0; i < numLMPC; ++i) {
    if(lmpc[i]->getSource() == mpc::ContactSurfaces) {
      if(count < contactSurfElems.size()) { // replace
        //cerr << "replacing element " << contactSurfElems[count] << " with lmpc " << i << endl;
        packedEset.deleteElem(contactSurfElems[count]);
        packedEset.mpcelemadd(contactSurfElems[count], lmpc[i]); // replace 
        packedEset[contactSurfElems[count]]->setProp(p);
        if(packedEset[contactSurfElems[count]]->numInternalNodes() == 1) { // i.e. lagrange multiplier
          int in[1] = { nNode++ };
          packedEset[contactSurfElems[count]]->setInternalNodes(in);
          geomState->addMultiplierNode(lmpc[i]->id, mu[i]);
        }
        count1++;
      }
      else { // new
        //cerr << "adding lmpc " << i << " to elemset at index " << nEle << endl;
        packedEset.mpcelemadd(nEle, lmpc[i]); // new
        packedEset[nEle]->setProp(p);
        if(packedEset[nEle]->numInternalNodes() == 1) {
          int in[1] = { nNode++ };
          packedEset[nEle]->setInternalNodes(in);
          geomState->addMultiplierNode(lmpc[i]->id, mu[i]);
        }
        contactSurfElems.push_back(nEle);
        nEle++;
      }
      count++;
    }
  }
  int count2 = 0;
  while(count < contactSurfElems.size()) {
    //cerr << "deleting elemset " << contactSurfElems.back() << endl;
    packedEset.deleteElem(contactSurfElems.back());
    contactSurfElems.pop_back();
    count2++;
  }
  packedEset.setEmax(nEle-count2); // because element set is packed
  //cerr << "replaced " << count1 << " and added " << count-count1 << " new elements while removing " << count2 << endl;
  numele = packedEset.last(); 
  numnodes = geomState->numNodes();
}

