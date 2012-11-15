#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Domain.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <HelmAxi.d/AxiHElem.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/Memory.h>
#include <Utils.d/pstress.h>
#include <Timers.d/GetTime.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>

#include <Driver.d/GeoSource.h>
#include <list>

// const double defaultTemp = -10000000;
extern int verboseFlag;

void
Domain::buildPrescDisp(Vector &pDisp, double lambda)
{
 int i;
 pDisp.zero();
 // Set the dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
  int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) {
     //fprintf(stderr," *** WARNING: dof does not exist: node %d dof %d ***\n",
     //                dbc[i].nnum+1,dbc[i].dofnum+1);
     continue;
   }

  int cdof = c_dsa->invRCN(dof);
  pDisp[cdof] = lambda*dbc[i].val;
 }
}

void
Domain::buildPrescDisp(Vector &pDisp, double t, double)
{
 int i;
 pDisp.zero();
 // Set the dirichlet boundary conditions
 for(i=0; i<numDirichlet; ++i) {
  int dof  = dsa->locate(dbc[i].nnum, 1 << dbc[i].dofnum);
   if(dof < 0) {
     //fprintf(stderr," *** WARNING: dof does not exist: node %d dof %d ***\n",
     //                dbc[i].nnum+1,dbc[i].dofnum+1);
     continue;
   }

  int cdof = c_dsa->invRCN(dof);
  pDisp[cdof] = t*dbc[i].val;
 }
}

double *
Domain::getNodalTemperatures()
{
  return temprcvd;
}

void
Domain::initNodalTemperatures()
{
  if(sinfo.thermalLoadFlag && sinfo.thermoeFlag == -1) {

    temprcvd = new double[numnodes];

    int i;
    for(i = 0; i < numnodes; ++i)
      temprcvd[i] = defaultTemp;

    for(i = 0; i < numDirichlet; ++i)
      if((1 << dbc[i].dofnum) == DofSet::Temp) {
        temprcvd[dbc[i].nnum] = dbc[i].val;
        //fprintf(stderr," Temp %f\n", temprcvd[dbc[i].nnum]);
      }
  }
}

FILE *
Domain::openFile(char *fileName, const char *extension)
{
 // Open decomposition file
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 int len1 = strlen(cinfo->checkfile);
 int len2 = strlen(extension);
 char *file = (char *) dbg_alloca(sizeof(char)*(len1+len2+1));
 strcpy(file, cinfo->checkfile);
 strcat(file, extension);

 FILE *filePtr;
 if((filePtr= fopen(file,"w")) == (FILE *) 0 )
   fprintf(stderr," *** ERROR: Cannot open %s ***\n",file );

 return filePtr;
}

void
Domain::printStatistics()
{
   filePrint(stderr,"\n ---------- PROBLEM PARAMETERS --------");
   filePrint(stderr,"\n ... # Nodes              = %7d ...",numnodes);
   filePrint(stderr,"\n ... # Elements           = %7d ...",numele);
   filePrint(stderr,"\n ... # Unconstrained dofs = %7d ...",numUncon());
   filePrint(stderr,"\n ... # Constrained dofs   = %7d ...",
           numDirichlet+numComplexDirichlet);
   filePrint(stderr,"\n ... Total # dofs         = %7d ...",numdof());
   filePrint(stderr,"\n ... # Loaded dofs        = %7d ...",
           numNeuman+numComplexNeuman);
   if(gravityFlag())
     filePrint(stderr,"\n ... Gravity Load is Applied        ...");
   filePrint(stderr,"\n ... # Output Files       = %7d ...",geoSource->getNumOutInfo());
   filePrint(stderr,"\n --------------------------------------\n");
}

double
Domain::computeStructureMass(bool printFlag)
{
  // Compute total mass and mass center of gravity
  // Added calculation for moments of inertia

  double totmas = 0.0; // total mass of structural system (including fluid elemets and discrete masses)
  double xc     = 0.0;
  double yc     = 0.0;
  double zc     = 0.0;
  double Mx     = 0.0;
  double My     = 0.0;
  double Mz     = 0.0;
  double Ixx    = 0.0;
  double Iyy    = 0.0;
  double Izz    = 0.0;
  double Ixy    = 0.0;
  double Iyz    = 0.0;
  double Ixz    = 0.0;
  maxNumNodes = 0;

  // compute max number nodes per element
  int iele, i;
  for(iele=0; iele < numele; ++iele) {
    int numNodesPerElement = packedEset[iele]->numNodes();
    maxNumNodes = myMax(maxNumNodes, numNodesPerElement);
  }

  // allocate one array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  for(iele=0; iele < numele; ++iele) {
    double elementMass = packedEset[iele]->getMass(nodes);
    totmas += elementMass;
    int numNodesPerElement = packedEset[iele]->numNodes();
    packedEset[iele]->nodes(nodeNumbers);
    int numRealNodes = 0;
    for(i=0; i<numNodesPerElement; ++i)
      if(nodes[nodeNumbers[i]]) numRealNodes++;

    double massPerNode = elementMass/numRealNodes;

    for(i=0; i<numNodesPerElement; ++i) {
       if(nodes[nodeNumbers[i]] == 0) continue;
       Node &node = nodes.getNode(nodeNumbers[i]);

       double x = node.x;
       double y = node.y;
       double z = node.z;

       xc  += massPerNode*x;
       yc  += massPerNode*y;
       zc  += massPerNode*z;

       Ixx += massPerNode*(y*y + z*z);
       Iyy += massPerNode*(x*x + z*z);
       Izz += massPerNode*(x*x + y*y);

       Ixy -= massPerNode*x*y;
       Ixz -= massPerNode*x*z;
       Iyz -= massPerNode*y*z;
    }
  }

  delete [] nodeNumbers;

  // now add the mass from the fluid elements (these are in a different element set)
  if(geoSource->numElemFluid() > 0) {
    for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
      int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
      maxNumNodes = myMax(maxNumNodes, numNodesPerElement);
    }

    nodeNumbers = new int[maxNumNodes];
  
    double fluidmas = 0.0;
    for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
      double elementMass = (*(geoSource->getPackedEsetFluid()))[iele]->getMass(nodes);
      fluidmas += elementMass;
      int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
      (*(geoSource->getPackedEsetFluid()))[iele]->nodes(nodeNumbers);
      int numRealNodes = 0;
      for(i=0; i<numNodesPerElement; ++i)
        if(nodes[nodeNumbers[i]]) numRealNodes++;

      double massPerNode = elementMass/numRealNodes;

      for(i=0; i<numNodesPerElement; ++i) {
         if(nodes[nodeNumbers[i]] == 0) continue;
         Node &node = nodes.getNode(nodeNumbers[i]);
  
         double x = node.x;
         double y = node.y;
         double z = node.z;

         xc  += massPerNode*x;
         yc  += massPerNode*y;
         zc  += massPerNode*z;

         Ixx += massPerNode*(y*y + z*z);
         Iyy += massPerNode*(x*x + z*z);
         Izz += massPerNode*(x*x + y*y);

         Ixy -= massPerNode*x*y;
         Ixz -= massPerNode*x*z;
         Iyz -= massPerNode*y*z;
      }
    }

    delete [] nodeNumbers;

    if(printFlag) {
      filePrint(stderr," Fluid Mass = %f\n",fluidmas);
      filePrint(stderr," --------------------------------------\n");
    }

    totmas += fluidmas;
  
  }

  double Iyx = Ixy;
  double Izx = Ixz;
  double Izy = Iyz;

  // Add discrete masses
  DMassData *current = firstDiMass;
  while(current != 0) {
    int n = current->node;
    // check if the node is in the model
    //if(nodes.exist(n)) {
      Node &node = nodes.getNode(n);
      switch(current->dof) {
        case 0: {
          Mx += current->diMass;
          xc += current->diMass*node.x;
          Iyy += current->diMass*(node.z*node.z);
          Izz += current->diMass*(node.y*node.y);
          Ixy -= current->diMass*(node.x*node.y);
          Ixz -= current->diMass*(node.x*node.z);
        } break;
        case 1: {
          My += current->diMass;
          yc += current->diMass*node.y;
          Ixx += current->diMass*(node.z*node.z);
          Izz += current->diMass*(node.x*node.x);
          Iyx -= current->diMass*(node.x*node.y);
          Iyz -= current->diMass*(node.y*node.z);
        } break;
        case 2: {
          Mz += current->diMass;
          zc += current->diMass*node.z;
          Ixx += current->diMass*(node.y*node.y);
          Iyy += current->diMass*(node.x*node.x);
          Izx -= current->diMass*(node.x*node.z);
          Izy -= current->diMass*(node.y*node.z);
        } break;
        case 3: {
          if(current->jdof == -1 || current->jdof == 3) Ixx += current->diMass;
          else if(current->jdof == 4) Ixy += current->diMass;
          else if(current->jdof == 5) Ixz += current->diMass;
        } break;
        case 4: {
          if(current->jdof == -1 || current->jdof == 4) Iyy += current->diMass;
          else if(current->jdof == 3) Ixy += current->diMass;
          else if(current->jdof == 5) Iyz += current->diMass;
        } break;
        case 5: {
          if(current->jdof == -1 || current->jdof == 5) Izz += current->diMass;
          else if(current->jdof == 3) Ixz += current->diMass;
          else if(current->jdof == 4) Iyz += current->diMass;
        } break;
      }
    //}
    current = current->next;
  }

  Mx += totmas;
  My += totmas;
  Mz += totmas;

  if(Mx != 0.0) xc /= Mx;
  if(My != 0.0) yc /= My;
  if(Mz != 0.0) zc /= Mz;
  // change moments of inertia to centroidal axes using parallel axes theorem: I_z = I_cm + m*d^2
  Ixx -= (My*(yc*yc)+Mz*(zc*zc));
  Iyy -= (Mz*(xc*xc)+Mz*(zc*zc));
  Izz -= (Mx*(xc*xc)+My*(yc*yc));

  Ixy += Mx*xc*yc;
  Iyx += My*xc*yc;
  Ixz += Mx*xc*zc;
  Izx += Mz*xc*zc;
  Iyz += My*yc*zc;
  Izy += Mz*yc*zc;

  totmas = (Mx+My+Mz)/3.0;
  if(printFlag) {
    if(Mx != My || Mx != Mz || My != Mz) {
      filePrint(stderr," Directional Mass\n");
      filePrint(stderr," Mx = %f My = %f Mz = %f\n",Mx,My,Mz);
      filePrint(stderr," --------------------------------------\n");
    }

    filePrint(stderr," Moments of Inertia\n");
    filePrint(stderr," Ixx = %e Iyy = %e Izz = %e\n",Ixx,Iyy,Izz);
    filePrint(stderr," --------------------------------------\n");

    filePrint(stderr," Products of Inertia\n");
    filePrint(stderr," Ixy = %e Iyz = %e Ixz = %e\n",Ixy,Iyz,Ixz);
    filePrint(stderr," --------------------------------------\n");

    if (Iyx != Ixy || Ixz != Izx || Iyz != Izy)  {
      filePrint(stderr," WARNING: Non-Symmetric Products of Inertia (Check DiMASS) \n");
      filePrint(stderr," Iyx = %e Izy = %e Izx = %e\n",Iyx,Izy,Izx);
      filePrint(stderr," --------------------------------------\n");
    }

    filePrint(stderr," Center of Gravity\n");
    filePrint(stderr," x = %f y = %f z = %f\n",xc,yc,zc);
    filePrint(stderr," --------------------------------------\n");
  }

  // Computation to find node closest to the center of gravity
  int nodeMarker = 0;
  double minDistance = 0.0;
  for(i=0; i<numnodes; ++i) {
    Node *node = nodes[i];
    if(node) {
      double x = node->x;
      double y = node->y;
      double z = node->z;

      double dx = xc - x;
      double dy = yc - y;
      double dz = zc - z;

      double distance = sqrt(dx*dx+dy*dy+dz*dz);

      if(i == 0) minDistance = distance;

      if(distance < minDistance) {
        minDistance = distance;
        nodeMarker  = i;
      }
    }
  }

  if(printFlag) {
    filePrint(stderr," Node %d is closest to the Center of Gravity\n",nodeMarker+1);
    Node *thisNode = nodes[nodeMarker];
    if(thisNode) {
      filePrint(stderr," Node %d has coordinates: %e %e %e \n",
                nodeMarker + 1, thisNode->x, thisNode->y, thisNode->z);
      filePrint(stderr," It is %e from the center of gravity\n",minDistance);
      filePrint(stderr," --------------------------------------\n");
    }
  }

  // Compute Geometric center of gravity of structure
  double xmax = 0.0;
  double ymax = 0.0;
  double zmax = 0.0;

  int nComponents = renumb.numComp + renumbFluid.numComp;

  if(printFlag) filePrint(stderr," Number of Components = %d\n",renumb.numComp);
  for(int n=0; n<renumb.numComp; ++n) {
    xc = 0.0, yc = 0.0, zc = 0.0;
    int realNodeCnt = 0;
    for(i = renumb.xcomp[n]; i<renumb.xcomp[n+1]; ++i) {
      int inode = renumb.order[i];

      if(dsa->firstdof(inode) == -1 || nodes[inode] == 0) continue;
      realNodeCnt++;
      Node &nd = nodes.getNode(inode);
      xc += nd.x;
      yc += nd.y;
      zc += nd.z;
      if(nd.x > xmax) xmax = nd.x;
      if(nd.y > ymax) ymax = nd.y;
      if(nd.z > zmax) zmax = nd.z;
    }

    double xg = xc/realNodeCnt;
    double yg = yc/realNodeCnt;
    double zg = zc/realNodeCnt;

    if(printFlag) {
      filePrint(stderr," Component %d: Centroid\n", n+1);
      filePrint(stderr," x = %f y = %f z = %f\n",xg,yg,zg);
    }
  }
  for(int n=0; n<renumbFluid.numComp; ++n) {
    xc = 0.0, yc = 0.0, zc = 0.0;
    int realNodeCnt = 0;
    for(i = renumbFluid.xcomp[n]; i<renumbFluid.xcomp[n+1]; ++i) {
      int inode = renumbFluid.order[i];
  
      if(dsa->firstdof(inode) == -1 || nodes[inode] == 0) continue;
      realNodeCnt++;
      Node &nd = nodes.getNode(inode);
      xc += nd.x;
      yc += nd.y;
      zc += nd.z;
      if(nd.x > xmax) xmax = nd.x;
      if(nd.y > ymax) ymax = nd.y;
      if(nd.z > zmax) zmax = nd.z;
    }
  
    double xg = xc/realNodeCnt;
    double yg = yc/realNodeCnt;
    double zg = zc/realNodeCnt;
    
    if(printFlag) {
      filePrint(stderr," Component %d (Fluid): Centroid\n", renumb.numComp+n+1);
      filePrint(stderr," x = %f y = %f z = %f\n",xg,yg,zg);
    }
  }
 
  if(printFlag) {
    filePrint(stderr," --------------------------------------\n");

    filePrint(stderr," Maximum x dimension = %f\n",xmax);
    filePrint(stderr," Maximum y dimension = %f\n",ymax);
    filePrint(stderr," Maximum z dimension = %f\n",zmax);
    filePrint(stderr," --------------------------------------\n");
  }

  return totmas;
}

double
Domain::computeFluidMass()
{
  // Compute total mass and mass center of gravity
  // Added calculation for moments of inertia

  double totmas = 0.0; // total mass of fluid
  double xc     = 0.0;
  double yc     = 0.0;
  double zc     = 0.0;
  double Mx     = 0.0;
  double My     = 0.0;
  double Mz     = 0.0;
  double Ixx    = 0.0;
  double Iyy    = 0.0;
  double Izz    = 0.0;
  double Ixy    = 0.0;
  double Iyz    = 0.0;
  double Ixz    = 0.0;
  maxNumNodes = 0;

  // compute max number nodes per element
  int iele, i;
  for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
    int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
    maxNumNodes = myMax(maxNumNodes, numNodesPerElement);
  }

  // allocate one array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  for(iele=0; iele < geoSource->numElemFluid(); ++iele) {
    double elementMass = (*(geoSource->getPackedEsetFluid()))[iele]->getMass(nodes);
    totmas += elementMass;
    int numNodesPerElement = (*(geoSource->getPackedEsetFluid()))[iele]->numNodes();
    (*(geoSource->getPackedEsetFluid()))[iele]->nodes(nodeNumbers);
    int numRealNodes = 0;
    for(i=0; i<numNodesPerElement; ++i)
      if(nodes[nodeNumbers[i]]) numRealNodes++;

    double massPerNode = elementMass/numRealNodes;

    for(i=0; i<numNodesPerElement; ++i) {
       if(nodes[nodeNumbers[i]] == 0) continue;
       Node &node = nodes.getNode(nodeNumbers[i]);

       double x = node.x;
       double y = node.y;
       double z = node.z;

       xc  += massPerNode*x;
       yc  += massPerNode*y;
       zc  += massPerNode*z;

       Ixx += massPerNode*(y*y + z*z);
       Iyy += massPerNode*(x*x + z*z);
       Izz += massPerNode*(x*x + y*y);

       Ixy -= massPerNode*x*y;
       Ixz -= massPerNode*x*z;
       Iyz -= massPerNode*y*z;
    }
  }

  delete [] nodeNumbers;

  if(totmas != 0.0) xc /= totmas;
  if(totmas != 0.0) yc /= totmas;
  if(totmas != 0.0) zc /= totmas;
  // change moments of inertia to centroidal axes using parallel axes theorem: I_z = I_cm + m*d^2
  Ixx -= (totmas*(yc*yc)+totmas*(zc*zc));
  Iyy -= (totmas*(xc*xc)+totmas*(zc*zc));
  Izz -= (totmas*(xc*xc)+totmas*(yc*yc));

  Ixy += totmas*xc*yc;
  Ixz += totmas*xc*zc;
  Iyz += totmas*yc*zc;

  filePrint(stderr," Moments of Inertia\n");
  filePrint(stderr," Ixx = %e Iyy = %e Izz = %e\n",Ixx,Iyy,Izz);
  filePrint(stderr," --------------------------------------\n");

  filePrint(stderr," Products of Inertia\n");
  filePrint(stderr," Ixy = %e Iyz = %e Ixz = %e\n",Ixy,Iyz,Ixz);
  filePrint(stderr," --------------------------------------\n");

  filePrint(stderr," Center of Gravity\n");
  filePrint(stderr," x = %f y = %f z = %f\n",xc,yc,zc);
  filePrint(stderr," --------------------------------------\n");

  // Computation to find node closest to the center of gravity
  int nodeMarker = 0;
  double minDistance = 0.0;
  for(i=0; i<numnodes; ++i) {
    Node *node = nodes[i];
    if(node) {
      double x = node->x;
      double y = node->y;
      double z = node->z;

      double dx = xc - x;
      double dy = yc - y;
      double dz = zc - z;

      double distance = sqrt(dx*dx+dy*dy+dz*dz);

      if(i == 0) minDistance = distance;

      if(distance < minDistance) {
        minDistance = distance;
        nodeMarker  = i;
      }
    }
  }

  filePrint(stderr," Node %d is closest to the Center of Gravity\n",nodeMarker+1);
  Node *thisNode = nodes[nodeMarker];
  if(thisNode) {
    filePrint(stderr," Node %d has coordinates: %e %e %e \n",
              nodeMarker + 1, thisNode->x, thisNode->y, thisNode->z);
    filePrint(stderr," It is %e from the center of gravity\n",minDistance);
    filePrint(stderr," --------------------------------------\n");
  }

  // Compute Geometric center of gravity of fluid
  double xmax = 0.0;
  double ymax = 0.0;
  double zmax = 0.0;

  int nComponents = renumbFluid.numComp;

  filePrint(stderr," Number of Components = %d\n",renumbFluid.numComp);
  for(int n=0; n<nComponents; ++n) {
    xc = 0.0, yc = 0.0, zc = 0.0;
    int realNodeCnt = 0;
    for(i = renumbFluid.xcomp[n]; i<renumb.xcomp[n+1]; ++i) {
      int inode = renumbFluid.order[i];

      if(dsa->firstdof(inode) == -1 || nodes[inode] == 0) continue;
      realNodeCnt++;
      Node &nd = nodes.getNode(inode);
      xc += nd.x;
      yc += nd.y;
      zc += nd.z;
      if(nd.x > xmax) xmax = nd.x;
      if(nd.y > ymax) ymax = nd.y;
      if(nd.z > zmax) zmax = nd.z;
    }

    double xg = xc/realNodeCnt;
    double yg = yc/realNodeCnt;
    double zg = zc/realNodeCnt;

    filePrint(stderr," Component (%d,%d): Centroid\n",
              renumbFluid.xcomp[n], renumbFluid.xcomp[n+1]);
    filePrint(stderr," x = %f y = %f z = %f\n",xg,yg,zg);
  }
  filePrint(stderr," --------------------------------------\n");

  filePrint(stderr," Maximum x dimension = %f\n",xmax);
  filePrint(stderr," Maximum y dimension = %f\n",ymax);
  filePrint(stderr," Maximum z dimension = %f\n",zmax);
  filePrint(stderr," --------------------------------------\n");

  return totmas;
}

Rbm *
Domain::constructHzem(bool printFlag)
{
  Rbm *rbm = new Rbm(dsa, c_dsa);
  if(printFlag)
    cerr << " ... GRBM algorithm detected " << rbm->numRBM() << " rigid body or zero energy modes ...\n";
  return rbm;
}

// Function to construct zero energy modes, ADDED FOR SLOSHING PROBLEM, EC, 20070723
Rbm *
Domain::constructSlzem(bool printFlag)
{
  Rbm *rbm = new Rbm(dsa, c_dsa);
  if(printFlag)
    cerr << " ... GRBM algorithm detected " << rbm->numRBM() << " rigid body or zero energy modes ...\n";
  return rbm;
}

Rbm *
Domain::constructRbm(bool printFlag)
{
  Rbm *rbm = 0;
  if(numLMPC) {
    if(renumb_nompc.numComp == 0)
      rbm = new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd, renumb,  numLMPC, lmpc);
    else
      rbm = new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd, renumb_nompc,  numLMPC, lmpc);
  }
  else 
    rbm = new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd, renumb);
  if(printFlag)
    cerr << " ... GRBM algorithm detected " << rbm->numRBM() << " rigid body or zero energy modes ...\n";
  return rbm;
}

/*
// Function to construct rigid body modes
Rbm *
Domain::constructRbm(IntFullM *fm)
{
 return new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd, renumb, fm);
}
*/

/*
Static, Dynamic, Modal, NonLinStatic, NonLinDynam,
ArcLength, ConditionNumber, TempDynamic, Top,
AxiHelm, MatNonLinStatic, MatNonLinDynam,
Helmholtz, HelmholtzFreqSweep, HelmholtzDirSweep, HelmholtzMF, HelmholtzSO
*/

const char* problemTypeMessage[] = {
" ... Linear Static Analysis         ... \n",
" ... Linear Dynamics/Quasistatics   ... \n",
" ... Modal Analysis                 ... \n",
" ... Nonlinear Static Analysis      ... \n",
" ... Nonlinear Dynamics/Quasistatics... \n",
" ... Nonlin. Stat. Anal. Arc Length ... \n",
" ... Condition Number Analysis      ... \n",
" ... Thermal Analysis               ... \n",
" ... Outputing TOPDOMDEC File       ... \n",
" ... Axisymmetric Acoustic Analysis ... \n",
" ... Material Nonlin. Static Analysis ... \n",
" ... Material Nonlin. Dynamic Analysis ... \n",
" ... Acoustic Scattering Analysis   ... \n",
" ... Acoustic Frequency Sweep Analysis   ... \n",
" ... Helmholtz Direction Sweep Analysis ... \n",
" ... HelmholtzMF Analysis           ... \n",
" ... HelmholtzSO Analysis           ... \n",
" ... Computing Decomposition        ... \n",
" ... Nonlin. Themal Dynamic Analysis... \n",
" ... Discontinuous Enrichment Method... \n",
" ... POD ROM Offline Computations   ... \n",
""
};

const char* solverTypeMessage[] = {
" ... Skyline Solver is Selected     ... \n",
" ... Sparse Solver is Selected      ... \n",
" ... BlockSky Solver is Selected    ...\n",
#ifdef USE_EIGEN3
" ... SimplicalLLT Solver is Selected... \n",
" ... SimplicalLDLT Solver is Selec'd... \n",
#ifdef EIGEN_CHOLMOD_SUPPORT
" ... Cholmod Solver is Selected     ... \n",
#else
" ... Sparse Solver is Selected      ... \n",
#endif
#ifdef EIGEN_UMFPACK_SUPPORT
" ... UmfPack Solver is Selected     ... \n",
#else
" ... Sparse Solver is Selected      ... \n",
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
" ... SuperLU Solver is Selected     ... \n",
#else
" ... Sparse Solver is Selected      ... \n",
#endif
#else
" ... Sparse Solver is Selected      ... \n",
" ... Sparse Solver is Selected      ... \n",
" ... Sparse Solver is Selected      ... \n",
" ... Sparse Solver is Selected      ... \n",
" ... Sparse Solver is Selected      ... \n",
#endif
#ifdef USE_SPOOLES
" ... Spooles Solver is Selected     ... \n",
#else
" ... Sparse Solver is Selected      ... \n",
#endif
#ifdef USE_MUMPS
" ... Mumps Solver is Selected       ... \n",
#else
" ... Sparse Solver is Selected      ... \n",
#endif
"",
" ... POD-GN Solver is Selected      ... \n",
" ... POD-Galerkin Solver is Selected... \n",
" ... POD-Galerkin Solver is Selected... \n",
" ... Goldfarb-Idnani Solver is Sel'd... \n",
#ifdef EIGEN_SPARSELU_SUPPORT
" ... SparseLU Solver is Selected    ... \n",
#else
" ... Sparse Solver is Selected      ... \n",
#endif
#ifdef EIGEN_SPQR_SUPPORT
" ... SparseQR Solver is Selected    ... \n",
#else
" ... Sparse Solver is Selected      ... \n"
#endif
};

void
Domain::preProcessing()
{
 // ... CONSTRUCT DOMAIN ELEMENT TO NODE CONNECTIVITY
 matrixTimers->makeConnectivity -= getTime();
 if(elemToNode) delete elemToNode;
 elemToNode = new Connectivity(&packedEset);
 if(sinfo.HEV) {
   if(elemToNodeFluid) delete elemToNodeFluid;
   elemToNodeFluid = new Connectivity(geoSource->getPackedEsetFluid());
 }
 matrixTimers->makeConnectivity += getTime();

 // ... RENUMBER IF NECESSARY AND ASKED FOR
 matrixTimers->renumbering -= getTime();

/*
 // Check if we are renumbering or not
 if(sinfo.renum || sinfo.sparse_renum)
   filePrint(stderr," ... Renumbering as Specified       ...\n");
 else {
   if(sinfo.subtype != 1 && sinfo.subtype != 8 && sinfo.subtype != 9) // sparse and spooles are always renumbered internally
     filePrint(stderr," *** WARNING: Renumbering Turned Off ***\n");
 }

 // Check if we are using a Solver that does not require renumbering
 if(sinfo.renum && (sinfo.subtype == 3 || sinfo.subtype == 4)) {
   filePrint(stderr," *** WARNING: Renumbering is NOT necessary \n");
   filePrint(stderr,"     for the selected solver; continuing \n");
   filePrint(stderr,"     with renumbering OFF! ***\n");

   // Turn off renumbering (GRBM needs renumbering due to components!)
   sinfo.renum = 1;
 }
*/

 // Perform renumbering if necessary
 Renumber rnum = getRenumbering();
 Renumber* rnumFluid = 0;
 if(sinfo.HEV) rnumFluid = getRenumberingFluid();
 matrixTimers->renumbering += getTime();

 // ... CONSTRUCT DOF SET ARRAY
 matrixTimers->createDofs -= getTime();
 if(dsa) delete dsa;
 dsa = new DofSetArray(numnodes, packedEset, rnum.renumb);
 if(sinfo.HEV) {
   dsaFluid = new DofSetArray(numnodesFluid, *(geoSource->getPackedEsetFluid()), rnumFluid->renumb);
 }
 matrixTimers->createDofs += getTime();

// TEMP CODE TO OUTPUT ACTIVE DOFS IN PROFILE
// dsa = dof set array

/*
  ofstream dof_file;
  dof_file.open("profile.dat");
  if(dof_file.fail()) {
    cout << "ERROR: opening profile.dat" << endl;
  }
  dof_file << setw(6) << numnodes << endl;

  int inode,d;
  int flag[6];
  for(inode = 0; inode < numnodes; ++inode) {
    flag[0] = dsa->locate(inode, DofSet::Xdisp);
    flag[1] = dsa->locate(inode, DofSet::Ydisp);
    flag[2] = dsa->locate(inode, DofSet::Zdisp);
    flag[3] = dsa->locate(inode, DofSet::Xrot);
    flag[4] = dsa->locate(inode, DofSet::Yrot);
    flag[5] = dsa->locate(inode, DofSet::Zrot);

    for(d=0;d<6; ++d) {
      if(flag[d] > 0) // active
        flag[d] = 1;
      else // inactive
        flag[d] = 0;
    }
    dof_file << setw(6) << (inode+1);
    for(d=0;d<6; ++d) {
      dof_file << " " << setw(1) << flag[d];
    }
    dof_file << endl;
  }
  dof_file.close();
  exit(1);
*/
}

void
Domain::make_constrainedDSA()
{
 if(c_dsa) delete c_dsa;
 c_dsa = new ConstrainedDSA(*dsa, numDirichlet, dbc);
 //ADDED FOR HEV PROBLEM, EC, 20070820
 if (solInfo().HEV)  {
   c_dsaFluid = new ConstrainedDSA(*dsaFluid, numDirichletFluid, dbcFluid);
 }
}

void
Domain::make_constrainedDSA(int *bc)
{
 // c_dsa = constrained dof set array
 // dsa   = dof set array
 // dbc   = dirichlet boundary conditions

 // numDirichlet = number of dirichlet boundary conditions
 // cdbc         = complex dirichlet boundary conditions

 // numComplexDirichlet = number of complex
 //                       dirichlet boundary conditions

 // bc = integer array marking dofs that are
 //      constrained or have forces applied

 //ADDED FOR HEV PROBLEM, EC, 20070820
 if (solInfo().HEV)  {
   //int numdofFluid = dsaFluid->size();
   //int numdof = dsa->size();
   //for (int i=0; i<num
   //c_dsaFluid = new ConstrainedDSA(*dsaFluid, dbcFluid, numDirichletFluid, cdbc, 0, bc);
                              //numComplexDirichlet, bc);
   fprintf(stderr," *** c_dsaFluid NOT built in this version of constructor! ***\n");
 }
 if(c_dsa) delete c_dsa;
 c_dsa = new ConstrainedDSA(*dsa, dbc, numDirichlet, cdbc,
                            numComplexDirichlet, bc);
}

void
Domain::make_constrainedDSA(int fake)
{ // make a fake constrainedDSA if fake !=0; ie lie to the code by telling
  //   it that there are noconstraints
  if(fake){
    if(c_dsa) delete c_dsa;
    c_dsa = new ConstrainedDSA(*dsa, 0, dbc);
  }
  else{
    make_constrainedDSA();
  }
}

static const char* topMes[] = {
" ... Writing TOPDOMDEC File         ... \n",
" ... Writing renumbered TOP File    ... \n",
" ... Writing material TOP File      ... \n",
"",
"",
" ... Writing axisymetric TOP File with 2D planes ... \n",
" ... Writing axisymetric 3D TOP File ... \n",
" ... Writing renumbered material TOP File      ... \n"
};

void Domain::writeTopFileElementSets(ControlInfo *cinfo, int * nodeTable, int* nodeNumber, int topFlag)
{
  // ... WRITE ELEMENT CONNECTIVITY
 if(packedEset.last() > 0)  // PJSA: replaced numele with packedEset.last() to include phantoms
   fprintf(cinfo->checkfileptr,"Elements %s using %s\n",
           cinfo->elemSetName, cinfo->nodeSetName);

 int inode, iele;
 int nEls = packedEset.last();

 list<int> phantoms;
 list<int> constraints;

 for(iele=0; iele<nEls; ++iele) {
   if(!packedEset[iele]->isPhantomElement() && !packedEset[iele]->isConstraintElement()) {
     int numNodesPerElement = packedEset[iele]->numTopNodes();
     if(numNodesPerElement <= 1) continue;
     int eletype = packedEset[iele]->getTopNumber();
     int eid = (topFlag == 1 || topFlag == 7) ? iele+1 : packedEset[iele]->getGlNum()+1; // only renumber for -T and -M
     fprintf(cinfo->checkfileptr,"%6d  %4d ",eid,eletype);
     packedEset[iele]->nodes(nodeNumber);
     for(inode=0; inode<numNodesPerElement; ++inode)
       // Avoid to print nodes that are internally created
       if(nodeNumber[inode] < numnodes && nodes[nodeNumber[inode]] != 0)
	 fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber[inode]]);
     fprintf(cinfo->checkfileptr,"\n");
   }
   else
     {
       if(packedEset[iele]->isPhantomElement())
	 phantoms.push_back(iele); //TG create a list of phantom elements, shouldn't be too big
       else
	 constraints.push_back(iele);
     }
 }

 std::cerr << " ... " << phantoms.size() << " phantoms " << constraints.size() << " constraints elements. " << endl;

 //TG output dimasses in a separate element set if there are any
 // as of now (7/5/06) a dimass will be represented by a bar element in xpost between the node the
 // dimass is attached to and itself.
 if(nDimass > 0)
   {
     fprintf(stderr, " ... Putting %d Dimasses as bars in a separate ElemSet\n", nDimass);
     DMassData* curMass = firstDiMass;
     fprintf(cinfo->checkfileptr,"Elements %s_dimass using %s\n",
	     cinfo->elemSetName, cinfo->nodeSetName);
     int iele = 0;
     int eletype = 506; //101; // BAR element //TG now 506 to detect dimasses in Xpost
     while(curMass != NULL)
       {
	 int node = curMass->node;
	 if(node < numnodes && nodes[node] != 0)
	   {
	     fprintf(cinfo->checkfileptr,"%6d  %4d ",iele+1,eletype);

	     fprintf(cinfo->checkfileptr,"%6d %6d\n",nodeTable[node], nodeTable[node]);
	     curMass = curMass->next;
	     iele ++;
	   }
	 else
	   {
	     cout << " Warning : virtual dimass" << endl;
             curMass = curMass->next;
	   }
       }
   }

 // output phantom elements in a separate elementset if there are any
 if (phantoms.size() > 0) {
     fprintf(cinfo->checkfileptr,"Elements %s_phantom using %s\n",
	     cinfo->elemSetName, cinfo->nodeSetName);
     int m_phantoms = phantoms.size();
     for(int i=0; i<m_phantoms; ++i){
       iele = phantoms.front();
       phantoms.pop_front();

       //** same as main element writing in loop above for non phantom elements
       int numNodesPerElement = packedEset[iele]->numTopNodes();
       if(numNodesPerElement <= 1) continue;
       int eletype = packedEset[iele]->getTopNumber();
       fprintf(cinfo->checkfileptr,"%6d  %4d ",iele+1,eletype);
       packedEset[iele]->nodes(nodeNumber);
       for(inode=0; inode<numNodesPerElement; ++inode)
	 // Avoid to print nodes that are internally created
	 if(nodeNumber[inode] < numnodes && nodes[nodeNumber[inode]] != 0)
	   fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber[inode]]);
       fprintf(cinfo->checkfileptr,"\n");
       //**
     }
   }

 // output constraint elements in a separate elementset if there are any
 if (constraints.size() > 0) {
     fprintf(cinfo->checkfileptr,"Elements %s_constraints using %s\n",
	     cinfo->elemSetName, cinfo->nodeSetName);
     int m_constraints = constraints.size();
     for(int i=0; i<m_constraints; ++i){
       iele = constraints.front();
       constraints.pop_front();

       //** same as main element writing in loop above for non constraints elements
       int numNodesPerElement = packedEset[iele]->numTopNodes();
       if(numNodesPerElement <= 1) continue;
       int eletype = packedEset[iele]->getTopNumber();
       fprintf(cinfo->checkfileptr,"%6d  %4d ",iele+1,eletype);
       packedEset[iele]->nodes(nodeNumber);
       for(inode=0; inode<numNodesPerElement; ++inode)
	 // Avoid to print nodes that are internally created
	 if(nodeNumber[inode] < numnodes && nodes[nodeNumber[inode]] != 0)
	   fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber[inode]]);
       fprintf(cinfo->checkfileptr,"\n");
       //**
     }
   }

 // output surface elements, each in a separate element set if there are any
 // for now, just output the vertices
 for(int iSurf=0; iSurf<nSurfEntity; iSurf++) {
   if(SurfEntities[iSurf]->GetId() == 0) continue;
   fprintf(cinfo->checkfileptr,"Elements surface_%d using %s\n",
           SurfEntities[iSurf]->GetId(), cinfo->nodeSetName);
   FaceElemSet &faceElemSet = SurfEntities[iSurf]->GetFaceElemSet();
   for(iele=0; iele<faceElemSet.last(); ++iele) {
     if(SurfEntities[iSurf]->GetIsShellFace() && iele%2==1) continue;
     int nVertices = faceElemSet[iele]->nVertices();
     int eletype;
     if(nVertices == 3) eletype = 104;
     else if(nVertices == 4) eletype = 2;
     else { cerr << "don't know xpost eletype for surface " << SurfEntities[iSurf]->GetId() << " element " << iele << " nVertices = " << nVertices << endl; continue; }
     int eleID = (SurfEntities[iSurf]->GetIsShellFace()) ? iele/2+1 : iele+1;
     fprintf(cinfo->checkfileptr,"%6d  %4d ",eleID,eletype);
     for(inode=0; inode<nVertices; ++inode) {
       int nodeNumber = faceElemSet[iele]->GetVertex(inode);
       fprintf(cinfo->checkfileptr,"%6d ",nodeTable[nodeNumber]);
     }
     fprintf(cinfo->checkfileptr,"\n");
   }
 }

 // PJSA 3-24-05 write default pattern for ElemScalar output
  if(nEls > 0)
   fprintf(cinfo->checkfileptr,"Pattern default using %s\n", cinfo->elemSetName);

}

void
Domain::makeNodeTable(int topFlag)
{
 long m4 = - memoryUsed();
 nodeTable = new int[numnodes];
 m4 += memoryUsed();
 //fprintf(stderr," ... Node Table %14.3f Mb\n",m4/(1024.0*1024.0));

 int inode;
 for(inode=0; inode<numnodes; ++inode)
   nodeTable[inode] = -1;

 // if topFlag = 0, output a topdomdec file with all nodes
 // including nodes that are not defined or used.

 // if topFlag = 1, output a topdomdec file without gaps
 // (nodes renumbered sequentially)

 // if topFlag = 2, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and with all nodes including nodes that are not defined or used.

 // if topFlag = 7, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and without gaps (nodes renumbered sequentially)

 exactNumNodes = 0;
 for(inode=0; inode<numnodes; ++inode) {
   if(nodes[inode] == 0) continue;
   exactNumNodes += 1;
   nodeTable[inode] = ((topFlag == 1) || (topFlag == 7)) ? exactNumNodes : inode+1;
 }
}

void
Domain::makeTopFile(int topFlag)
{
 //MODIFIED FOR HEV PROBLEM, EC, 20070820
 if (solInfo().HEV)  {
   solInfo().HEV = 0;
   packedEset.deleteElems();
   numele = geoSource->getElems(packedEset);
 }

 fprintf(stderr," ... Memory Used so far    %14.3f Mb\n",
                 memoryUsed()/(1024.0*1024.0));

 if(topFlag >= 0)
   fprintf(stderr,"%s",topMes[topFlag]);

 // ... CONSTRUCT DOMAIN ELEMENT TO NODE CONNECTIVITY
 long m1 = - memoryUsed();
 elemToNode = new Connectivity( &packedEset );
 m1 += memoryUsed();

 // ... CONSTRUCT DOF SET ARRAY
 long m2 = - memoryUsed();
 dsa = new DofSetArray( numnodes, packedEset );
 m2 += memoryUsed();

 // Renumber rnum = getRenumbering();
 // nodeToNode->findMaxDist(rnum.renumb);
 // nodeToNode->findProfileSize(dsa);

 fprintf(stderr," ... Elem. to Node Connectivity %14.3f Mb\n",
         m1/(1024.0*1024.0));
 fprintf(stderr," ... DOF set array              %14.3f Mb\n",
         m2/(1024.0*1024.0));

 // ... WRITE INPUT FILE FOR TOP/DOMDEC
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 cinfo->checkfileptr = openFile(cinfo->checkfile, ".top");

 // Allocate memory for a node table to compact node numbers
 // into sequential format
   fprintf(stderr," ... Total Memory         = %14.3f Mb\n",
                memoryUsed()/(1024.0*1024.0));
 long m4 = - memoryUsed();

 int *nodeTable = new int[numnodes];
 m4 += memoryUsed();
 fprintf(stderr," ... Node Table %14.3f Mb\n",m4/(1024.0*1024.0));

 int inode;
 for(inode=0; inode<numnodes; ++inode)
   nodeTable[inode] = -1;

 // if topFlag = 0, output a topdomdec file with all nodes
 // including nodes that are not defined or used.

 // if topFlag = 1, output a topdomdec file without gaps
 // (nodes renumbered sequentially)

 // if topFlag = 2, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and with all nodes including nodes that are not defined or used.

 // if topFlag = 7, output a topdomdec file with each group of
 // elements with the same material output as a separate element set.
 // and without gaps (nodes renumbered sequentially)

 // ... WRITE NODE COORDINATES
 fprintf(cinfo->checkfileptr,"Nodes %s\n",cinfo->nodeSetName);
 int exactNumNodes = 0;
 for(inode=0; inode<numnodes; ++inode) {
   if(nodes[inode] == 0) continue;
   Node &nd = nodes.getNode(inode);
   exactNumNodes += 1;
   nodeTable[inode] = ((topFlag == 1) || (topFlag == 7)) ? exactNumNodes : inode+1;
   double x = nd.x;
   double y = nd.y;
   double z = nd.z;
   fprintf(cinfo->checkfileptr,"%d\t % 14.6f\t% 14.6f\t % 14.6f\n",
                              nodeTable[inode],x,y,z);
 }

 // First find the maximum number of nodes per Element
 // this used to be done in the write element connectivity loop
 maxNumNodes=0;
 int iele;
 int nEls = packedEset.last(); //HB: to avoid calling packedEset.last() at each loop
 for(iele=0; iele<nEls; ++iele) {
   int numNodesPerElement = packedEset[iele]->numNodes();
   maxNumNodes = myMax(numNodesPerElement, maxNumNodes);
 }
 // allocate integer array to store node numbers
 int *nodeNumber = new int[maxNumNodes];

 // ... WRITE ELEMENT CONNECTIVITY
 // TG moved to this function as will need to be outputted in .top file for -m and -M as well.
 writeTopFileElementSets(cinfo, nodeTable, nodeNumber, topFlag);
 //
 // ElemSet creation happens in there
 // ElemSet_phantom, ElemSet_dimass and ElementSet_constraints aswell

 // write structure displacement boundary conditions (dirichlet)
 // if these type of boundary conditions exist in the input file.
 int i;
 if(numDirichlet + numComplexDirichlet > 0)
   fprintf(cinfo->checkfileptr,"SDBoundary %s using %s\n",
           cinfo->bcondSetName,cinfo->nodeSetName);
 for(i=0; i<numDirichlet; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[dbc[i].nnum],dbc[i].dofnum+1,dbc[i].val);
 }

/*
 // ... ADD CONVECTIVE FLUXES
 // KHP: get the correct header for the TOPDOMDEC file!
 if(numConvBC > 0)
   fprintf(cinfo->checkfileptr,"SDTemperature %s using %s\n",
           cinfo->bcondSetName,cinfo->nodeSetName);
 for(i=0; i<numConvBC; ++i) {
  fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[cvbc[i].nnum], cvbc[i].dofnum+1,cvbc[i].val);
  }
*/


 // also print complex dirichlet boundary conditions, for Helmholtz problem
 for(i=0; i<numComplexDirichlet; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[cdbc[i].nnum],cdbc[i].dofnum+1,cdbc[i].reval);
 }

 // write structure force boundary conditions (Neumann)
 // if these type of boundary conditions exist in the input file.
 if(numNeuman + numComplexNeuman > 0)
   fprintf(cinfo->checkfileptr,"SFBoundary %s using %s\n",
           cinfo->bcondSetName,cinfo->nodeSetName);

 for(i=0; i<numNeuman; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[nbc[i].nnum],nbc[i].dofnum+1,nbc[i].val);
 }

 // also write complex Neumann boundary conditions, for Helmholtz problem
 // note we are only outputing the real value of the complex bc as topdomdec
 // will not be able to understand a complex entry (i.e. 2 double values
 // instead of 1)
 for(i=0; i<numComplexNeuman; ++i) {
   fprintf(cinfo->checkfileptr,"%d\t%d\t%f\n",
           nodeTable[cnbc[i].nnum],cnbc[i].dofnum+1,cnbc[i].reval);
 }

 // If we have more than one component, output a decomposition file
 // specifing the elements in each component. Very useful for debugging
 // a model.

 int *hasAppeared = new int[packedEset.last()];
 for(i=0; i<nEls; ++i)
   hasAppeared[i] = 0;

 // number of components in the structure input file
 //int nComponents = renumb.numComp;
 int nComponents = 1;

 // always output component file if there is more than one component.

 if(nComponents > 1) {

   FILE *componentFile = fopen("component.dec","w");

   int count = 0;

   fprintf(componentFile,"Decomposition COMPONENTS for %s\n%d\n",
           cinfo->elemSetName,nComponents);
   int n;
   for(n = 0; n<nComponents; ++n) {
     for(i = renumb.xcomp[n]; i<renumb.xcomp[n+1]; ++i) {
       inode = renumb.order[i];
       int numElemPerNode = nodeToElem->num(inode);
       int j;
       for(j=0; j<numElemPerNode; ++j) {
         if(hasAppeared[(*nodeToElem)[inode][j]] == 0) {
           hasAppeared[(*nodeToElem)[inode][j]] = (*nodeToElem)[inode][j]+1;
           count++;
         }
       }
     }
     fprintf(componentFile,"%d\n",count);
     for(i=0; i<nEls; ++i) {
       if(hasAppeared[i] != 0) {
         fprintf(componentFile," %d\n",hasAppeared[i]);
         hasAppeared[i] = 0;
       }
     }
     count = 0;
   }
   fflush(componentFile);
 }

 delete [] hasAppeared;

// FOR DANIEL - to output element lists for each attribute

 if((topFlag == 2) || (topFlag == 7)) {
   FILE *matList = fopen("material.top","w");

   // ... WRITE NODE COORDINATES
   fprintf(matList,"Nodes %s\n", cinfo->nodeSetName);
   for(inode=0; inode<numnodes; ++inode) {
     if(nodes[inode] == 0) continue;
     Node &nd = nodes.getNode(inode);
     // if(topFlag && nodes[inode] != 0) continue;
     double x = 0.0;
     double y = 0.0;
     double z = 0.0;
     if(nodes[inode] != 0) {
       x = nd.x;
       y = nd.y;
       z = nd.z;
     }
     fprintf(matList,"%d\t % 14.6f\t% 14.6f\t % 14.6f\n",
             nodeTable[inode],x,y,z);
   }

   int n;
   map<int, Attrib> &attrib = geoSource->getAttributes();
   int na = geoSource->getNumAttributes(); // this is actually the number of elements !!
   SPropContainer &sProps = geoSource->getStructProps();
   for(map<int, StructProp>::iterator it = sProps.begin(); it != sProps.end(); ++it) {
       n = it->first;
       bool first = true; // PJSA to deal with case of empty EleSet
       for(i=0; i<na; ++i) {
         if(attrib[i].attr == n) {
           int e = attrib[i].nele; 
           Element *elem = geoSource->getElem(e);
           if(elem) {
             if(first) {
               fprintf(matList,"Elements EleSet%d using %s\n", n+1, cinfo->nodeSetName);
               first = false;
             }
             int eletype = elem->getTopNumber();
             fprintf(matList,"%d %d ",e+1,eletype);
             int numNodesPerElement = elem->numTopNodes();
             if(numNodesPerElement <= 1) continue;
             for(inode=0; inode<numNodesPerElement; ++inode) {
               elem->nodes(nodeNumber);
               fprintf(matList,"%d ",nodeTable[nodeNumber[inode]]);
             }
             fprintf(matList,"\n");
           }
         }
       }
       if (output_match_in_top)// Match case for ElemScalar
	 fprintf(matList,"Match EleSet%d with default\n", n+1);
   }

   // PJSA 9-2-03 write boundary conditions to material.top
   if(numDirichlet + numComplexDirichlet > 0)
     fprintf(matList,"SDBoundary %s using %s\n",
             cinfo->bcondSetName, cinfo->nodeSetName);
   for(i=0; i<numDirichlet; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[dbc[i].nnum],dbc[i].dofnum+1,dbc[i].val);
   }

   for(i=0; i<numComplexDirichlet; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[cdbc[i].nnum],cdbc[i].dofnum+1,cdbc[i].reval);
   }
   if(numNeuman + numComplexNeuman > 0)
     fprintf(matList,"SFBoundary %s using %s\n",
             cinfo->bcondSetName,cinfo->nodeSetName);
   for(i=0; i<numNeuman; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[nbc[i].nnum],nbc[i].dofnum+1,nbc[i].val);
   }
   for(i=0; i<numComplexNeuman; ++i) {
     fprintf(matList,"%d\t%d\t%f\n",
             nodeTable[cnbc[i].nnum],cnbc[i].dofnum+1,cnbc[i].reval);
   }

 }

 // print model statistics to screen
 //printStatistics();

 // ... CALCULATE STRUCTURE MASS IF REQUESTED
 if(sinfo.massFlag)  {
   Renumber rnum = getRenumbering();
   double mass = computeStructureMass();
   fprintf(stderr," --------------------------------------\n");
   fprintf(stderr," ... Structure mass = %e  ...\n",mass);
   fprintf(stderr," --------------------------------------\n");
 }

 // Compute Total memory requested for constructing the TOP/DOMDEC file
 long m3 = memoryUsed();
 fprintf(stderr," ... Total Memory         = %14.3f Mb\n",m3/(1024.0*1024.0));
}

void
Domain::makeAxiTopFile(int topFlag, int numSlices) {

 double pi=4*atan(1.0);
 double angle;
 int j, inode, iele;
 //int maxNode = numdof();

 fprintf(stderr," ... Memory Used so far    %14.3f Mb\n",
                 memoryUsed()/(1024.0*1024.0));

 fprintf(stderr,"%s",topMes[topFlag]);

 // ... CONSTRUCT DOMAIN ELEMENT TO NODE CONNECTIVITY
 long m1 = - memoryUsed();
 elemToNode = new Connectivity( &packedEset );
 nodeToElem = elemToNode->reverse();
 nodeToNode = nodeToElem->transcon(elemToNode);
 m1 += memoryUsed();

 // ... CONSTRUCT DOF SET ARRAY
 long m2 = - memoryUsed();
 dsa = new DofSetArray( numNodes(), packedEset );
 m2 += memoryUsed();

 fprintf(stderr," ... Elem. to Node Connectivity %14.3f Mb\n",
         m1/(1024.0*1024.0));
 fprintf(stderr," ... DOF set array              %14.3f Mb\n",
         m2/(1024.0*1024.0));

 // ... WRITE INPUT FILE FOR TOP/DOMDEC
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 cinfo->checkfileptr = openFile(cinfo->checkfile, ".top");

  fprintf(stderr," ... Total Memory         = %14.3f Mb\n",
                memoryUsed()/(1024.0*1024.0));

 // ... WRITE NODE COORDINATES
 fprintf(cinfo->checkfileptr,"Nodes %s\n",cinfo->nodeSetName);
 int exactNumNodes = 0;
 for (j=0; j<numSlices; ++j) {

   angle = j*2*pi/numSlices;

   for(inode=0; inode<numNodes(); ++inode) {
      Node &nd = nodes.getNode(inode);
      if (nodes[inode] == 0) continue;
      exactNumNodes += 1;
      double x = nd.x*cos(angle);
      double y = nd.x*sin(angle);
      double z = nd.y;
      fprintf(cinfo->checkfileptr,"%d\t % 14.6f\t% 14.6f\t % 14.6f\n",
                              exactNumNodes,x,y,z);
   }
 }
 // ... WRITE ELEMENT CONNECTIVITY
 if (numele > 0)
   fprintf(cinfo->checkfileptr,"Elements %s using %s\n",
           cinfo->elemSetName,cinfo->nodeSetName);
 int exactNumEle=0;

/*
 fprintf(stderr,"\n ---------- PROBLEM PARAMETERS --------");
 fprintf(stderr,"\n ... # nodes              = %7d ...",numNodes());
 fprintf(stderr,"\n ... # elements           = %7d ...",numElements());
 fprintf(stderr,"\n ... # unconstrained dofs = %7d ...",numUncon());
 fprintf(stderr,"\n ... total # dofs         = %7d ...",numdof());
 fprintf(stderr,"\n --------------------------------------\n");
*/
 for (iele=0; iele<numele; ++iele) {
   AxiHElement *elem = dynamic_cast<AxiHElement *>(packedEset[iele]);
   if (elem == 0)  {
     int one=1;
     fprintf(stderr,"Element chosen non axisymmetric for TOP file. Aborting \n");
     exit(one);
   }
   if (topFlag == 5)
       elem->buildMesh2D(exactNumEle,cinfo->checkfileptr,numdof(),numSlices);
   if (topFlag == 6)
       elem->buildMesh3D(exactNumEle,cinfo->checkfileptr,numdof(),numSlices);
 }

 // Compute Total memory requested for constructing the TOP/DOMDEC file
 long m3 = memoryUsed();
  fprintf(stderr," ... Total Memory         = %14.3f Mb\n",m3/(1024.0*1024.0));
}


void
Domain::setsizeSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1)  sizeSfemStress = geoSource->numNode();  // node-based output
  else if(avgnum == 0) {  // element-based output
   sizeSfemStress = 0;
   Connectivity *elemToNode = new Connectivity(domain->getEset());
   int numele = geoSource->getNumAttributes();  // number of elements; another option domain->numElements();
   for(int iele=0; iele<numele; ++iele)   {
//     cerr << "number of nodes in this element  = " << elemToNode->num(iele) << endl;
     sizeSfemStress = sizeSfemStress + elemToNode->num(iele); // add number of nodes for each element
   }
  }
  else {
   cerr << "avgnum = " << avgnum << " not implemented in Domain::setsizeSfemStress()" << endl;
   sizeSfemStress = 0;
  }
}

double*
Domain::getSfemStress(int fileNumber, double* dummystress)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  if(avgnum == 1)  return stress->data();
  else if(avgnum == 0) return stressAllElems->data();
  else {cerr << "avgnum = " << avgnum << " not implemented in Domain::getSfemStress()" << endl; return 0;}
}


void
Domain::updateSfemStress(double* str, int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  if(avgnum == 1)  for (int i=0;i<stress->size();++i) (*stress)[i] = str[i];
  else if(avgnum == 0) for (int i=0;i<stressAllElems->size();++i) (*stressAllElems) = str[i]; // YYY DG
  else {cerr << "avgnum = " << avgnum << " not implemented in Domain::updateSfemStress()" << endl;}
}

void
Domain::getStressStrain(Vector &sol, double *bcx, int fileNumber,
                        int stressIndex, double time, int printFlag)
{
  int numNodes = (outFlag) ? exactNumNodes : geoSource->numNode();

  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS, WHERE
  // ... ONLY THE AXIAL STRAIN (EXX) AND AXIAL STRESS (SXX) ARE CALCULATED

  OutputInfo *oinfo = geoSource->getOutputInfo();

  // avgnum = 2 --> do not include stress/strain of bar/beam element in averaging
  int avgnum = oinfo[fileNumber].averageFlg;

  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;
  int surface = oinfo[fileNumber].surface;
    // upper  surface = 1
    // median surface = 2
    // lower  surface = 3

  int k;
  int iele;

  double *nodalTemperatures = 0;
  // Either get the nodal temperatures from the input file or
  // from the thermal model
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  if(printFlag != 2) {
    // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
    if(avgnum != 0) {
      if(stress == 0) stress = new Vector(numNodes,0.0);
      if(weight == 0) weight = new Vector(numNodes,0.0);
    }
    else if(stressAllElems == 0) stressAllElems = new Vector(sizeSfemStress,0.0);
    if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);


    if((elstress == 0)||(elweight == 0)) {
      int NodesPerElement, maxNodesPerElement=0;
      for(iele=0; iele<numele; ++iele) {
        NodesPerElement = elemToNode->num(iele);
        maxNodesPerElement = myMax(maxNodesPerElement, NodesPerElement);
      }
      if(elstress == 0) elstress = new Vector(maxNodesPerElement, 0.0);
      if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
    }

    if(avgnum != 0) {
    // zero the vectors
      stress->zero();
      weight->zero();
    }
    else if (printFlag == 1) stressAllElems->zero();
  }

  int count = 0;
  for(iele = 0; iele < numele; ++iele) {

    int NodesPerElement = elemToNode->num(iele);
    packedEset[iele]->nodes(nodeNumbers);

    if(printFlag != 2) {

      // Don't do anything if element is a phantom or constraint
      if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isConstraintElement()) continue;

      // Don't include beams or bars in the averaging if nodalpartial (avgnum = 2) is requested
      if ((avgnum == 2 && packedEset[iele]->getElementType() == 6) ||
          (avgnum == 2 && packedEset[iele]->getElementType() == 7) ||
          (avgnum == 2 && packedEset[iele]->getElementType() == 1)) continue;

      elDisp->zero();
      elstress->zero();
      elweight->zero();

      // DETERMINE ELEMENT DISPLACEMENT VECTOR
      for (k=0; k < allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if (cn >= 0)
          (*elDisp)[k] = sol[cn];
        else
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
      }

      int iNode;
      if (sinfo.thermalLoadFlag || (sinfo.thermoeFlag>=0))
        for (iNode = 0; iNode < NodesPerElement; ++iNode) {
          if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
            elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
          else
            elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
        }

      // CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
      packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                    *elDisp, stressIndex, surface,
                                    elemNodeTemps.data(), ylayer, zlayer, avgnum);
      if(avgnum != 0) {
        // ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
        for(k = 0; k < NodesPerElement; ++k) {
          int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
          (*stress)[node] += (*elstress)[k];
          (*weight)[node] += (*elweight)[k];
        }
      }

    } // end of (printFlag != 2)

    // PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
    if(avgnum == 0) {
      int offset[2];
      offset[0] = 0;
      //offset[1] = NodesPerElement;
      offset[1] = packedEset[iele]->numTopNodes(); //HB 06-25-05: avoid the internal nodes for MpcElement

      if(printFlag == 0) {
        if(iele == 0)
          geoSource->outputElemStress(fileNumber, (double *) 0, 0, offset, time); // print time
        geoSource->outputElemStress(fileNumber, elstress->data(), 1, offset); // print stresses
      }
      if(printFlag == 1) {
        for(k = 0; k < NodesPerElement; ++k) {
          stressAllElems[count] = (*elstress)[k];
          count++;
        }
      }
      if(printFlag == 2) {
        if(iele == 0)
          geoSource->outputElemStress(fileNumber, (double *) 0, 0, offset, time); // print time
//        geoSource->outputElemStress(fileNumber, stressAllElems[count], 1, offset); // print stresses YYY DG
        count=count+NodesPerElement;
      }
    }


  } // end of the iele loop

  // AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
  // ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.
  if(avgnum == 1 || avgnum == 2) {

    if(printFlag != 2) {
    // assemble stress vector
     for(k = 0; k < numNodes; ++k)  {
       if((*weight)[k] == 0.0)
         (*stress)[k] = 0.0;
       else
         (*stress)[k] /= (*weight)[k];
     }
    }

    if(printFlag != 1) {
     if(oinfo[fileNumber].nodeNumber == -1)
       geoSource->outputNodeScalars(fileNumber, stress->data(), numNodes, time);
     else
       geoSource->outputNodeScalars(fileNumber, stress->data()+oinfo[fileNumber].nodeNumber, 1, time);
    }

  }
}


void
Domain::getPrincipalStress(Vector &sol, double *bcx, int fileNumber,
                         int stressIndex, double time)
{
  int numNodes = (outFlag) ? exactNumNodes : geoSource->numNode();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // set stress VS. strain for element subroutines
  int strInd;
  int strDir;
  if ((stressIndex==0)||(stressIndex==1)||(stressIndex==2)) {
    strInd = 0;
    strDir = stressIndex+1;
  }
  else if ((stressIndex==3)||(stressIndex==4)||(stressIndex==5)) {
    strInd = 1;
    strDir = stressIndex-2;
  }
  else {
    fprintf(stderr,"*** ERROR: Bad Principal Stress Direction ***\n");
    exit(-1);
  }

  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS

  int avgnum = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  int j,k;
  // ... OUTPUT NODE NUMBER
  int n = oinfo[fileNumber].nodeNumber;

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(p_stress == 0) p_stress = new FullM(numNodes,6);
  if(weight == 0) weight = new Vector(numNodes,0.0);
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  int iele;
  if((p_elstress == 0)||(elweight == 0)) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = myMax(maxNodesPerElement, NodesPerElement);
    }
    if(p_elstress == 0) p_elstress = new FullM(maxNodesPerElement,9);
    if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
  }

  // zero the vectors
  p_stress->zero();
  weight->zero();

  for(iele=0; iele<numele; ++iele) {

    elDisp->zero();
    p_elstress->zero();
    elweight->zero();

    int NodesPerElement = elemToNode->num(iele);

    packedEset[iele]->nodes(nodeNumbers);

// ... DETERMINE ELEMENT DISPLACEMENT VECTOR

    for(k=0; k<allDOFs->num(iele); ++k) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
      if(cn >= 0)
        (*elDisp)[k] = sol[cn];
      else
        (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
    }

    int iNode;
    if(sinfo.thermalLoadFlag)
      for(iNode=0; iNode<NodesPerElement; ++iNode) {
        if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }

// ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
    if(packedEset[iele]->getProperty())
      packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                     *elDisp, strInd, surface,
                                     elemNodeTemps.data());

// ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

    for(k = 0; k < NodesPerElement; ++k) {
      int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
      for(j = 0; j < 6; ++j) {
        (*p_stress)[node][j] += (*p_elstress)[k][j];
      }
      (*weight)[node] += (*elweight)[k];
    }

// ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
//     THIS WRITES THE CHOSEN PRINCIPAL STRESS FOR EACH ELEMENT
    if(avgnum == 0) {
      geoSource->outputNodeScalars(fileNumber, (*p_elstress)[0]+(5+strDir), 1);

    }
  }

// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

  if(avgnum == 1 || avgnum == 2) {

    if(n == -1) {
      for(k=0; k<numNodes; ++k) {
        if((*weight)[k] == 0.0) {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] = 0.0;
          }
        } else  {
          for(j=0; j<6; ++j) {
            (*p_stress)[k][j] /= (*weight)[k];
          }
        }
      }
    } else {
      if((*weight)[n] == 0.0) {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] = 0.0;
        }
      } else {
        for(j=0; j<6; ++j) {
          (*p_stress)[n][j] /= (*weight)[n];
        }
      }
    }

    // CALCULATE PRINCIPALS AT EACH NODE
    double svec[6], pvec[3];
    if(n == -1) {
      double *globalPVec = new double[numNodes];
      for(k=0; k<numNodes; ++k) {
        for (j=0; j<6; ++j)
          svec[j] = (*p_stress)[k][j];

        // Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        globalPVec[k] = pvec[strDir-1];
      }
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodes, time);
    }
    else {
      for (j=0; j<6; ++j)
        svec[j] = (*p_stress)[n][j];

      // Convert Engineering to Tensor Strains
      if(strInd != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      geoSource->outputNodeScalars(fileNumber, pvec+strDir-1, 1, time);
    }

  }

}

double
Domain::getNodalStressStrain(Vector &sol, double *bcx,
                             int node, int stressIndex, int surface)
{

  // ... Number of Connected Elements
  int nconele = nodeToElem->num(node);

  // Watch out for this, as the maximum number of nodes may increase!
  int maxNodesPerElement = 20;

  Vector elstress(maxNodesPerElement,0.0);
  Vector elweight(maxNodesPerElement,0.0);

  Vector elDisp(maxNumDOFs,0.0);

  double stress=0.0;
  double weight=0.0;

  int icon;
  for(icon=0; icon<nconele; ++icon) {

     int k;
     int iele = (*nodeToElem)[node][icon];
     int NodesPerElement = elemToNode->num(iele);

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          elDisp[k] = sol[cn];
        else
          elDisp[k] = bcx[(*allDOFs)[iele][k]];
     }

     // ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT

     packedEset[iele]->getVonMises(elstress, elweight, nodes,
                                   elDisp, stressIndex, surface);

     // ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

     for(k=0; k<NodesPerElement; ++k) {
        int actnod = (*elemToNode)[iele][k];
        if (actnod == node) {
           stress += elstress[k]*elweight[k];
           weight += elweight[k];
        }
     }
  }
  return (stress/weight);
}

void
Domain::getElementForces(Vector &sol, double *bcx, int fileNumber,
                         int forceIndex, double time)
{
  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();

  // ...REMARK:  FORCES ARE CALCULATED FOR BARS AND BEAMS CURRENTLY
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  // ... Watch out for this as the max number of nodes may increase
  int NodesPerElement = 2;

  FullM forces(numele, NodesPerElement);

  int i, iele;
  // note, we are reusing elstress here for an elemental force vector
  if(elstress == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = myMax(maxNodesPerElement, NodesPerElement);
    }
    elstress = new Vector(maxNodesPerElement, 0.0);
  }

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  for(iele=0; iele<numele; ++iele) {
    packedEset[iele]->nodes(nodeNumbers);

    // ... DETERMINE ELEMENT DISPLACEMENT VECTOR
    int k;
    for(k=0; k<allDOFs->num(iele); ++k) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
      if(cn >= 0)
        (*elDisp)[k] = sol[cn];
      else
        (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
    }

    // ... CALCULATE INTERNAL FORCE VALUE FOR EACH ELEMENT
    // ... taking into account the new temperature in case of thermal coupling
    if((sinfo.thermalLoadFlag || (sinfo.thermoeFlag>=0)) && packedEset[iele]->getProperty()) {
      int iNode;
      for(iNode=0; iNode<NodesPerElement; ++iNode) {
        if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
    }
    packedEset[iele]->getIntrnForce(*elstress,nodes,elDisp->data(),
                                    forceIndex,elemNodeTemps.data());

    // ... COPY ELEMENT'S NODAL FORCES INTO A TOTAL FORCE MATRIX
    for(i=0; i<NodesPerElement; ++i)
      forces[iele][i] = (*elstress)[i];

  }

   // ... PRINT THE ELEMENT FORCES TO A FILE
   geoSource->outputElemVectors(fileNumber, forces.data(), numele, time);
}

void
Domain::getElasticStiff(FullSquareMatrix *kelArray)
{
 int iele;
 for(iele=0; iele < numele; ++iele) {
    kelArray[iele] = packedEset[iele]->stiffness(nodes,kelArray[iele].data());
  }
}

void
Domain::getElasticForces(Vector &dsp, double *bcx, Vector &ext_f, double eta,
                         FullSquareMatrix *kelArray)
{

  int size    = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *kel = (double *) dbg_alloca(size);

  int iele;
  for(iele=0; iele<numele; ++iele) {

     //int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);

     Vector elDisp   (numEleDOFs,0.0);
     Vector elForce  (numEleDOFs,0.0);

     FullSquareMatrix karray(numEleDOFs,kel);

     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          elDisp[k] = dsp[cn]; }
        else {
          elDisp[k] = bcx[(*allDOFs)[iele][k]]; }
     }

     if (kelArray) {
       karray=kelArray[iele];
     }
     else {
       karray=packedEset[iele]->stiffness(nodes,karray.data());
     }

     int i,j;
     for (i=0;i<numEleDOFs;i++) {
       for (j=0;j<numEleDOFs;j++) {
         elForce[i] += karray[i][j]*elDisp[j] ;
       }
     }

     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          ext_f[cn] += eta * elForce[k]; }
     }
  }
}

void
Domain::getKtimesU(Vector &dsp, double *bcx, Vector &ext_f, double eta,
                   FullSquareMatrix *kelArray)
{
  int size = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *karray = (double *) dbg_alloca(size);

  for(int iele=0; iele<numele; ++iele) {

     int numEleDOFs = allDOFs->num(iele);

     Vector elDisp(numEleDOFs,0.0);
     Vector elForce(numEleDOFs,0.0);

     FullSquareMatrix kel(numEleDOFs,karray);

     for(int k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          elDisp[k] = dsp[cn];
        }
        else {
          //XXXXcerr << "bcx not added!!!\n";
          elDisp[k] = 0.0;
        }
     }

     if(kelArray) {
       //kel=kelArray[iele];
       kel.copy(kelArray[iele]); // PJSA 4-1-08
     }
     else {
       kel=packedEset[iele]->stiffness(nodes,karray);
     }

     for(int i=0;i<numEleDOFs;i++) {
       for(int j=0;j<numEleDOFs;j++) {
         elForce[i] += kel[i][j]*elDisp[j] ;
       }
     }

     for(int k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          ext_f[cn] += eta * elForce[k]; }
     }
  }
}

double
Domain::getStructureMass()
{
  double totmas = 0.0; // total mass of structure

  int iele;
  for(iele=0; iele < numele; ++iele) {
    double elementMass = packedEset[iele]->getMass(nodes);
    totmas += elementMass;
  }

  double *nodeCount = (double *) dbg_alloca(sizeof(double)*numnodes);
  double *nodeMass  = (double *) dbg_alloca(sizeof(double)*numnodes);

  int n;
  for(n=0; n<numnodes; ++n)
    nodeCount[n] = nodeMass[n] = 0.0;

  // ADD DISCRETE MASS
  DMassData *current = firstDiMass;
  while(current != 0) {
    int n = current->node;
    nodeMass[n]  += current->diMass;
    nodeCount[n] += 1;
    current = current->next;
  }

  for(n=0; n<numnodes; ++n) {
    if(nodeCount[n] == 0.0) continue;
    if(nodeCount[n] > 6) nodeCount[n] = 6;
    totmas += nodeMass[n]/nodeCount[n];
  }

   return totmas;
}

double
Domain::getStrainEnergy(Vector &sol, double *bcx, SparseMatrix * gStiff,
                        FullSquareMatrix *kelArray)
{
  double energy=0.0;

  // Evaluation of Strain Energy if global Stiffness Matrix is given

  if (gStiff) {

    Vector tmpVec(c_dsa->size());

    gStiff->mult(sol,tmpVec);

    energy = 0.5 * (sol * tmpVec);

    return energy;
  }

  // Evaluation of Strain Energy if local Stiffness Matices are given
  // or has to be determined

  int size    = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *kel = (double *) dbg_alloca(size);

  Vector elDisp (maxNumDOFs,0.0);
  Vector elKelD (maxNumDOFs,0.0);

  int iele;
  for(iele=0; iele<numele; ++iele) {

     //int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);

     elDisp.zero();
     elKelD.zero();

     FullSquareMatrix karray(numEleDOFs,kel);

     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          elDisp[k] = sol[cn];
        else
          elDisp[k] = bcx[(*allDOFs)[iele][k]];
     }

     if (kelArray) {
       karray=kelArray[iele];
     }
     else {
       karray=packedEset[iele]->stiffness(nodes,karray.data());
     }

     int i,j;
     for (i=0;i<numEleDOFs;i++) {
       for (j=0;j<numEleDOFs;j++) {
         elKelD[i]+=karray[i][j]*elDisp[j];
       }
     }

     energy += 0.5*(elKelD*elDisp);
  }
  return energy;
}

double
Domain::getNodalDisp(Vector &sol, double *bcx, int node, int dispTyp)
{
    int loc,loc1;
    double val;

    switch (dispTyp) {

    case 0:
      loc  = c_dsa->locate( node, DofSet::Xdisp);
      loc1 =   dsa->locate( node, DofSet::Xdisp);
      if      ( loc >= 0  ) { val = sol[loc];  }
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 1:
      loc  = c_dsa->locate( node, DofSet::Ydisp);
      loc1 =   dsa->locate( node, DofSet::Ydisp);
      if      ( loc >= 0  ) { val = sol[loc];  }
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 2:
      loc  = c_dsa->locate( node, DofSet::Zdisp);
      loc1 =   dsa->locate( node, DofSet::Zdisp);
      if      ( loc >= 0  ) { val = sol[loc];  }
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 3:
      loc  = c_dsa->locate( node, DofSet::Xrot);
      loc1 =   dsa->locate( node, DofSet::Xrot);
      if      ( loc >= 0  ) { val = sol[loc];  }
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 4:
      loc  = c_dsa->locate( node, DofSet::Yrot);
      loc1 =   dsa->locate( node, DofSet::Yrot);
      if      ( loc >= 0  ) { val = sol[loc];  }
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 5:
      loc  = c_dsa->locate( node, DofSet::Zrot);
      loc1 =   dsa->locate( node, DofSet::Zrot);
      if      ( loc >= 0  ) { val = sol[loc];  }
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 6:
      double dummy;
      loc  = c_dsa->locate( node, DofSet::Xdisp);
      loc1 =   dsa->locate( node, DofSet::Xdisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  }
      else if ( loc1 >= 0 ) { dummy = bcx[loc1]; }
      else                  { dummy = 0.0; }
      val += dummy*dummy;
      loc  = c_dsa->locate( node, DofSet::Ydisp);
      loc1 =   dsa->locate( node, DofSet::Ydisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  }
      else if ( loc1 >= 0 ) { dummy = bcx[loc1]; }
      else                  { dummy = 0.0; }
      val += dummy*dummy;
      loc  = c_dsa->locate( node, DofSet::Zdisp);
      loc1 =   dsa->locate( node, DofSet::Zdisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  }
      else if ( loc1 >= 0 ) { dummy = bcx[loc1]; }
      else                  { dummy = 0.0; }
      val += dummy*dummy;
      val =  sqrt(val);
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
    }
    return val;
}

