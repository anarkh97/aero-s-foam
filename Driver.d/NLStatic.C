#include <stdio.h>
#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>

// New include files for Restart file
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#include <Corotational.d/Corotator.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/mathUtility.h>
#include <Math.d/matrix.h>
#include <Control.d/ControlInterface.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>

#include <Driver.d/GeoSource.h>

void
Domain::getStiffAndForce(GeomState &geomState, Vector& elementForce,
		         Corotator **corotators, FullSquareMatrix *kel,
                         Vector &residual, double lambda, double time)
/*******************************************************************
 *
 * Purpose :
 *
 *  Compute element tangential stiffness and element internal force
 *  and assemble element internal force into global internal force.
 *  Also compute follower external force contribution to both 
 *  tangential stiffness matrix and residual
 *
 * Input :
 *
 *  corotators : element corotator array
 *  geomState  : current node geometric state
 *  nodes      : undeformed nodal coordinate set
 *
 * Output :
 *
 *  residual   : residual vector = external force - internal force
 *  kel        : array of element tangential stiffness matrices
 *               in current configuration
 *
 *****************************************************************/

{
  for(int iele = 0; iele < numele; ++iele) {

    elementForce.zero();

    // Get updated tangent stiffness matrix and element internal force
    if(corotators[iele]) {
      corotators[iele]->getStiffAndForce(geomState, nodes, kel[iele],
                                         elementForce.data(), sinfo.getTimeStep(), time);
    }
    // Compute k and internal force for an element with x translation (or temperature) dofs
    else if(solInfo().soltyp == 2) {
      kel[iele].zero();
      Vector temp(packedEset[iele]->numNodes());
      for(int i=0; i<packedEset[iele]->numNodes(); ++i) {
        temp[i] = geomState[packedEset[iele]->nodes()[i]].x;
      }
      kel[iele] = packedEset[iele]->stiffness(nodes, kel[iele].data());
      kel[iele].multiply(temp, elementForce, 1.0); // elementForce = kel*temp
    }
    // Assemble element internal force into residual force vector
    for(int idof = 0; idof < kel[iele].dim(); ++idof) {
      int dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
      if(dofNum >= 0)
        residual[dofNum] -= elementForce[idof];
    }
  }

  for(int iele = numele; iele < packedEset.size(); ++iele) {
    Corotator *c = dynamic_cast<Corotator*>(packedEset[iele]);
    if(c) {
      // TODO implicit
      FullSquareMatrix kelTmp(packedEset[iele]->numDofs());
      Vector elementForceTmp(packedEset[iele]->numDofs());  
      kelTmp.zero();
      elementForceTmp.zero();
      c->getStiffAndForce(geomState, nodes, kelTmp, elementForceTmp.getData(), sinfo.getTimeStep(), time);
      int *p = new int[packedEset[iele]->numDofs()];
      packedEset[iele]->dofs(*c_dsa, p);
      //cerr << "iele = " << iele << ", force = "; elementForceTmp.print();
      for(int idof = 0; idof < packedEset[iele]->numDofs(); ++idof) {
        if(p[idof] > -1) residual[p[idof]] -= elementForceTmp[idof];
      }
      delete [] p;
    }
  }

  if(domain->pressureFlag()) {
    double cflg = (sinfo.newmarkBeta == 0.0) ? 0.0 : 1.0;
    for(int iele = 0; iele < numele;  ++iele) {
      // If there is a zero pressure defined, skip the element
      if(packedEset[iele]->getPressure() == 0) continue;

      // Compute element pressure force in the local coordinates
      elementForce.zero();
      packedEset[iele]->computePressureForce(nodes, elementForce, &geomState, 1);
      elementForce *= lambda;
//#define PRESSURE_MFTT
#ifdef PRESSURE_MFTT
      double mfttFactor = (domain->mftval) ? domain->mftval->getVal(time) : 1.0;
      elementForce *= mfttFactor; // TODO consider
#endif
      // Include the "load stiffness matrix" in kel[iele]
      if(sinfo.newmarkBeta != 0.0)
        corotators[iele]->getDExternalForceDu(geomState, nodes, kel[iele],
                                              elementForce.data());

      // Determine the elemental force for the corrotated system
      corotators[iele]->getExternalForce(geomState, nodes, elementForce.data());

      // Assemble element pressure forces into residual force vector
      for(int idof = 0; idof < kel[iele].dim(); ++idof) {
        int dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
        if(dofNum >= 0)
          residual[dofNum] += elementForce[idof];
      }
    }
  }

  if(domain->thermalFlag()) {
    if(!temprcvd) initNodalTemperatures(); // XXXX to be moved
    Vector elementTemp(maxNumNodes);
    for(int iele = 0; iele < numele;  ++iele) {
      // By convention phantom elements do not have thermal load
      if(packedEset[iele]->getProperty() == 0) continue;

      // Extract the element nodal temperatures from temprcvd and/or element property ambient temperature
      for(int inod = 0; inod < elemToNode->num(iele); ++inod) {
        double t = temprcvd[(*elemToNode)[iele][inod]];
        elementTemp[inod] = (t == defaultTemp) ? packedEset[iele]->getProperty()->Ta : t;
      }

      // Compute element thermal force in the local coordinates
      elementForce.zero();
      packedEset[iele]->getThermalForce(nodes, elementTemp, elementForce, 1);
      elementForce *= lambda;

      // Include the "load stiffness matrix" in kel[iele]
      if(sinfo.newmarkBeta != 0.0)
        corotators[iele]->getDExternalForceDu(geomState, nodes, kel[iele],
                                              elementForce.data());

      // Determine the elemental force for the corrotated system
      corotators[iele]->getExternalForce(geomState, nodes, elementForce.data());

      // Assemble element thermal forces into residual force vector
      for(int idof = 0; idof < kel[iele].dim(); ++idof) {
        int dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
        if(dofNum >= 0)
          residual[dofNum] += elementForce[idof];
      }
    }
  }

  if(claw && claw->numActuator) {
    for(int i = 0; i < numNeuman; ++i) {
      int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
      if(dof < 0) continue;
      if(nbc[i].type == BCond::Actuators) residual[dof] += lambda*nbc[i].val; // XXXX need to multiply by weight for multidomain
    }
  }
 
  if(!solInfo().getNLInfo().unsymmetric)
    for(int iele = 0; iele < numele;  ++iele) 
      kel[iele].symmetrize();
    
}

// used in nonlinear statics
void
Domain::createKelArray(FullSquareMatrix *&kArray)
{
  // Allocate array of pointers to FullSquareMatrix to store the
  // element stiffness matrices
  kArray = new FullSquareMatrix[numele];

  // Allocate the correct size for each elements stiffness matrix
  int iele;
  for(iele = 0; iele<numele; ++iele)
    kArray[iele].setSize(packedEset[iele]->numDofs());
}

// used in nonlinear dynamics

void
Domain::createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray)
{

 // Allocate array of pointers to FullSquareMatrix to store
 // the element stiffness matrices and element mass matrices
 kArray = new FullSquareMatrix[numele];
 mArray = new FullSquareMatrix[numele];

 // Allocate the correct size for each elements stiffness & mass matrix
 int iele;
 for(iele = 0; iele<numele; ++iele) {
   int dimension = packedEset[iele]->numDofs();
   kArray[iele].setSize(dimension);
   mArray[iele].setSize(dimension);
 }

 // Form and store element mass matrices into an array
 for(iele=0; iele<numele; ++iele)
   mArray[iele] = packedEset[iele]->massMatrix(nodes, mArray[iele].data());

 // zero rotational degrees of freedom within element mass matrices
 int *dofType = dsa->makeDofTypeArray();

 int i,j;
 for(iele=0; iele<numele; ++iele) {
   for(i=0; i<mArray[iele].dim(); ++i)
     for(j=0; j<mArray[iele].dim(); ++j)
        if( dofType[ (*allDOFs)[iele][i] ] == 1 ||
            dofType[ (*allDOFs)[iele][j] ] == 1) {
           mArray[iele][i][j] = 0.0;
        }
 }

}

// used in nonlinear dynamics
void
Domain::createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray, FullSquareMatrix *&cArray)
{
 // Allocate array of pointers to FullSquareMatrix to store
 // the element stiffness matrices, element damping matrices and element mass matrices
 kArray = new FullSquareMatrix[numele];
 mArray = new FullSquareMatrix[numele];
 cArray = new FullSquareMatrix[numele];

 // Allocate the correct size for each elements stiffness & mass matrix
 int iele;
 for(iele = 0; iele<numele; ++iele) {
   int dimension = packedEset[iele]->numDofs();
   kArray[iele].setSize(dimension);
   mArray[iele].setSize(dimension);
   cArray[iele].setSize(dimension);
 }

 // Form and store element damping matrices and mass matrices into arrays
 for(iele=0; iele<numele; ++iele) {
   mArray[iele] = packedEset[iele]->massMatrix(nodes, mArray[iele].data());
   cArray[iele] = packedEset[iele]->dampingMatrix(nodes, cArray[iele].data());
 }

 // add Rayleigh damping
 int i,j;
 FullSquareMatrix kel;
 for(iele=0; iele<numele; ++iele) {
   kel.setSize(packedEset[iele]->numDofs());
   kel = packedEset[iele]->stiffness(nodes, kel.data());
   for(i=0; i<cArray[iele].dim(); ++i)
     for(j=0; j<cArray[iele].dim(); ++j)
       cArray[iele][i][j] += (sinfo.alphaDamp*mArray[iele][i][j] + sinfo.betaDamp*kel[i][j]);
 }

 // zero rotational degrees of freedom within element mass matrices and damping matrices
 int *dofType = dsa->makeDofTypeArray();
 for(iele=0; iele<numele; ++iele) {
   for(i=0; i<mArray[iele].dim(); ++i)
     for(j=0; j<mArray[iele].dim(); ++j)
        if( dofType[ (*allDOFs)[iele][i] ] == 1 ||
            dofType[ (*allDOFs)[iele][j] ] == 1) {
           mArray[iele][i][j] = 0.0;
           cArray[iele][i][j] = 0.0;
        }
 }

}

void
Domain::postProcessing(GeomState *geomState, Vector& force, Vector &aeroForce,
                       double time, int step, double* velocity, double *vcx,
                       Corotator **allCorot, FullSquareMatrix *mel, double *acceleration, double *acx)
{

  if(time == sinfo.initialTime) {
    geoSource->openOutputFiles();
    //printStatistics();
  }

  if( sinfo.nRestart > 0 && velocity !=0) {
    StackVector v_n(velocity, numUncon());
    writeRestartFile(time, step, v_n, geomState);
  }

  int numOutInfo = geoSource->getNumOutInfo();
  for(int iInfo = 0; iInfo < numOutInfo; ++iInfo)
  {
    postProcessingImpl(iInfo, geomState, force, aeroForce, time, step, velocity, vcx, allCorot, mel, acceleration, acx);
  }

}

void
Domain::postProcessingImpl(int iInfo, GeomState *geomState, Vector& force, Vector &aeroForce,
                           double time, int step, double* velocity, double *vcx,
                           Corotator **allCorot, FullSquareMatrix *mel, double *acceleration, double *acx)
{
 int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes

 //fprintf(stderr,"Running postPro in NLSTATIC\n");
 enum {SXX=0,SYY=1,SZZ=2,SXY= 3,SYZ= 4,SXZ= 5,VON=6,
       EXX=7,EYY=8,EZZ=9,EXY=10,EYZ=11,EXZ=12,STRAINVON=13,
       VONTOP=14,VONBOT=15};

 enum {INX,INY,INZ,AXM,AYM,AZM};

 enum {PSTRESS1=0,PSTRESS2=1,PSTRESS3=2,
       PSTRAIN1=3,PSTRAIN2=4,PSTRAIN3=5};

 /*if(time == sinfo.initialTime) {
   geoSource->openOutputFiles();
   printStatistics();
 }

// call to writeRestartFile for Nonlinear!
  if( sinfo.nRestart > 0 && velocity !=0) {
    StackVector v_n(velocity, numUncon());
    writeRestartFile(time, step, v_n, geomState);
  }*/

  int i; //, iInfo;
  //int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  //for(iInfo=0; iInfo < numOutInfo; ++iInfo) {

     // Check output interval
    if(step%oinfo[iInfo].interval != 0 && time != 0.0) return;

    int w = oinfo[iInfo].width;
    int p = oinfo[iInfo].precision;

    switch(oinfo[iInfo].type) {
       case OutputInfo::Displacement:
	 if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i = 0; i < numNodes; ++i) {
             //cerr << "i = " << i << ", nodes[i]     = " << nodes[i]->x << ", " << nodes[i]->y << ", " << nodes[i]->z << endl;
             //cerr << "i = " << i << ", geomState[i] = " << (*geomState)[i].x << ", " << (*geomState)[i].y << ", " << (*geomState)[i].z << endl;
             double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
             double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
             double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
             fprintf(oinfo[iInfo].filptr," % *.*E % *.*E % *.*E\n"
                                 ,w,p,x,w,p,y,w,p,z);
           }
         } else {
           // if only one node was requested for output
      	   // print 1 set of data per line
	   i = oinfo[iInfo].nodeNumber;
           double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
           double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
           double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
      	   fprintf(oinfo[iInfo].filptr,"  % *.*E   % *.*E % *.*E % *.*E\n"
                                 ,w,p,time,w,p,x,w,p,y,w,p,z);
         }
         fflush(oinfo[iInfo].filptr);
         break;
      case OutputInfo::Temperature:
         if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double x = (nodes[i]&& i < geomState->numNodes()) ? (*geomState)[i].x : 0;
             fprintf(oinfo[iInfo].filptr," % *.*E\n"
                                 ,w,p,x);
           }
         } else {
           // if only one node was requested for output
           // print 1 set of data per line
           i = oinfo[iInfo].nodeNumber;
           double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x : 0;
           fprintf(oinfo[iInfo].filptr,"  % *.*E   % *.*E\n"
                                 ,w,p,time,w,p,x);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::Disp6DOF:
         // 6 dof output should include node number
         if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
             double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
             double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
             double rot[3];
             mat_to_vec((*geomState)[i].R,rot);
             fprintf(oinfo[iInfo].filptr,
               "%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n", i+1
               ,w,p,x,w,p,y,w,p,z,w,p,rot[0],w,p,rot[1],w,p,rot[2]);
           }
         } else {
           // if only one node was requested for output
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
           double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
           double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
           double rot[3];
           mat_to_vec((*geomState)[i].R,rot);
           fprintf(oinfo[iInfo].filptr,
             "%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n", i+1, w, p, x, w, p, y,
              w,p,z,w,p,rot[0],w,p,rot[1],w,p,rot[2]);
         }
         fflush(oinfo[iInfo].filptr);
         break;

       case OutputInfo::Velocity6:               // all nodes or only one
       case OutputInfo::Velocity:
       {
         StackVector v_n(velocity, numUncon());
         int first_node, last_node;
         if (oinfo[iInfo].nodeNumber == -1) { fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time); first_node=0; last_node=numNodes; }
         else { fprintf(oinfo[iInfo].filptr,"  % *.*E",w,p,time); first_node=oinfo[iInfo].nodeNumber; last_node=first_node+1; }

         for (int iNode=first_node; iNode<last_node; ++iNode)  {

           double x,y,z;
           double xr,yr,zr;
           getOrAddDofForPrint(false, v_n, vcx, iNode, &x, &DofSet::Xdisp, &y, &DofSet::Ydisp, &z, &DofSet::Zdisp);

           if (oinfo[iInfo].type == OutputInfo::Velocity6)
             getOrAddDofForPrint(false, v_n, vcx, iNode, &xr, &DofSet::Xrot, &yr, &DofSet::Yrot, &zr, &DofSet::Zrot);

           if (oinfo[iInfo].type == OutputInfo::Velocity6)  {
             fprintf(oinfo[iInfo].filptr, "%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                                          iNode+1,w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
           } else
             fprintf(oinfo[iInfo].filptr," % *.*E % *.*E % *.*E\n",w,p,x,w,p,y,w,p,z);
         }
         fflush(oinfo[iInfo].filptr);
       } break;

       case OutputInfo::TemperatureFirstTimeDerivative: {
         StackVector v_n(velocity, numUncon());
         int first_node, last_node;
         if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n", w, p, time);
           first_node = 0;
           last_node = numNodes;
         }
         else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E", w, p, time);
           first_node = oinfo[iInfo].nodeNumber;
           last_node = first_node+1;
         }
         for(int iNode = first_node; iNode < last_node; ++iNode) {
           double x;
           getOrAddDofForPrint(false, v_n, vcx, iNode, &x, &DofSet::Temp);
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,x);
         }
         fflush(oinfo[iInfo].filptr);
       } break;

       case OutputInfo::Accel6:
       case OutputInfo::Acceleration:
       {
         StackVector a_n(acceleration, numUncon()); // XXXX acceleration not passed
         int first_node, last_node;
         if (oinfo[iInfo].nodeNumber == -1) { fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time); first_node=0; last_node=numNodes; }
         else { fprintf(oinfo[iInfo].filptr,"  % *.*E",w,p,time); first_node=oinfo[iInfo].nodeNumber; last_node=first_node+1; }

         for (int iNode=first_node; iNode<last_node; ++iNode) {

           double x,y,z;
           double xr,yr,zr;
           getOrAddDofForPrint(false, a_n, acx, iNode, &x, &DofSet::Xdisp, &y, &DofSet::Ydisp, &z, &DofSet::Zdisp);

           if (oinfo[iInfo].type == OutputInfo::Accel6)
             getOrAddDofForPrint(false, a_n, acx, iNode, &xr, &DofSet::Xrot, &yr, &DofSet::Yrot, &zr, &DofSet::Zrot);

           if (oinfo[iInfo].type == OutputInfo::Accel6)  {
             fprintf(oinfo[iInfo].filptr, "%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                                          iNode+1,w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
           } else
             fprintf(oinfo[iInfo].filptr," % *.*E % *.*E % *.*E\n",w,p,x,w,p,y,w,p,z);
         }
         fflush(oinfo[iInfo].filptr);
       }  break;

       case OutputInfo::DispX:
         if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,x);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,x);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::DispY:
         if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,y);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,y);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::DispZ:
         if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,z);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,z);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::RotX:
          if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double rot[3];
             mat_to_vec((*geomState)[i].R,rot);
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,rot[0]);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double rot[3];
           mat_to_vec((*geomState)[i].R,rot);
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,rot[0]);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::RotY:
          if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double rot[3];
             mat_to_vec((*geomState)[i].R,rot);
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,rot[1]);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double rot[3];
           mat_to_vec((*geomState)[i].R,rot);
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,rot[1]);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::RotZ:
          if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double rot[3];
             mat_to_vec((*geomState)[i].R,rot);
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,rot[2]);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double rot[3];
           mat_to_vec((*geomState)[i].R,rot);
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,rot[2]);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::DispMod:
        if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
             double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
             double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
             double mod = sqrt(x*x+y*y+z*z);
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,mod);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double x = (nodes[i]) ? (*geomState)[i].x - nodes[i]->x : 0;
           double y = (nodes[i]) ? (*geomState)[i].y - nodes[i]->y : 0;
           double z = (nodes[i]) ? (*geomState)[i].z - nodes[i]->z : 0;
           double mod = sqrt(x*x+y*y+z*z);
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,mod);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::RotMod:
         if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double rot[3];
             mat_to_vec((*geomState)[i].R,rot);
             double mod = sqrt(rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,mod);
           }
         } else {
           // if only one node was requested for output
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double rot[3];
           mat_to_vec((*geomState)[i].R,rot);
           double mod = sqrt(rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,mod);
         }
         fflush(oinfo[iInfo].filptr);
         break;
      case OutputInfo::TotMod:
        if(oinfo[iInfo].nodeNumber == -1) {
           fprintf(oinfo[iInfo].filptr,"  % *.*E\n",w,p,time);
           for(i=0; i < numNodes; ++i) {
             double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
             double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
             double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
             double rot[3];
             mat_to_vec((*geomState)[i].R,rot);
             double mod = sqrt(x*x+y*y+z*z+rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
             fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,mod);
           }
         } else {
           fprintf(oinfo[iInfo].filptr,"  % *.*E  ",w,p,time);
           i = oinfo[iInfo].nodeNumber;
           double x = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].x - nodes[i]->x : 0;
           double y = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].y - nodes[i]->y : 0;
           double z = (nodes[i] && i < geomState->numNodes()) ? (*geomState)[i].z - nodes[i]->z : 0;
           double rot[3];
           mat_to_vec((*geomState)[i].R,rot);
           double mod = sqrt(x*x+y*y+z*z+rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);
           fprintf(oinfo[iInfo].filptr," % *.*E\n",w,p,mod);
         }
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::Rigid:
/*
         Vector rigid(maxNumDOFs);
         Vector rDisp(numUncon(),0.0);
         for(i=0; i<numele; ++i) {
           allCorot[i]->extractRigidMotion(geomState, nodes, rigid.data());
           int NodesPerElement = packedEset[iele]->numNodes();
           for(k=0; k<NodesPerElement; ++k) {
             rDisp[(*elemToNode)[iele][k]] += rigid[k];
             (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
           }
         }
*/
       case OutputInfo::StressXX:
         getStressStrain( *geomState, allCorot,  iInfo, SXX, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressYY:
         getStressStrain( *geomState, allCorot,  iInfo, SYY, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressZZ:
         getStressStrain( *geomState, allCorot,  iInfo, SZZ, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressXY:
         getStressStrain( *geomState, allCorot,  iInfo, SXY, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressYZ:
         getStressStrain( *geomState, allCorot,  iInfo, SYZ, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressXZ:
         getStressStrain( *geomState, allCorot,  iInfo, SXZ, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainXX:
         getStressStrain( *geomState, allCorot,  iInfo, EXX, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainYY:
         getStressStrain( *geomState, allCorot,  iInfo, EYY, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainZZ:
         getStressStrain( *geomState, allCorot,  iInfo, EZZ, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainXY:
         getStressStrain( *geomState, allCorot,  iInfo, EXY, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainYZ:
         getStressStrain( *geomState, allCorot,  iInfo, EYZ, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainXZ:
         getStressStrain( *geomState, allCorot,  iInfo, EXZ, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressVM:
         getStressStrain( *geomState, allCorot,  iInfo, VON, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainVM:
         getStressStrain( *geomState, allCorot,  iInfo, STRAINVON, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressPR1:
         getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS1, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressPR2:
         getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS2, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StressPR3:
         getPrincipalStress(*geomState,allCorot,iInfo,PSTRESS3, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainPR1:
         getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN1, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainPR2:
         getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN2, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::StrainPR3:
         getPrincipalStress(*geomState,allCorot,iInfo,PSTRAIN3, time);
         fflush(oinfo[iInfo].filptr);
         break;
       case OutputInfo::InXForce:
         getElementForces(*geomState, allCorot, iInfo, INX, time);
         break;
       case OutputInfo::InYForce:
         getElementForces(*geomState, allCorot, iInfo, INY, time);
         break;
       case OutputInfo::InZForce:
         getElementForces(*geomState, allCorot, iInfo, INZ, time);
         break;
       case OutputInfo::AXMoment:
         getElementForces(*geomState, allCorot, iInfo, AXM, time);
         break;
       case OutputInfo::AYMoment:
         getElementForces(*geomState, allCorot, iInfo, AYM, time);
         break;
       case OutputInfo::AZMoment:
         getElementForces(*geomState, allCorot, iInfo, AZM, time);
         break;
       case OutputInfo::Energies:
         {
           // Since this file is used for both Non-linear Statics
           // and Non-linear Dynamics
           // we need both included in this routine.

           // Build displacement Vector
           Vector sol(numUncon(), 0.0 );
           int i;
           for(i=0; i<numNodes; ++i) {
             int xloc  = c_dsa->locate(i, DofSet::Xdisp);
             if(xloc >= 0)
               sol[xloc]  = ( (*geomState)[i].x - nodes[i]->x);
             int yloc  = c_dsa->locate(i, DofSet::Ydisp);
             if(yloc >= 0)
               sol[yloc]  = ( (*geomState)[i].y - nodes[i]->y);
             int zloc  = c_dsa->locate(i, DofSet::Zdisp);
             if(zloc >= 0)
               sol[zloc]  = ( (*geomState)[i].z - nodes[i]->z);
             double rot[3];
             mat_to_vec((*geomState)[i].R,rot);
             int xrot  = c_dsa->locate(i, DofSet::Xrot);
             if(xrot >= 0)
               sol[xrot]  = rot[0];
             int yrot  = c_dsa->locate(i, DofSet::Yrot);
             if(yrot >= 0)
               sol[yrot]  = rot[1];
             int zrot  = c_dsa->locate(i, DofSet::Zrot);
             if(zrot >= 0)
               sol[zrot]  = rot[2];
           }

           // Non-Linear Dynamics
           if (sinfo.probType == SolverInfo::NonLinDynam) {
             double dW=0.0, dWaero = 0.0, Wela=0.0, Wkin=0.0;

             if(time==sinfo.initialTime) {
               Wext=0.0;
               Waero=0.0;
               Wdmp=0.0;
               pWela=0.0;
               pWkin=0.0;
               previousExtForce = new Vector(force);
               previousAeroForce = new Vector(aeroForce);
               dW = force*sol;
               dWaero = aeroForce*sol;
             } else {
               double c = solInfo().newmarkGamma;
/*
     NOTE: The formula for dW should be (1-c)*f^n + c*f^{n+1}.
           However, force stores f^{n+1/2}
           and previousExtForce stores f^{n-1/2}. For this reason,
           the formula below looksk
           different. Furthermore, we had to assume
           f^n = (f^{n-1/2} + f^{n+1/2})/2 to minimize
           code changes. On the other hand, f^{n+1/2} = (f^n + f^{n+1})/2
           holds without any assumption.
*/

               dW = (c*force + (1.0-c)*(*previousExtForce)) * (sol - (*previousDisp));
               dWaero = (c*aeroForce + (1.0-c)*(*previousAeroForce)) *(sol - (*previousDisp));

               if(step==sinfo.initialTimeIndex) { dW*=2.0; dWaero *= 2.0; }
             }

             Wext += dW;
             Waero += dWaero;

             Vector vtmp(velocity,numUncon());
             Vector tmpVec(numUncon(),0.0);
             int iele, idof, jdof, dofn1, dofn2;

             // Compute Kinetic Energy
             // Compute vtmp^t M vtmp
             for(iele = 0; iele < numele; ++iele) {
               for(idof = 0; idof <  mel[iele].dim(); ++idof) {
                 dofn1 = c_dsa->getRCN((*allDOFs)[iele][idof]);
                 if(dofn1 >= 0) {
                   for(jdof = 0; jdof <  mel[iele].dim(); ++jdof) {
                     dofn2 = c_dsa->getRCN((*allDOFs)[iele][jdof]);
                     if(dofn2 >= 0)
                       tmpVec[dofn1] += vtmp[dofn2]*mel[iele][idof][jdof];
                   }
                 }
               }
             }
             Wkin = 0.5 * (vtmp * tmpVec);

             // Compute Internal Energy
             // This is done at the element level.
             double EleWela;
             for(iele = 0; iele < numele; ++iele) {
               EleWela = 0.0;
               EleWela = allCorot[iele]->getElementEnergy(*geomState,nodes);
               Wela += EleWela;
             }

             double error = (time==sinfo.initialTime) ? 0.0 :
                            (Wela+Wkin)-(pWela+pWkin)-dW;
             geoSource->outputEnergies(iInfo, time, Wext, Waero, Wela, Wkin, 0.0, error);

             pWela=Wela;
             pWkin=Wkin;

             if(time==sinfo.initialTime) {
               previousDisp = new Vector(sol);
             } else {
               (*previousExtForce) = force;
               (*previousAeroForce) = aeroForce;
               (*previousDisp)     = sol;
             }

           // Non-Linear Statics
           } else {
             double lambda = time;
             if (time == 0.0) {
               double deltaLambda = solInfo().getNLInfo().dlambda;
               double maxLambda = solInfo().getNLInfo().maxLambda;
               if(deltaLambda == maxLambda) lambda = 1.0;
             }
             Wext=0.0;
             Waero=0.0;
             Wdmp=0.0;
             pWela=0.0;
             pWkin=0.0;

             // Wkin = kinetic energy
             // Total Energy = Wext+Wela+Wkin
             double Wkin=0.0;

             // Wext = external energy
             Wext = lambda*force * sol;
             Waero = lambda*aeroForce * sol;

             // Compute Internal Energy
             // This is done at the element level.
             // Wela = elastic energy
             int iele;
             double Wela = 0.0;
             double EleWela;
             for(iele = 0; iele < numele; ++iele) {
               EleWela = 0.0;
               EleWela = allCorot[iele]->getElementEnergy(*geomState,nodes);
               Wela += EleWela;
             }

             double error = Wext+Wela+Wkin;
             geoSource->outputEnergies(iInfo,time,Wext, Waero, Wela,Wkin,0.0,error);
           }
         }
         break;
       default:
         fprintf(stderr," *** WARNING: Output case %d not implemented for non-linear direct solver \n", iInfo);
	 break;
    }

  //}

}

void
Domain::createCorotators(Corotator **allCorot)
{
 // Allocate memory for element stiffness matrix
 double *kmatrix = new double[maxNumDOFs*maxNumDOFs];

 // Loop over elements to get their Corotators
 for(int iele = 0; iele < numele; ++iele)
   allCorot[iele] = packedEset[iele]->getCorotator(nodes, kmatrix,
                                      sinfo.getNLInfo().fitAlgShell,
                                      sinfo.getNLInfo().fitAlgBeam);

 delete [] kmatrix;
}

void
Domain::getGeometricStiffness(GeomState &geomState,Vector& elementInternalForce,
                          Corotator **allCorot, FullSquareMatrix *&geomKelArray)
{

   // Get Geometric Stiffness

   int iele;
   for(iele = 0; iele < numele; ++iele) {

      // Get updated tangent stiffness matrix and element internal force
      allCorot[iele]->formGeometricStiffness(geomState, nodes,
                                             geomKelArray[iele],
                                             elementInternalForce.data());
      if(!solInfo().getNLInfo().unsymmetric) geomKelArray[iele].symmetrize();
   }

}

void
Domain::computeGeometricPreStress(Corotator **&allCorot, GeomState *&geomState,
                                  FullSquareMatrix *&kelArray, StaticTimers *times,
                                  FullSquareMatrix *&geomKelArray)
{
   // ... ALLOCATE MEMORY FOR THE ARRAY OF COROTATORS
   times->corotatorTime -= getTime();
   allCorot = new Corotator *[numElements()];

   // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
   createCorotators(allCorot);
   times->corotatorTime += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Create Element Corotators %19.5f s\n", times->corotatorTime/1000.0);
#endif

   // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
   times->kelArrayTime -= getTime();
   createKelArray(kelArray);
   times->kelArrayTime += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Create Element Stiffness Array %14.5f s\n", times->kelArrayTime/1000.0);
#endif

   times->timeGeom -= getTime();
   geomState = new GeomState( *getDSA(), *getCDSA(), getNodes());
   times->timeGeom += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Build GeomState %29.5f s\n", times->timeGeom/1000.0);
#endif

   // geometric prestress for linear: update geomState with initial displacements and compute element stiffness matrices
   if(!sinfo.isNonLin() && (numInitDisp6() > 0) && (sinfo.gepsFlg == 1)) {
     times->timePresc -= getTime();
     geomState->updatePrescribedDisplacement(getInitDisp6(), numInitDisp6(), 1.0);
     times->timePresc += getTime();
#ifdef PRINT_NLTIMERS
     fprintf(stderr," ... Update Geometry to initial state %12.5f s\n", times->timePresc/1000.0);
#endif
   }

   times->buildStiffAndForce -= getTime();
   Vector elementInternalForce(maxNumDOF(), 0.0);
   Vector residual(numUncon(), 0.0);
   getStiffAndForce(*geomState, elementInternalForce, allCorot,
                    kelArray, residual);
   times->buildStiffAndForce += getTime();
#ifdef PRINT_NLTIMERS
   fprintf(stderr," ... Build Element Tangent Stiffness %13.5f s\n",
           times->buildStiffAndForce/1000.0);
#endif

   // Buckling analysis if requested
   if(geomKelArray) {
     Vector elementInternalForce(maxNumDOF(), 0.0);
     getGeometricStiffness(*geomState, elementInternalForce, allCorot, geomKelArray);
   }
}


void
Domain::getStressStrain(GeomState &geomState, Corotator **allCorot,
                        int fileNumber, int stressIndex, double time)
{

  OutputInfo *oinfo = geoSource->getOutputInfo();

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS, WHERE
  // ... ONLY THE AXIAL STRAIN (EXX) AND AXIAL STRESS (SXX) ARE CALCULATED

  int avgnum  = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  double ylayer = oinfo[fileNumber].ylayer;
  double zlayer = oinfo[fileNumber].zlayer;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  int k;

  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;

  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;

  // ... WRITE CURRENT TIME VALUE
  if(oinfo[fileNumber].nodeNumber == -1)
    fprintf(oinfo[fileNumber].filptr,"  % *.*E\n",w,p,time);

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(stress == 0)
    stress = new Vector(numnodes,0.0);
  else
    stress->zero();

  if(weight == 0)
    weight = new Vector(numnodes,0.0);
  else
    weight->zero();

  if(elDisp == 0)
    elDisp = new Vector(maxNumDOFs,0.0);

  int iele;
  if((elstress == 0)||(elweight == 0)) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = myMax(maxNodesPerElement, NodesPerElement);
    }
    if(elstress == 0) elstress = new Vector(maxNodesPerElement, 0.0);
    if(elweight == 0) elweight = new Vector(maxNodesPerElement, 0.0);
  }


  int flag;
  for(iele = 0; iele < numElements(); ++iele) {

    elDisp->zero();
    elstress->zero();
    elweight->zero();

     // extract deformations from current Geometry State of structure
     allCorot[iele]->extractDeformations( geomState, nodes,
                                          elDisp->data(), flag);
//---------------------------------------------------------------------------

   int *nodeNumbers = new int[maxNumNodes];
   Vector elemNodeTemps(maxNumNodes);
   elemNodeTemps.zero();

   packedEset[iele]->nodes(nodeNumbers);

   int NodesPerElement = packedEset[iele]->numNodes();
   double *nodalTemperatures = 0;
   // Either get the nodal temperatures from the input file or
   // from the thermal model
   if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
   if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

   int iNode;
   if (sinfo.thermalLoadFlag || (sinfo.thermoeFlag>=0))
     for (iNode = 0; iNode < NodesPerElement; ++iNode) {
       if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
         elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
       else
         elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
     }

//----------------------------------------------------------------------------

     if (flag == 1) {
// USE LINEAR STRESS ROUTINE
// ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
     packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                   *elDisp, stressIndex, surface,
				   elemNodeTemps.data(),ylayer,zlayer,avgnum);

     } else if (flag == 2) {
// USE NON-LINEAR STRESS ROUTINE
     allCorot[iele]->getNLVonMises(*elstress, *elweight, geomState,
                                   nodes, stressIndex);

     } else {
// NO STRESS RECOVERY
     }

// ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

//     int NodesPerElement = packedEset[iele]->numNodes();

     for(k=0; k<NodesPerElement; ++k) {
       (*stress)[(*elemToNode)[iele][k]] += (*elstress)[k];
       (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
     }

// ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
     if(avgnum == 0) {
       for(k=0; k<NodesPerElement; ++k)
         fprintf(oinfo[fileNumber].filptr," % *.*E",w,p,(*elstress)[k]);
       fprintf(oinfo[fileNumber].filptr,"\n");
     }
  }


// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

   int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes
   if(avgnum == 1 || avgnum == 2) {

     if(oinfo[fileNumber].nodeNumber == -1) {
       for(k=0; k<numNodes; ++k) {
          if((*weight)[k] == 0.0)
            fprintf(oinfo[fileNumber].filptr," % *.*E\n",w,p,0.0);
          else
            fprintf(oinfo[fileNumber].filptr," % *.*E\n",w,p,
                    (*stress)[k]/=(*weight)[k]);
       }
     } else {
       if((*weight)[oinfo[fileNumber].nodeNumber] == 0.0)
         fprintf(oinfo[fileNumber].filptr," %*.*E % *.*E\n",w,p,time,w,p,0.0);
       else
         fprintf(oinfo[fileNumber].filptr," %*.*E % *.*E\n",w,p,time,w,p,(*stress)[oinfo[fileNumber].nodeNumber]/=(*weight)[oinfo[fileNumber].nodeNumber]);
     }
   }
   fflush(oinfo[fileNumber].filptr);
}

void
Domain::getPrincipalStress(GeomState &geomState, Corotator **allCorot,
                         int fileNumber, int stressIndex, double time)
{
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
    fprintf(stderr,"*** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // ... STRESSES ARE CALCULATED FOR EVERYTHING EXCEPT BARS & BEAMS

  int avgnum  = oinfo[fileNumber].averageFlg;
  int surface = oinfo[fileNumber].surface;

  // upper  surface = 1
  // median surface = 2
  // lower  surface = 3

  int j,k;

  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;
  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;

  // ... OUTPUT NODE NUMBER
  int n = oinfo[fileNumber].nodeNumber;

  // ... ALLOCATE VECTORS STRESS AND WEIGHT AND INITIALIZE TO ZERO
  if(p_stress == 0) p_stress = new FullM(numnodes,6);
  if(weight == 0) weight = new Vector(numnodes,0.0);
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

  // ... WRITE CURRENT TIME VALUE
  if(n == -1) {
    fprintf(oinfo[fileNumber].filptr,"  % *.*E\n",w,p,time);
  }

  int flag;

  for(iele = 0; iele < numElements(); ++iele) {

    elDisp->zero();
    p_elstress->zero();
    elweight->zero();

     // extract deformations from current Geometry State of structure
     allCorot[iele]->extractDeformations( geomState, nodes,
                                          elDisp->data(), flag);

     if (flag == 1) {
// USE LINEAR STRESS ROUTINE
// ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT
     packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
                                    *elDisp, strInd, surface);

     } else if (flag == 2) {
// USE NON-LINEAR STRESS ROUTINE
     allCorot[iele]->getNLAllStress(*p_elstress, *elweight, geomState,
                                    nodes, strInd);

     } else {
// NO STRESS RECOVERY
     }

// ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

    int NodesPerElement = elemToNode->num(iele);

    for(k=0; k<NodesPerElement; ++k) {
      for(j=0; j<6; ++j) {
        (*p_stress)[(*elemToNode)[iele][k]][j] += (*p_elstress)[k][j];
      }
      (*weight)[(*elemToNode)[iele][k]] += (*elweight)[k];
    }

// ... PRINT NON-AVERAGED STRESS VALUES IF REQUESTED
//     THIS WRITES THE CHOSEN PRINCIPAL STRESS FOR EACH ELEMENT
    if(avgnum == 0) {
      fprintf(oinfo[fileNumber].filptr," % *.*E\n",
        w,p,(*p_elstress)[0][5+strDir]);
    }
  }

// ... AVERAGE STRESS/STRAIN VALUE AT EACH NODE BY THE NUMBER OF
// ... ELEMENTS ATTACHED TO EACH NODE IF REQUESTED.

  if(avgnum == 1 || avgnum == 2) {

    if(n == -1) {
      for(k=0; k<numnodes; ++k) {
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

// ... CALCULATE PRINCIPALS AT EACH NODE

    double svec[6], pvec[3];
    if(n == -1) {
      for(k=0; k<numnodes; ++k) {
        for (j=0; j<6; ++j) {
          svec[j] = (*p_stress)[k][j];
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        fprintf(oinfo[fileNumber].filptr," % *.*E\n",w,p,pvec[strDir-1]);
      }
    }
    else {
      for (j=0; j<6; ++j) {
        svec[j] = (*p_stress)[n][j];
      }
// Convert Engineering to Tensor Strains
      if(strInd != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      fprintf(oinfo[fileNumber].filptr," % *.*E\n",w,p,pvec[strDir-1]);
    }

    fflush(oinfo[fileNumber].filptr);
  }
  else {
    fflush(oinfo[fileNumber].filptr);
  }
}


// Nonlinear version of getElementForces

void
Domain::getElementForces( GeomState &geomState, Corotator **allCorot,
                          int fileNumber, int forceIndex, double time)
{

  // allocate integer array to store node numbers
  // make sure that maxNumNodes is computed
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes);
  elemNodeTemps.zero();
  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag)  nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  // ... FORCES ARE CALCULATED FOR BARS AND BEAMS CURRENTLY

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

  int flag;

  for(iele=0; iele<numele; ++iele) {

  packedEset[iele]->nodes(nodeNumbers);

// ... DETERMINE ELEMENT DISPLACEMENT VECTOR

    allCorot[iele]->extractDeformations(geomState, nodes,
                                        elDisp->data(), flag);

// ... CALCULATE INTERNAL FORCE VALUE FOR EACH ELEMENT
     if(sinfo.thermalLoadFlag || (sinfo.thermoeFlag >=0)) {
       int iNode;
       for(iNode=0; iNode<NodesPerElement; ++iNode) {
        if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
     }

    packedEset[iele]->getIntrnForce(*elstress,nodes,elDisp->data(),forceIndex,elemNodeTemps.data());

// ... COPY ELEMENT'S NODAL FORCES INTO A TOTAL FORCE MATRIX

    for(i=0; i<NodesPerElement; ++i)
      forces[iele][i] = (*elstress)[i];

  }

// ... PRINT THE ELEMENT FORCES TO A FILE
   geoSource->outputElemVectors(fileNumber, forces.data(), numele, time);
/*
   fprintf(oinfo[fileNumber].filptr,"%*.*e\n",w,p,time);

   int k;
   for(k=0; k<numele; ++k)
     fprintf(oinfo[fileNumber].filptr,"% *.*e\t% *.*e\n",
             w,p,forces[k][0],w,p,forces[k][1]);

   fflush(oinfo[fileNumber].filptr);
*/

}



// Nonlinear restart file

void
Domain::writeRestartFile(double time, int timeIndex, Vector &v_n,
                         GeomState *geomState)
{
// either test for pointer or frequency > 0

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if((timeIndex % sinfo.nRestart == 0) || (time >= sinfo.tmax-0.1*sinfo.getTimeStep())) {
   int fn = open(cinfo->currentRestartFile, O_WRONLY | O_CREAT, 0666);
   if(fn >= 0) {
     cerr << "here in NLStatic.C writeRestartFile\n";
     int writeSize;
     writeSize = write(fn, &timeIndex, sizeof(int));
     if(writeSize != sizeof(int))
       fprintf(stderr," *** ERROR: Writing restart file time index\n");

     writeSize = write(fn, &time, sizeof(double));
     if(writeSize != sizeof(double))
       fprintf(stderr," *** ERROR: Writing restart file time\n");

      writeSize = write(fn, v_n.data(), v_n.size()*sizeof(double));
      if(int(writeSize) != int(v_n.size()*sizeof(double)))
        fprintf(stderr," *** ERROR: Writing restart file velocity\n");

     double *positions = (double *) dbg_alloca(sizeof(double)*3*numnodes);
     geomState->getPositions(positions);
     writeSize = write(fn, positions,numnodes*3*sizeof(double));
     if(int(writeSize) != int(numnodes*3*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file geometry_state\n");

     double *rotations = (double *) dbg_alloca(sizeof(double)*9*numnodes);
     geomState->getRotations(rotations);
     writeSize = write(fn, rotations,numnodes*9*sizeof(double));
     if(int(writeSize) != int(numnodes*9*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file geometry_state\n");

     // PJSA 9-17-2010
     int numEle = packedEset.last();
     for(int i = 0; i < numEle; ++i)
       packedEset[i]->writeHistory(fn);

     close(fn);
   } else {
      perror(" *** ERROR: Restart file could not be opened: ");
      exit(-1);
   }
 }
}

void
Domain::readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                        Vector &v_p, double *bcx, double *vcx,
                        GeomState &geomState)
{
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(cinfo->lastRestartFile) {
   fprintf(stderr, " ... Restart From a Prev. Nonl. Run ...\n");
   int fn = open(cinfo->lastRestartFile,O_RDONLY );
   if(fn >= 0) {
     int restartTIndex;
     double restartT;
     int readSize = read(fn, &restartTIndex, sizeof(int));
     if(readSize != sizeof(int))
       fprintf(stderr," *** ERROR: Inconsistent restart file 1\n");
     sinfo.initialTimeIndex = restartTIndex;

     readSize = read(fn, &restartT, sizeof(double));
     if(readSize != sizeof(double))
       fprintf(stderr," *** ERROR: Inconsistent restart file 2\n");
     fprintf(stderr,"Initial Time = %f\n\n",restartT);
     sinfo.initialTime = restartT;

     v_n.zero();
     readSize = read(fn, v_n.data(), sizeof(double)*v_n.size());
     if(int(readSize) != int(sizeof(double)*v_n.size()))
       fprintf(stderr," *** ERROR: Inconsistent restart file 2.5\n");

     double *positions = (double *)dbg_alloca(sizeof(double)*3*numnodes);
     readSize = read(fn, positions, numnodes*3*sizeof(double));
     if(int(readSize) != int(numnodes*3*sizeof(double)))
       fprintf(stderr," *** ERROR: Inconsistent restart file 3\n");
     geomState.setPositions(positions);

     double *rotations = (double *)dbg_alloca(sizeof(double)*9*numnodes);
     readSize = read(fn, rotations, numnodes*9*sizeof(double));
     if(int(readSize) != int(numnodes*9*sizeof(double)))
       fprintf(stderr," *** ERROR: Inconsistent restart file 4\n");
     geomState.setRotations(rotations);

     // PJSA 9-17-2010
     int numEle = packedEset.last();
     for(int i = 0; i < numEle; ++i) 
       packedEset[i]->readHistory(fn);

     close(fn);

     d_n.zero();
     a_n.zero();
     v_p.zero();
     for(int i = 0; i < numNodes(); ++i) {

       int xloc  = c_dsa->locate(i, DofSet::Xdisp );
       int xloc1 =   dsa->locate(i, DofSet::Xdisp );

       if(xloc >= 0)
         d_n[xloc]  = ( (geomState)[i].x - nodes[i]->x);
       else if (xloc1 >= 0)
         bcx[xloc1] = ( (geomState)[i].x - nodes[i]->x);

       int yloc  = c_dsa->locate(i, DofSet::Ydisp );
       int yloc1 =   dsa->locate(i, DofSet::Ydisp );

       if(yloc >= 0)
         d_n[yloc]  = ( (geomState)[i].y - nodes[i]->y);
       else if (yloc1 >= 0)
         bcx[yloc1] = ( (geomState)[i].y - nodes[i]->y);

       int zloc  = c_dsa->locate(i, DofSet::Zdisp);
       int zloc1 =   dsa->locate(i, DofSet::Zdisp);

       if(zloc >= 0)
         d_n[zloc]  = ( (geomState)[i].z - nodes[i]->z);
       else if (zloc1 >= 0)
         bcx[zloc1] = ( (geomState)[i].z - nodes[i]->z);
     }

     if(solInfo().aeroFlag >= 0)
       aeroPreProcess( d_n, v_n, a_n, v_p, bcx, vcx );

   } else {
      perror(" *** ERROR: Restart file could not be opened: ");
      //exit(-1);
   }

 }
 //else {
 //    fprintf(stderr, " ... No restart                     ...\n");
 //}


}
