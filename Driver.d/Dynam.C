#include <Utils.d/dbg_alloca.h>
#include <stdio.h>
//#include <stdlib.h>

#ifndef TFLOP
#ifndef WINDOWS 
#include <dlfcn.h>
#endif
#endif

// New include files for Restart file
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <Math.d/DBSparseMatrix.h>
#include <Driver.d/Dynam.h>
#include <Driver.d/Domain.h>
#include <Hetero.d/FlExchange.h>
#include <Utils.d/ModeData.h>
#include <Element.d/State.h>
#include <Driver.d/GeoSource.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>

extern ModeData modeData; 

void
Domain::initDispVeloc(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p)
{
 // ... INITIALIZE DISPLACEMENT, VELOCITY, AND ACCELERATION 
 // ... VECTORS TO ZERO
 d_n = v_n = a_n = v_p = 0.0;

 // ... SET INITIAL VELOCITY
 int i;
 for(i = 0; i < numIVel; ++i) {
   int dof = c_dsa->locate(iVel[i].nnum, 1 << iVel[i].dofnum);
   if(dof >= 0)
     v_n[dof] = iVel[i].val;
 }

 // If zeroInitialDisp is set, then return as the displacement is
 // already set to zero above.
 // note: geps is always zero for nonlinear
 if(sinfo.zeroInitialDisp == 0) {
   // ... SET INITIAL DISPLACEMENT FROM IDISP IF IDISP6 DOES NOT EXIST
   // ... OR IF WE ARE USING GEOMETRIC PRE-STRESS (GEPS)
   if(domain->numInitDisp6() == 0 || sinfo.gepsFlg == 1 ) { // note: always use global num to do this check
     if(domain->solInfo().modalIDisp) { //HB: BCond array iDis contains y0
       filePrint(stderr, " ... Compute initial displacement from given modal basis (u0=X.y0)\n"); //HB
       modeData.addMultY(numIDis, iDis, d_n, c_dsa);
     } 
     else { //HB: BCond array iDis contains directly u0
       //if(domain->numInitDisp()) cerr << "In Domain::initDispVeloc: adding IDISP to d_n\n"; //HB
       for(i = 0; i < numIDis; ++i) {
         int dof = c_dsa->locate(iDis[i].nnum, 1 << iDis[i].dofnum);
         if(dof >= 0) d_n[dof] = iDis[i].val;
       } 
     }
   }

   // ... SET INITIAL DISPLACEMENT FROM IDISP6
   if(domain->numInitDisp() == 0 || sinfo.gepsFlg == 0) {
     //if(domain->numInitDisp6()) cerr << "In Domain::initDispVeloc: add IDISP6 to d_n\n"; //HB
     for(i = 0; i < numIDis6; ++i) {
       int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
       if(dof >= 0)
         d_n[dof] = iDis6[i].val;
     }   
   }
 }

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(probType() != SolverInfo::NonLinDynam) {
   if(cinfo->lastRestartFile) {
     fprintf(stderr, " ... Restarting From a Previous Run ...\n");
     int fn = open(cinfo->lastRestartFile,O_RDONLY );
     if(fn >= 0) {
       int vsize, restartTIndex;
       double restartT;
       int readSize = read(fn, &vsize, sizeof(int));
       if(vsize != d_n.size() || readSize != (int) sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       readSize = read(fn, &restartTIndex, sizeof(int));
       if(readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0) 
         sinfo.initialTimeIndex = 0;
       else
         sinfo.initialTimeIndex = restartTIndex;
       readSize = read(fn, &restartT, sizeof(double));
       if(readSize != sizeof(double))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTime = 0.0;
       else
         sinfo.initialTime = restartT;
       readSize = read(fn, d_n.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
     
       readSize = read(fn, v_n.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, v_p.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double))) {
         fprintf(stderr," *** WARNING: Older version of restart file"
                        " -- Missing velocity field is set to zero\n");
         v_p.zero();
       }
       readSize = read(fn, &(sinfo.initExtForceNorm), sizeof(double));
       if(readSize != sizeof(double))
         fprintf(stderr," *** ERROR: Inconsistent restart file: External Force Norm\n");

       close(fn);
     } else {
        perror(" *** ERROR: Restart file could not be opened: ");
        //exit(-1);
     }
   
   } 
   //else {
   //    fprintf(stderr, " ... No restart                     ...\n");
   //}
 }
}

void
Domain::writeRestartFile(double time, int timeIndex, Vector &d_n, 
                         Vector &v_n, Vector &v_p, double Fref)
{
// either test for pointer or frequency > 0
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(timeIndex % sinfo.nRestart == 0 || time >= sinfo.tmax-0.1*sinfo.dt){
   int fn = open(cinfo->currentRestartFile, O_WRONLY | O_CREAT, 0666);
   if(fn >= 0) {

     int vsize = d_n.size();
     int writeSize = write(fn, &vsize, sizeof(int));
     if(writeSize != sizeof(int))
       fprintf(stderr," *** ERROR: Writing restart file vector size\n");

     writeSize = write(fn, &timeIndex, sizeof(int));
     if(writeSize != sizeof(int))
       fprintf(stderr," *** ERROR: Writing restart file time index\n");

     writeSize = write(fn, &time, sizeof(double));
     if(writeSize != sizeof(double))
       fprintf(stderr," *** ERROR: Writing restart file time\n");

     writeSize = write(fn, d_n.data(), d_n.size()*sizeof(double));
     if(int(writeSize) != int(d_n.size()*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file displacement \n");

     writeSize = write(fn, v_n.data(), v_n.size()*sizeof(double));
     if(int(writeSize) != int(v_n.size()*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file velocity\n");

     writeSize = write(fn, v_p.data(), v_p.size()*sizeof(double));
     if(int(writeSize) != int(v_p.size()*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file prev. velocity\n");

     writeSize = write(fn, &Fref, sizeof(double));
     if(writeSize != sizeof(double))
       fprintf(stderr," *** ERROR: Writing restart file external force norm\n");
     close(fn);
   } else {
      perror(" *** ERROR: Restart file could not be opened: ");
      exit(-1);
   }
 }
}
//---------------------------------------------------------------------------------------------

void
Domain::getOrAddDofForPrint(bool ad, Vector& d_n, double* bcx, int iNode, double *xdata, 
                        int *dofx, double *ydata, int *dofy, double *zdata, int *dofz)
{
    // ad==true -> add :   d_n[xloc]+=data  
    // ad==false -> get:  data=d_n[xloc]
    
    //beginning of dynamOutput with: xdata=vx=variation, dn=dn_aero, bcx=bcx_aero 

    // c_dsa mean: dofsetarray without constrained dof ie free dof
    //   dsa mean: constrained dofsetarray

    // xloc>=0  : dof free
    // xloc1>=0 : dof exist but not free
    // xloc<0 && xloc1<0 : dof doesn't exist
  
    if (dofx){                                                                                              
      int xloc  = c_dsa->locate( iNode, *dofx);
      int xloc1 =   dsa->locate( iNode, *dofx);
      if(xloc >= 0)       { if(!ad) *xdata=d_n[xloc];  else d_n[xloc]+=*xdata; }
      else if(xloc1 >= 0 && bcx) { if(!ad) *xdata=bcx[xloc1]; else bcx[xloc1]+=*xdata;}
      else                { if(!ad) *xdata=0.0;}
    }
    
    if(dofy){     
     int yloc  = c_dsa->locate( iNode, *dofy);
     int yloc1 =   dsa->locate( iNode, *dofy);
     if(yloc >= 0)       { if(!ad) *ydata=d_n[yloc];  else d_n[yloc]+=*ydata; }
     else if(yloc1 >= 0 && bcx) { if(!ad) *ydata=bcx[yloc1]; else bcx[yloc1]+=*ydata;}
     else                { if(!ad) *ydata=0.0;}
    }
   
    if(dofz){                                                                                         
      int zloc  = c_dsa->locate( iNode, *dofz);
      int zloc1 =   dsa->locate( iNode, *dofz);
      if(zloc >= 0)       { if(!ad) *zdata=d_n[zloc];  else d_n[zloc]+=*zdata; }
      else if(zloc1 >= 0 && bcx) { if(!ad) *zdata=bcx[zloc1]; else bcx[zloc1]+=*zdata;}
      else                { if(!ad) *zdata=0.0;}
    }
}
//-----------------------------------------------------------------------------------------------

void
Domain::addVariationOfShape_StructOpt(int iNode, CoordSet *nodescopy, double &x, double &y, double &z)
{
  // ... considering the variation of shape
  Node oldnod;      
  oldnod=nodescopy->getNode(iNode);
        
  double vx = nodes[iNode]->x - oldnod.x ;
  double vy = nodes[iNode]->y - oldnod.y ;
  double vz = nodes[iNode]->z - oldnod.z ;
       
  x=x+vx;   y=y+vy;    z=z+vz;
}
//------------------------------------------------------------------------------------------------

void
Domain::aeroSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double *bcx, double* vcx)
{
  // Send u + IDISP6 to fluid code.
  // IDISP6 is used to compute pre-stress effects.
  if(sinfo.aeroFlag < 0) return;

  getTimers().sendFluidTime -= getTime();
  Vector d_n_aero(d_n);
  double * bcx_aero = bcx;

  if(sinfo.gepsFlg == 1) {

    for(int i = 0; i < numIDis6; ++i) {
      int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
      if(dof >= 0)
        d_n_aero[dof] += iDis6[i].val;
    }
  }

  State state(c_dsa, dsa, bcx_aero, vcx, d_n_aero, v_n, a_n, v_p);

  flExchanger->sendDisplacements(nodes, state);
  if(verboseFlag) fprintf(stderr, " ... [E] Sent displacements ...\n");

  getTimers().sendFluidTime += getTime();
}

void
Domain::dynamOutput(int tIndex, double *bcx, DynamMat& dMat, Vector& ext_f, Vector &aeroForce, 
                    Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* vcx) {
  // Current time stamp
  double time = tIndex*sinfo.dt;
   
  if (sinfo.nRestart > 0 && !sinfo.modal) {
    writeRestartFile(time, tIndex, d_n, v_n, v_p, sinfo.initExtForceNorm);
  }
  
  // Send to fluid code (except C0)
  if(sinfo.aeroFlag >= 0 && !sinfo.lastIt && tIndex != sinfo.initialTimeIndex && sinfo.aeroFlag != 20)
    aeroSend(d_n, v_n, a_n, v_p, bcx, vcx);

  int numOutInfo = geoSource->getNumOutInfo();

  // Open the file and write the header in the first time step
  if (tIndex == sinfo.initialTimeIndex) {
    if(numOutInfo > 0)
      //fprintf(stderr," ... Outputting Dynam./Quasi. Sol.  ...\n");
      geoSource->openOutputFiles();
  }

  this->dynamOutputImpl(tIndex, bcx, dMat, ext_f, aeroForce, d_n, v_n, a_n, v_p, vcx, time, 0, numOutInfo);
} 

// Output for time-parallel linear dynamics
// There is a full set of files for each time-slice on the CPU
// Preconditions:
// 1) The function GeoSource::duplicateFilesForPita must have been called during the initialization
// 2) The corresponding output files must already be open
void
Domain::pitaDynamOutput(int tIndex, double *bcx, DynamMat& dMat, Vector& ext_f, Vector &aeroForce, 
                    Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* vcx, 
                    int sliceRank, double time)
{
  std::pair<int, int> requestIndices = geoSource->getTimeSliceOutputFileIndices(sliceRank);
  this->dynamOutputImpl(tIndex, bcx, dMat, ext_f, aeroForce, d_n, v_n, a_n, v_p, vcx, time, requestIndices.first, requestIndices.second);
} 

// Perform output for files with indices in [firstRequest, lastRequest)
// Precondition: The corresponding output files must already be open
void
Domain::dynamOutputImpl(int tIndex, double *bcx, DynamMat& dMat, Vector& ext_f, Vector &aeroForce,
                         Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* vcx,
                         double time, int firstRequest, int lastRequest)
{
  // Print out the displacement info
  // XXXX numnodes is no longer reliable since it comes from the reduced node to node for direct mpcs
  double (*glDisp)[11] = new double[numnodes][11];//DofSet::max_known_nonL_dof
  for (int i = 0; i < numnodes; ++i)
    for (int j = 0 ; j < 11 ; j++)
      glDisp[i][j] = 0.0;
  /*int exactNumNodes =*/ mergeDistributedDisp(glDisp, d_n.data(), bcx);

  for (int i = firstRequest; i < lastRequest; ++i) {
    enum {YOUNG,MDENS,THICK};
    int iNode;
    int first_node, last_node;
    
    int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes
    OutputInfo *oinfo = geoSource->getOutputInfo();
    
    //CD: ad and get will be used in  addVariationOfShape_StructOpt and getOrAddDofForPrint which were 
    //    added to "clean" dynamOutput 
    bool get=false;
     
    if (oinfo[i].nodeNumber == -1){first_node=0; last_node=numNodes;}           
    else { first_node=oinfo[i].nodeNumber; last_node=first_node+1;}

    if((oinfo[i].interval != 0) && (tIndex % oinfo[i].interval == 0)) {
      int dof=-1;
      int w = oinfo[i].width;
      int p = oinfo[i].precision;

      switch(oinfo[i].type) {

        case OutputInfo::Displacement:

          if (oinfo[i].nodeNumber == -1) {                             // all nodes
            geoSource->outputNodeVectors(i, glDisp, numNodes, time);
          } else {                                                      // one node
            iNode = oinfo[i].nodeNumber;
            
            double x,y,z;
            getOrAddDofForPrint(get, d_n, bcx, iNode, &x, &DofSet::Xdisp, &y, &DofSet::Ydisp, &z, &DofSet::Zdisp);
            
            fprintf(oinfo[i].filptr,"  % *.*E   % *.*E % *.*E % *.*E\n",
                    w, p, time, w, p, x, w, p, y, w, p, z);
          }
          fflush(oinfo[i].filptr);
          break;
         
        case OutputInfo::Helmholtz: {

          if (oinfo[i].nodeNumber == -1) 
            fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          double x; 
          for (iNode=first_node; iNode<last_node; iNode++){
            x=0.0;
            getOrAddDofForPrint(get, d_n, bcx, iNode, &x, &DofSet::Helm, 0, 0, 0, 0);                    
            if (oinfo[i].nodeNumber == -1) 
              fprintf(oinfo[i].filptr,"  % *.*E\n",w, p, x);
            else 
              fprintf(oinfo[i].filptr,"  % *.*E   % *.*E\n",w, p, time, w, p, x);
          }
          fflush(oinfo[i].filptr);
          } break;

        case OutputInfo::Temperature: {

          if (oinfo[i].nodeNumber == -1)
            fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          double x;
          for (iNode=first_node; iNode<last_node; iNode++){
            x=0.0;
            getOrAddDofForPrint(get, d_n, bcx, iNode, &x, &DofSet::Temp, 0, 0, 0, 0);
            if (oinfo[i].nodeNumber == -1)
              fprintf(oinfo[i].filptr,"  % *.*E\n",w, p, x);
            else
              fprintf(oinfo[i].filptr,"  % *.*E   % *.*E\n",w, p, time, w, p, x);
          }
          fflush(oinfo[i].filptr);
          } break;
 
        case OutputInfo::Disp6DOF:
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          for (iNode=first_node; iNode<last_node; iNode++) {                   // all node or one node

            double x,y,z,xr,yr,zr;
            getOrAddDofForPrint(get, d_n, bcx, iNode, &x, &DofSet::Xdisp, &y, &DofSet::Ydisp, &z, &DofSet::Zdisp);
            getOrAddDofForPrint(get, d_n, bcx, iNode, &xr, &DofSet::Xrot, &yr, &DofSet::Yrot, &zr, &DofSet::Zrot);

            if (oinfo[i].nodeNumber == -1) {
              fprintf(oinfo[i].filptr,"%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                                       iNode+1,w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
            } else {
              fprintf(oinfo[i].filptr, "  % *.*E   %d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                                        w,p,time,iNode+1,w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
            }
            fflush(oinfo[i].filptr);
          }
          break;
                 
                
        case OutputInfo::Velocity6:               // all nodes or only one
        case OutputInfo::Velocity:

          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          else fprintf(oinfo[i].filptr,"  % *.*E",w,p,time);

          for (iNode=first_node; iNode<last_node; ++iNode)  {
               
            double x,y,z;
            double xr,yr,zr;
            getOrAddDofForPrint(get, v_n, vcx, iNode, &x, &DofSet::Xdisp, &y, &DofSet::Ydisp, &z, &DofSet::Zdisp);
              
            if (oinfo[i].type == OutputInfo::Velocity6)
              getOrAddDofForPrint(get, v_n, vcx, iNode, &xr, &DofSet::Xrot, &yr, &DofSet::Yrot, &zr, &DofSet::Zrot);
             
            if (oinfo[i].type == OutputInfo::Velocity6)  {
              fprintf(oinfo[i].filptr, "%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                                           iNode+1,w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
            } else
              fprintf(oinfo[i].filptr," % *.*E % *.*E % *.*E\n",w,p,x,w,p,y,w,p,z);
          }
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::PressureFirstTimeDerivative: {

          if (oinfo[i].nodeNumber == -1)
            fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          double x;
          for (iNode=first_node; iNode<last_node; iNode++){
            x=0.0;
            getOrAddDofForPrint(get, v_n, vcx, iNode, &x, &DofSet::Helm, 0, 0, 0, 0);
            if (oinfo[i].nodeNumber == -1)
              fprintf(oinfo[i].filptr,"  % *.*E\n",w, p, x);
            else
              fprintf(oinfo[i].filptr,"  % *.*E   % *.*E\n",w, p, time, w, p, x);
          }
          fflush(oinfo[i].filptr);
          } break;

        case OutputInfo::TemperatureFirstTimeDerivative: {

          if (oinfo[i].nodeNumber == -1)
            fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          double x;
          for (iNode=first_node; iNode<last_node; iNode++){
            x=0.0;
            getOrAddDofForPrint(get, v_n, vcx, iNode, &x, &DofSet::Temp, 0, 0, 0, 0);
            if (oinfo[i].nodeNumber == -1)
              fprintf(oinfo[i].filptr,"  % *.*E\n",w, p, x);
            else
              fprintf(oinfo[i].filptr,"  % *.*E   % *.*E\n",w, p, time, w, p, x);
          }
          fflush(oinfo[i].filptr);
          } break;
            
        case OutputInfo::Accel6:
        case OutputInfo::Acceleration:
          // XXXX acx (prescribed accelerations not implemented)
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          else fprintf(oinfo[i].filptr,"  % *.*E",w,p,time);

          for (iNode=first_node; iNode<last_node; ++iNode) {

            double x,y,z;
            double xr,yr,zr;
            getOrAddDofForPrint(get, a_n, (double *) 0 /*acx*/, iNode, &x, &DofSet::Xdisp, &y, &DofSet::Ydisp, &z, &DofSet::Zdisp);

            if (oinfo[i].type == OutputInfo::Accel6)
              getOrAddDofForPrint(get, a_n, (double *) 0 /*acx*/, iNode, &xr, &DofSet::Xrot, &yr, &DofSet::Yrot, &zr, &DofSet::Zrot);

            if (oinfo[i].type == OutputInfo::Accel6)  {
              fprintf(oinfo[i].filptr, "%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                                           iNode+1,w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
            } else
              fprintf(oinfo[i].filptr," % *.*E % *.*E % *.*E\n",w,p,x,w,p,y,w,p,z);
          }
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::PressureSecondTimeDerivative: {

          if (oinfo[i].nodeNumber == -1)
            fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          double x;
          for (iNode=first_node; iNode<last_node; iNode++){
            x=0.0;
            getOrAddDofForPrint(get, a_n, (double *) 0, iNode, &x, &DofSet::Helm, 0, 0, 0, 0);
            if (oinfo[i].nodeNumber == -1)
              fprintf(oinfo[i].filptr,"  % *.*E\n",w, p, x);
            else
              fprintf(oinfo[i].filptr,"  % *.*E   % *.*E\n",w, p, time, w, p, x);
          }
          fflush(oinfo[i].filptr);
          } break;

        case OutputInfo::StressXX:
          getStressStrain(d_n,bcx,i,SXX, time);
          break;
        case OutputInfo::StressYY:
          getStressStrain(d_n,bcx,i,SYY, time);
          break;
        case OutputInfo::StressZZ:
          getStressStrain(d_n,bcx,i,SZZ, time);
          break;
        case OutputInfo::StressXY:
          getStressStrain(d_n,bcx,i,SXY, time);
          break;
        case OutputInfo::StressYZ:
          getStressStrain(d_n,bcx,i,SYZ, time);
          break;
        case OutputInfo::StressXZ:
          getStressStrain(d_n,bcx,i,SXZ, time);
          break;
        case OutputInfo::StrainXX:
          getStressStrain(d_n,bcx,i,EXX, time);
          break;
        case OutputInfo::StrainYY:
          getStressStrain(d_n,bcx,i,EYY, time);
          break;
        case OutputInfo::StrainZZ:
          getStressStrain(d_n,bcx,i,EZZ, time);
          break;
        case OutputInfo::StrainXY:
          getStressStrain(d_n,bcx,i,EXY, time);
          break;
        case OutputInfo::StrainYZ:
          getStressStrain(d_n,bcx,i,EYZ, time);
          break;
        case OutputInfo::StrainXZ:
          getStressStrain(d_n,bcx,i,EXZ, time);
          break;
        case OutputInfo::StressVM:
          getStressStrain(d_n,bcx,i,VON, time);
          break;
        case OutputInfo::Damage:
          getStressStrain(d_n,bcx,i,DAMAGE, time);
          break;
        case OutputInfo::StressPR1:
          getPrincipalStress(d_n,bcx,i,PSTRESS1,time);
          break;
        case OutputInfo::StressPR2:
          getPrincipalStress(d_n,bcx,i,PSTRESS2,time);
          break;
        case OutputInfo::StressPR3:
          getPrincipalStress(d_n,bcx,i,PSTRESS3,time);
          break;
        case OutputInfo::StrainPR1:
          getPrincipalStress(d_n,bcx,i,PSTRAIN1,time);
          break;
        case OutputInfo::StrainPR2:
          getPrincipalStress(d_n,bcx,i,PSTRAIN2,time);
          break;
        case OutputInfo::StrainPR3:
          getPrincipalStress(d_n,bcx,i,PSTRAIN3,time);
          break;
        case OutputInfo::InXForce:
          getElementForces(d_n, bcx, i, INX, time);
          break;
        case OutputInfo::InYForce:
          getElementForces(d_n, bcx, i, INY, time);
          break;
        case OutputInfo::InZForce:
          getElementForces(d_n, bcx, i, INZ, time);
          break;
        case OutputInfo::AXMoment:
          getElementForces(d_n, bcx, i, AXM, time);
          break;
        case OutputInfo::AYMoment:
          getElementForces(d_n, bcx, i, AYM, time);
          break;
        case OutputInfo::AZMoment:
          getElementForces(d_n, bcx, i, AZM, time);
          break;
        case OutputInfo::Energies: {

          double dW=0.0, dWaero = 0.0, dWdmp=0.0, Wela=0.0, Wkin=0.0;
          int numUncon = c_dsa->size();

          if (time==sinfo.initialTime) {
            Wext=0.0;
            Waero=0.0;
            Wdmp=0.0;
            pWela=0.0;
            pWkin=0.0;
            Vector F(numUncon);
            buildRHSForce(F);
            previousExtForce = new Vector(F);
            if(sinfo.aeroFlag >= 0) {
              previousAeroForce = new Vector(aeroForce);
              dWaero = aeroForce*d_n;
            }
            dW = F*d_n;
          } else {
            double c = solInfo().newmarkGamma;

//NOTE:   The formula for dW should be (1-c)*f^n + c*f^{n+1}.
//        However, ext_f stores f^{n+1/2}
//        and previousExtForce stores f^{n-1/2}. For this reason,
//        the formula below looks
//        different. Furthermore, we had to
//        assume f^n = (f^{n-1/2} + f^{n+1/2})/2 to minimize
//        code changes. On the other hand, f^{n+1/2} = (f^n + f^{n+1})/2
//        holds without any assumption.
//
            dW = (c*ext_f + (1.0-c)*(*previousExtForce)) * (d_n - (*previousDisp));
            if(sinfo.aeroFlag >= 0) {
              dWaero = (c*aeroForce + (1.0-c)*(*previousAeroForce)) *(d_n - (*previousDisp));
            }
            if (tIndex==sinfo.initialTimeIndex) { dW *= 2.0; dWaero *= 2.0; }
          }
          Wext += dW;
          Waero += dWaero;
          Vector tmpVec(numUncon);
          dMat.M->mult(v_n,tmpVec);
          Wkin = 0.5 * (v_n * tmpVec);
          dMat.K->mult(d_n,tmpVec);
          Wela = 0.5 * (d_n * tmpVec);
          if (dMat.C && (time > sinfo.initialTime)) {                                //??????????????????????
            dMat.C->mult(v_n, tmpVec);
            double c = solInfo().newmarkGamma;
            dWdmp = (c*tmpVec + (1.0-c)*(*previousCq))*(d_n - (*previousDisp));
            Wdmp += dWdmp;
          }
          double error = (time==sinfo.initialTime) ? 0.0 : (Wela+Wkin)-(pWela+pWkin)+dWdmp-dW; //??????????????????????
          geoSource->outputEnergies(i, time, Wext, Waero, Wela, Wkin, -Wdmp, error);
          pWela=Wela;
          pWkin=Wkin;
          if (time==sinfo.initialTime) {                                             //?????????????????????
            previousDisp = new Vector(d_n);
            previousCq   = new Vector(tmpVec);
          } else {
            (*previousExtForce) = ext_f;
            if(sinfo.aeroFlag >= 0) (*previousAeroForce) = aeroForce;
            (*previousDisp)     = d_n;
            (*previousCq)       = tmpVec;
          }
        }
        break;

        case OutputInfo::StrainVM:
          getStressStrain(d_n,bcx,i,STRAINVON,time);
          break;
        case OutputInfo::AeroXForce:
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode=first_node; iNode<last_node; ++iNode) {
            int xloc  = c_dsa->locate( iNode, DofSet::Xdisp);
            double fx  = (xloc >= 0) ? ext_f[xloc] : 0.0;
            fprintf(oinfo[i].filptr," % *.*E\n",w,p,fx);
          }
          fflush(oinfo[i].filptr);
          break;
        case OutputInfo::AeroYForce:
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode=first_node; iNode<last_node; ++iNode) {
            int yloc  = c_dsa->locate( iNode, DofSet::Ydisp);
            double fy  = (yloc >= 0) ? ext_f[yloc] : 0.0;
            fprintf(oinfo[i].filptr," % *.*E\n",w,p,fy);
          }
          fflush(oinfo[i].filptr);
          break;
        case OutputInfo::AeroZForce:
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode=first_node; iNode<last_node; ++iNode) {
            int zloc  = c_dsa->locate( iNode, DofSet::Zdisp);
            double fz = (zloc >= 0) ? ext_f[zloc] : 0.0;
            fprintf(oinfo[i].filptr," % *.*E\n",w,p,fz);
          }
          fflush(oinfo[i].filptr);
          break;
        case OutputInfo::AeroXMom:
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if(oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode = first_node; iNode<last_node; ++iNode) {
            int xrot  = c_dsa->locate( iNode, DofSet::Xrot);
            double mx = (xrot >= 0) ? ext_f[xrot] : 0.0;
            fprintf(oinfo[i].filptr," % *.*E\n",w,p,mx);
          }
          fflush(oinfo[i].filptr);
          break;
        case OutputInfo::AeroYMom:
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode=first_node; iNode<last_node; ++iNode) {
            int yrot  = c_dsa->locate( iNode, DofSet::Yrot);
            double my = (yrot >= 0) ? ext_f[yrot] : 0.0;
            fprintf(oinfo[i].filptr," % *.*E\n",w,p,my);
          }
          fflush(oinfo[i].filptr);
          break;
        case OutputInfo::AeroZMom:
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode=first_node; iNode<last_node; ++iNode) {
            int zrot  = c_dsa->locate( iNode, DofSet::Zrot);
            double mz = (zrot >= 0) ? ext_f[zrot] : 0.0;
            fprintf(oinfo[i].filptr," % *.*E\n",w,p,mz);
          }
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::DispX:
          if(dof==-1) dof = DofSet::Xdisp;
        case OutputInfo::DispY:
          if(dof==-1) dof = DofSet::Ydisp;
        case OutputInfo::DispZ:
          if(dof==-1) dof = DofSet::Zdisp;
        case OutputInfo::RotX:
          if(dof==-1) dof = DofSet::Xrot;
        case OutputInfo::RotY:
          if(dof==-1) dof = DofSet::Yrot;
        case OutputInfo::RotZ:
          if(dof==-1) dof = DofSet::Zrot;
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode=first_node; iNode<last_node; iNode++){
            if(nodes[iNode] == 0) continue;
            double data;
            getOrAddDofForPrint(get, d_n, bcx, iNode, &data, &dof);
            fprintf(oinfo[i].filptr," % *.*E\n",w,p,data);
          }
          break;
                                                                                                                                       
        case OutputInfo::DispMod:
        case OutputInfo::RotMod:
        case OutputInfo::TotMod:
          fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          if(oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"\n");
          for (iNode=first_node; iNode<last_node; iNode++){
            if(nodes[iNode] == 0) continue;
            double x,y,z,xr,yr,zr;
            int dofx,dofy,dofz;
            if (oinfo[i].type == OutputInfo::DispMod || oinfo[i].type == OutputInfo::TotMod) {
              dofx=DofSet::Xdisp; dofy=DofSet::Ydisp; dofz=DofSet::Zdisp;
              getOrAddDofForPrint(get, d_n, bcx, iNode, &x, &dofx, &y, &dofy, &z, &dofz);
            }
            if (oinfo[i].type == OutputInfo::RotMod || oinfo[i].type == OutputInfo::TotMod){
              dofx=DofSet::Xrot; dofy=DofSet::Yrot; dofz=DofSet::Zrot;
              getOrAddDofForPrint(get, d_n, bcx, iNode, &xr, &dofx, &yr, &dofy, &zr, &dofz);
            }
            if (oinfo[i].type == OutputInfo::DispMod) fprintf(oinfo[i].filptr," % *.*E\n",w,p,sqrt(x*x+y*y+z*z));
            if (oinfo[i].type == OutputInfo::RotMod)  fprintf(oinfo[i].filptr," % *.*E\n",w,p,sqrt(xr*xr+yr*yr+zr*zr));
            if (oinfo[i].type == OutputInfo::TotMod)
              fprintf(oinfo[i].filptr," % *.*E\n",w,p,sqrt(x*x+y*y+z*z+xr*xr+yr*yr+zr*zr));
          }
          break;
/*
        case OutputInfo::Accel6:
        case OutputInfo::Acceleration:
            //CD: as it is not quite correct, it wasn't modified by using
            //    getOrAddDofForPrint 
          if (oinfo[i].nodeNumber == -1) fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
          else fprintf(oinfo[i].filptr,"  % *.*E  ",w,p,time);
          for (iNode=first_node; iNode<last_node; iNode++){

            double x,y,z;
            int xloc  = c_dsa->locate( iNode, DofSet::Xdisp);
            int xloc1 =   dsa->locate( iNode, DofSet::Xdisp);
            if(xloc >= 0) x = a_n[xloc];
            else if(xloc1 >= 0)  x = 0.0; // not quite correct, need to pass acceleration for constrained points  x=bcx[xloc1];
            else x = 0.0;
           
            int yloc  = c_dsa->locate( iNode, DofSet::Ydisp);
            int yloc1 =   dsa->locate( iNode, DofSet::Ydisp);
            if (yloc >= 0) y = a_n[yloc]; 
            else if(yloc1 >= 0) y = 0.0; 
            else y = 0.0; //same problem

            int zloc  = c_dsa->locate( iNode, DofSet::Zdisp);
            int zloc1 =   dsa->locate( iNode, DofSet::Zdisp);
            if (zloc >= 0) z = a_n[zloc]; 
            else if(zloc1 >= 0) z = 0.0; 
            else z = 0.0;
          
            if (oinfo[i].type == OutputInfo::Accel6) {
         
              double xr,yr,zr;
              int xloc  = c_dsa->locate( iNode, DofSet::Xrot);
              int xloc1 =   dsa->locate( iNode, DofSet::Xrot);
              if(xloc >= 0) xr = a_n[xloc];
              else if(xloc1 >= 0)  xr = 0.0; // not quite correct, need to pass acceleration for constrained points  x=bcx[xloc1];
              else xr = 0.0;
     
              int yloc  = c_dsa->locate( iNode, DofSet::Yrot);
              int yloc1 =   dsa->locate( iNode, DofSet::Yrot);
              if(yloc >= 0) yr = a_n[yloc]; else if(yloc1 >= 0) yr = 0.0; else yr = 0.0; //same problem
 
              int zloc  = c_dsa->locate( iNode, DofSet::Zrot);
              int zloc1 =   dsa->locate( iNode, DofSet::Zrot);
              if(zloc >= 0) zr = a_n[zloc]; else if(zloc1 >= 0) zr = 0.0; else zr = 0.0;

              fprintf(oinfo[i].filptr,"%d  % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n"
                                        ,iNode+1,w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
            }
            else if (oinfo[i].nodeNumber == -1){
              fprintf(oinfo[i].filptr," % *.*E % *.*E % *.*E\n",iNode+1,w,p,x,w,p,y,w,p,z);
            }
            else{ fprintf(oinfo[i].filptr," % *.*E % *.*E % *.*E\n",w,p,x,w,p,y,w,p,z); }
          }
          break;
*/            
        case OutputInfo::YModulus:
          getElementAttr(i,YOUNG);
          break;
        case OutputInfo::MDensity:
          getElementAttr(i,MDENS);
          break;
        case OutputInfo::Thicknes:
          getElementAttr(i,THICK);
          break;
        case OutputInfo::Composit:
          getCompositeData(i,time);
          break;          
        case OutputInfo::TDEnforcement: {
          double *plot_data = new double[numNodes];
          for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.0;
          for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
            MortarConds[iMortar]->get_plot_variable(oinfo[i].tdenforc_var,plot_data);
          }
          if(oinfo[i].nodeNumber == -1)
            geoSource->outputNodeScalars(i, plot_data, numNodes, time);
          else
            geoSource->outputNodeScalars(i, &plot_data[oinfo[i].nodeNumber], 1, time);
          delete [] plot_data;
        } break;
        default:
          break;
        
      }
    }
  }
       
  if (glDisp) delete [] glDisp;

}

//----------------------------------------------------------------------------------------------

/*
void
Domain::outputHeader(int fileNumber)
{
  // If only one node is requested for output,
  // then do not write the header to the file as this
  // output file will be used as gnuplot input
  if(oinfo[fileNumber].nodeNumber != -1) return;

  char prbType[20];

  if (sinfo.probType == SolverInfo::Static )       strcpy(prbType,"Static");
  if (sinfo.probType == SolverInfo::Dynamic)       strcpy(prbType,"Dynam");
  if (sinfo.probType == SolverInfo::Helmholtz)     strcpy(prbType,"FAcoustic");
  if (sinfo.probType == SolverInfo::NonLinStatic)  strcpy(prbType,"NLStatic");
  if (sinfo.probType == SolverInfo::NonLinDynam)   strcpy(prbType,"NLDynamic");
  if (sinfo.probType == SolverInfo::ArcLength)     strcpy(prbType,"Arclength");
  if (sinfo.probType == SolverInfo::TempDynamic)   strcpy(prbType,"Temp");
  if (sinfo.probType == SolverInfo::AxiHelm)       strcpy(prbType,"AxiHelm");

  if(oinfo[fileNumber].type == OutputInfo::EigenPair) {
      strcpy(prbType,"modal");
  }
  if(oinfo[fileNumber].type == OutputInfo::YModulus) {
     strcpy(prbType,"Attributes");
  }
  if(oinfo[fileNumber].type == OutputInfo::MDensity) {
     strcpy(prbType,"Attributes");
  }
  if(oinfo[fileNumber].type == OutputInfo::Thicknes) {
     strcpy(prbType,"Attributes");
  }

  if(oinfo[fileNumber].type == OutputInfo::ShapeAtt) {
     strcpy(prbType,"Attributes");
  }

  if(oinfo[fileNumber].type == OutputInfo::ShapeStc) {
     strcpy(prbType,"Static");
  }
  
  int avgnum = oinfo[fileNumber].averageFlg;
  int type   = oinfo[fileNumber].type;
  //int node   = oinfo[fileNumber].nodeNumber;

  // No header for AeroForce
  if(oinfo[fileNumber].type == OutputInfo::AeroForce) return;

  // No header for CompositeData
  if(oinfo[fileNumber].type == OutputInfo::Composit) return;

  if(oinfo[fileNumber].type == OutputInfo::InXForce || 
     oinfo[fileNumber].type == OutputInfo::InYForce ||
     oinfo[fileNumber].type == OutputInfo::InZForce || 
     oinfo[fileNumber].type == OutputInfo::AXMoment ||
     oinfo[fileNumber].type == OutputInfo::AYMoment || 
     oinfo[fileNumber].type == OutputInfo::AZMoment ) {
    fprintf(oinfo[fileNumber].filptr,header[type],prbType,
            cinfo->nodeSetName,numele);
    return;
  }
*/
  // Compute the last node number
/*
  int exactNumNodes = 0;
  int inode;
  for(inode=0; inode<numnodes; ++inode) {
//    if(nodes[inode] && ((*dsa)[inode].list() != 0 ))
    if(nodes[inode])
       exactNumNodes = (exactNumNodes+1);
  }
*/
/*
  int exactNumNodes = numnodes;

  if(avgnum == 1)
    fprintf(oinfo[fileNumber].filptr,header[type],prbType,
                                     cinfo->nodeSetName,exactNumNodes);
  else
    fprintf(oinfo[fileNumber].filptr,header[type],prbType,
                                     cinfo->nodeSetName,numele);

}
*/

ControlInterface*
Domain::getUserSuppliedFunction()
{
  ControlInterface *userSupFunc = 0;

#ifndef TFLOP
#ifndef WINDOWS
  if( claw ) {
    void *handle;


    dlerror(); // forget about the last error
    handle = dlopen(claw->fileName, RTLD_NOW);
    char *errorMsg;
    if((errorMsg = dlerror() ) != 0) {
      fprintf(stderr," *** ERROR: in dynamic loading of %s: %s\n",
             claw->fileName,errorMsg);
      exit(-1);
    }

    ControlInterface ** fcp =
        (ControlInterface **) dlsym(handle, claw->routineName);

    if(fcp ==0) {
      fprintf(stderr," *** ERROR: in dynamic loading of %s: "
                     "control function not found\n",
                     claw->routineName);
      exit(-1);
    }

    userSupFunc = *fcp;
 }
#endif
#else
// Should have something here
#endif

 return userSupFunc;

}

void
Domain::aeroPreProcess(Vector& d_n, Vector& v_n, Vector& a_n,
                       Vector& v_p, double *bcx, double *vcx)
{
  if(sinfo.aeroFlag >= 0) {

    Vector d_n_aero(d_n);
 
    // Send u + IDISP6 to fluid code.
    // IDISP6 is used to compute pre-stress effects.

    if(sinfo.gepsFlg == 1) {
      // If we are in the first time step, and we initialized with
      // IDISP6, do not send IDISP6
      if(numIDis == 0 && sinfo.zeroInitialDisp != 1) {
        fprintf(stderr," ... DO NOT SEND IDISP6 0\n"); //HB
      } else {
        fprintf(stderr," ... SENDING IDISP6 0\n"); //HB
        int i;
        for(i = 0; i < numIDis6; ++i) {
          int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
          if(dof >= 0)
            d_n_aero[dof] += iDis6[i].val;
        }
      }
    }

    State curState(c_dsa, dsa, bcx, vcx, d_n_aero, v_n, a_n, v_p);

    int numOutInfo = geoSource->getNumOutInfo();
    OutputInfo *oinfo = geoSource->getOutputInfo();

    int flag = 0;

    // Check if aero forces are requested for output
    int iInfo;
    for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
      if(oinfo[iInfo].type == OutputInfo::AeroForce) { 
        flag = 1;
        break;
      }
    }

    if(flag)
      flExchanger = new FlExchanger(packedEset, c_dsa, oinfo+iInfo);
    else
      flExchanger = new FlExchanger(packedEset, c_dsa );

    char *matchFile = geoSource->getMatchFileName();
    if(matchFile == 0)
      matchFile = "MATCHER";
    flExchanger->read(0, matchFile);

    //XML New step of negotiation with fluid code
    flExchanger->negotiate();

    int restartinc = (solInfo().nRestart >= 0) ? (solInfo().nRestart) : 0;


    if(sinfo.aeroFlag == 8) {
      flExchanger->sendParam(sinfo.aeroFlag, sinfo.dt, sinfo.mppFactor,
                             restartinc, sinfo.isCollocated, sinfo.alphas);
      flExchanger->sendModeFreq(modeData.frequencies, modeData.numModes);
      flExchanger->sendModeShapes(modeData.numModes, modeData.numNodes,
                   modeData.modes, nodes, curState, sinfo.mppFactor);
      return;
    }
    else
      flExchanger->sendParam(sinfo.aeroFlag, sinfo.dt, sinfo.tmax, restartinc,
                             sinfo.isCollocated, sinfo.alphas);
    if(verboseFlag) fprintf(stderr,"... [E] Sent parameters and modes ...\n");

    if(sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
      flExchanger->initRcvParity(1);
      flExchanger->initSndParity(1);
    } else {
      flExchanger->initRcvParity(-1);
      flExchanger->initSndParity(-1);
    }

    flExchanger->sendDisplacements(nodes, curState);
    if(verboseFlag) fprintf(stderr,"... [E] Sent initial displacements ...\n");

    if(sinfo.aeroFlag == 1) { // Ping pong only
      fprintf(stderr, "Ping Pong Only requested. Structure code exiting\n");
      return;
    }
  }
}

void
Domain::thermoePreProcess()
{
  if(sinfo.thermoeFlag >=0) {

    int buffLen = numnodes;

    temprcvd = new double[numnodes]; // Initialize received temperature

    // if sinfo.aeroFlag >= 0, flExchanger has already been initialize before,
    // thus, only when sinfo.aeroFlag < 0 is necessary.
    if(sinfo.aeroFlag < 0)
      flExchanger = new FlExchanger(packedEset, c_dsa );

    flExchanger->thermoread(buffLen);
 
    flExchanger->getStrucTemp(temprcvd) ;
    //if(verboseFlag) fprintf(stderr,"... [E] Received initial temperatures ...\n");
  }
}

/*
** This subroutine computes the maximum stability time step,
** which is equal to 2.0 divided by the maximum natural frequency
** of the structure.  The maxium natural frequency is computed
** by power method.
*/

double
Domain::computeStabilityTimeStep(DynamMat& dMat)
{
      if(outFile) 
        fprintf(stderr, " ... Checking Newmark Stability     ...\n");

      double eigmax;
      double relTol    = sinfo.stable_tol; // stable_tol default is 1.0e-3
      double preeigmax = 0.0;

      int numdofs = dMat.K->dim();
      int maxIte  = sinfo.stable_maxit; // stable_maxit default is 100

      Vector v(numdofs);
      Vector z(numdofs);

// Starts from an arbitrary array.
      int i,j;
      for (i=0; i<numdofs; ++i)
        v[i] = (double) (i+1) / (double) numdofs;

// Power iteration loop

      for (i=0; i<maxIte; ++i) {
        dMat.K->mult(v,z);

        for (j=0; j< numdofs; ++j)
          z[j] /= dMat.M->diag(j);

// Normalize

        double zmax = z[0];
        for (j=1; j< numdofs; ++j)
          if (abs(z[j])>zmax) zmax = abs(z[j]);

        eigmax = zmax;

        v = (1.0/zmax)*z;

        if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;

        preeigmax = eigmax;
      }

      // compute stability maximum time step
      double sdt = 2.0 / sqrt(eigmax);

      return sdt;
}

//------------------------------------------------------------------------------


double Domain::getKineticEnergy( Vector& vel, SparseMatrix * gMass )
{
  double energy=0.0;
  
  Vector tmpVec(c_dsa->size());
    
  gMass->mult(vel,tmpVec);

  energy = 0.5 * (vel * tmpVec);
  
  return energy;
}

//------------------------------------------------------------------------------


double Domain::getDampingEnergy( Vector& vel, SparseMatrix * gDamp )
{
  double energy=0.0;
  
  Vector tmpVec(c_dsa->size());
    
  gDamp->mult(vel,tmpVec);

  energy = 0.5 * (vel * tmpVec);
  
  return energy;
}
