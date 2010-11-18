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
 if(numIVelModal) {
   filePrint(stderr, " ... Compute initial velocity from given modal basis (v0=X.y0) ... \n"); //HB
   modeData.addMultY(numIVelModal, iVelModal, v_n, c_dsa);
 }
 for(int i = 0; i < numIVel; ++i) {
   int dof = c_dsa->locate(iVel[i].nnum, 1 << iVel[i].dofnum);
   if(dof >= 0)
     v_n[dof] += iVel[i].val;
 }

 // If zeroInitialDisp is set, then return as the displacement is
 // already set to zero above.
 // note: geps is always zero for nonlinear
 if(sinfo.zeroInitialDisp == 0) {
   // ... SET INITIAL DISPLACEMENT FROM IDISP IF IDISP6 DOES NOT EXIST
   // ... OR IF WE ARE USING GEOMETRIC PRE-STRESS (GEPS)
   //if(domain->numInitDisp6() == 0 || sinfo.gepsFlg == 1) { // note: always use global num to do this check
   if(sinfo.gepsFlg == 1) { // note: always use global num to do this check
     if(numIDisModal) {
       filePrint(stderr, " ... Compute initial displacement from given modal basis (u0=X.y0) ... \n"); //HB
       modeData.addMultY(numIDisModal, iDisModal, d_n, c_dsa);
     }
     for(int i = 0; i < numIDis; ++i) {
       int dof = c_dsa->locate(iDis[i].nnum, 1 << iDis[i].dofnum);
       if(dof >= 0) d_n[dof] += iDis[i].val;
     }
   }

   // ... SET INITIAL DISPLACEMENT FROM IDISP6
   if((domain->numInitDisp() == 0 && domain->numInitDispModal() == 0) || sinfo.gepsFlg == 0) {
     for(int i = 0; i < numIDis6; ++i) {
       int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
       if(dof >= 0)
         d_n[dof] = iDis6[i].val;
     }   
     // also add any modal idisps set under IDISP
     if(numIDisModal) {
       filePrint(stderr, " ... Compute initial displacement from given modal basis (u0=X.y0) ... \n"); //HB
       modeData.addMultY(numIDisModal, iDisModal, d_n, c_dsa);
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
 if(timeIndex % sinfo.nRestart == 0 || time >= sinfo.tmax-0.1*sinfo.getTimeStep()) {
   int fn = open(cinfo->currentRestartFile, O_WRONLY | O_CREAT, 0666);
   if(fn >= 0) {
     cerr << "here in Dynam.C writeRestartFile\n";
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
Domain::aeroSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState)
{
  // Send u + IDISP6 to fluid code.
  // IDISP6 is used to compute pre-stress effects.
  getTimers().sendFluidTime -= getTime();
  Vector d_n_aero(d_n);

  if(sinfo.gepsFlg == 1) {

    for(int i = 0; i < numIDis6; ++i) {
      int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
      if(dof >= 0)
        d_n_aero[dof] += iDis6[i].val;
    }
  }

  State state(c_dsa, dsa, bcx, vcx, d_n_aero, v_n, a_n, v_p);

  flExchanger->sendDisplacements(state, -1, geomState);
  if(verboseFlag) fprintf(stderr, " ... [E] Sent displacements ...\n");

  getTimers().sendFluidTime += getTime();
}

void
Domain::aeroheatSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState)
{
  State tempState(c_dsa, dsa, bcx, d_n, v_n, v_p);

  flExchanger->sendTemperature(tempState);
  if(verboseFlag) fprintf(stderr," ... [T] Sent temperatures ...\n");
}

void
Domain::thermohSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState)
{
/* we have to send the vector of temperatures in NODAL order, not
int DOF order (in which is d_n)!
print(), zero() can be put behind VECTORS!
numUncon() is a function and sets the vector in its correspond. size.
d_n.print("comment"); 
*/

  int iNode;
  Vector tempsent(numnodes);

  for(iNode=0; iNode<numnodes; ++iNode) {
    int tloc  = c_dsa->locate( iNode, DofSet::Temp);
    int tloc1 =   dsa->locate( iNode, DofSet::Temp);
    double temp  = (tloc >= 0) ? d_n[tloc] : bcx[tloc1];
    if(tloc1 < 0) temp = 0.0;
    tempsent[iNode] = temp;
  }

  flExchanger->sendStrucTemp(tempsent);
  if(verboseFlag) filePrint(stderr," ... [T] Sent temperatures ...\n");
}

void
Domain::buildAeroelasticForce(Vector& aero_f, PrevFrc& prevFrc, int tIndex, double t, double gamma, double alphaf, GeomState* geomState)
{
  // ... COMPUTE AEROELASTIC FORCE 
  getTimers().receiveFluidTime -= getTime();

  // ... Temporary variable for inter(extra)polated force
  double *tmpFmem = new double[numUncon()];
  StackVector tmpF(tmpFmem, numUncon());
  tmpF.zero();

  int iscollocated;
  double tFluid = flExchanger->getFluidLoad(tmpF, tIndex, t,
                                            alphaf, iscollocated, geomState);
  if(verboseFlag) filePrint(stderr," ... [E] Received fluid load ...\n");

  if(sinfo.aeroFlag == 20) {
    if(prevFrc.lastTIndex >= 0)
      aero_f.linC(0.5,tmpF,0.5,prevFrc.lastFluidLoad);
    else
      aero_f = tmpF;
  }
  else {
    if(iscollocated == 0) {
      if(prevFrc.lastTIndex >= 0) {
        tmpF *= (1/gamma);
        tmpF.linAdd(((gamma-1.0)/gamma), prevFrc.lastFluidLoad);
      }
    }

    double alpha = (prevFrc.lastTIndex < 0) ? 1.0 : 1.0-alphaf;
    aero_f.linC(alpha, tmpF, (1.0-alpha), prevFrc.lastFluidLoad);
  }
  prevFrc.lastFluidLoad = tmpF;
  prevFrc.lastFluidTime = tFluid;
  prevFrc.lastTIndex = tIndex;

  //fprintf(stderr,"... alpha = %e, gamma = %e, isCollocated = %d ...\n", alpha, gamma, iscollocated);

  delete [] tmpFmem;
  getTimers().receiveFluidTime += getTime();
}

void
Domain::buildAeroheatFlux(Vector &f, Vector &prev_f, int tIndex, double t)
{
  // ... ADD FLUID FLUX
  getTimers().receiveFluidTime -= getTime();

  // ... Temporary variable for inter(extra)polated force
  double *tmpFmem = (double *) dbg_alloca(sizeof(double)*numUncon());
  StackVector tmpF(tmpFmem, numUncon());

  double sflux = 0;
  double bfl ;

/*  linAdd is in Maths.d/Vector.C , f.linAdd(a,g,b,h) 
    means f = f+a*g+b*h 
    f.linAdd(a,g)= f+a*g  */

  tmpF.zero();
  flExchanger->getFluidFlux(tmpF, t, bfl);

  if(verboseFlag) filePrint(stderr, " ... [T] Received fluid fluxes ...\n");

  int vectlen = tmpF.size();

/*  Compute fluid flux at n+1/2, since we use midpoint rule in thermal */

  int useProjector = domain->solInfo().filterFlags;

  if(tIndex == 0)
    f += tmpF;
  else {
    if(useProjector) f = tmpF;
    else
      f.linAdd(0.5, tmpF, 0.5, prev_f);
  }

  prev_f = tmpF;

/*  Print out sum of fluxes received at n+1 and at n+1/2 in a file
    specified by raerotfl in entry */

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  for(int i=0; i < numOutInfo; i++) {
    if(oinfo[i].interval != 0 && tIndex % oinfo[i].interval == 0) {
      switch(oinfo[i].type) {
        case OutputInfo::AeroFlux:
          fprintf(oinfo[i].filptr,"%e   ",t);
          fprintf(oinfo[i].filptr,"%e   ",bfl);
          for(int j=0; j<vectlen; j++)  {
            sflux += 0.5*(tmpF[j]+prev_f[j]);
          }
          fprintf(oinfo[i].filptr,"%e\n",sflux);
          fflush(oinfo[i].filptr);
        break;
        default: /* do nothing */
        break;
      }
    }
  }
}

void
Domain::thermoeComm()
{
  flExchanger->getStrucTemp(temprcvd);
  if(verboseFlag) fprintf(stderr," ... [E] Received temperatures ...\n");
  //buildThermalForce(temprcvd, f, geomState);
}

void
Domain::dynamOutput(int tIndex, double *bcx, DynamMat& dMat, Vector& ext_f, Vector &aeroForce, 
                    Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* vcx) {
  // Current time stamp
  double time = tIndex*sinfo.getTimeStep();
   
  if(sinfo.nRestart > 0 && !sinfo.modal && probType() != SolverInfo::NonLinDynam) {
    writeRestartFile(time, tIndex, d_n, v_n, v_p, sinfo.initExtForceNorm);
  }
  
  // Send to fluid code (except C0)
  if(sinfo.aeroFlag >= 0 && !sinfo.lastIt && tIndex != sinfo.initialTimeIndex && sinfo.aeroFlag != 20)
    aeroSend(d_n, v_n, a_n, v_p, bcx, vcx);

  if(sinfo.aeroheatFlag >= 0 && tIndex != 0 )
    aeroheatSend(d_n, v_n, a_n, v_p, bcx, vcx);

  if(sinfo.thermohFlag >=0 && tIndex != sinfo.initialTimeIndex) 
    thermohSend(d_n, v_n, a_n, v_p, bcx, vcx);

  int numOutInfo = geoSource->getNumOutInfo();

  // Open the file and write the header in the first time step
  if (tIndex == sinfo.initialTimeIndex) {
    if(numOutInfo > 0)
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
  int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes
  int numNodeLim = myMax(numNodes,numnodes);
  double (*glDisp)[11] = new double[numNodeLim][11];//DofSet::max_known_nonL_dof
  for (int i = 0; i < numNodeLim; ++i)
    for (int j = 0 ; j < 11 ; j++)
      glDisp[i][j] = 0.0;
  mergeDistributedDisp(glDisp, d_n.data(), bcx);

  for (int i = firstRequest; i < lastRequest; ++i) {
    enum {YOUNG,MDENS,THICK};
    int iNode;
    int first_node, last_node;
    
    OutputInfo *oinfo = geoSource->getOutputInfo();
    
    //CD: ad and get will be used in  addVariationOfShape_StructOpt and getOrAddDofForPrint which were 
    //    added to "clean" dynamOutput 
    bool get=false;
     
    if (oinfo[i].nodeNumber == -1){first_node=0; last_node=numNodes;}           
    else { first_node=oinfo[i].nodeNumber; last_node=first_node+1;}
    //first_node=0; last_node=numNodes;     

    if ((oinfo[i].interval != 0) && (tIndex % oinfo[i].interval == 0)) {
      int dof=-1;
      int w = oinfo[i].width;
      int p = oinfo[i].precision;

      int success = processDispTypeOutputs(oinfo[i], glDisp, numNodes, i, time);
      if (success) continue;
      success = processOutput(oinfo[i].type, d_n, bcx, i, time);
      if (success) continue;
      success = 1;


      int nNodes = last_node-first_node;
      switch(oinfo[i].type) {

        case OutputInfo::Velocity6: {

          double (*data)[6] = new double[nNodes][6];

          for (iNode = 0; iNode < nNodes; ++iNode)  {

            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data[iNode], &DofSet::Xdisp,
                                data[iNode]+1, &DofSet::Ydisp, data[iNode]+2, &DofSet::Zdisp);

            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data[iNode]+3, &DofSet::Xrot,
                                data[iNode]+4, &DofSet::Yrot, data[iNode]+5, &DofSet::Zrot);
          }
          geoSource->outputNodeVectors6(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::Velocity:  {

          double (*data)[3] = new double[nNodes][3];

          for (iNode = 0; iNode < nNodes; ++iNode)  {
            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data[iNode], &DofSet::Xdisp,
                                data[iNode]+1, &DofSet::Ydisp, data[iNode]+2, &DofSet::Zdisp);

          }
          geoSource->outputNodeVectors(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::PressureFirstTimeDerivative: {

          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)
            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data+iNode, &DofSet::Helm, 0, 0, 0, 0);
          
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;

        } 
          break;
        case OutputInfo::TemperatureFirstTimeDerivative: {
          double *data = new double[nNodes]; 
          for (iNode = 0; iNode < nNodes; ++iNode)
            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data+iNode, &DofSet::Temp, 0,0,0,0);
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::Accel6:  {
          double (*data)[6] = new double[nNodes][6];
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            getOrAddDofForPrint(get, a_n, (double *) 0 /*acx*/, first_node+iNode, data[iNode],
                                &DofSet::Xdisp, data[iNode]+1, &DofSet::Ydisp,
                                data[iNode]+2, &DofSet::Zdisp);
            getOrAddDofForPrint(get, a_n, (double *) 0 /*acx*/, first_node+iNode, data[iNode]+3, 
                                &DofSet::Xrot, data[iNode]+4, &DofSet::Yrot, data[iNode]+5, 
                                &DofSet::Zrot);
          }
          geoSource->outputNodeVectors6(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::Acceleration:  {
          // XXXX acx (prescribed accelerations not implemented)
          double (*data)[3] = new double[nNodes][3];

          for (iNode = 0; iNode < nNodes; ++iNode)  {
            getOrAddDofForPrint(get, a_n, (double *) 0 /*acx*/, first_node+iNode, data[iNode], 
            &DofSet::Xdisp, data[iNode]+1, &DofSet::Ydisp, data[iNode]+2, &DofSet::Zdisp);
          }
          geoSource->outputNodeVectors(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::PressureSecondTimeDerivative: {
          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)
            getOrAddDofForPrint(get, a_n, (double *) 0, first_node+iNode, data+iNode, &DofSet::Helm, 
                                0, 0, 0, 0);
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
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
          if (dMat.C && (time > sinfo.initialTime)) {  //??????????????????????
            dMat.C->mult(v_n, tmpVec);
            double c = solInfo().newmarkGamma;
            dWdmp = (c*tmpVec + (1.0-c)*(*previousCq))*(d_n - (*previousDisp));
            Wdmp += dWdmp;
          }
          //??????????????????????
          double error = (time==sinfo.initialTime) ? 0.0 : (Wela+Wkin)-(pWela+pWkin)+dWdmp-dW; 
          geoSource->outputEnergies(i, time, Wext, Waero, Wela, Wkin, -Wdmp, error);
          pWela=Wela;
          pWkin=Wkin;
          if (time==sinfo.initialTime) {  //?????????????????????
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
       
        case OutputInfo::AeroXForce:  {
          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
            data[iNode]  = (xloc >= 0) ? aeroForce[xloc] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::AeroYForce:  {
          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
            data[iNode]  = (yloc >= 0) ? aeroForce[yloc] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::AeroZForce:  {
          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
            data[iNode] = (zloc >= 0) ? aeroForce[zloc] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::AeroXMom:  {
          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
            data[iNode] = (xrot >= 0) ? aeroForce[xrot] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::AeroYMom:  {
          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
            data[iNode] = (yrot >= 0) ? aeroForce[yrot] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
          break;
        case OutputInfo::AeroZMom:  {
          double *data = new double[nNodes];
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
            data[iNode] = (zrot >= 0) ? aeroForce[zrot] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodes, time);
          delete [] data;
        }
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
        case OutputInfo::ModeError: // don't print warning message since these are
        case OutputInfo::ModeAlpha: // output in SingleDomainDynamic::modeDecomp
          break;

        default:
          success = 0;
          break;
        
      }
      if (success == 0)
        fprintf(stderr, " *** AS.WRN: output %d is not supported \n", i);
    }
  }
       
  
  if (glDisp) delete [] glDisp;

}

//----------------------------------------------------------------------------------------------

ControlInterface* Domain::getUserSuppliedFunction() {

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

    if(aeroEmbeddedSurfaceId.size()!=0) {
      int iSurf = -1;
      for(int i=0; i<nSurfEntity; i++)
        if(aeroEmbeddedSurfaceId.find(SurfEntities[i]->ID())!=aeroEmbeddedSurfaceId.end()) {
          iSurf = i; 
          break; //only allows one Surface.
        }
      if(iSurf<0) {
        fprintf(stderr,"ERROR: Embedded wet surface not found! Aborting...\n");
        exit(-1);
      }
      flExchanger = new FlExchanger(nodes, packedEset, SurfEntities[iSurf], c_dsa); //packedEset is not used, but flExchanger needs
                                                                                    // to have a reference of it at construction.
    } else {
      if(flag)
        flExchanger = new FlExchanger(nodes, packedEset, c_dsa, oinfo+iInfo);
      else
        flExchanger = new FlExchanger(nodes, packedEset, c_dsa );
    }

    char *matchFile = geoSource->getMatchFileName();
    if(matchFile == 0)
      matchFile = (char*) "MATCHER";

    if(aeroEmbeddedSurfaceId.size()!=0) 
      flExchanger->matchup();
    else
      flExchanger->read(0, matchFile);

    //KW: send the embedded wet surface to fluid 
    if(aeroEmbeddedSurfaceId.size()!=0) {
      flExchanger->sendEmbeddedWetSurface();
      if(verboseFlag) fprintf(stderr,"... [E] Sent embedded wet surface ...\n");
    }

    //XML New step of negotiation with fluid code
    flExchanger->negotiate();

    int restartinc = (solInfo().nRestart >= 0) ? (solInfo().nRestart) : 0;

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

    if(sinfo.aeroFlag == 8) {
      flExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.mppFactor,
                             restartinc, sinfo.isCollocated, sinfo.alphas);
      flExchanger->sendModeFreq(modeData.frequencies, modeData.numModes);
      if(verboseFlag) fprintf(stderr," ... [E] Sent parameters and mode frequencies ...\n");
      flExchanger->sendModeShapes(modeData.numModes, modeData.numNodes,
                   modeData.modes, curState, sinfo.mppFactor);
      if(verboseFlag) fprintf(stderr," ... [E] Sent mode shapes ...\n");
    }
    else {
      flExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.tmax, restartinc,
                             sinfo.isCollocated, sinfo.alphas);
      if(verboseFlag) fprintf(stderr," ... [E] Sent parameters ...\n");

      if(sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
        flExchanger->initRcvParity(1);
        flExchanger->initSndParity(1);
      } else {
        flExchanger->initRcvParity(-1);
        flExchanger->initSndParity(-1);
      }

      flExchanger->sendDisplacements(curState);
      if(verboseFlag) fprintf(stderr," ... [E] Sent initial displacements ...\n");

      if(sinfo.aeroFlag == 1) { // Ping pong only
        fprintf(stderr, "Ping Pong Only requested. Structure code exiting\n");
      }
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
      flExchanger = new FlExchanger(nodes, packedEset, c_dsa );

    flExchanger->thermoread(buffLen);
 
    flExchanger->getStrucTemp(temprcvd) ;
    if(verboseFlag) fprintf(stderr,"... [E] Received initial temperatures ...\n");
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

void
Domain::updateUsddInDbc(double* userDefineDisp, int* map)
{
  int j = 0;
  for(int i = 0; i < numDirichlet; ++i)
    if(dbc[i].type == BCond::Usdd) {
      int k = (map) ? map[j] : j;
      dbc[i].val = userDefineDisp[k];
      j++;
    }
}

void
Domain::updateUsdfInNbc(double* userDefineForce, int* map, double* weight)
{
  int j = 0;
  for(int i = 0; i < numNeuman; ++i) 
    if(nbc[i].type == BCond::Usdf) {
      int k = (map) ? map[j] : j;
      nbc[i].val = userDefineForce[k];
      if(weight) {
        int dof = c_dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
        nbc[i].val *= weight[dof];
      }
      j++;
    }
}

void
Domain::updateActuatorsInNbc(double* actuatorsForce, int* map, double* weight)
{
  int j = 0;
  for(int i = 0; i < numNeuman; ++i)   
    if(nbc[i].type == BCond::Actuators) {
      int k = (map) ? map[j] : j;
      nbc[i].val = actuatorsForce[k];
      if(weight) {
        int dof = c_dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
        nbc[i].val *= weight[dof];
      }
      j++;
    }
}

