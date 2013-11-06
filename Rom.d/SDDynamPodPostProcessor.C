#include <Driver.d/SysState.h>
#include "PodProjectionNonLinDynamic.h"

namespace Rom {

SDDynamPodPostProcessor::SDDynamPodPostProcessor(Domain *d, double *bcx, double *vcx, double *acx,
                                                 StaticTimers* times, GeomState *geomState, Corotator **corot) :
    SDDynamPostProcessor(d, bcx, vcx, acx, times, geomState, allCorot),
    DispSensorValues(NULL),
    AccSensorValues(NULL),
    VelSensorValues(NULL)
{
  oinfo = geoSource->getOutputInfo();
  numOutInfo = geoSource->getNumOutInfo();

  DispSensor = false;
  AccSensor  = false;
  VelSensor  = false;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber != -1) {
          if(oinfo[iOut].type == OutputInfo::Acceleration || (oinfo[iOut].angularouttype == OutputInfo::total && oinfo[iOut].rescaling == false)) {
            AccSensor = true;
          }
          else {
            filePrint(stderr, " *** WARNING: unsupported probe output type. Use OUTPUT instead of OUTPUT6, or\n");
            filePrint(stderr, "     OUTPUT6 with ROTVECOUTTYPE=euler, ANGULAROUTTYPE=total, and RESCALING=off.\n");
          }
        } else {
          oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          filePrint(oinfo[iOut].filptr, "0\n"); 
        }
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement :
        if(oinfo[iOut].nodeNumber != -1) {
          if(oinfo[iOut].type == OutputInfo::Displacement || (oinfo[iOut].rotvecouttype == OutputInfo::Euler && oinfo[iOut].rescaling == false)) {
            DispSensor = true;
          }
          else {
            filePrint(stderr, " *** WARNING: unsupported probe output type. Use OUTPUT instead of OUTPUT6, or\n");
            filePrint(stderr, "     OUTPUT6 with ROTVECOUTTYPE=euler, ANGULAROUTTYPE=total, and RESCALING=off.\n");
          }
        } else {
          oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          filePrint(oinfo[iOut].filptr, "1\n");
        }
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity :
        if(oinfo[iOut].nodeNumber != -1) {
          if(oinfo[iOut].type == OutputInfo::Velocity || (oinfo[iOut].angularouttype == OutputInfo::total && oinfo[iOut].rescaling == false)) {
            VelSensor = true;
          }
          else {
            filePrint(stderr, " *** WARNING: unsupported probe output type. Use OUTPUT instead of OUTPUT6, or\n");
            filePrint(stderr, "     OUTPUT6 with ROTVECOUTTYPE=euler, ANGULAROUTTYPE=total, and RESCALING=off.\n");
          }
        } else {
          oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          filePrint(oinfo[iOut].filptr, "2\n");
        }
        break; 
      default:
        filePrint(stderr, " *** WARNING: Online ROM output only supports GDISPLAC, GVELOCIT, and GACCELER.\n");
        filePrint(stderr, "     output type selected is %d \n", oinfo[iOut].type);
    }
  }

  if(DispSensor || VelSensor || AccSensor) {
    geoSource->openSensorOutputFiles();
  }

  buildSensorNodeVector();
}

SDDynamPodPostProcessor::~SDDynamPodPostProcessor() {

  if(DispSensorValues) delete DispSensorValues;
  if(AccSensorValues) delete AccSensorValues;
  if(VelSensorValues) delete VelSensorValues;
}

void
SDDynamPodPostProcessor::printPODSize(int PODsize) {

  podSize = PODsize;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      default:
        break;
    }
  }
}

void
SDDynamPodPostProcessor::makeSensorBasis(VecBasis *fullBasis) {

  SensorBasis = fullBasis;
  SensorBasis->makeSparseBasis2(nodeVector, domain->getCDSA());

  // allocate space for containers to hold sensor values
  if(DispSensor) {
    DispSensorValues = new GenVector<double>(SensorBasis->vectorInfo());
    DispSensorValues->zero();
  }
  if(AccSensor) {
    AccSensorValues = new GenVector<double>(SensorBasis->vectorInfo());
    AccSensorValues->zero();
  }
  if(VelSensor) {
    VelSensorValues = new GenVector<double>(SensorBasis->vectorInfo());
    VelSensorValues->zero();
  }
}
  
void
SDDynamPodPostProcessor::buildSensorNodeVector() {

  //load vector of sensor nodes converted to local numbering
  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    if(oinfo[iOut].nodeNumber != -1) {
      nodeVector.push_back(oinfo[iOut].nodeNumber);
    }
  }

  //if multiple outputs are requested for a single node, cull redundancies from node vector
  std::sort(nodeVector.begin(), nodeVector.end());
  std::vector<int>::iterator packedNodeIt = std::unique(nodeVector.begin(), nodeVector.end());
  nodeVector.resize(packedNodeIt-nodeVector.begin());
}

void
SDDynamPodPostProcessor::printSensorValues(GenVector<double> &SensorData, OutputInfo *OINFO, double *time) {

  int locNode = OINFO->nodeNumber;
  int ndofs;
  int dofs[6];
  if(OINFO->type == OutputInfo::Disp6DOF || OINFO->type == OutputInfo::Velocity6 || OINFO->type == OutputInfo::Accel6) {
    ndofs = domain->getCDSA()->number(locNode, DofSet::XYZdisp | DofSet::XYZrot, dofs);
  }
  else {
    ndofs = domain->getCDSA()->number(locNode, DofSet::XYZdisp, dofs);
  }
  int w = OINFO->width;
  int p = OINFO->precision; 
  fprintf(OINFO->filptr, "  % *.*E  ", w, p, *time);
  for(int j = 0; j < ndofs; ++j) {
    if(dofs[j] != -1) {
      fprintf(OINFO->filptr, " % *.*E", w, p, SensorData[dofs[j]]);
    }
    else {
      fprintf(OINFO->filptr, " % *.*E", w, p, 0.0); // TODO constrained dofs
    }
  }
  fprintf(OINFO->filptr, "\n");
  fflush(OINFO->filptr);
}

void
SDDynamPodPostProcessor::dynamOutput(int tIndex, double t, DynamMat &dynOps, Vector &externalForce,
                                     Vector *AeroF, SysState<Vector>& systemState) {

  //all MPI processes have a full copy of reduced coordinates, only master processes needs to print
  int p = std::numeric_limits<double>::digits10+1;

  bool DispProjected = false, AccProjected = false, VelProjected = false;

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(tIndex % oinfo[iOut].interval == 0) {

      switch(oinfo[iOut].type) {
         case OutputInfo::Accel6 : case OutputInfo::Acceleration :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, systemState.getAccel()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else if(oinfo[iOut].type == OutputInfo::Acceleration || (oinfo[iOut].angularouttype == OutputInfo::total && oinfo[iOut].rescaling == false)) {
               if(!AccProjected) {
                 SensorBasis->expand2(systemState.getAccel(), *AccSensorValues);
                 AccProjected = true;
               }
               printSensorValues(*AccSensorValues, &oinfo[iOut], &t);
             }
           }
           break;
         case OutputInfo::Disp6DOF : case OutputInfo::Displacement : 
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, systemState.getDisp()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else if(oinfo[iOut].type == OutputInfo::Displacement || (oinfo[iOut].rotvecouttype == OutputInfo::Euler && oinfo[iOut].rescaling == false)) {
               if(!DispProjected) {
                 SensorBasis->expand2(systemState.getDisp(), *DispSensorValues);
                 DispProjected = true;
               }
               printSensorValues(*DispSensorValues, &oinfo[iOut], &t);
             }
           }
           break;
         case OutputInfo::Velocity6 : case OutputInfo::Velocity :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, systemState.getVeloc()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else if(oinfo[iOut].type == OutputInfo::Velocity || (oinfo[iOut].angularouttype == OutputInfo::total && oinfo[iOut].rescaling == false)) {
               if(!VelProjected) {
                 SensorBasis->expand2(systemState.getVeloc(), *VelSensorValues);
                 VelProjected = true;
               }
               printSensorValues(*VelSensorValues, &oinfo[iOut], &t);
             }
           }
           break;
         default:
           filePrint(stderr, " ... ROM output only supports Acceleration, Displacement, and Velocity ... \n");
      }
    }
  }
}

} /* end namespace Rom */
