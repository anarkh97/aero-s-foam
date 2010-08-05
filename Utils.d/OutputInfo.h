#ifndef _OUTPUTINFO_H_
#define _OUTPUTINFO_H_

// To add an output type
//
// 1. add the type to the enum Type
// 2. add this new type to the lexer.l file
// 3. add an output message in Driver.d/Domain.C
// 4. add a header message in Driver.d/StructProp.h
// 

#include <stdio.h>

struct OutputInfo {
   enum Type
        { Displacement, Temperature, StressXX,    StressYY,    StressZZ,
          StressXY,     StressYZ,    StressXZ,    StrainXX,    StrainYY,
          StrainZZ,     StrainXY,    StrainYZ,    StrainXZ,    HeatFlXX,
          HeatFlXY,     HeatFlXZ,    GrdTempX,    GrdTempY,    GrdTempZ,
          StressVM,     StressPR1,   StressPR2,   StressPR3,   StrainPR1,
          StrainPR2,    StrainPR3,   InXForce,    InYForce,    InZForce,
          AXMoment,     AYMoment,    AZMoment,    Energies,    AeroForce,
          EigenPair,    StrainVM,    Helmholtz,   Disp6DOF,    AeroXForce,
          AeroYForce,   AeroZForce,  AeroXMom,    AeroYMom,    AeroZMom,
          Velocity,     Acceleration,YModulus,    MDensity,    Thicknes,
          ShapeAtt,     ShapeStc,    Composit,    DispX,       DispY,
          DispZ,        RotX,        RotY,        RotZ,        DispMod,
          RotMod,       TotMod,      Rigid,       ElemToNode,  NodeToElem,
          NodeToNode,   AeroFlux,    HeatFlX,     GrdTemp,     Velocity6,
          Accel6,       ModeAlpha,   ModeError,   Reactions,   ModalDsp,
          ModalExF,     
          ContactPressure, NodalMassMat, NodalStifMat, Farfield, EigenPressure,
          FreqRespModes, HelmholtzModes,
          StressPR1Direc, StressPR2Direc, StressPR3Direc,
          StrainPR1Direc, StrainPR2Direc, StrainPR3Direc,
          EigenSlosh, SloshDispX, SloshDispY, SloshDispZ, SloshDisplacement,
          TDEnforcement, Damage, EffPStrn, HardVar,
          TemperatureFirstTimeDerivative, PressureFirstTimeDerivative, PressureSecondTimeDerivative
         };

   Type  type;
   int   interval;
   char* filename;
   FILE* filptr;
   int   averageFlg;
   int   surface;
   int   width;
   int   precision;
   int   nodeNumber;	// To output just one node's information to a file.
   int   dim;           // dimension of output data
   int   dataType;      // 1 for nodal (nodal_full or nodal_partial), 2 for elemental
   int   loadcase;      // for loadcase =-1 file is assigned to all load cases
   double ylayer,zlayer;
   int timeSliceRank; // Global rank of the corresponding time-slice (used by PITA)
   int ndtype; // 0 = deterministic, 1 = mean, 2 = stddev, 3 = pdf
   enum { realimag, modulusphase, animate };
   int complexouttype;
   int ncomplexout;   
   int tdenforc_var; // CONFACE=1, NORMAL_FORCE_MAG, NORMAL_TRACTION_MAG, TANGENTIAL_FORCE_MAG, TANGENTIAL_TRACTION_MAG,
                     // CDIRNORX, CDIRNORY, CDIRNORZ, CDIRTANX, CDIRTANY, CDIRTANZ, SLIP_MAG, NODAL_DISSIPATION,
                     // CONTACT_AREA, GAP_CUR, GAP_OLD
   void initialize() {
     width = 10; 
     precision = 4;
     interval = 0;
     filptr = 0;
     filename = 0;
     averageFlg = 1;
     surface = 2;
     nodeNumber = -1;
     ylayer = 0.0;
     zlayer = 0.0;
     timeSliceRank = 0;
     ndtype = 0;
     complexouttype = OutputInfo::realimag;
     ncomplexout = 16;
     tdenforc_var = 3;
   }

   void finalize(int numColumns) {
     if(numColumns == 6 && (type == Displacement || type == EigenPair)) type = Disp6DOF;
     if(numColumns == 6 && (type == Velocity)) type = Velocity6;
     if(numColumns == 6 && (type == Acceleration)) type = Accel6;

     if (averageFlg == 1 || averageFlg == 3)
       dataType = 1;
     else if (averageFlg == 0)
       dataType = 2;
 
     if (type == InXForce || type == InYForce ||
       type == InZForce || type == AXMoment ||
       type == AYMoment || type == AZMoment)
     dataType = 2;

     switch(type) {

       case EigenPair:
       case Displacement:
         dim = 3;
         break;
  
       case AeroForce:
         dim = 3;
         break;

       case Disp6DOF:
         dim = 6;
         break;

       case Energies:
         dim = 5;
         break;

       case Velocity:
       case Acceleration:
         dim = 3;
         break;
 
       case Velocity6:
       case Accel6:
         dim = 6;
         break;

       case StressPR1Direc:
       case StressPR2Direc:
       case StressPR3Direc:
       case StrainPR1Direc:
       case StrainPR2Direc:
       case StrainPR3Direc:
         dim = 3;
         break;

       default:
         dim = 1;
         break;
     }
   };

 
 void copyParam(const OutputInfo& oI){
  
    filptr = 0;
    
    filename       = oI.filename;
    type           = oI.type;
    interval       = oI.interval;  
    averageFlg     = oI.averageFlg;
    surface        = oI.surface;
    width          = oI.width;
    precision      = oI.precision;
    nodeNumber     = oI.nodeNumber;
    dim            = oI.dim;
    dataType       = oI.dataType;
    ylayer         = oI.ylayer;
    zlayer         = oI.zlayer;
    complexouttype = oI.complexouttype;
    ncomplexout    = oI.ncomplexout;     
    timeSliceRank = oI.timeSliceRank; 
 
 }
  
};


#endif
