#include <stdio.h>
#include <Element.d/Beam.d/EulerBeam.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <math.h>

#include <stdlib.h>
#include <Corotational.d/BeamCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>
#include <iostream>

// Define FORTRAN routines as external function
extern "C" {
 void _FORTRAN(modmstif6)(double&, double&, double*,
                          double*, double&, double&, double&,
                          double&, double*, double*, double*,
                          int&);

 void _FORTRAN(e3dmas)(double&, double*,
                       double&, double*, double*, double*,
                       double*, double*, const int&, double&, const int&);

 void _FORTRAN(mass6)(double&, double&, double&, double&, double*,
                      double&, double*, double*, double*,
                      double*, double*, const int&, double&, 
                      const int&, double*);

 void _FORTRAN(sands6)(double&, double&, const int&, double*, const int&,
                       const int&, const int&, double*, double&, double&,
                       double&, double&, double*, double*, double*, double*,
                       double&, double&, double*);

 void _FORTRAN(transform)(double*, double*, double*, double*, double*, double*, double*);
}

EulerBeam::EulerBeam(int* nodenums)
{
  //  cout <<  nodenums[0] << "," << nodenums[1] << "," <<  nodenums[2] << endl;
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
        elemframe = 0;
	offset = 0;
}

void
EulerBeam::buildFrame(CoordSet& cs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  if(nn[2] < 0) {
    if(elemframe != 0) {
      EFrame &theFrame = *elemframe;
      theFrame[0][0] = nd2.x-nd1.x;
      theFrame[0][1] = nd2.y-nd1.y;
      theFrame[0][2] = nd2.z-nd1.z;
      normalize(theFrame[0]);
      crossprod(theFrame[0],theFrame[1],theFrame[2]);
      normalize(theFrame[2]);
      crossprod(theFrame[2],theFrame[0],theFrame[1]);
    }
    else {
      fprintf(stderr," *** WARNING: No element frame exists for beam. Constructing default frame.\n",
              nn[0]+1,nn[1]+1);
      c0[0][0] = nd2.x-nd1.x;
      c0[0][1] = nd2.y-nd1.y;
      c0[0][2] = nd2.z-nd1.z;
      normalize(c0[0]);
      double N1 = sqrt( c0[0][0]*c0[0][0] + c0[0][1]*c0[0][1] );
      double N2 = sqrt( c0[0][0]*c0[0][0] + c0[0][2]*c0[0][2] );

      if (N1 > N2) {
        c0[1][0] = -c0[0][1]/N1;
        c0[1][1] = c0[0][0]/N1;
        c0[1][2] = 0.0;
      }
      else {
        c0[1][0] = c0[0][2]/N2;
        c0[1][1] = 0.0;
        c0[1][2] = -c0[0][0]/N2;
      }

      c0[2][0] = c0[0][1] * c0[1][2] - c0[0][2] * c0[1][1];
      c0[2][1] = c0[0][2] * c0[1][0] - c0[0][0] * c0[1][2];
      c0[2][2] = c0[0][0] * c0[1][1] - c0[0][1] * c0[1][0];

      elemframe = &c0;
    }
    return;
  }
  else {
    Node &nd3 = cs.getNode(nn[2]);

    elemframe = new EFrame[1];
    EFrame &theFrame = *elemframe;
    theFrame[0][0] = nd2.x-nd1.x;
    theFrame[0][1] = nd2.y-nd1.y;
    theFrame[0][2] = nd2.z-nd1.z;
    normalize(theFrame[0]);

    double xz[3] = {nd3.x-nd2.x, nd3.y-nd2.y, nd3.z-nd2.z};
    crossprod(xz,theFrame[0],theFrame[1]);
    normalize(theFrame[1]);
    crossprod(theFrame[0],theFrame[1],theFrame[2]);
    return;
  }
}

Element *
EulerBeam::clone()
{
 return new EulerBeam(*this);
}

void
EulerBeam::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

void
EulerBeam::getIntrnForce(Vector& elForce, CoordSet& cs,
			 double *elDisp, int forceIndex, double *ndTemps)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        int maxsze = 1;
        int maxstr = 7;
        int maxgus = 2;
        int elm    = 1;

        double elStress[2][7];
        //fprintf(stderr," TEMPERATURE[0] = %f\n", ndTemps[0]);
        //fprintf(stderr," TEMPERATURE[1] = %f\n", ndTemps[1]);

    _FORTRAN(sands6)(prop->A, prop->E, elm, (double*)elStress, maxsze, maxgus,
                     maxstr, (double*)*elemframe, prop->Ixx, prop->Iyy,
                     prop->Izz, prop->nu,x,y,z,elDisp, prop->W, prop->Ta, 
                     ndTemps);


// forceIndex    = 0 =  Nodal stresses along longitudinal axis of the beam
//               = 1 =  Nodal  strains along longitudinal axis of the beam 
//               = 2 =  Nodal curvatures in local x-y plane
//               = 3 =  Nodal moments in local x-y plane
//               = 4 =  Nodal curvatures in local x-z plane
//               = 5 =  Nodal moments in local x-z plane

        elForce[0] = elStress[0][forceIndex];
        elForce[1] = elStress[1][forceIndex];
}

double
EulerBeam::getMass(CoordSet& cs)
{
        double length;
    
        getLength(cs, length);

        double mass = length*(prop->rho)*(prop->A);
        return mass;
}

void
EulerBeam::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                           Vector& gravityForce, int gravflg, GeomState *geomState)
{
        double massPerNode = 0.5*getMass(cs);
	
        double t0n[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
	double length;
	  
	if (geomState) {
	    updTransMatrix(cs, geomState, t0n, length);
		 
        }  else  {
           getLength(cs, length);
	     
           for(int i=0; i<3; ++i) {   
              for(int j=0; j<3; ++j) {
                  t0n[i][j] = (*elemframe)[i][j] ;
             
	      }
           }
        }         

        int i;
        double localg[3];

        for(i=0; i<3; ++i)
          localg[i] = 0.0;

        for(i=0; i<3; ++i) {
          localg[0] += t0n[0][i]*gravityAcceleration[i];
          localg[1] += t0n[1][i]*gravityAcceleration[i];
          localg[2] += t0n[2][i]*gravityAcceleration[i];
        }
        double localf[3], localm[3];
        double globalf[3], globalm[3];
        localf[0] =  massPerNode*localg[0];
        localf[1] =  massPerNode*localg[1];
        localf[2] =  massPerNode*localg[2];
        if (gravflg == 2) { // consistent
          localm[0] =  0.0;
          localm[1] = -massPerNode*localg[2]*length/6.0;
          localm[2] =  massPerNode*localg[1]*length/6.0;
        }
        else if (gravflg == 1) { // lumped with fixed-end moments
          localm[0] =  0.0;
          localm[1] = -massPerNode*localg[2]*length/8.0;
          localm[2] =  massPerNode*localg[1]*length/8.0;
        }
        else localm[0] = localm[1] = localm[2] = 0.0; // lumped without fixed-end moments

        for(i=0; i<3; ++i) {
          globalf[i] = (t0n[0][i]*localf[0]) + (t0n[1][i]*localf[1]) + (t0n[2][i]*localf[2]);
          globalm[i] = (t0n[1][i]*localm[1]) + (t0n[2][i]*localm[2]);
        }

	gravityForce[0]  =  globalf[0];
	gravityForce[1]  =  globalf[1];
	gravityForce[2]  =  globalf[2];
	gravityForce[3]  =  globalm[0];
	gravityForce[4]  =  globalm[1];
	gravityForce[5]  =  globalm[2];
        gravityForce[6]  =  globalf[0];
	gravityForce[7]  =  globalf[1];
	gravityForce[8]  =  globalf[2];
	gravityForce[9]  = -globalm[0];
	gravityForce[10] = -globalm[1];
	gravityForce[11] = -globalm[2];

}

FullSquareMatrix
EulerBeam::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
       // Check for phantom element, which has no stiffness
        if(prop == NULL) {
           FullSquareMatrix ret(12,mel);
           ret.zero();
           return ret;
        }

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double *gravityAcceleration = 0, *grvfor = 0, totmas = 0.0;

        int grvflg = 0, masflg = 0;

// Lumped Mass Matrix
        if (cmflg == 0) {
          _FORTRAN(e3dmas)(prop->rho,(double*)mel,prop->A,
	          x,y,z,gravityAcceleration,grvfor,grvflg,totmas,masflg);

// Consistent (Full) Mass Matrix
        } else {
          _FORTRAN(mass6)(prop->rho,prop->Ixx, prop->Iyy, prop->Izz,
                  (double*)mel,prop->A,x,y,z,gravityAcceleration,
                   grvfor,grvflg,totmas,masflg,(double *) *elemframe);
        }

        FullSquareMatrix ret(12, mel);

        //fprintf(stderr,"mass matrix BEFORE offset\n");
        //ret.print();
	
	if(offset) offsetAxis(ret);

        //fprintf(stderr,"mass matrix AFTER offset\n");
        //ret.print();

        return ret;
}

FullSquareMatrix
EulerBeam::stiffness(CoordSet &cs, double *d, int flg)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

//      EFrame is stored in following order.
//                      X_x, X_y, X_z,
//                      Y_x, Y_y, Y_z,
//                      Z_x, Z_y, Z_z


	// Check for zero area
        if(prop->A <= 0.0) {
          fprintf(stderr," *** WARNING: Euler beam has zero area. nodes %d %d\n",
          nn[0]+1, nn[1]+1);
        }

        if(sqrt((x[0]-x[1])*(x[0]-x[1])+(y[0]-y[1])*(y[0]-y[1])+(z[0]-z[1])*(z[0]-z[1]))==0.0)
          fprintf(stderr," *** WARNING: Euler beam has zero lenth. nodes %d %d\n",
                 nn[0]+1,nn[1]+1);

	// Check for zero young's modulus
        if(prop->E <= 0.0) {
          fprintf(stderr," *** WARNING: Euler beam element %d has zero Young's Modulus. nodes %d %d\n", glNum, nn[0]+1, nn[1]+1);
          fprintf(stderr,"E = %e  A = %e  nu = %e  rho = %e\n", prop->E, prop->A, prop->nu, prop->rho);
        }

       // Check for the frame
       if(elemframe == 0) {
         fprintf(stderr," ************************************************\n");
         fprintf(stderr," *** ERROR: Euler beam %d lacks a frame (Nodes %d %d)\n",                          glNum, nn[0]+1,nn[1]+1);
         fprintf(stderr," ************************************************\n");
         exit(-1);
       }
       //cerr << "elemframe = " << (*elemframe)[0][0] << " " << (*elemframe)[0][1] << " " << (*elemframe)[0][2] << "   "
       //                       << (*elemframe)[1][0] << " " << (*elemframe)[1][1] << " " << (*elemframe)[1][2] << "   "
       //                       << (*elemframe)[2][0] << " " << (*elemframe)[2][1] << " " << (*elemframe)[2][2] << endl;

	_FORTRAN(modmstif6)(prop->A, prop->E, (double *)d,
		             (double *) *elemframe, 
			     prop->Ixx, prop->Iyy, prop->Izz,
			     prop->nu, x, y, z, flg);

	FullSquareMatrix ret(12,d);

	if(offset) offsetAxis(ret);

        return ret;
}

int
EulerBeam::numNodes()
{
        return 2;
}

int *
EulerBeam::nodes(int *p)
{
        if(p == 0) p = new int[2];
        p[0] = nn[0];
        p[1] = nn[1];
        return p;
}

int
EulerBeam::numDofs()
{
        return 12;
}

int *
EulerBeam::dofs(DofSetArray &dsa, int *p)
{
        if(p == 0) p = new int[12];

	dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);

        return p;
}

void
EulerBeam::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
}

Corotator *
EulerBeam::getCorotator(CoordSet &cs, double *kel, int , int fitAlgBeam)
{
  int flag = 0; 
  FullSquareMatrix myStiff = stiffness(cs, kel, flag);
  return new BeamCorotator(nn[0], nn[1], (*elemframe)[2], myStiff, fitAlgBeam);
}

int
EulerBeam::getTopNumber()
{
  return 106;
}

void
EulerBeam::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                GeomState *geomState, int cflg)
{
  double normal[3], normal2[3];
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double length;
  
  if(geomState) {
    double t0n[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    updTransMatrix(cs, geomState, t0n, length);
    normal[0] = t0n[1][0];
    normal[1] = t0n[1][1];
    normal[2] = t0n[1][2];
    normal2[0] = t0n[2][0];
    normal2[1] = t0n[2][1];
    normal2[2] = t0n[2][2];
  } 
  else {
     // Obtain normal to beam (second vector in element frame)
     normal[0] = (*elemframe)[1][0];
     normal[1] = (*elemframe)[1][1];
     normal[2] = (*elemframe)[1][2];
     normal2[0] = (*elemframe)[2][0];
     normal2[1] = (*elemframe)[2][1];
     normal2[2] = (*elemframe)[2][2]; 
     getLength(cs, length);
  } 
  double pressureForce = 0.5*pressure*length;
  px = pressureForce*normal[0];
  py = pressureForce*normal[1];
  pz = pressureForce*normal[2]; 
  
  // Consistent
  double localMz = pressureForce*length/6.0;
  double mx = localMz*normal2[0];
  double my = localMz*normal2[1];
  double mz = localMz*normal2[2]; 
  elPressureForce[0]  = px;
  elPressureForce[1]  = py;
  elPressureForce[2]  = pz;
  elPressureForce[3]  = mx;
  elPressureForce[4]  = my;
  elPressureForce[5]  = mz;
  elPressureForce[6]  = px;
  elPressureForce[7]  = py;
  elPressureForce[8]  = pz;
  elPressureForce[9]  = -mx;
  elPressureForce[10] = -my;
  elPressureForce[11] = -mz;
}

void
EulerBeam::getThermalForce(CoordSet &cs, Vector &ndTemps,
                           Vector &elementThermalForce, int glflag, GeomState *geomState)
{
// Called from Dynam.C
// Computes the thermal-mechanical coupling force C*theta
// A = cross section, W= dilatation coeff

    double localThF[12];
    int i, j;

    double Tref  = prop->Ta;
    double coeff = prop->E*prop->W*prop->A;
    double length;

    // Local Thermal Forces : There are only axial forces (see notes)
    // indices 0-5 --> node1, 6-11 -->node2

    double deltaT1 = ndTemps[0]-Tref;
    double deltaT2 = ndTemps[1]-Tref;

    //fprintf(stderr," deltaT1 = %f\n", deltaT1);
    //fprintf(stderr," deltaT2 = %f\n", deltaT2);

    for(i=0; i<12; i++) {localThF[i] = 0. ;}

    localThF[0] = -coeff * (0.5 * deltaT1 + 0.5 * deltaT2);
    localThF[6] =  coeff * (0.5 * deltaT1 + 0.5 * deltaT2);

    // Compute Global Thermal Forces:: From local to global -->
    //                                transpose(Transform. Matrix)

   double t0n[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

   if (geomState) {
//   fprintf(stderr," *** UPDATING EFRAMES for Non-Linear [T forces]\n");
    updTransMatrix(cs, geomState, t0n, length);
   }  else  {
//    fprintf(stderr," *** USING initial FRAME\n");
    for(i=0; i<3; ++i) {   
     for(j=0; j<3; ++j) {
      t0n[i][j] = (*elemframe)[i][j] ;
     }
    }
   }

//Node 1:
    for(i=0; i<3; ++i) {
     elementThermalForce[i] = 0.;
     for(j=0; j<3; ++j) {
      elementThermalForce[i] += t0n[j][i]*localThF[j];
//       fprintf(stderr,"Thermal Frame = %f\n", t0n[j][i]);
     }
    }
    for(i=0; i<3; ++i) {
     elementThermalForce[i+3] = 0.;
     for(j=0; j<3; ++j) {
      elementThermalForce[i+3] += t0n[j][i]*localThF[j+3];
     }
    }
//Node 2:
    for(i=0; i<3; ++i) {
     elementThermalForce[i+6] = 0.;
     for(j=0; j<3; ++j) {
      elementThermalForce[i+6] += t0n[j][i]*localThF[j+6];
     }
    }
 
    for(i=0; i<3; ++i) {
     elementThermalForce[i+9] = 0.;
     for(j=0; j<3; ++j) {
      elementThermalForce[i+9] += t0n[j][i]*localThF[j+9];
     }
    }                                

//     ndTemps.print();
//     elementThermalForce.print("Element Thermal Force");
}

void
EulerBeam::computeDisp(CoordSet&cs,
      State &state, const InterpPoint &ip, double*res, GeomState *geomState)
{
//REM: No Rotation needs to be sent to fluid2D
//12=3translations, 3rotations + their respect. velocities
//2=number of nodes in element

    double xyz[2][12];
    double locxyz[2][12];
    double locres[6];

    state.getDVRot(nn[0], xyz[0], xyz[0]+6);
    state.getDVRot(nn[1], xyz[1], xyz[1]+6);

    double t0n[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    double L;
    int i, j;

    if (geomState) {
      updTransMatrix(cs, geomState, t0n, L);
    }  else  {
                getLength(cs, L);
                for(i=0; i<3; ++i) {
                  for(j=0; j<3; ++j) {
                   t0n[i][j] = (*elemframe)[i][j] ;
                  }
                }
    }


// Rotate displacements from global frame to local frame

//Node 1 (displacements):
  for(i=0; i<3; ++i) {
   locxyz[0][i] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[0][i] += t0n[i][j]*xyz[0][j];
   }
  }

  for(i=0; i<3; ++i) {
   locxyz[0][i+3] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[0][i+3] += t0n[i][j]*xyz[0][j+3];
   }
  }

//Node 1 (velocities):

  for(i=0; i<3; ++i) {
   locxyz[0][i+6] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[0][i+6] += t0n[i][j]*xyz[0][j+6];
   }
  }

  for(i=0; i<3; ++i) {
   locxyz[0][i+9] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[0][i+9] += t0n[i][j]*xyz[0][j+9];
   }
  }

//Node 2 (displacements):

  for(i=0; i<3; ++i) {
   locxyz[1][i] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[1][i] += t0n[i][j]*xyz[1][j];
   }
  }

  for(i=0; i<3; ++i) {
   locxyz[1][i+3] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[1][i+3] += t0n[i][j]*xyz[1][j+3];
   }
  }

//Node 2 (velocities):

  for(i=0; i<3; ++i) {
   locxyz[1][i+6] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[1][i+6] += t0n[i][j]*xyz[1][j+6];
   }
  }

  for(i=0; i<3; ++i) {
   locxyz[1][i+9] = 0.;
   for(j=0; j<3; ++j) {
    locxyz[1][i+9] += t0n[i][j]*xyz[1][j+9];
   }
  }

  double locGap[3];
  // Rotate the gap; note the Z component should be zero
  for(i=0; i<3; ++i) {
    locGap[i] = 0;
    for(j=0; j<3; ++j)
      locGap[i] += t0n[i][j]*ip.gap[j];
  }

// Multiplying by beam shape functions

  const double *gp = ip.xy;
  double gp2 = gp[0]*gp[0] ;
  double gp3 = gp2*gp[0] ;
  double Nwi = 1-3*gp2+2*gp3;
  double Nthi= L*(gp[0]-2*gp2+gp3);
  double Nwj = 3*gp2-2*gp3;
  double Nthj= L*(-gp2+gp3);
  double theta = gp[0]*locxyz[1][3] + (1-gp[0])*locxyz[0][3];
  double thetadot = gp[0]*locxyz[1][9] + (1-gp[0])*locxyz[0][9];

//x-direction
 locres[0] = (1-gp[0])*locxyz[0][0] + gp[0]*locxyz[1][0];

//y-direction
 locres[1] = Nwi*locxyz[0][1] + Nthi*locxyz[0][5] +
             Nwj*locxyz[1][1] + Nthj*locxyz[1][5]
             - theta*locGap[2] ;

//z-direction (BE CAREFUL OF MOMENT SIGN)
 locres[2] = Nwi*locxyz[0][2] - Nthi*locxyz[0][4] +
             Nwj*locxyz[1][2] - Nthj*locxyz[1][4]
             + theta*locGap[1] ;

//x-velocity
  locres[3] = (1-gp[0])*locxyz[0][6] + gp[0]*locxyz[1][6];

//y-velocity
  locres[4] = Nwi*locxyz[0][7] + Nthi*locxyz[0][11] +
              Nwj*locxyz[1][7] + Nthj*locxyz[1][11]
             - thetadot*locGap[2] ;

//z-velocity
  locres[5] = Nwi*locxyz[0][8] + Nthi*locxyz[0][10] +
              Nwj*locxyz[1][8] + Nthj*locxyz[1][10]
             + thetadot*locGap[1] ;

// Transform back to global (Transpose of matrix)

  for(i=0; i<3; ++i) {
   res[i] = 0.;
   for(j=0; j<3; ++j) {
    res[i] += t0n[j][i]*locres[j];
//       fprintf(stderr,"Displacement Frame = %f\n", t0n[j][i]);
   }
  }

  for(i=0; i<3; ++i) {
   res[i+3] = 0.;
   for(j=0; j<3; ++j) {
    res[i+3] += t0n[j][i]*locres[j+3];
   }
  }

/* int j;
 for(j=3; j<6; ++j)
    res[j] = (1-gp[0])*xyz[0][j] + gp[0]*xyz[1][j]; */

/* double xyz[2][6];
 state.getDV(nn[0], xyz[0], xyz[0]+3);
 state.getDV(nn[1], xyz[1], xyz[1]+3);

 int j;
 for(j=0; j<6; ++j)
    res[j] = (1-gp[0])*xyz[0][j] + gp[0]*xyz[1][j]; */

}

void
EulerBeam::getFlLoad(CoordSet&cs, const InterpPoint &ip, double *flF,
                     double *resF, GeomState *geomState)
{
// fprintf(stderr," Gauss Points: %f %f\n", gp[0], gp[1]);
// 2DCode sends Fx, Fy, Fz

// Transform forces from global to local

  double locload[3];
  double locresF[12];
  int i, j;

  double t0n[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  double L;

  if (geomState) {
    updTransMatrix(cs, geomState, t0n, L);
   }  else  {
              getLength(cs, L);
              for(i=0; i<3; ++i) {
                for(j=0; j<3; ++j) {
                  t0n[i][j] = (*elemframe)[i][j] ;
                }
              }
   }


   for(i=0; i<3; ++i) {
    locload[i] = 0.0;
     for(j=0; j<3; j++) {
      locload[i] += t0n[i][j]*flF[j];
//      fprintf(stderr,"Send Forces %f ", flF[j]);
     }
   }

  double locGap[3];
  // Rotate the gap; note the Z component should be zero
  for(i=0; i<3; ++i) {
    locGap[i] = 0;
    for(j=0; j<3; ++j)
      locGap[i] += t0n[i][j]*ip.gap[j];
  }


// Multiply by shape functions

    const double *gp = ip.xy;
    double gp2 = gp[0]*gp[0] ;
    double gp3 = gp2*gp[0] ;
    double Nwi = 1-3*gp2+2*gp3;
    double Nthi= L*(gp[0]-2*gp2+gp3);
    double Nwj = 3*gp2-2*gp3;
    double Nthj= L*(-gp2+gp3);

    for(i=0; i<3; ++i) {
     locresF[i]   = Nwi*locload[i];  //Node 1: Fx, Fy, Fz
     locresF[i+3] = Nthi*locload[i]; //        Mx, My, Mz
     locresF[i+6] = Nwj*locload[i];  //Node 2
     locresF[i+9] = Nthj*locload[i];
//     fprintf(stderr,"Moments %f %f\n ",locresF[i+3], locresF[i+9]);
//     fprintf(stderr,"Forces %f %f\n ",locresF[i], locresF[i+6]);
    }

   double mLx = locGap[1] * locload[2] - locGap[2] * locload[1];
   //fprintf(stderr, "Moment in local x: %e will go to %e %e %e\n", mLx,
          //t0n[0][0], t0n[0][1], t0n[0][2]);
   locresF[3] += mLx*(1-ip.xy[0]);
   locresF[9] += mLx*(ip.xy[0]);
// Transform back to global load

//Node 1:
    for(i=0; i<3; ++i) {
     resF[i] = 0.;
     for(j=0; j<3; ++j) {
      resF[i] += t0n[j][i]*locresF[j];
//      fprintf(stderr,"Load Frame = %f\n", t0n[j][i]);
     }
//     fprintf(stderr,"Forces 1: %f", resF[i]);
    }
     for(i=0; i<3; ++i) {
     resF[i+3] = 0.;
     for(j=0; j<3; ++j) {
      resF[i+3] += t0n[j][i]*locresF[j+3];
     }
    }
//Node 2:
    for(i=0; i<3; ++i) {
     resF[i+6] = 0.;
     for(j=0; j<3; ++j) {
      resF[i+6] += t0n[j][i]*locresF[j+6];
     }
    }

     for(i=0; i<3; ++i) {
     resF[i+9] = 0.;
     for(j=0; j<3; ++j) {
      resF[i+9] += t0n[j][i]*locresF[j+9];
     }
    }
}

void
EulerBeam::updTransMatrix(CoordSet& cs, GeomState *geomState, double t0n[3][3], double &length)
{
// Returns t0n[3][3] and length

   double  xn[2][3];

   double  zVecL[2][3];
   double (* rot[2])[3][3];

   // Get Nodes current coordinates
   NodeState &ns1 = (*geomState)[nn[0]];
   NodeState &ns2 = (*geomState)[nn[1]];

   xn[0][0]  = ns1.x; // x coordinate of node state 1
   xn[0][1]  = ns1.y; // y coordinate of node state 1
   xn[0][2]  = ns1.z; // z coordinate of node state 1

   xn[1][0]  = ns2.x; // x coordinate of node state 2
   xn[1][1]  = ns2.y; // y coordinate of node state 2
   xn[1][2]  = ns2.z; // z coordinate of node state 2

   double dx = xn[1][0] - xn[0][0];
   double dy = xn[1][1] - xn[0][1];
   double dz = xn[1][2] - xn[0][2];

   length = sqrt(dx*dx + dy*dy + dz*dz);

   rot[0]    = &(ns1.R); // rotation tensor of node state 1
   rot[1]    = &(ns2.R); // rotation tensor of node state 2

// Compute nodal rotated Z-axis in global coordinate system

   int i, nod;
   for(nod=0; nod<2; ++nod ) {
      for(i=0; i<3; ++i ) {
         zVecL[nod][i] = (*rot[nod])[i][0]*(*elemframe)[2][0]
                        +(*rot[nod])[i][1]*(*elemframe)[2][1]
                        +(*rot[nod])[i][2]*(*elemframe)[2][2];
      }
   }

/* Fitalg 1: Z-axis from node 1 */
   // We are setting fit Alg. to 2 to average z vectors.
   int fitAlg = 2;

   if (fitAlg == 1) {
      t0n[2][0] = zVecL[0][0];
      t0n[2][1] = zVecL[0][1];
      t0n[2][2] = zVecL[0][2];
   }

/* Fitalg .ne. 1: Z-axis as sum of nodal z-axis */
   else {
      t0n[2][0] = zVecL[0][0] + zVecL[1][0];
      t0n[2][1] = zVecL[0][1] + zVecL[1][1];
      t0n[2][2] = zVecL[0][2] + zVecL[1][2];
   }


   t0n[0][0]  = xn[1][0] - xn[0][0];
   t0n[0][1]  = xn[1][1] - xn[0][1];
   t0n[0][2]  = xn[1][2] - xn[0][2];

/* X-axis along element in Cn */
   normalize( t0n[0] );

/* Y-axis as cross product between z and x */
    crossprod( t0n[2], t0n[0], t0n[1]);
    normalize( t0n[1] );

/* Z-axis as cross product between x and y */
    crossprod( t0n[0], t0n[1], t0n[2]);
    normalize( t0n[2] );
}

void
EulerBeam::getLength(CoordSet& cs, double &length)
{
// Returns length of element

     Node &nd1 = cs.getNode(nn[0]);
     Node &nd2 = cs.getNode(nn[1]);

     double x[2], y[2], z[2];

     x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
     x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

     double dx = x[1] - x[0];
     double dy = y[1] - y[0];
     double dz = z[1] - z[0];

     length = sqrt(dx*dx + dy*dy + dz*dz);

}

void EulerBeam::offsetAxis(FullSquareMatrix& mat)
{
  double cross[3][3] = { 
    { 0.0, -offset[2], offset[1] },
    { offset[2], 0.0, -offset[0] },
    { -offset[1], offset[0], 0.0 }
  };
  double tmp[12][12];
  int i, j;

  for(i = 0; i < 12; ++i){
    for(j = 0; j < 12; ++j){
      tmp[i][j] = mat[i][j];
    }
    for(j = 0; j < 3; ++j){
      tmp[i][j+3] -= mat[i][0]*cross[0][j] + mat[i][1]*cross[1][j]
                     + mat[i][2]*cross[2][j];
      tmp[i][j+9] -= mat[i][6]*cross[0][j] + mat[i][7]*cross[1][j]
                     + mat[i][8]*cross[2][j];
    }
  }

  for(j = 0; j < 12; ++j){
    for(i = 0; i < 12; ++i){
      mat[i][j] = tmp[i][j];
    }
    for(i = 0; i < 3; ++i){
      mat[i+3][j] -= cross[0][i]*tmp[0][j] + cross[1][i]*tmp[1][j]
                     + cross[2][i]*tmp[2][j];
      mat[i+9][j] -= cross[0][i]*tmp[6][j] + cross[1][i]*tmp[7][j]
                     + cross[2][i]*tmp[8][j];
    }
  }
}


void
EulerBeam::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                       Vector& elDisp, int strInd, int surface, 
                       double *ndTemps, double ylayer, double zlayer, int avgnum)
{
   // Calculates the axial strain and stress for the beam element.

   weight = 1.0;
	
   Node &nd1 = cs.getNode(nn[0]);
   Node &nd2 = cs.getNode(nn[1]);

   double x[2], y[2], z[2];

   x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
   x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

   // Find the internal forces and moments

   int maxsze = 1;
   int maxstr = 7;
   int maxgus = 2;
   int elm    = 1;  
	
   double elStress[2][7];
   double elForce[3][2]={{0.0,0.0},{0.0,0.0},{0.0,0.0}};


  _FORTRAN(sands6)(prop->A, prop->E, elm, (double*)elStress, maxsze, maxgus,
                   maxstr, (double*)*elemframe, prop->Ixx, prop->Iyy,
                   prop->Izz, prop->nu,x,y,z,elDisp.data(), prop->W, prop->Ta, 
                   ndTemps);

   // elForce[0] -> Axial Force (x-direction)
   // elForce[1] -> Moment around the y-axis (My)
   // elForce[2] -> Moment around the z-axis (Mz) 

   elForce[0][0] = -elStress[0][0];
   elForce[1][0] = -elStress[0][4];
   elForce[2][0] = -elStress[0][5];
   elForce[0][1] =  elStress[1][0];
   elForce[1][1] =  elStress[1][4];
   elForce[2][1] =  elStress[1][5];

   // Y is the position in the cross section where the stress is to be 
   // calculated and is defined as Y = -ylayer * ymin if ylayer is < 0.0 (ymin 
   // is negative in the input file) and Y = ylayer * ymax if ylayer is > 0.0 
   // (ymax is positive), where ylayer is the percentage of the position ymin 
   // or ymax.
   // Same thing for Z

   double Y = 0.0;
   double Z = 0.0;
	
   if (ylayer < 0.0) {
      Y = -ylayer*prop->ymin;
   } else if (ylayer > 0.0) {
      Y =  ylayer*prop->ymax;
   }

   if (zlayer < 0.0) {
      Z = -zlayer*prop->zmin;
   } else if (zlayer > 0.0) {
      Z =  zlayer*prop->zmax;
   }
	

   switch (avgnum) {
    
      case 0:
      {
        if (strInd == 0) {
	
	  // Axial Stress
	
 	  double IY = prop->Iyy;
          double IZ = prop->Izz;
          double cA = prop->A;
	      
	  stress[0] = elForce[0][0]/cA - elForce[2][0]*Y/IZ + elForce[1][0]*Z/IY;
  	  stress[1] = elForce[0][1]/cA - elForce[2][1]*Y/IZ + elForce[1][1]*Z/IY;
	
	} else if (strInd == 7) {
	
	  // Axial Strain
        
	  double EA  = prop->E*prop->A;
	  double EIZ = prop->E*prop->Izz;
	  double EIY = prop->E*prop->Iyy;
	
	  double Tref  = prop->Ta;
	  double alpha = prop->W;
	
	  double dT1 = ndTemps[0]-Tref;
	  double dT2 = ndTemps[1]-Tref;
	
	  double localThS;
          localThS = alpha * (0.5 * dT1 + 0.5 * dT2);
	
          stress[0] = elForce[0][0]/EA - elForce[2][0]*Y/EIZ + elForce[1][0]*Z/EIY + localThS;
  	  stress[1] = elForce[0][1]/EA - elForce[2][1]*Y/EIZ + elForce[1][1]*Z/EIY + localThS;
	 
	} else {
	  stress[0] = 0.0;
	  stress[1] = 0.0;
	}	
        break;
      }

      case 1:
      {
       if (strInd < 6) {

          // Axial Stress

          double IY = prop->Iyy;
          double IZ = prop->Izz;
          double cA = prop->A;

          double xg[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
          double tmpStr1[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
          double tmpStr2[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

          tmpStr1[0] = elForce[0][0]/cA - elForce[2][0]*Y/IZ + elForce[1][0]*Z/IY;
          tmpStr2[0] = elForce[0][1]/cA - elForce[2][1]*Y/IZ + elForce[1][0]*Z/IY;

         _FORTRAN(transform)((double*)*elemframe[0],(double*)*elemframe[1],(double*)*elemframe[2],
                              xg[0], xg[1], xg[2], tmpStr1);

         _FORTRAN(transform)((double*)*elemframe[0],(double*)*elemframe[1],(double*)*elemframe[2],
                              xg[0], xg[1], xg[2], tmpStr2);

          stress[0] = tmpStr1[strInd];
          stress[1] = tmpStr2[strInd];

        } else if (strInd > 6 && strInd < 13) {

          // Axial Strain

          double EA  = prop->E*prop->A;
          double EIZ = prop->E*prop->Izz;
          double EIY = prop->E*prop->Iyy;

          double Tref  = prop->Ta;
          double alpha = prop->W;

          double dT1 = ndTemps[0]-Tref;
          double dT2 = ndTemps[1]-Tref;

          double localThS;
          localThS = alpha * (0.5 * dT1 + 0.5 * dT2);

          double xg[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
          double tmpStr1[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
          double tmpStr2[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

          tmpStr1[0] = elForce[0][0]/EA - elForce[2][0]*Y/EIZ + elForce[1][0]*Z/EIY + localThS;
          tmpStr2[0] = elForce[0][1]/EA - elForce[2][1]*Y/EIZ + elForce[1][1]*Z/EIY + localThS;

         _FORTRAN(transform)((double*)*elemframe[0],(double*)*elemframe[1],(double*)*elemframe[2],
                              xg[0], xg[1], xg[2], tmpStr1);

         _FORTRAN(transform)((double*)*elemframe[0],(double*)*elemframe[1],(double*)*elemframe[2],
                              xg[0], xg[1], xg[2], tmpStr2);

          stress[0] = tmpStr1[strInd-7];
          stress[1] = tmpStr2[strInd-7];

        }
        else {
          stress[0] = 0.0;
          stress[1] = 0.0;
        }
        break;
      }

      case 2:
      {
        weight = 0.0;
        stress[0] = 0.0;
        stress[1] = 0.0;
        break;
      }

      default:
        cerr << "avgnum = " << avgnum << " is not a valid number\n";
    }
}

/* this is for a 2D beam!!!
void
EulerBeam::subtractLoadStiffness(GeomState &geomState, CoordSet &cs, FullSquareMatrix &elK)
{
   // PJSA 1/15/2010 subtract load stiffness due to pressure
 double dx = ns2.x - ns1.x;
 double dy = ns2.y - ns1.y;
 double dz = ns2.z - ns1.z;
 double length = sqrt(dx*dx + dy*dy + dz*dz);
 
 double p = pressure/2, q = pressure*length/12;
 elK[0][1] -= p;
 elK[0][2] -= q; 
 elK[0][4] += p;
 elK[0][5] += q;
 elK[1][0] += p;
 elK[1][3] -= p;
 elK[2][0] -= q;
 elK[2][3] += q;
 elK[3][1] -= p;
 elK[3][2] += q;
 elK[3][4] += p;
 elK[3][5] -= q;
 elK[4][0] += p;
 elK[4][3] -= p;
 elK[5][0] += q;
 elK[5][3] -= q;
}
*/
