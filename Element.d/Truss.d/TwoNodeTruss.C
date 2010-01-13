#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <Element.d/Truss.d/TwoNodeTruss.h>
#include <Corotational.d/BarCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Corotational.d/utilities.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

extern "C" {
 void _FORTRAN(transform)(double*, double*, double*, double*, double*, double*, double*);
}

TwoNodeTruss::TwoNodeTruss(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        preload = 0.0;
}

void 
TwoNodeTruss::setPreLoad(double load, int &flg)
{ 
        preload = load;
        flg = 1;
}

Element *
TwoNodeTruss::clone()
{
	return new TwoNodeTruss(*this);
}

void
TwoNodeTruss::renum(int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
TwoNodeTruss::getIntrnForce(Vector& elForce, CoordSet& cs,
			    double *elDisp, int forceIndex, double *ndTemps)
{
        // ... BARS ONLY CARRY AXIAL FORCES
        if(forceIndex > 0) {
          elForce[0] = 0.0;
          elForce[1] = 0.0;
          return;
        }

        Node &nd1 = cs.getNode( nn[0] );
        Node &nd2 = cs.getNode( nn[1] );

	double dx = nd2.x - nd1.x;
	double dy = nd2.y - nd1.y;
	double dz = nd2.z - nd1.z;

	double length = sqrt(dx*dx + dy*dy + dz*dz);

	// scale dx, dy, and dz by the length
	dx /= length;
	dy /= length;
	dz /= length;

	// Compute the change in length of the element
	double dq = dx*(elDisp[3]-elDisp[0]) 
                  + dy*(elDisp[4]-elDisp[1]) 
                  + dz*(elDisp[5]-elDisp[2]);

	// Compute axial strain
	double exx = dq/length;

	// Compute force
	double f = prop->A*prop->E*exx;

        // Add Preload
        f  += preload;

        // Compute thermal force
        double coefficient = prop->A*prop->E*prop->W;
        double Tref = prop->Ta;

        double fth1 = coefficient*(ndTemps[0]-Tref);
        double fth2 = coefficient*(ndTemps[1]-Tref);

        // If We take the mean temperature:
//      double fth = coefficient*((ndTemps[0]+ndTemps[1])/2-Tref);

	// return nodal forces
        elForce[0] = -f + fth1;
        elForce[1] =  f - fth2;
}

double
TwoNodeTruss::getMass(CoordSet& cs)
{
        Node &nd1 = cs.getNode( nn[0] );
        Node &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );
        double   mass = (length*prop->A*prop->rho);

        return mass;
}

void
TwoNodeTruss::getGravityForce(CoordSet& cs, double *gravityAcceleration, Vector &gravityForce,
                              int gravflg, GeomState *geomState)
{
        double fx,fy,fz;
        double massPerNode = 0.5*getMass(cs);

        fx = massPerNode*gravityAcceleration[0];
        fy = massPerNode*gravityAcceleration[1];
        fz = massPerNode*gravityAcceleration[2];

        gravityForce[0] = fx;
        gravityForce[1] = fy;
        gravityForce[2] = fz;
        gravityForce[3] = fx;
        gravityForce[4] = fy;
        gravityForce[5] = fz;

}

FullSquareMatrix
TwoNodeTruss::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
        double mass = getMass(cs);
        double massPerNode = 0.5*mass;

        FullSquareMatrix elementMassMatrix(6,mel);

// zero the element mass matrix
	elementMassMatrix.zero();

// set the diagonal elements 
	int i;
        for(i=0; i<6; ++i)
           elementMassMatrix[i][i] = massPerNode;

//        fprintf(stderr," ... Element Mass Matrix for Truss ...\n");
//        elementMassMatrix.print();
        return elementMassMatrix;
}

FullSquareMatrix
TwoNodeTruss::stiffness(CoordSet &cs, double *k, int flg)
{
        Node &nd1 = cs.getNode( nn[0] );
        Node &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt( dx*dx + dy*dy + dz*dz );

        if(length == 0.0) {
          fprintf(stderr,"ERROR: truss has zero length %d %d\n",nn[0],nn[1]);
          fprintf(stderr," ... exiting fem program ...\n");
          exit(-1);
        }

        double c1[3];
        
        c1[0] = dx/length;
        c1[1] = dy/length;
        c1[2] = dz/length;

	FullSquareMatrix ret(6,k);

        // Check for negative or zero area and zero modulus
        if(prop->A <= 0.0)
          fprintf(stderr,"ERROR: truss has zero area %d %d\n",nn[0],nn[1]);
        if(prop->E <= 0.0)
          fprintf(stderr,"ERROR: truss has zero modulus %d %d\n",nn[0],nn[1]);

        double elementK = prop->E*prop->A/length;

//        fprintf(stderr,"element K = %g \n", elementK);

	int i,j;
        for(i=0; i < 3; ++i)
          for(j=0; j < 3; ++j) {
             ret[i][j]     = elementK*c1[i]*c1[j];
             ret[i+3][j+3] = elementK*c1[i]*c1[j];
             ret[i+3][j]   = -ret[i][j];
             ret[i][j+3]   = -ret[i][j];
          }

        if (preload != 0.0)  {
//        fprintf(stderr," ... Adding Preload for Truss ...\n");
          for(i=0; i < 3; ++i) {
             ret[i][i]     += preload/length;
             ret[i+3][i+3] += preload/length;
             ret[i+3][i]   = -ret[i][i];
             ret[i][i+3]   = -ret[i][i];
          }
        }
//        fprintf(stderr," ... Element Stiffness Matrix for Truss ...\n");
//        ret.print();
        return ret;
}

int
TwoNodeTruss::numNodes()
{
        return 2;
}

int *
TwoNodeTruss::nodes(int *p)
{
        if(p == 0) p = new int[2];
        p[0] = nn[0];
        p[1] = nn[1];
        return p;
}

int
TwoNodeTruss::numDofs()
{
        return 6;
}

int *
TwoNodeTruss::dofs(DofSetArray &dsa, int *p)
{
        if(p == 0) p = new int[6];

        dsa.number(nn[0],DofSet::XYZdisp, p  );
        dsa.number(nn[1],DofSet::XYZdisp, p+3);

        return p;
}

void
TwoNodeTruss::markDofs( DofSetArray &dsa )
{
        dsa.mark(nn, 2, DofSet::XYZdisp);
}

// For Corotational Non Linear problems

Corotator *
TwoNodeTruss::getCorotator(CoordSet &cs, double *kel, int, int)
{
 return new BarCorotator(nn[0], nn[1], prop->E, prop->A, preload, cs);
}

int
TwoNodeTruss::getTopNumber()
{
  return 101;
}

void
TwoNodeTruss::getThermalForce(CoordSet &cs, Vector &ndTemps,
                              Vector &elementThermalForce, int glflag, GeomState *geomState)
{
// Called from Dynam.C
// Computes the thermal-mechanical coupling force C*theta
// See cupb3d.f of RCFEM. C is 6x2 matrix
// A = cross section, W= dilatation coeff


    double x[2], y[2], z[2], elC[6][2];
    double dx,dy,dz,length;
    int i, j;
    
    if (geomState) {

       // Update the transformation matrix for nonlinear
         
        NodeState &nd1 = (*geomState)[nn[0]];
	NodeState &nd2 = (*geomState)[nn[1]];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;	
	
    } 
     
    else {
    
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
	
        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;	
     
     }
     
    dx = x[1] - x[0];
    dy = y[1] - y[0];
    dz = z[1] - z[0];

    length = sqrt( dx*dx + dy*dy + dz*dz ); 
     
    double Tref  = prop->Ta;
    double coeff = prop->E*prop->W*prop->A/(2*length);

//  Coupling matrix:

    elC[0][0] = -coeff*dx;
    elC[0][1] = -coeff*dx;
    elC[1][0] = -coeff*dy;
    elC[1][1] = -coeff*dy;
    elC[2][0] = -coeff*dz;
    elC[2][1] = -coeff*dz;

    elC[3][0] = coeff*dx;
    elC[3][1] = coeff*dx;
    elC[4][0] = coeff*dy;
    elC[4][1] = coeff*dy;
    elC[5][0] = coeff*dz;
    elC[5][1] = coeff*dz;

    for(i=0; i<6; ++i) {
      elementThermalForce[i] = 0.0;
        for(j=0; j<2; ++j)
          elementThermalForce[i] += elC[i][j]*(ndTemps[j]-Tref);
    }
//     elementThermalForce.print("Thermal Forces");
}


void
TwoNodeTruss::getVonMises(Vector& stress, Vector& weight, CoordSet& cs,
			  Vector& elDisp, int strInd, int surface, 
			  double *ndTemps, double ylayer, double zlayer, int avgnum)
{
#ifndef SALINAS 
   weight = 1.0;

   Node &nd1 = cs.getNode( nn[0] );
   Node &nd2 = cs.getNode( nn[1] );

   double dx = nd2.x - nd1.x;
   double dy = nd2.y - nd1.y;
   double dz = nd2.z - nd1.z;

   double length = sqrt(dx*dx + dy*dy + dz*dz);

   // scale dx, dy, and dz by the length
   dx /= length;
   dy /= length;
   dz /= length;

   // Compute the change in length of the element
   double dq = dx*(elDisp[3]-elDisp[0])
             + dy*(elDisp[4]-elDisp[1])
             + dz*(elDisp[5]-elDisp[2]);

   // Compute axial strain
   double exx = dq/length;


    switch (avgnum) {

      case 0:
      {
	if (strInd == 0) {
	   // Compute axial force
           double f = prop->A*prop->E*exx;
	
           // Add Preload
           f  += preload;

           // Compute thermal force
           double coefficient = prop->E*prop->A*prop->W;
           double Tref = prop->Ta;
 
           double fth1 = coefficient*(ndTemps[0]-Tref);
           double fth2 = coefficient*(ndTemps[1]-Tref);

           // compute stresses 
	   double elForce[2]={0.0,0.0};
           elForce[0] = -f + fth1;
           elForce[1] =  f - fth2;
	   stress[0] = -elForce[0]/prop->A;
	   stress[1] =  elForce[1]/prop->A;
	}
	else if (strInd == 7) {
	   stress[0] =  exx;
	   stress[1] =  exx;
	}
	else {
	   stress[0] = 0.0;
	   stress[1] = 0.0;
	}
        break;
      }

      case 1:
      { 
        double xl[3][3];
        double xg[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
        buildBarFrame(cs, xg, xl);

        if (strInd < 6) {
           // Compute axial force
           double f = prop->A*prop->E*exx;

           // Add Preload
           f  += preload;

           // Compute thermal force
           double coefficient = prop->E*prop->A*prop->W;
           double Tref = prop->Ta;
           double fth1 = coefficient*(ndTemps[0]-Tref);
           double fth2 = coefficient*(ndTemps[1]-Tref);

           // return stresses
           double elForce[2]={0.0,0.0};
           elForce[0] = -f + fth1;
           elForce[1] =  f - fth2;
           double tmpStr1[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
           double tmpStr2[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
           tmpStr1[0] = -elForce[0]/prop->A;
           tmpStr2[0] =  elForce[1]/prop->A;
          _FORTRAN(transform)(xl[0], xl[1], xl[2], xg[0], xg[1], xg[2], tmpStr1);
          _FORTRAN(transform)(xl[0], xl[1], xl[2], xg[0], xg[1], xg[2], tmpStr2);
           stress[0] = tmpStr1[strInd];
           stress[1] = tmpStr2[strInd];
        }
        else if (strInd > 6 && strInd < 13) {
           double tmpStr[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
           tmpStr[0] = exx;
           _FORTRAN(transform)(xl[0], xl[1], xl[2], xg[0], xg[1], xg[2], tmpStr);
           stress[0] = tmpStr[strInd-7];
           stress[1] = tmpStr[strInd-7];
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
#endif
}

void
TwoNodeTruss::buildBarFrame(CoordSet& cs, double xg[3][3], double xl[3][3])
{

        Node &nd1 = cs.getNode( nn[0] );
        Node &nd2 = cs.getNode( nn[1] );

        xl[0][0] = nd2.x - nd1.x;
        xl[0][1] = nd2.y - nd1.y;
        xl[0][2] = nd2.z - nd1.z;

        // local x-axis
        normalize( xl[0] );

        // local y-axis
        double tmp1[3],tmp2[3];

        crossprod(xg[1], xl[0], tmp1); // cross product of global y-axis with local x-axis
        crossprod(xg[2], xl[0], tmp2); // cross product of global z-axis with local x-axis

        if(magnitude(tmp1) > magnitude(tmp2)) {
           xl[1][0] = tmp1[0];
           xl[1][1] = tmp1[1];
           xl[1][2] = tmp1[2];
        } else {
           xl[1][0] = tmp2[0];
           xl[1][1] = tmp2[1];
           xl[1][2] = tmp2[2];
        }
        normalize( xl[1]);

        // local z-axis
        crossprod(xl[0], xl[1], xl[2]);
        normalize( xl[2] );

}
