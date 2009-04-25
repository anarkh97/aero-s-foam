#include        <stdlib.h>

#include	<Element.d/Shell.d/Therm3NoShell.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Corotational.d/Shell3Corotator.h>
#include        <Corotational.d/utilities.h>
#include	<Element.d/State.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>


// localxy - Computes local coordinates and area

extern "C"      {
void _FORTRAN(localxy)(double*, double*, double*, double*, double*, int&, double&);
}

// THREE NODE SHELL THERMAL ELEMENT

Therm3NoShell::Therm3NoShell(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}

Element *
Therm3NoShell::clone()
{
	return new Therm3NoShell(*this);
}

void
Therm3NoShell::renum(int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
}

double
Therm3NoShell::getMass(CoordSet& cs)
{
        return 0.0;
}

FullSquareMatrix
Therm3NoShell::massMatrix(CoordSet &cs, double *mel,int cmflg)
{

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], area;
        double xl[3], yl[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        double capacitance = prop->rho*prop->Q*prop->eh;

// Construct Area

        int localxyflag = 0;
       _FORTRAN(localxy)(x, y, z, xl, yl, localxyflag, area);

// Mass elements

        FullSquareMatrix ret(3,mel);

        ret.zero();

        for (int i =0; i<3; ++i)
           ret[i][i] = capacitance*area/3;

        return ret;

}

FullSquareMatrix
Therm3NoShell::stiffness(CoordSet &cs, double *d, int flg)
{

	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3];
	double xl[3], yl[3], area;

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        double h = prop->eh;

        // Check for a zero thickness
	if(h <= 0.0) {
          fprintf(stderr," *** ERROR: Zero shell thickness "
                         "(Therm3NoShell.C) %d %d %d\n",
                         nn[0]+1, nn[1]+1, nn[2]+1);
          fprintf(stderr," ... exiting fem program ...\n");
          exit(-1);
        }

        // Construct Local Coordinates and Area

        int localxyflag = 1;

       _FORTRAN(localxy)(x, y, z, xl, yl, localxyflag, area);
  
        // Stiffness Elements based on local coordinates
   
        double area2 = area*2;

        double m21 = (yl[1]-yl[2])/area2;
        double m22 = (yl[2]-yl[0])/area2;
        double m23 = (yl[0]-yl[1])/area2;

        double m31 = (xl[2]-xl[1])/area2;
        double m32 = (xl[0]-xl[2])/area2;
        double m33 = (xl[1]-xl[0])/area2;

        FullSquareMatrix K(3,d);

        double ke = prop->k*h*area;

        K[0][0] = ke*(m21*m21 + m31*m31);
        K[0][1] = ke*(m21*m22 + m31*m32);
        K[0][2] = ke*(m23*m21 + m33*m31);
        K[1][0] = ke*(m21*m22 + m31*m32);
        K[1][1] = ke*(m22*m22 + m32*m32);
        K[1][2] = ke*(m22*m23 + m32*m33);
        K[2][0] = ke*(m23*m21 + m33*m31);
        K[2][1] = ke*(m22*m23 + m32*m33);
        K[2][2] = ke*(m23*m23 + m33*m33);

        return K;
}

int
Therm3NoShell::numNodes()
{
 	return 3;
}

int*
Therm3NoShell::nodes(int *p)
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
Therm3NoShell::numDofs()
{
 	return 3;
}

int*
Therm3NoShell::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[3];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);

	return p;
}

void
Therm3NoShell::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 3, DofSet::Temp);
}

void
Therm3NoShell::computeTemp(CoordSet&cs, State &state, double gp[2],
                            double *tres)
{
 double Temp[3][2];

 state.getTemp(nn[0], Temp[0], Temp[0]+1);
 state.getTemp(nn[1], Temp[1], Temp[1]+1);
 state.getTemp(nn[2], Temp[2], Temp[2]+1);

// tres[0] = temperature
// tres[1] = d(Temperature)/dt

 int j;
 for(j=0; j<2; ++j)
    tres[j] = (1-gp[0]-gp[1])* Temp[0][j] +
                       gp[0] * Temp[1][j] +
                       gp[1] * Temp[2][j] ;

}

void
Therm3NoShell::getFlFlux(double gp[2], double *flF, double *tresF)
{
// Projects a fluid flux contained in flF[0] to all 3 nodes of triangle
// Returns tresF
// fprintf(stderr, "Gauss Points %f %f\n ", gp[0], gp[1]);

   tresF[0]  = (1-gp[0]-gp[1])* flF[0];
   tresF[1]  = gp[0] * flF[0];
   tresF[2]  = gp[1] * flF[0];

}

int
Therm3NoShell::getTopNumber()
{
  return 146;//4;
}
