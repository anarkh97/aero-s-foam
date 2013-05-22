#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Element.d/Radiation.d/BarRadiation.h>
#include <Corotational.d/BarThermalCorotator.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>

BarRadiation::BarRadiation(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

Element *
BarRadiation::clone()
{
	return new BarRadiation(*this);
}

void
BarRadiation::renum(int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
BarRadiation::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

double
BarRadiation::getMass(CoordSet& cs)
{
        return 0.0;
}

FullSquareMatrix
BarRadiation::massMatrix(CoordSet &cs, double *mel, int cmflg)
{

        FullSquareMatrix elementMassMatrix(2,mel);

// zero the element mass matrix
	elementMassMatrix.zero();

        return elementMassMatrix;
}

FullSquareMatrix
BarRadiation::stiffness(CoordSet &cs, double *Kcv, int flg)
{
// This is the additional matrix when radiation is present.
// It is added into the conductance matrix.

        Node &nd1 = cs.getNode( nn[0] );
        Node &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	//double dx = x[1] - x[0];
	//double dy = y[1] - y[0];
	//double dz = z[1] - z[0];

	//double length = sqrt( dx*dx + dy*dy + dz*dz );

//... Construct radiative matrix ...

        FullSquareMatrix ret(2,Kcv);

        ret[0][0] = 0.0;
        ret[1][1] = 0.0;
        ret[1][0] = 0.0;
        ret[0][1] = 0.0;
        
        return ret;
}


Corotator *
BarRadiation::getCorotator(CoordSet &cs, double* kel, int, int)
{
 return new BarThermalCorotator(nn[0], nn[1], prop->P, prop->eps, prop->sigma, prop->Tr, cs); 
}

int
BarRadiation::numNodes()
{
        return 2;
}

int *
BarRadiation::nodes(int *p)
{
        if(p == 0) p = new int[2];
        p[0] = nn[0];
        p[1] = nn[1];
        return p;
}

int
BarRadiation::numDofs()
{
        return 2;
}

int *
BarRadiation::dofs(DofSetArray &dsa, int *p)
{
        if(p == 0) p = new int[2];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);

        return p;
}

void
BarRadiation::markDofs( DofSetArray &dsa )
{
        dsa.mark( nn, 2, DofSet::Temp);
}

int
BarRadiation::getTopNumber()
{
  return 147;
}
