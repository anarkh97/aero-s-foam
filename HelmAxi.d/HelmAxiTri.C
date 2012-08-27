#include <cmath>
#include <Driver.d/PolygonSet.h>
#include <Math.d/matrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Math.d/FullSquareMatrix.h>
#include <HelmAxi.d/HelmAxiTri.h>

extern "C" {
void _FORTRAN(tri3aximas)(double*, double*, const int&, double*, const int&);
};
extern "C" {
void _FORTRAN(tri3axistif1)(double*, double*, const int&, double*, const int&);
};
extern "C" {
void _FORTRAN(tri3axistif2)(double*, double*, const int&, double*, const int&);
};


HelmAxiTri::HelmAxiTri(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}


Element *
HelmAxiTri::clone()
{
 return new HelmAxiTri(*this);
}


void
HelmAxiTri::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}


FullSquareMatrix
HelmAxiTri::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        FullSquareMatrix ret(3,mel);

	ret.zero();

        return ret;
}


FullSquareMatrix
HelmAxiTri::stiffness(CoordSet &cs, double *K, int flg) {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	int i;
        double x[3], y[3], Kstiff[9];
	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 

        // Calculate a part of the stiffness matrix

        // First get the stiffness
        _FORTRAN(tri3axistif1)(x, y, 3, K, 3);

        // Put the stiffness temporarely in Kstiff
        for (i=0; i<9; ++i)
          Kstiff[i] = K[i];

        // Then get the mass
        // For while the wave number = el. Area
        _FORTRAN(tri3aximas)(x, y, 3, K, 3);

        double kappa = prop->kappaHelm;

        // Now compute [K] - kappa^2 [M]
        for (i=0; i<9; ++i)
            K[i] = Kstiff[i] - (kappa*kappa)*K[i];

        FullSquareMatrix ret(3, K);

        return ret;
}


FullSquareMatrix
HelmAxiTri::stiffteta(CoordSet &cs,double *K)
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);

        double x[3], y[3];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;

        // Calculate stiffness matrix here:
        _FORTRAN(tri3axistif2)(x, y, 3, K, 3);

        FullSquareMatrix ret(3, K);

        return ret;
}


int
HelmAxiTri::numNodes() {

 	return 3;
}


int*
HelmAxiTri::nodes(int *p)
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}


int
HelmAxiTri::numDofs()
{
 	return 3;
}


int*
HelmAxiTri::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[3];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p+1);
	dsa.number(nn[2], DofSet::Helm, p+2);

	return p;
}


void
HelmAxiTri::markDofs(DofSetArray &dsa)
{
 	dsa.mark(nn[0], DofSet::Helm);
 	dsa.mark(nn[1], DofSet::Helm);
 	dsa.mark(nn[2], DofSet::Helm);
}


void
HelmAxiTri::addFaces(PolygonSet *pset)
{

        pset->addLine(nn[0], nn[1]);
        pset->addLine(nn[1], nn[2]);
        pset->addLine(nn[2], nn[0]);

}


void HelmAxiTri::buildMesh3D(int &elemNum, FILE *outF,
                              int nodeInc, int numSlices)
{
 for(int iSlice = 0; iSlice < numSlices; ++iSlice) {
   fprintf(outF, "%d  3  %d %d %d %d %d %d %d %d\n",
     ++elemNum,
     nn[0] + iSlice*nodeInc +1,
     nn[1] + iSlice*nodeInc +1,
     nn[2] + iSlice*nodeInc +1,
     nn[0] + iSlice*nodeInc +1,
     nn[0] + ((iSlice+1)%numSlices)*nodeInc +1,
     nn[1] + ((iSlice+1)%numSlices)*nodeInc +1,
     nn[2] + ((iSlice+1)%numSlices)*nodeInc +1,
     nn[0] + ((iSlice+1)%numSlices)*nodeInc +1 );
 }
}


void HelmAxiTri::buildMesh2D(int &elemNum, FILE *outF,
                              int nodeInc, int numSlices)
{
 for(int iSlice = 0; iSlice < numSlices; ++iSlice) {
   fprintf(outF, "%d  4  %d %d %d \n",
     ++elemNum,
     nn[0] + iSlice*nodeInc +1,
     nn[1] + iSlice*nodeInc +1,
     nn[2] + iSlice*nodeInc +1);
 }
}
