#include        <cmath>
#include <Driver.d/PolygonSet.h>
#include        <Math.d/matrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <HelmAxi.d/HelmAxiQuad.h>

extern "C"      {
void    _FORTRAN(quadaximas)(double*, double*, const int&, double*, const int&);
};
extern "C"	{
void    _FORTRAN(quadaxistif1)(double*, double*, const int&, double*, const int&);
};
extern "C"      {
void    _FORTRAN(quadaxistif2)(double*, double*, const int&, double*, const int&);
};

HelmAxiQuad::HelmAxiQuad(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

//AxiHElement * HelmAxiQuad::clone()
Element * HelmAxiQuad::clone()
{
 return new HelmAxiQuad(*this);
}

void HelmAxiQuad::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}


void HelmAxiQuad::renum(EleRenumMap &table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

FullSquareMatrix
HelmAxiQuad::massMatrix(CoordSet &cs, double *d, int cmflg)
{
        FullSquareMatrix mass(4,d);

	mass.zero();

        return mass;
}


FullSquareMatrix
HelmAxiQuad::stiffness(CoordSet &cs, double *K, int flg) {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);

        int i;
	double x[4], y[4], Kstiff[16];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;

        // Calculate a part of the stiffness matrix 

        // First get the stiffness  
        _FORTRAN(quadaxistif1)(x, y, 2, K, 4);

        // Put the stiffness temporarely in Kstiff
        for (i=0;i<16;++i)
          Kstiff[i] = K[i];

        // Then get the mass 
        // For while the wave number = el. Area
        _FORTRAN(quadaximas)(x, y, 2, K, 4);

        double kappa = prop->kappaHelm;

        // Now compute [K] - kappa^2 [M] 
        for (i=0;i<16;++i)
            K[i] = Kstiff[i] - (kappa*kappa)*K[i];
       
        FullSquareMatrix ret(4, K);

        return ret;
}


FullSquareMatrix
HelmAxiQuad::stiffteta(CoordSet &cs,double *K)
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);

        double x[4], y[4];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        // Calculate stiffness matrix here:
        _FORTRAN(quadaxistif2)(x, y, 2, K, 4);

        FullSquareMatrix ret(4, K);

        return ret;
}

int
HelmAxiQuad::numNodes()
{
 	return 4;
}

int*
HelmAxiQuad::nodes(int *p)
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
HelmAxiQuad::numDofs()
{
 	return 4;
}

int*
HelmAxiQuad::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[4];

        dsa.number(nn[0],DofSet::Helm , p);
        dsa.number(nn[1],DofSet::Helm , p+1);
        dsa.number(nn[2],DofSet::Helm , p+2);
        dsa.number(nn[3],DofSet::Helm , p+3);

	return p;
}

void
HelmAxiQuad::markDofs(DofSetArray &dsa)
{
 	dsa.mark(nn[0],DofSet::Helm);
 	dsa.mark(nn[1],DofSet::Helm);
 	dsa.mark(nn[2],DofSet::Helm);
 	dsa.mark(nn[3],DofSet::Helm);
}

void
HelmAxiQuad::addFaces(PolygonSet *pset)
{

	pset->addLine(nn[0], nn[1]);	
	pset->addLine(nn[1], nn[2]);	
	pset->addLine(nn[2], nn[3]);	
	pset->addLine(nn[3], nn[0]);	

}


void HelmAxiQuad::buildMesh3D(int &elemNum, FILE *outF, 
                              int nodeInc, int numSlices)
{
 int iSlice;
 for(iSlice = 0; iSlice < numSlices; ++iSlice) {
   fprintf(outF, "%d  3  %d %d %d %d %d %d %d %d\n",
     ++elemNum,
     nn[0] + iSlice*nodeInc +1,
     nn[1] + iSlice*nodeInc +1,
     nn[2] + iSlice*nodeInc +1,
     nn[3] + iSlice*nodeInc +1,
     nn[0] + ((iSlice+1)%numSlices)*nodeInc +1,
     nn[1] + ((iSlice+1)%numSlices)*nodeInc +1,
     nn[2] + ((iSlice+1)%numSlices)*nodeInc +1,
     nn[3] + ((iSlice+1)%numSlices)*nodeInc +1 );
 }
}


void HelmAxiQuad::buildMesh2D(int &elemNum, FILE *outF,
                              int nodeInc, int numSlices)
{
 int iSlice;
 for(iSlice = 0; iSlice < numSlices; ++iSlice) {
   fprintf(outF, "%d  2  %d %d %d %d \n",
     ++elemNum,
     nn[0] + iSlice*nodeInc +1,
     nn[1] + iSlice*nodeInc +1,
     nn[2] + iSlice*nodeInc +1,
     nn[3] + iSlice*nodeInc +1);
 }
}

