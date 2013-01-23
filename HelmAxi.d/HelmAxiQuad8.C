#include <cmath>
#include <Driver.d/PolygonSet.h>
#include <Math.d/matrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Math.d/FullSquareMatrix.h>
#include <HelmAxi.d/HelmAxiQuad8.h>

extern "C"      {
void _FORTRAN(quad8aximas)(double*,double*,const int&,double*,const int&);
};
extern "C"	{
void _FORTRAN(quad8axistif1)(double*,double*, const int&, double*, const int&);
};
extern "C"      {
void _FORTRAN(quad8axistif2)(double*,double*, const int&, double*, const int&);
};

HelmAxiQuad8::HelmAxiQuad8(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
	nn[4] = nodenums[4];
	nn[5] = nodenums[5];
	nn[6] = nodenums[6];
	nn[7] = nodenums[7];
}

//AxiHElement * HelmAxiQuad8::clone()
Element * HelmAxiQuad8::clone()
{
 return new HelmAxiQuad8(*this);
}

void HelmAxiQuad8::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
  nn[6] = table[nn[6]];
  nn[7] = table[nn[7]];
}

FullSquareMatrix
HelmAxiQuad8::massMatrix(CoordSet &cs, double *d, int cmflg)
{
        FullSquareMatrix mass(8,d);

	mass.zero();

        return mass;
}


FullSquareMatrix
HelmAxiQuad8::stiffness(CoordSet &cs, double *K, int flg) {

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
	Node nd5 = cs.getNode(nn[4]);
	Node nd6 = cs.getNode(nn[5]);
	Node nd7 = cs.getNode(nn[6]);
	Node nd8 = cs.getNode(nn[7]);

        int i;
	double x[8], y[8], Kstiff[64];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;
	x[4] = nd5.x; y[4] = nd5.y; 
	x[5] = nd6.x; y[5] = nd6.y;
	x[6] = nd7.x; y[6] = nd7.y; 
	x[7] = nd8.x; y[7] = nd8.y;

        // Calculate a part of the stiffness matrix 

        // First get the stiffness  
        _FORTRAN(quad8axistif1)(x, y, 3, K, 8);

        // Put the stiffness temporarely in Kstiff
        for (i=0;i<64;++i)
          Kstiff[i] = K[i];

        // Then get the mass 
        // For while the wave number = el. Area
        _FORTRAN(quad8aximas)(x, y, 3, K, 8);

        double kappa = prop->kappaHelm;

        // Now compute [K] - kappa^2 [M] 
        for (i=0;i<64;++i)
            K[i] = Kstiff[i] - (kappa*kappa)*K[i];
       
        FullSquareMatrix ret(8, K);

        return ret;
}


FullSquareMatrix
HelmAxiQuad8::stiffteta(CoordSet &cs,double *K)
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
	Node nd5 = cs.getNode(nn[4]);
	Node nd6 = cs.getNode(nn[5]);
	Node nd7 = cs.getNode(nn[6]);
	Node nd8 = cs.getNode(nn[7]);

	double x[8], y[8];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;
	x[4] = nd5.x; y[4] = nd5.y; 
	x[5] = nd6.x; y[5] = nd6.y;
	x[6] = nd7.x; y[6] = nd7.y; 
	x[7] = nd8.x; y[7] = nd8.y;

        // Calculate stiffness matrix here:
        _FORTRAN(quad8axistif2)(x, y, 3, K, 8);

        FullSquareMatrix ret(8, K);

        return ret;
}

int
HelmAxiQuad8::numNodes()
{
 	return 8;
}

int*
HelmAxiQuad8::nodes(int *p)
{
 	if(p == 0) p = new int[8];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
 	p[4] = nn[4];
 	p[5] = nn[5];
 	p[6] = nn[6];
 	p[7] = nn[7];
	return p;
}

int
HelmAxiQuad8::numDofs()
{
 	return 8;
}

int*
HelmAxiQuad8::dofs(DofSetArray &dsa, int *p) {

 if(p == 0) p = new int[8];

 dsa.number(nn[0],DofSet::Helm , p);
 dsa.number(nn[1],DofSet::Helm , p+1);
 dsa.number(nn[2],DofSet::Helm , p+2);
 dsa.number(nn[3],DofSet::Helm , p+3);
 dsa.number(nn[4],DofSet::Helm , p+4);
 dsa.number(nn[5],DofSet::Helm , p+5);
 dsa.number(nn[6],DofSet::Helm , p+6);
 dsa.number(nn[7],DofSet::Helm , p+7);

 return p;

}

void
HelmAxiQuad8::markDofs(DofSetArray &dsa) {

 dsa.mark(nn[0],DofSet::Helm);
 dsa.mark(nn[1],DofSet::Helm);
 dsa.mark(nn[2],DofSet::Helm);
 dsa.mark(nn[3],DofSet::Helm);
 dsa.mark(nn[4],DofSet::Helm);
 dsa.mark(nn[5],DofSet::Helm);
 dsa.mark(nn[6],DofSet::Helm);
 dsa.mark(nn[7],DofSet::Helm);

}

void
HelmAxiQuad8::addFaces(PolygonSet *pset) {

 pset->addLine(nn[0], nn[4]);
 pset->addLine(nn[4], nn[1]);
 pset->addLine(nn[1], nn[5]);
 pset->addLine(nn[5], nn[2]);
 pset->addLine(nn[2], nn[6]);
 pset->addLine(nn[6], nn[3]);
 pset->addLine(nn[3], nn[7]);
 pset->addLine(nn[7], nn[0]); 

}


void HelmAxiQuad8::buildMesh3D(int &elemNum, FILE *outF, 
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


void HelmAxiQuad8::buildMesh2D(int &elemNum, FILE *outF,
                              int nodeInc, int numSlices)
{
 int iSlice;
 for(iSlice = 0; iSlice < numSlices; ++iSlice) {
   fprintf(outF, "%d  12  %d %d %d %d %d %d %d %d \n",
     ++elemNum,
     nn[0] + iSlice*nodeInc +1, nn[1] + iSlice*nodeInc +1,
     nn[2] + iSlice*nodeInc +1, nn[3] + iSlice*nodeInc +1,
     nn[4] + iSlice*nodeInc +1, nn[5] + iSlice*nodeInc +1,
     nn[6] + iSlice*nodeInc +1, nn[7] + iSlice*nodeInc +1);
 }
}

