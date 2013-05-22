#include <cmath>
#include <Driver.d/PolygonSet.h>
#include <Math.d/matrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Math.d/FullSquareMatrix.h>
#include <HelmAxi.d/HelmAxiTri6.h>

extern "C" {
void _FORTRAN(tri6aximas)(double*, double*, const int&, double*, const int&);
};
extern "C" {
void _FORTRAN(tri6axistif1)(double*, double*, const int&, double*, const int&);
};
extern "C" {
void _FORTRAN(tri6axistif2)(double*, double*, const int&, double*, const int&);
};


HelmAxiTri6::HelmAxiTri6(int* nodenums) {

 nn[0] = nodenums[0];
 nn[1] = nodenums[1];
 nn[2] = nodenums[2];
 nn[3] = nodenums[3];
 nn[4] = nodenums[4];
 nn[5] = nodenums[5];

}


Element *
HelmAxiTri6::clone() {

 return new HelmAxiTri6(*this);

}


void
HelmAxiTri6::renum(int *table) {

 nn[0] = table[nn[0]];
 nn[1] = table[nn[1]];
 nn[2] = table[nn[2]];
 nn[3] = table[nn[3]];
 nn[4] = table[nn[4]];
 nn[5] = table[nn[5]];

}

void
HelmAxiTri6::renum(EleRenumMap &table) {

 nn[0] = table[nn[0]];
 nn[1] = table[nn[1]];
 nn[2] = table[nn[2]];
 nn[3] = table[nn[3]];
 nn[4] = table[nn[4]];
 nn[5] = table[nn[5]];

}


FullSquareMatrix
HelmAxiTri6::massMatrix(CoordSet &cs,double *mel,int cmflg) {

 FullSquareMatrix ret(6,mel);
 ret.zero();

 return ret;

}


FullSquareMatrix
HelmAxiTri6::stiffness(CoordSet &cs, double *K, int flg) {

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);
 Node nd3 = cs.getNode(nn[2]);
 Node nd4 = cs.getNode(nn[3]);
 Node nd5 = cs.getNode(nn[4]);
 Node nd6 = cs.getNode(nn[5]);

 int i;
 double x[6], y[6], Kstiff[36];
 x[0] = nd1.x; y[0] = nd1.y; 
 x[1] = nd2.x; y[1] = nd2.y;
 x[2] = nd3.x; y[2] = nd3.y; 
 x[3] = nd4.x; y[3] = nd4.y; 
 x[4] = nd5.x; y[4] = nd5.y;
 x[5] = nd6.x; y[5] = nd6.y; 

 // Calculate a part of the stiffness matrix

 // First get the stiffness
 _FORTRAN(tri6axistif1)(x, y, 3, K, 6);

 // Put the stiffness temporarely in Kstiff
 for (i=0; i<36; ++i)
   Kstiff[i] = K[i];

 // Then get the mass
 // For while the wave number = el. Area
 _FORTRAN(tri6aximas)(x, y, 3, K, 6);

 double kappa = prop->kappaHelm;

 // Now compute [K] - kappa^2 [M]
 for (i=0; i<36; ++i)
   K[i] = Kstiff[i] - (kappa*kappa)*K[i];

 FullSquareMatrix ret(6, K);

 return ret;

}


FullSquareMatrix
HelmAxiTri6::stiffteta(CoordSet &cs,double *K) {

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);
 Node nd3 = cs.getNode(nn[2]);
 Node nd4 = cs.getNode(nn[3]);
 Node nd5 = cs.getNode(nn[4]);
 Node nd6 = cs.getNode(nn[5]);

 double x[6], y[6];
 x[0] = nd1.x; y[0] = nd1.y; 
 x[1] = nd2.x; y[1] = nd2.y;
 x[2] = nd3.x; y[2] = nd3.y; 
 x[3] = nd4.x; y[3] = nd4.y; 
 x[4] = nd5.x; y[4] = nd5.y;
 x[5] = nd6.x; y[5] = nd6.y; 

 // Calculate stiffness matrix here:
 _FORTRAN(tri6axistif2)(x, y, 3, K, 6);

 FullSquareMatrix ret(6, K);

 return ret;

}


int
HelmAxiTri6::numNodes() {

 return 6;

}


int*
HelmAxiTri6::nodes(int *p) {

 if(p == 0) p = new int[6];

 p[0] = nn[0];
 p[1] = nn[1];
 p[2] = nn[2];
 p[3] = nn[3];
 p[4] = nn[4];
 p[5] = nn[5];

 return p;

}


int
HelmAxiTri6::numDofs() {

 return 6;

}


int*
HelmAxiTri6::dofs(DofSetArray &dsa, int *p) {

 if (p == 0) p = new int[6];

 dsa.number(nn[0], DofSet::Helm, p);
 dsa.number(nn[1], DofSet::Helm, p+1);
 dsa.number(nn[2], DofSet::Helm, p+2);
 dsa.number(nn[3], DofSet::Helm, p+3);
 dsa.number(nn[4], DofSet::Helm, p+4);
 dsa.number(nn[5], DofSet::Helm, p+5);

 return p;

}


void
HelmAxiTri6::markDofs(DofSetArray &dsa) {

 dsa.mark(nn[0], DofSet::Helm);
 dsa.mark(nn[1], DofSet::Helm);
 dsa.mark(nn[2], DofSet::Helm);
 dsa.mark(nn[3], DofSet::Helm);
 dsa.mark(nn[4], DofSet::Helm);
 dsa.mark(nn[5], DofSet::Helm);

}


void
HelmAxiTri6::addFaces(PolygonSet *pset) {

 pset->addLine(nn[0], nn[3]);
 pset->addLine(nn[3], nn[1]);
 pset->addLine(nn[1], nn[4]);
 pset->addLine(nn[4], nn[2]);
 pset->addLine(nn[2], nn[5]);
 pset->addLine(nn[5], nn[0]);

}


void HelmAxiTri6::buildMesh3D(int &elemNum, FILE *outF,
                              int nodeInc, int numSlices) {

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


void HelmAxiTri6::buildMesh2D(int &elemNum, FILE *outF,
                              int nodeInc, int numSlices) {

 for(int iSlice = 0; iSlice < numSlices; ++iSlice) {
   fprintf(outF, "%d  8  %d %d %d %d %d %d\n",
     ++elemNum,
     nn[0] + iSlice*nodeInc +1, nn[1] + iSlice*nodeInc +1,
     nn[2] + iSlice*nodeInc +1, nn[3] + iSlice*nodeInc +1,
     nn[4] + iSlice*nodeInc +1, nn[5] + iSlice*nodeInc +1);
 }

}
