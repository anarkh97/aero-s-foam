#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/LEIsoParamHexa.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>

#define O3 4
#define O4 7

LEIsoParamHexa::LEIsoParamHexa(int o, int* nodenums) {
 int i;
 order = int(rint(pow(double(o),1.0/3.0)));
 int orderc = order*order*order;
 nn = new int[orderc];
 for(i=0;i<orderc;i++) nn[i] = nodenums[i];
}


LEIsoParamHexa::LEIsoParamHexa(const LEIsoParamHexa& e) {
 order = e.order;
 int orderc = order*order*order;
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = e.nn[i];
}


Element * LEIsoParamHexa::clone() {
 return new LEIsoParamHexa(*this);
}


void LEIsoParamHexa::renum(int *table) {
 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}

extern bool useFull;

int* LEIsoParamHexa::nodes(int *p) {
if (useFull) {

 int orderc = order*order*order;
 if(p == 0) p = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) p[i] = nn[i];
} else {
 if(p == 0) p = new int[8];
 p[0] = nn[0];
 p[1] = nn[order-1];
 p[2] = nn[order*order-1];
 p[3] = nn[order*order-order];                                                                    p[4] = nn[order*order*order-order*order];
 p[5] = nn[order*order*order-order*order+order-1];
 p[6] = nn[order*order*order-1];
 p[7] = nn[order*order*order-order];

}
 return p;
}


int* LEIsoParamHexa::dofs(DofSetArray &dsa, int *p) {

 int i;
 int orderc = order*order*order;
 if(p == 0) p = new int[3*orderc];
 for(i=0;i<orderc;i++) dsa.number(nn[i],DofSet::XYZdisp,p+3*i);
 return p;
}


void LEIsoParamHexa::markDofs(DofSetArray &dsa) {

 int i;
 int orderc = order*order*order;
 for(i=0;i<orderc;i++) dsa.mark(nn[i],DofSet::XYZdisp);
}


double LEIsoParamHexa::getMass(CoordSet&) {
 fprintf(stderr,"LEIsoParamHexa::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix LEIsoParamHexa::massMatrix(CoordSet &cs, double *K, int fl) {

 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 LEMassFunction f(3*orderc,prop->rho,K);
 ipu.zeroOut<double> (9*orderc*orderc,K);
 int gorder = 7;
 if (order<=3) gorder = O3;
 else if (order<=4) gorder = O4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(3*orderc,K);

 FullSquareMatrix ret(3*orderc,K);
 return ret;
}


FullSquareMatrix LEIsoParamHexa::stiffness(CoordSet &cs, double *K, int flg ) {
 IsoParamUtils ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 LEStiffFunction f(3*orderc,prop->E,prop->nu,K);
 ipu.zeroOut<double> (9*orderc*orderc,K);
  int gorder = 7;
 if (order<=3) gorder = O3;
 else if (order<=4) gorder = O4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(3*orderc,K);

 FullSquareMatrix ret(3*orderc,K);
 return ret;
}

void   LEIsoParamHexa::getGravityForce(CoordSet &cs,double *gravity,
                           Vector &force, int gravflag, GeomState *gs) {
 if (gravflag != 2) {
   fprintf(stderr,"LEIsoParamHexa::getGravityForce: lumped is not supported\n");
   return;
 } else {
   IsoParamUtils ipu(order);
   int orderc = ipu.getorderc();
   double *xyz=(double*)alloca(sizeof(double)*3*orderc);
   cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

   LEGravityFunction f(3*orderc,prop->rho, gravity, force.data());
   ipu.zeroOut<double> (3*orderc,force.data());
   int gorder = 7;
   if (order<=3) gorder = O3;
   else if (order<=4) gorder = O4;
   ipu.volumeInt3d(xyz, f, gorder);
 }
}



extern bool useFull;

int
LEIsoParamHexa::numNodes() {
  if(useFull)
    return order*order*order;
  else
    return(8);   // to ignore effect of mid-size nodes in dec
}

