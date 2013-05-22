#include <cstdio>
#include <cstdlib>
#include <alloca.h>

#include <Element.d/Helm.d/LEIsoParamTetra.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>


LEIsoParamTetra::LEIsoParamTetra(int o, int* nodenums) {

 if (o==4) order = 2;
 else if (o==10) order = 3;
 else if (o==20) order = 4;
 else if (o==35) order = 5;
 else if (o==56) order = 6;
 else {
   fprintf(stderr,"Order too high in LEIsoParamTetra::LEIsoParamTetra.\n");
   exit(-1);
 }
 
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = nodenums[i];
}


LEIsoParamTetra::LEIsoParamTetra(const LEIsoParamTetra& e) {
 order = e.order;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 nn = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) nn[i] = e.nn[i];
}


Element * LEIsoParamTetra::clone() {
 return new LEIsoParamTetra(*this);
}


void LEIsoParamTetra::renum(int *table) {
 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}


void LEIsoParamTetra::renum(EleRenumMap& table) {
 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) nn[i] = table[nn[i]];
}


int* LEIsoParamTetra::nodes(int *p) {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 if(p == 0) p = new int[orderc];
 int i;
 for(i=0;i<orderc;i++) p[i] = nn[i];
 return p;
}


int* LEIsoParamTetra::dofs(DofSetArray &dsa, int *p) {

 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 if(p == 0) p = new int[3*orderc];
 for(i=0;i<orderc;i++) dsa.number(nn[i],DofSet::XYZdisp,p+3*i);
 return p;
}


void LEIsoParamTetra::markDofs(DofSetArray &dsa) {

 int i;
 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 for(i=0;i<orderc;i++) dsa.mark(nn[i],DofSet::XYZdisp);
}


double LEIsoParamTetra::getMass(CoordSet&) {
 fprintf(stderr,"LEIsoParamTetra::getMass not implemented.\n");
 return 0.0;
}


FullSquareMatrix LEIsoParamTetra::massMatrix(CoordSet &cs, double *K, int fl) {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 LEMassFunction f(3*orderc,prop->rho,K);
 ipu.zeroOut<double> (9*orderc*orderc,K);
 int gorder = 7*7*7;
 if (order<=3) gorder = 4*4*4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(3*orderc,K);

 FullSquareMatrix ret(3*orderc,K);
 return ret;
}


FullSquareMatrix LEIsoParamTetra::stiffness(CoordSet &cs, double *K, int flg ) {

 IsoParamUtilsTetra ipu(order);
 int orderc = ipu.getorderc();
 double *xyz=(double*)alloca(sizeof(double)*3*orderc);
 cs.getCoordinates(nn,orderc,xyz,xyz+orderc,xyz+2*orderc);

 LEStiffFunction f(3*orderc,prop->E, prop->nu,K);
 ipu.zeroOut<double> (9*orderc*orderc,K);
 int gorder = 7*7*7;
 if (order<=3) gorder = 4*4*4;
 ipu.volumeInt3d(xyz, f, gorder);
 ipu.symmetrize(3*orderc,K);

 FullSquareMatrix ret(3*orderc,K);
 return ret;
}


extern bool useFull;

int
LEIsoParamTetra::numNodes() {
  //Not tested -JF
  if(useFull)
    return (order*(order+1)*(order+2))/6;
  else
    return(4);   // to ignore effect of mid-size nodes in dec
}

