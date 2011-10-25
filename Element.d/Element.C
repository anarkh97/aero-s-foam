#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <Element.d/Element.h>

void
Element::setCompositeData(int, int, double*, double*, double*)
{
  fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for non composite element type %d\n", elementType);
}

double *
Element::setCompositeData2(int, int, double*, double*, CoordSet&, double)
{
  fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for non composite element type %d\n", elementType);
  return 0;
}

void
Element::getVonMisesInt(CoordSet &,Vector &,double &,double &, int,
			double &,double &, double* dT)
{
  assert(0);
}

void
Element::getVonMises(Vector &stress, Vector &weight,
		     CoordSet &, Vector &, int, int, double*,
		     double,double,int)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: getVonMises not implemented for element type %d\n", elementType);
  stress.zero();
  weight.zero();
}

extern "C" {
  void _FORTRAN(vmelmvc)(DComplex*, const int &, const int &, const int &, const int &, const int &);
  void _FORTRAN(strainvmc)(DComplex*, const int &, const int &, const int &, const int &);
}

void
Element::getVonMises(ComplexVector &stress, Vector &weight,
                     CoordSet &cs, ComplexVector &disp, int strInd, int surface,
                     double *ndTemps, double ylayer, double zlayer, int avgnum)
{
 // PJSA 6-23-04
 // split displacement into real & imaginary vectors to compute independently the real
 // and imaginary components of the normal & shear stresses/strains
 // these can then be combined into a complex vector & used to compute the von mises stress or strain
 Vector realDisp(disp.size());
 for(int i=0; i<disp.size(); ++i) realDisp[i] = disp[i].real();
 Vector imagDisp(disp.size());
 for(int i=0; i<disp.size(); ++i) imagDisp[i] = disp[i].imag();

 if( ((strInd >= 0) && (strInd <=5)) || ((strInd >= 7) && (strInd <=12)) ) { // not von mises
   Vector realStress(stress.size(), 0.0);
   getVonMises(realStress, weight, cs, realDisp, strInd, surface, ndTemps, ylayer, zlayer, avgnum);
   Vector imagStress(stress.size(), 0.0);
   getVonMises(imagStress, weight, cs, imagDisp, strInd, surface, ndTemps, ylayer, zlayer, avgnum);
   for(int i=0; i<stress.size(); ++i) stress[i] = DComplex(realStress[i], imagStress[i]);
 }
 else if((strInd==6) || (strInd==13)) { // von mises
   if((ylayer!=0) || (zlayer!=0)) {
     fprintf(stderr," *** WARNING: complex getVonMises not implemented for ylayer %f zlayer %f \n",ylayer,zlayer);
     stress.zero();
     weight.zero();
     return;
   }
   int maxgus = numNodes();
   int si = (strInd == 6) ? 0 : 1;
   FullM realStress(maxgus,9); realStress.zero();
   getAllStress(realStress, weight, cs, realDisp, si, surface, ndTemps);
   FullM imagStress(maxgus,9); imagStress.zero();
   getAllStress(imagStress, weight, cs, imagDisp, si, surface, ndTemps);
   GenFullM<DComplex> complexStress(maxgus,7);  complexStress.zero();
   for(int i=0; i<maxgus; ++i)
     for(int j=0; j<6; ++j) complexStress[i][j] = DComplex(realStress[i][j], imagStress[i][j]);

   if(strInd == 6) _FORTRAN(vmelmvc)(complexStress.data(),maxgus,7,1,1,maxgus);
   else _FORTRAN(strainvmc)(complexStress.data(),maxgus,7,1,maxgus);

   for(int i=0; i<maxgus; ++i) stress[i] = complexStress[i][6];
 }
 else {
   if(!isConstraintElement() && !isSpring())
     fprintf(stderr," *** WARNING: complex getVonMises not implemented for strInd %d\n",strInd);
   stress.zero();
   weight.zero();
 }
}

void
Element::getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                      Vector &elDisp, int strInd, int surface,
                      double *ndTemps)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: getAllStress not implemented for element type %d\n", elementType);
  stress.zero();
  weight.zero();
}

PrioInfo Element::examine(int sub, MultiFront *mf)
{
  fprintf(stderr," *** ERROR: Element type %d cannot be decomposed since examine function is not implemented \n", elementType);
  PrioInfo p;
  p.isReady = false;
  return p;
}

void
Element::getGravityForce(CoordSet&, double *, Vector &force, int, GeomState *)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: Gravity force not implemented for element (%6d), type %3d\n", getGlNum()+1, elementType);
  force.zero();
}

void
Element::getThermalForce(CoordSet&, Vector &, Vector &force, int glflag,
                         GeomState *)
{
  if(!isConstraintElement() && !isSpring())
   // fprintf(stderr," *** WARNING: Thermal force not implemented for element type %d\n", elementType);
  force.zero();
}

void
Element::computeHeatFluxes(Vector &heatflux, CoordSet&, Vector &, int)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: Heat Fluxes not implemented for element type %d\n", elementType);
  heatflux.zero();
}


void
Element::computeSloshDisp(Vector &fluidDispSlosh, CoordSet&, Vector &, int)
{
  if (elementType != 302)  {
      fprintf(stderr," *** WARNING: Fluid Displacements not implemented for element type %d\n", elementType);
  fluidDispSlosh.zero();
  }
}

void
Element::computeSloshDispAll(Vector &fluidDispSlosh, CoordSet&, Vector &)
{
  if (elementType != 302)  {
      fprintf(stderr," *** WARNING: Fluid Displacements not implemented for element type %d\n", elementType);
  fluidDispSlosh.zero();
  }
}

void
Element::getIntrnForce(Vector &elForce, CoordSet&, double *, int,double *)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: Internal force not implemented for element type %d\n", elementType);
  elForce.zero();
}


void
Element::computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *, int cflg)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: Pressure force not implemented for element type %d\n", elementType);
  elPressureForce.zero();
}

int
Element::getTopNumber()
{
 return -1;
}

void
Element::computePressureForce(Node *, Vector& elPressureForce,
                              double *, int cflg)
{
 if(!isConstraintElement() && !isSpring())
   fprintf(stderr," *** WARNING: Pressure force 2 not implemented for element %d\n", glNum);
 elPressureForce.zero();
}

void
Element::addFaces(PolygonSet *)
{
  fprintf (stderr, "\n\n ERROR: addFaces not implemented for this element !!!"
                   " \n   exiting... \n\n");
  exit (1);
}

void
Element::setMaterial(NLMaterial *)
{
  fprintf(stderr, "WARNING: Trying to use Material on unsupported element!\n");
}

void
Element::getCG(CoordSet &cset, double &xcg, double &ycg, double &zcg)
{
  int *myNodes;
  int preAllocNodes[125];
  int nnd = numNodes()-numInternalNodes();
  if(numNodes() > 125)
    myNodes = new int[numNodes()];
  else
    myNodes = preAllocNodes;

  nodes(myNodes);

  xcg = ycg = zcg = 0.0;

  int i;
  for(i = 0; i < nnd; ++i) {
    Node *n = cset[myNodes[i]];
    xcg += n->x;
    ycg += n->y;
    zcg += n->z;
  }
  xcg /= nnd;
  ycg /= nnd;
  zcg /= nnd;
  if(numNodes() > 125)
    delete[] myNodes;
}

// END DEC

FullSquareMatrix Element::stiffness(CoordSet&, double* kel, int cmflg)
{
  FullSquareMatrix result(1, kel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrix Element::massMatrix(CoordSet&, double* mel, int cmflg)
{
  FullSquareMatrix result(1, mel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrix
Element::dampingMatrix(CoordSet& cs, double *m, int cmflg)
{
  FullSquareMatrix ret(1,m);
  ret.setSize(numDofs());
  ret.zero();
  return ret;
}

FullSquareMatrix
Element::imStiffness(CoordSet& cs, double *m, int cmflg)
{
  fprintf(stderr,"  ElementC imStiffness not implmented for this element\n");
  FullSquareMatrix ret(4,m);
  ret.zero();
  return ret;
}

void Element::lumpMatrix(FullSquareMatrix& m)
{
  double MM = 0.0; // total mass of element
  double MD = 0.0; // mass of the diagonal
  const int dim = m.dim();
  for(int i = 0; i < dim; ++i) {
    MD += m[i][i];
    for(int j = 0; j < i; ++j) { 
      MM += m[i][j]; 
    }
  }
  MM = MD + 2.0*MM;
  if (MD > 0) {
    const double factor = MM/MD;
    for(int i = 0; i < dim; ++i) {
      m[i][i] *= factor;
      //cerr << "i = " << i << ", m[i][i] = " << m[i][i] << endl;
      factors.push_back(m[i][i]/MM); // PJSA store this for getGravityForce
      for(int j = 0; j < i; ++j) {
        m[i][j] = m[j][i] = 0.0;
      }
    }
  }
  return;
}

// mratio = 1.0 for consistent
// mratio = 0.0 for lumped
FullSquareMatrix Element::massMatrix(CoordSet& cs, double* m, double mratio)
{
  if(mratio == 0.0 && getMassType() == 1) { // in this case, get the consistent mass matrix and lump it using diagonal scaling
    FullSquareMatrix result = massMatrix(cs, m, 1);
    lumpMatrix(result);
    return result;
  }
  else return massMatrix(cs, m, int(mratio));
}

FullSquareMatrixC Element::complexStiffness(CoordSet&, DComplex* kel, int cmflg)
{
  FullSquareMatrixC result(1, kel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrixC Element::complexDampingMatrix(CoordSet&, DComplex* cel, int cmflg)
{
  FullSquareMatrixC result(1, cel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrixC Element::complexMassMatrix(CoordSet&, DComplex* mel, double mratio)
{
  FullSquareMatrixC result(1, mel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

#include <Element.d/Helm.d/HelmElement.h>

bool Element::isFluidElement() { return dynamic_cast<HelmElement *>(this); }

double
Element::computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs, GeomState *gs,
                                  double stable_tol, int stable_maxit)
{
  if(prop) {

      using std::sqrt;
      double eigmax;
      double relTol    = stable_tol; // stable_tol default is 1.0e-3
      double preeigmax = 0.0;

      int numdofs = K.dim();
      int maxIte  = stable_maxit; // stable_maxit default is 100

      Vector v(numdofs);
      Vector z(numdofs);

// Starts from an arbitrary array.
      int i,j;
      for (i=0; i<numdofs; ++i)
        v[i] = (double) (i+1) / (double) numdofs;

// Power iteration loop

      for (i=0; i<maxIte; ++i) {
        z.zero();
        K.multiply(v,z,1);

        for (j=0; j< numdofs; ++j)
          z[j] /= M[j][j];

// Normalize

        double zmax = z[0];
        for (j=1; j< numdofs; ++j)
          if (abs(z[j])>zmax) zmax = abs(z[j]);

        eigmax = zmax;

        v = (1.0/zmax)*z;

        if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;

        preeigmax = eigmax;
      }

      // compute stability maximum time step
      double sdt = 2.0 / sqrt(eigmax);

      return sdt;
  }
  else { // phantom
      return std::numeric_limits<double>::infinity();
  }
}
