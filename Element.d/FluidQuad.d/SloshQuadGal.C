#include        <Element.d/FluidQuad.d/SloshQuadGal.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <Element.d/State.h>

extern "C"      {
// Overload thermquad3b for the stiffness operator of the sloshing problem
void    _FORTRAN(thermquad3b)(double*, double*, double&, double&, double&, 
                              double*, const int *, const int *);
// void    _FORTRAN(q4d1dofmas)(double*, double*, const int&, double*, const int&);
void    _FORTRAN(slsas2)(char*, double*, double*, double&, double*, double*, 
                         int&, const int&, const int&, const int&);
}


// Four Node Quad element (Galerkin)

SloshQuadGal::SloshQuadGal(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
SloshQuadGal::clone()
{
 return new SloshQuadGal(*this);
}

void
SloshQuadGal::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
SloshQuadGal::getMass(CoordSet&)
{
 return 0.0;
}

FullSquareMatrix
SloshQuadGal::massMatrix(CoordSet &cs, double *d, int cmflg)
{
/*      Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);
        Node &nd4 = cs.getNode(nn[3]);

        int i;
        double x[4], y[4], mm[16];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        const int numgauss = 2;
        const int numdof   = 4;

        double capacitance = prop->rho*prop->Q*prop->eh;

// Consistent mass

        _FORTRAN(q4d1dofmas)(x, y, numgauss, mm, numdof);

// Lumped mass

        _FORTRAN(q4maslumpheat)(x,y,numgauss, mm, numdof);

        for (i=0;i<16;i++)
         d[i] = capacitance*mm[i];
*/
        FullSquareMatrix mass(4,d);
        
        mass.zero();
        return mass;
}

FullSquareMatrix
SloshQuadGal::stiffness(CoordSet &cs,double *Ks, int)
{

	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);
	Node &nd4 = cs.getNode(nn[3]);

        int i;
	double x[4], y[4], Kstiff[16];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;

        // Calculate stiffness matrix 

	const int numgauss = 2;
	const int numdof   = 4;
        
        double k = 1.;
        double h = prop ->eh;
        double c = 0.;

        _FORTRAN(thermquad3b)(x, y, k, c, h, Kstiff, &numgauss, &numdof);

        for (i=0;i<16;i++)
         Ks[i] = Kstiff[i];

        FullSquareMatrix ret(4, Ks);
        //ret.print();

        return ret;
}

int
SloshQuadGal::numNodes()
{
 	return 4;
}

int*
SloshQuadGal::nodes(int *p)
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
SloshQuadGal::numDofs()
{
 	return 4;
}

int*
SloshQuadGal::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Potential);
        p[1] = dsa.locate(nn[1],DofSet::Potential);
        p[2] = dsa.locate(nn[2],DofSet::Potential);
        p[3] = dsa.locate(nn[3],DofSet::Potential);

	return p;
}

void
SloshQuadGal::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 4, DofSet::Potential);
}

int
SloshQuadGal::getTopNumber()
{
  return 110; // 2
}

/*
void
SloshQuadGal::computeTemp(CoordSet&cs,
      State &state, double gp[2], double*tres)
{
// 4 is for the number of nodes, 2 is for temp and its derivative
// with respect to time
 double Temp[4][2];

 state.getTemp(nn[0], Temp[0], Temp[0]+1);
 state.getTemp(nn[1], Temp[1], Temp[1]+1);
 state.getTemp(nn[2], Temp[2], Temp[2]+1);
 state.getTemp(nn[3], Temp[3], Temp[3]+1);

// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[0][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[1][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[2][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[3][0]);

// tres[0] = temperature
// tres[1] = d(Temperature)/dt

 int j;
 for(j=0; j<2; ++j)
    tres[j] = (1-gp[0])*(1-gp[1])* Temp[0][j] +
              gp[0]*(1-gp[1])    * Temp[1][j] +
              gp[0]*gp[1]        * Temp[2][j] +
              (1-gp[0])*gp[1]    * Temp[3][j]; 
//     fprintf(stderr, "TEMP1 : %14.5e\n",tres[0]);
//     fprintf(stderr, "DTEMP1: %14.5e\n",tres[1]);
}

void
ThermQuadGal::getFlFlux(double gp[2], double *flF, double *tresF)
{
// Projects a fluid flux contained in flF[0] to all 4 nodes of quad
// Returns tresF
// fprintf(stderr, "Gauss Points %f %f\n ", gp[0], gp[1]);

   tresF[0]  = (1-gp[0])*(1-gp[1])* flF[0];
   tresF[1]  = gp[0]*(1-gp[1])    * flF[0];
   tresF[2]  = gp[0]*gp[1]        * flF[0];
   tresF[3]  = (1-gp[0])*gp[1]    * flF[0];

//   fprintf(stderr, "Fluxes are node 1: %f\n", tresF[0]);
//   fprintf(stderr, "Fluxes are node 2: %f\n", tresF[1]);
//   fprintf(stderr, "Fluxes are node 3: %f\n", tresF[2]);
//   fprintf(stderr, "Fluxes are node 4: %f\n", tresF[3]);
//   fflush(stderr); 
}

*/
void
SloshQuadGal::computeSloshDisp(Vector& fluidDispSlosh, CoordSet &cs, Vector& elPotSlosh,
                                   int hgInd)

{
// Fluid Displacement
// ... For Z-Direction :
// hgInd are defined in Eigen.C

   if(hgInd ==2) {
      fluidDispSlosh = 0.;
      return;
   }
 
   Node &nd1 = cs.getNode(nn[0]);
   Node &nd2 = cs.getNode(nn[1]);
   Node &nd3 = cs.getNode(nn[2]);
   Node &nd4 = cs.getNode(nn[3]);

   double x[4], y[4];

   x[0] = nd1.x; y[0] = nd1.y;
   x[1] = nd2.x; y[1] = nd2.y;
   x[2] = nd3.x; y[2] = nd3.y;
   x[3] = nd4.x; y[3] = nd4.y;

   int maxgus = 4;
   int maxstr = 6;
   int elm    = 1;
   int numel  = 1;

   //double k = prop ->k;
   double k = 1.;

   double elFluidDispSlosh[4][6];

   char ESCM[7] = "DIRECT"; // ... DIRECT Heat Fluxes CALCULATION
   //char ESCM[7] = "EXTRAP";   // ... EXTRAPOLATION FROM GAUSS POINTS

//   elPotSlosh.print();

   _FORTRAN(slsas2)(ESCM, x, y, k, elPotSlosh.data(),
                    (double*)elFluidDispSlosh, maxgus, maxstr, elm, numel);

   //..hgInd=0 : fluidDispSlosh-x     
   //..hgInd=1 : fluidDispSlosh-y    
   //..hgInd=2 : fluidDispSlosh-z = 0.

     fluidDispSlosh[0] = elFluidDispSlosh[0][hgInd];
     fluidDispSlosh[1] = elFluidDispSlosh[1][hgInd];
     fluidDispSlosh[2] = elFluidDispSlosh[2][hgInd];
     fluidDispSlosh[3] = elFluidDispSlosh[3][hgInd];
}

void
SloshQuadGal::computeSloshDispAll(Vector& fluidDispSlosh, CoordSet &cs, Vector& elPotSlosh)

{
// Fluid Displacement
// ... For Z-Direction :
// hgInd are defined in Eigen.C

//   if(hgInd ==2) {
//      fluidDispSlosh = 0.;
//      return;
//   }
 
   Node &nd1 = cs.getNode(nn[0]);
   Node &nd2 = cs.getNode(nn[1]);
   Node &nd3 = cs.getNode(nn[2]);
   Node &nd4 = cs.getNode(nn[3]);

   double x[4], y[4];

   x[0] = nd1.x; y[0] = nd1.y;
   x[1] = nd2.x; y[1] = nd2.y;
   x[2] = nd3.x; y[2] = nd3.y;
   x[3] = nd4.x; y[3] = nd4.y;

   int maxgus = 4;
   int maxstr = 6;
   int elm    = 1;
   int numel  = 1;

   double k = 1.;

   double elFluidDispSlosh[4][6];

   char ESCM[7] = "DIRECT"; // ... DIRECT Heat Fluxes CALCULATION
   //char ESCM[7] = "EXTRAP";   // ... EXTRAPOLATION FROM GAUSS POINTS

   _FORTRAN(slsas2)(ESCM, x, y, k, elPotSlosh.data(),
                    (double*)elFluidDispSlosh, maxgus, maxstr, elm, numel);

   //..hgInd=0 : fluidDispSlosh-x     
   //..hgInd=1 : fluidDispSlosh-y    
   //..hgInd=2 : fluidDispSlosh-z = 0.

     fluidDispSlosh[0] = elFluidDispSlosh[0][0];
     fluidDispSlosh[1] = elFluidDispSlosh[0][1];
     fluidDispSlosh[2] = 0.;
     fluidDispSlosh[3] = elFluidDispSlosh[1][0];
     fluidDispSlosh[4] = elFluidDispSlosh[1][1];
     fluidDispSlosh[5] = 0.;
     fluidDispSlosh[6] = elFluidDispSlosh[2][0];
     fluidDispSlosh[7] = elFluidDispSlosh[2][1];
     fluidDispSlosh[8] = 0.;
     fluidDispSlosh[9] = elFluidDispSlosh[3][0];
     fluidDispSlosh[10] = elFluidDispSlosh[3][1];
     fluidDispSlosh[11] = 0.;
}
