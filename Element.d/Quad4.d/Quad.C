#include	<Utils.d/dbg_alloca.h>

#include	<Element.d/Quad4.d/Quad.h>
#include	<Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Utils.d/dofset.h>
#include        <Element.d/State.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void    _FORTRAN(getcmt)(double&, double&, double* );
void    _FORTRAN(quad4m)(double*, double*, double*, double*, const int&,
                         double*, const int&);
void    _FORTRAN(q4dmas)(const int&, double&, double*, const int&, double*,
                         double*, double*, double*,double*,int&,double&,int&);

}

Quad::Quad(int _numNodes, int* nodenums)
{
  numnod = _numNodes;
  nn = new int[numnod];
  int i;
  for(i=0; i<numnod; ++i)
    nn[i] = nodenums[i];
}

Element *
Quad::clone()
{
 return new Quad(*this);
}

void
Quad::renum(int *table)
{
  int i;
  for(i=0; i<numnod; ++i)
    nn[i] = table[nn[i]];
}

void
Quad::getVonMises(Vector& stress,Vector& weight,CoordSet &cs, 
                          Vector& elDisp, int strInd,int,double *,
			  double ylayer, double zlayer, int avgnum)
{
        
        // NOTE: SIGMAZZ, SIGMAYZ, SIGMAXZ, STRAINZZ, STRAINYZ & STRAINXZ
        // ARE EQUAL TO ZERO FOR A 4 NODE QUAD ELEMENT

        if(strInd == 2 || strInd == 4  || strInd == 5  ||
           strInd == 11 || strInd == 12   ) {

	   weight = 0.0; stress = 0.0;

           return;
        }

        // ELSE CALCULATE SIGMAXX, SIGMAYY, SIGMAXY, STRAINXX, STRAINYY, 
        // STRAINXY AND VONMISES STRESS

 	weight = 1.0;

        double c[9];
        double *x = new double[numnod];
        double *y = new double[numnod];

        int i;
        for(i=0; i<numnod; ++i) {
          x[i] = cs.getNode(nn[i]).x;
          y[i] = cs.getNode(nn[i]).y;
	}

       _FORTRAN(getcmt)(prop->E, prop->nu, c);
}

double
Quad::getMass(CoordSet& cs)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);
        Node &nd4 = cs.getNode(nn[3]);

        Vector r1(3), r2(3), r3(3), r4(3);

        r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = 0.0;
        r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = 0.0;
        r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = 0.0;
        r4[0] = nd4.x; r4[1] = nd4.y; r4[2] = 0.0;

        Vector v1(3), v2(3), v3(3), v4(3), v5(3);

        v1 = r2 - r1;
        v2 = r3 - r1;
        v3 = r4 - r1;

        v4 = v1.cross(v2);
        v5 = v2.cross(v3);

        double area = 0.5*(v4.magnitude() + v5.magnitude());
        double mass = area*prop->rho*prop->eh;
        return mass;
}

FullSquareMatrix
Quad::massMatrix(CoordSet &cs,double *mel,int cmflg)
{

        double *x = (double *) dbg_alloca(sizeof(double)*numnod);
        double *y = (double *) dbg_alloca(sizeof(double)*numnod);
        double *h = (double *) dbg_alloca(sizeof(double)*numnod);

	int i;
	for(i=0; i<numnod; ++i) {
	  x[i] = cs.getNode(nn[i]).x;
	  y[i] = cs.getNode(nn[i]).y;
	  h[i] = prop->eh;
	}

	double rho = prop->rho;
        double *gravityAcceleration = 0;
        double *grvfor = 0;
        double totmas = 0.0;

        int grvflg = 0;
        int masflg = 0;

	// will change with higher order elements
	const int numgauss = 2;
 	
       _FORTRAN(q4dmas)(numgauss,rho,(double*)mel,numDofs(),h,x,y,
                        gravityAcceleration,grvfor,grvflg,totmas,masflg);


        FullSquareMatrix ret(numDofs(),mel);

        return ret;
}

FullSquareMatrix
Quad::stiffness(CoordSet &cs, double *d, int flg)
{
	double *x = (double *) dbg_alloca(sizeof(double)*numnod);
	double *y = (double *) dbg_alloca(sizeof(double)*numnod);
	double *h = (double *) dbg_alloca(sizeof(double)*numnod);

	int i;
	for(i=0; i<numnod; ++i) {
          x[i] = cs.getNode(nn[i]).x;
          y[i] = cs.getNode(nn[i]).y;
	  h[i] = prop->eh;
	}

	double c[9];

	double e   = prop->E;
	double nu  = prop->rho;

	_FORTRAN(getcmt)(e, nu, c);

	// numgauss is dependent on element order i.e. linear, quadratic
	const int numgauss = 2;

	_FORTRAN(quad4m)(x, y, h, c, numgauss, (double *)d, numDofs());

        FullSquareMatrix ret(numDofs(),d);

        return ret;
}

int
Quad::numNodes()
{
 	return numnod;
}

int*
Quad::nodes(int *p)
{
 	if(p == 0) p = new int[numnod];
	int i;
        for(i=0; i<numnod; ++i)
 	  p[i] = nn[i];
	return p;
}

int
Quad::numDofs()
{
 	return (numnod*2);
}

int*
Quad::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[2*numnod];

	int i;
        for(i=0; i<numnod; ++i)
	  dsa.number(nn[i], DofSet::Xdisp | DofSet::Ydisp, p+i*2);

	return p;
}

void
Quad::markDofs(DofSetArray &dsa)
{
	int i;
  for(i=0; i<numnod; ++i)
    dsa.mark(nn[i], DofSet::Xdisp | DofSet::Ydisp);
}

int
Quad::getTopNumber()
{
 return 2;
}

