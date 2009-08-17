#include	<Element.d/FluidTriangle3.d/SloshTriangleFS.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Utils.d/dofset.h>
#include        <math.h>
#include        <Element.d/State.h>

SloshTriangleFS::SloshTriangleFS(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}


Element *
SloshTriangleFS::clone()
{
 return new SloshTriangleFS(*this);
}


void
SloshTriangleFS::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}


double
SloshTriangleFS::getMass(CoordSet& cs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);

  Vector r1(3), r2(3), r3(3);

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

  Vector v1(3), v2(3), v3(3), v4(3), v5(3);

  v1 = r2 - r1;
  v2 = r3 - r1;

  v3 = v1.cross(v2);

  double area = 0.5*v3.magnitude();
  double mass = area*prop->rho*prop->eh;

  return mass;
}

FullSquareMatrix
SloshTriangleFS::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
	double mass = getMass(cs);
	//double massPerNode = mass/3.0;

        FullSquareMatrix ret(3,mel);

	ret.zero();

//Calculate entries
	int i;
	int j;

// This is the LUMPED mass 

        //for(i=0; i<3; ++i)
        //  ret[i][i] = massPerNode;

// This is the CONSISTENT MASS

        for(i=0; i<3; ++i) {
           for(j=0; j<3; ++j)
              ret[i][j] = mass/12;
        }

        for(i=0; i<3; ++i)
           ret[i][i] = mass/6;
        
//      ret.print();
//      cout<<endl;

        return ret;
}

FullSquareMatrix
SloshTriangleFS::stiffness(CoordSet &cs, double *d, int flg)
{
/*
	double t  = prop->eh;
        double k = prop ->k;
*/
        FullSquareMatrix K(3,d);

        K.zero();

//      K.print();
//      cout<<endl;

        return K;
}

int
SloshTriangleFS::numNodes()
{
 	return 3;
}

int*
SloshTriangleFS::nodes(int *p)
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
SloshTriangleFS::numDofs()
{
 	return 3;
}

int*
SloshTriangleFS::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[3];

        p[0] = dsa.locate(nn[0],DofSet::Potential);
        p[1] = dsa.locate(nn[1],DofSet::Potential);
        p[2] = dsa.locate(nn[2],DofSet::Potential);

	return p;
}

void
SloshTriangleFS::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 3, DofSet::Potential);
}

int
SloshTriangleFS::getTopNumber()
{
  return 153;//4;
}

/*
void
SloshTriangleFS::computeTemp(CoordSet&cs,
      State &state, double gp[2], double*tres)
{
// 3 is for the number of nodes, 2 is for temp and its derivative
// with respect to time
 double Temp[3][2];

 state.getTemp(nn[0], Temp[0], Temp[0]+1);
 state.getTemp(nn[1], Temp[1], Temp[1]+1);
 state.getTemp(nn[2], Temp[2], Temp[2]+1);

//   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[0][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[1][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[2][0]);

// tres[0] = temperature
// tres[1] = d(Temperature)/dt

 int j;
 for(j=0; j<2; ++j)
    tres[j] = (1-gp[0]-gp[1])* Temp[0][j] +
                       gp[0] * Temp[1][j] +
                       gp[1] * Temp[2][j] ;

//     fprintf(stderr, "TEMP1 : %14.5e\n",tres[0]);
//     fprintf(stderr, "DTEMP1: %14.5e\n",tres[1]);
}

void
SloshTriangleFS::getFlFlux(double gp[2], double *flF, double *tresF)
{
// Projects a fluid flux contained in flF[0] to all 3 nodes of triangle
// Returns tresF
// fprintf(stderr, "Gauss Points %f %f\n ", gp[0], gp[1]);

   tresF[0]  = (1-gp[0]-gp[1])* flF[0];
   tresF[1]  = gp[0] * flF[0];
   tresF[2]  = gp[1] * flF[0];

//   fprintf(stderr, "Fluxes are node 1: %f\n", tresF[0]);
//   fprintf(stderr, "Fluxes are node 2: %f\n", tresF[1]);
//   fprintf(stderr, "Fluxes are node 3: %f\n", tresF[2]);
//   fflush(stderr);
}
*/
