#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <Element.d/Truss.d/TwoNodeTrussRigid.h>
#include <Corotational.d/BarCorotator.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

TwoNodeTrussRigid::TwoNodeTrussRigid(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}


void
TwoNodeTrussRigid::renum(int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
}


FullSquareMatrix
TwoNodeTrussRigid::massMatrix(CoordSet &cs, double *mel, int cmflg)
{

        FullSquareMatrix elementMassMatrix(7,mel);

// zero the element mass matrix

	elementMassMatrix.zero();

        return elementMassMatrix;
}


FullSquareMatrix
TwoNodeTrussRigid::stiffness(CoordSet &cs, double *k, int flg)
{
        Node &nd1 = cs.getNode( nn[0] );
        Node &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt( dx*dx + dy*dy + dz*dz );

        if(length == 0.0) {
          fprintf(stderr,"ERROR: truss has zero length %d %d\n",nn[0],nn[1]);
          fprintf(stderr," ... exiting fem program ...\n");
          exit(-1);
        }

        double c1[3];
        
        c1[0] = dx/length;
        c1[1] = dy/length;
        c1[2] = dz/length;

        FullSquareMatrix ret(7,k);
        ret.zero();

        ret[0][6] = ret[6][0] = c1[0];
        ret[1][6] = ret[6][1] = c1[1];
        ret[2][6] = ret[6][2] = c1[2]; 
        ret[3][6] = ret[6][3] = -c1[0];
        ret[4][6] = ret[6][4] = -c1[1];
        ret[5][6] = ret[6][5] = -c1[2];

        return ret;
}

int
TwoNodeTrussRigid::numNodes()
{
        return 3;
}

int *
TwoNodeTrussRigid::nodes(int *p)
{
        if(p == 0) p = new int[3];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        return p;
}

int
TwoNodeTrussRigid::numDofs()
{
        return 7;
}

int *
TwoNodeTrussRigid::dofs(DofSetArray &dsa, int *p)
{
        if(p == 0) p = new int[7];

        dsa.number(nn[0],DofSet::XYZdisp, p  );
        dsa.number(nn[1],DofSet::XYZdisp, p+3);
        dsa.number(nn[2],DofSet::Lagrange1, p+6);

        return p;
}

void
TwoNodeTrussRigid::markDofs( DofSetArray &dsa )
{
        dsa.mark(nn, 2, DofSet::XYZdisp);
        dsa.mark(nn[2], DofSet::Lagrange1);
}

// HB: (01/04) add geometric stiffness contribution
// WARNING: NOT VALIDATED !!!
void
TwoNodeTrussRigid::getStiffAndForce(GeomState &gState, CoordSet &cs,
                                  FullSquareMatrix &Ktan, double *f)
{
  //cerr << " ... Get in TwoNodeTrussRigid::getStiffAndForce() " << endl; 
  // Get Nodes original coordinates
  Node &node1 = cs.getNode(nn[0]);
  Node &node2 = cs.getNode(nn[1]);
  // Get Nodes current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];
  NodeState lagrangeNode = gState[nn[2]];
  
  // Compute the original length 
  double l0x = node1.x-node2.x;
  double l0y = node1.y-node2.y;
  double l0z = node1.z-node2.z; 
  double l0 = sqrt(l0x*l0x + l0y*l0y + l0z*l0z);

  // compute current length
  double lx = ns1.x-ns2.x;
  double ly = ns1.y-ns2.y;
  double lz = ns1.z-ns2.z;
  double l = sqrt(lx*lx + ly*ly + lz*lz);

  // the error:
  double error = l - l0;
  //cerr  << " l0 = " << l0 <<", l = " << l <<endl;
  double coef = 1.0/l;

  double lambda = lagrangeNode.x;
  //printf("lagrangeNode.x = %f \n", lagrangeNode.x);

  // (internal) force vector
  f[0] = coef*lx*lambda;
  f[1] = coef*ly*lambda;
  f[2] = coef*lz*lambda;
  f[3] = -f[0];
  f[4] = -f[1];
  f[5] = -f[2];
  f[6] = error;
 
  //for(int i=0;i<7;i++)
  //  cerr << " f[" << i << "] = " << f[i] << endl;
 
  // tangent stiffness
  Ktan.zero();
  // geometric stiffness constribution 
  bool AddGeoStiff = true; // for testing purpose only
  if(AddGeoStiff){
  double coef1 = lambda/(l*l*l);
  double l2 = l*l;
  Ktan[0][0] = Ktan[3][3] = coef1*(l2-lx*lx); 
  Ktan[1][1] = Ktan[4][4] = coef1*(l2-ly*ly); 
  Ktan[2][2] = Ktan[5][5] = coef1*(l2-lz*lz); 
 
  Ktan[0][1] = Ktan[1][0] = Ktan[3][4] = Ktan[4][3] = -coef1*lx*ly;
  Ktan[0][2] = Ktan[2][0] = Ktan[3][5] = Ktan[5][3] = -coef1*lx*lz;
  Ktan[1][2] = Ktan[2][1] = Ktan[4][5] = Ktan[5][4] = -coef1*ly*lz;

  Ktan[0][3] = Ktan[3][0] = -Ktan[0][0];
  Ktan[1][4] = Ktan[4][1] = -Ktan[1][1];
  Ktan[2][5] = Ktan[5][2] = -Ktan[2][2];

  Ktan[0][4] = Ktan[4][0] = Ktan[1][3] = Ktan[3][1] = -Ktan[0][1];
  Ktan[0][5] = Ktan[5][0] = Ktan[2][3] = Ktan[3][2] = -Ktan[0][2];
  Ktan[1][5] = Ktan[5][1] = Ktan[2][4] = Ktan[4][2] = -Ktan[1][2];
  } 
  // constrain (Lagrange multipler) constribution 
  Ktan[6][0] = Ktan[0][6] = coef*lx;
  Ktan[6][1] = Ktan[1][6] = coef*ly;
  Ktan[6][2] = Ktan[2][6] = coef*lz;
  Ktan[6][3] = Ktan[3][6] = -Ktan[0][6];
  Ktan[6][4] = Ktan[4][6] = -Ktan[1][6];
  Ktan[6][5] = Ktan[5][6] = -Ktan[2][6];
}
