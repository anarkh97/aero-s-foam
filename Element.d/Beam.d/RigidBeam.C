#include	<stdio.h>
#include        <stdlib.h>
#include        <Element.d/Beam.d/RigidBeam.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include	<math.h>

RigidBeam::RigidBeam(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

void
RigidBeam::renum(int *table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
}


FullSquareMatrix
RigidBeam::massMatrix(CoordSet &cs,double *mel,int cmflg)
{

        FullSquareMatrix ret(18, mel);
  
        ret.zero();

        return ret;
}

FullSquareMatrix
RigidBeam::stiffness(CoordSet &cs, double *d, int flg)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        double dx = nd2.x - nd1.x;
        double dy = nd2.y - nd1.y;
        double dz = nd2.z - nd1.z;

        double length = sqrt( dx*dx + dy*dy + dz*dz );

        FullSquareMatrix ret(18,d);
        ret.zero();

        if(length == 0.0) {
          fprintf(stderr," *** WARNING: Rigid beam has zero length, nodes: %d %d\n",nn[0]+1,nn[1]+1);
          //fprintf(stderr," ... exiting fem program ...\n");
          //exit(-1);
          for(int i = 0; i < 6; ++i) {
            ret[i][12+i] = ret[12+i][i] = 1.0;
            ret[6+i][12+i] = ret[12+i][6+i] = -1.0;
          }
        }
	else {
          double c1[3], c2[3], c3[3];

          c1[0] = dx/length;
          c1[1] = dy/length;
          c1[2] = dz/length;
        
          // Translation in X
          ret[0][12] = ret[12][0] =  c1[0];
          ret[6][12] = ret[12][6] = -c1[0];

          // Translation in Y
          ret[1][12] = ret[12][1] =  c1[1];
          ret[7][12] = ret[12][7] = -c1[1];

          // Translation in Z
          ret[2][12] = ret[12][2] =  c1[2];
          ret[8][12] = ret[12][8] = -c1[2];

          // Rotation in X
          ret[3][13] = ret[13][3] = 1;
          ret[9][13] = ret[13][9] = -1; 

          // Rotation in Y
          ret[4][14] = ret[14][4] = 1;
          ret[10][14] = ret[14][10] = -1;

          // Rotation in Z
          ret[5][15] = ret[15][5] = 1;
          ret[11][15] = ret[15][11] = -1;

          // Relative Rotation in Y and Z

          double N1 = sqrt( c1[0]*c1[0] + c1[1]*c1[1] );
          double N2 = sqrt( c1[0]*c1[0] + c1[2]*c1[2] );
        
          if (N1 > N2) {
            c2[0] = -c1[1]/N1;
            c2[1] = c1[0]/N1;
            c2[2] = 0.0; 
          }
          else {
            c2[0] = c1[2]/N2;
            c2[1] = 0.0;
            c2[2] = -c1[0]/N2;
          }

          c3[0] = c1[1] * c2[2] - c1[2] * c2[1];
          c3[1] = c1[2] * c2[0] - c1[0] * c2[2];
          c3[2] = c1[0] * c2[1] - c1[1] * c2[0];
 
          ret[0][16] = ret[16][0] = c2[0];
          ret[1][16] = ret[16][1] = c2[1];
          ret[2][16] = ret[16][2] = c2[2];
          ret[3][16] = ret[16][3] = -c2[1] * dz + c2[2] * dy;
          ret[4][16] = ret[16][4] =  c2[0] * dz - c2[2] * dx;
          ret[5][16] = ret[16][5] = -c2[0] * dy + c2[1] * dx;
          ret[6][16] = ret[16][6] = -c2[0];
          ret[7][16] = ret[16][7] = -c2[1];
          ret[8][16] = ret[16][8] = -c2[2]; 

          ret[0][17] = ret[17][0] = c3[0];
          ret[1][17] = ret[17][1] = c3[1];
          ret[2][17] = ret[17][2] = c3[2];
          ret[3][17] = ret[17][3] = -c3[1] * dz + c3[2] * dy;
          ret[4][17] = ret[17][4] =  c3[0] * dz - c3[2] * dx;
          ret[5][17] = ret[17][5] = -c3[0] * dy + c3[1] * dx;
          ret[6][17] = ret[17][6] = -c3[0];
          ret[7][17] = ret[17][7] = -c3[1];
          ret[8][17] = ret[17][8] = -c3[2];
        }

        //ret *= 1.0e3; // XXXX
        return ret;
}

int
RigidBeam::numNodes()
{
        return 3;
}

int *
RigidBeam::nodes(int *p)
{
        if(p == 0) p = new int[3];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        return p;
}

int
RigidBeam::numDofs()
{
        return 18;
}

int *
RigidBeam::dofs(DofSetArray &dsa, int *p)
{
        if(p == 0) p = new int[18];

	dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);
        dsa.number(nn[2],DofSet::Lagrange1, p+12);
        dsa.number(nn[2],DofSet::Lagrange2, p+13);
        dsa.number(nn[2],DofSet::Lagrange3, p+14);
        dsa.number(nn[2],DofSet::Lagrange4, p+15);
        dsa.number(nn[2],DofSet::Lagrange5, p+16);
        dsa.number(nn[2],DofSet::Lagrange6, p+17);

        return p;
}

void
RigidBeam::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
        dsa.mark(nn[2], DofSet::Lagrange1);
        dsa.mark(nn[2], DofSet::Lagrange2);
        dsa.mark(nn[2], DofSet::Lagrange3);
        dsa.mark(nn[2], DofSet::Lagrange4);
        dsa.mark(nn[2], DofSet::Lagrange5);
        dsa.mark(nn[2], DofSet::Lagrange6);
}
