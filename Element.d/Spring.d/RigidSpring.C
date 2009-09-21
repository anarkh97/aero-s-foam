#include	<stdio.h>
#include        <stdlib.h>
#include        <Element.d/Spring.d/RigidSpring.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include	<math.h>

RigidSpring::RigidSpring(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

void
RigidSpring::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}


FullSquareMatrix
RigidSpring::massMatrix(CoordSet &cs,double *mel,int cmflg)
{

        FullSquareMatrix ret(18, mel);
  
        ret.zero();

        return ret;
}

FullSquareMatrix
RigidSpring::stiffness(CoordSet &cs, double *d, int flg)
{

        FullSquareMatrix ret(18,d);
        ret.zero();
        
        // Translation in X
        ret[0][12] = ret[12][0] =  1.0;
        ret[6][12] = ret[12][6] = -1.0; 

        // Translation in Y
        ret[1][13] = ret[13][1] =  1.0;
        ret[7][13] = ret[13][7] = -1.0;

        // Translation in Z
        ret[2][14] = ret[14][2] =  1.0;
        ret[8][14] = ret[14][8] = -1.0;

        // Rotation in X
        ret[3][15] = ret[15][3] =  1.0;
        ret[9][15] = ret[15][9] = -1.0;

        // Rotation in Y
        ret[4][16] = ret[16][4] = 1.0;
        ret[10][16] = ret[16][10] = -1.0;

        // Rotation in Z
        ret[5][17] = ret[17][5] = 1.0;
        ret[11][17]= ret[17][11] = -1.0;

        return ret;
}

int
RigidSpring::numNodes()
{
        return 3;
}

int *
RigidSpring::nodes(int *p)
{
        if(p == 0) p = new int[3];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        return p;
}

int
RigidSpring::numDofs()
{
        return 18;
}

int *
RigidSpring::dofs(DofSetArray &dsa, int *p)
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
RigidSpring::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
        dsa.mark(nn[2], DofSet::Lagrange1);
        dsa.mark(nn[2], DofSet::Lagrange2);
        dsa.mark(nn[2], DofSet::Lagrange3);
        dsa.mark(nn[2], DofSet::Lagrange4);
        dsa.mark(nn[2], DofSet::Lagrange5);
        dsa.mark(nn[2], DofSet::Lagrange6);
}
