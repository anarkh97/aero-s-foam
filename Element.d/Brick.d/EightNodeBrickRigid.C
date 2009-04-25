#include	<Element.d/Brick.d/EightNodeBrickRigid.h>
#include	<Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <math.h>
#include        <stdio.h>
#include        <stdlib.h>


EightNodeBrickRigid::EightNodeBrickRigid(int* nodenums)
{
        for (int i=0; i<8; i++)
            nn[i] = nodenums[i];
}


void
EightNodeBrickRigid::renum(int *table)
{
        for (int i=0; i<9; i++)
            nn[i] = table[nn[i]];
}

FullSquareMatrix
EightNodeBrickRigid::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        FullSquareMatrix ret(42, mel);

        ret.zero();

        return ret;

}

FullSquareMatrix
EightNodeBrickRigid::stiffness(CoordSet &cs, double *d, int flg)
{
	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);
	Node &nd4 = cs.getNode(nn[3]);
	Node &nd5 = cs.getNode(nn[4]);
	Node &nd6 = cs.getNode(nn[5]);
	Node &nd7 = cs.getNode(nn[6]);
	Node &nd8 = cs.getNode(nn[7]);

	double x[8], y[8], z[8];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z; 
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
	x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
	x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
	x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

        int indices[18][2] = {{0,1},
                              {1,2},
                              {2,3}, 
                              {4,5}, 
                              {5,6}, 
                              {6,7}, 
                              {0,2}, 
                              {4,6}, 
                              {0,3}, 
                              {4,7}, 
                              {0,4}, 
                              {1,5}, 
                              {2,6}, 
                              {3,7}, 
                              {0,5}, 
                              {1,6}, 
                              {2,7}, 
                              {0,7}}; 

        double c1[3];

        double dx;
        double dy;
        double dz;
        double length;

        int n1;
        int n2;

        FullSquareMatrix ret(42,d);
        ret.zero();
        
        int i;
        for (i=0; i<18; i++){

            n1 = indices[i][0];
            n2 = indices[i][1];

            dx = x[n2] - x[n1];
            dy = y[n2] - y[n1];
            dz = z[n2] - z[n1]; 

            length = sqrt( dx*dx + dy*dy + dz*dz );

            if(length == 0.0) {
            fprintf(stderr,"ERROR: zero distance %d %d\n",nn[n1],nn[n2]);
            fprintf(stderr," ... exiting fem program ...\n");
            exit(-1);
            }

            c1[0] = dx/length;
            c1[1] = dy/length;
            c1[2] = dz/length;

            ret[3*n1 + 0][i + 24] = ret[i + 24][3*n1 + 0] =  c1[0];
            ret[3*n1 + 1][i + 24] = ret[i + 24][3*n1 + 1] =  c1[1];
            ret[3*n1 + 2][i + 24] = ret[i + 24][3*n1 + 2] =  c1[2];
            ret[3*n2 + 0][i + 24] = ret[i + 24][3*n2 + 0] =  -c1[0];
            ret[3*n2 + 1][i + 24] = ret[i + 24][3*n2 + 1] =  -c1[1];
            ret[3*n2 + 2][i + 24] = ret[i + 24][3*n2 + 2] =  -c1[2];

        } 
        return ret;
}

int
EightNodeBrickRigid::numNodes()
{
 	return 9;
}

int*
EightNodeBrickRigid::nodes(int *p)
{
 	if(p == 0) p = new int[9];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
 	p[4] = nn[4];
 	p[5] = nn[5];
 	p[6] = nn[6];
 	p[7] = nn[7];
 	p[8] = nn[8];
	return p;
}

int
EightNodeBrickRigid::numDofs()
{
 	return 42;
}

int*
EightNodeBrickRigid::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[42];

	dsa.number(nn[0],DofSet::XYZdisp, p  );
        dsa.number(nn[1],DofSet::XYZdisp, p+3);
        dsa.number(nn[2],DofSet::XYZdisp, p+6);
        dsa.number(nn[3],DofSet::XYZdisp, p+9);
        dsa.number(nn[4],DofSet::XYZdisp, p+12);
        dsa.number(nn[5],DofSet::XYZdisp, p+15);
        dsa.number(nn[6],DofSet::XYZdisp, p+18);
        dsa.number(nn[7],DofSet::XYZdisp, p+21);
        dsa.number(nn[8],DofSet::Lagrange1, p+24);
        dsa.number(nn[8],DofSet::Lagrange2, p+25);
        dsa.number(nn[8],DofSet::Lagrange3, p+26);
        dsa.number(nn[8],DofSet::Lagrange4, p+27);
        dsa.number(nn[8],DofSet::Lagrange5, p+28);
        dsa.number(nn[8],DofSet::Lagrange6, p+29);
        dsa.number(nn[8],DofSet::Lagrange7, p+30);
        dsa.number(nn[8],DofSet::Lagrange8, p+31);
        dsa.number(nn[8],DofSet::Lagrange9, p+32);
        dsa.number(nn[8],DofSet::Lagrange10, p+33);
        dsa.number(nn[8],DofSet::Lagrange11, p+34);
        dsa.number(nn[8],DofSet::Lagrange12, p+35);
        dsa.number(nn[8],DofSet::Lagrange13, p+36);
        dsa.number(nn[8],DofSet::Lagrange14, p+37);
        dsa.number(nn[8],DofSet::Lagrange15, p+38);
        dsa.number(nn[8],DofSet::Lagrange16, p+39);
        dsa.number(nn[8],DofSet::Lagrange17, p+40);
        dsa.number(nn[8],DofSet::Lagrange18, p+41);

	return p;
}

void
EightNodeBrickRigid::markDofs(DofSetArray &dsa)
{
	dsa.mark(nn, 8, DofSet::XYZdisp);
        for (int i=0; i<18; i++)
            dsa.mark(nn[8], DofSet::AllLagrange[i]);
}

int
EightNodeBrickRigid::getTopNumber()
{
 return 117;
}

