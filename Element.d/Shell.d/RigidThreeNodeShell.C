#include	<stdio.h>
#include        <math.h>
#include	<stdlib.h>
#include	<Element.d/Shell.d/RigidThreeNodeShell.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>


RigidThreeNodeShell::RigidThreeNodeShell(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}


void
RigidThreeNodeShell::renum(int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

FullSquareMatrix
RigidThreeNodeShell::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
           FullSquareMatrix ret(36,mel);

	   ret.zero();

           return ret;
}

FullSquareMatrix
RigidThreeNodeShell::stiffness(CoordSet &cs,double *d, int flg)
{

	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3], h[3];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

	h[0] = h[1] = h[2] = prop->eh;

        // Check for a zero thickness
        if(h[0] <= 0.0) {
           fprintf(stderr," *** ERROR: Zero shell thickness "
                          "(RigidThreeNodeShell.C) %d %d %d\n",
                          nn[0]+1, nn[1]+1, nn[2]+1);
           fprintf(stderr," ... exiting fem program ...\n");
           exit(-1);
        }

        FullSquareMatrix ret(30,d);
        ret.zero();

        // The Rigid Shell is considered as being a set of 2 rigid beams connected together
	// the first beam being between node 0 and node 1 and the second beam being between
	// node 0 and node 2

        double c1[3];
        double c2[3];
        double c3[3];
        int index = 18;


        for (int i = 1; i < 3; i++){

            double dx = x[i] - x[0];
            double dy = y[i] - y[0];
            double dz = z[i] - z[0];

            double length = sqrt( dx*dx + dy*dy + dz*dz );

            if (length == 0.0){
               fprintf(stderr, "ERROR: Zero distance between node %d and node %d\n", nn[0], nn[i]);
               fprintf(stderr, " ... exiting fem program ...\n");
               exit(-1);
            }

            c1[0] = dx/length;
            c1[1] = dy/length;
            c1[2] = dz/length;

            
            // Translation in X
            ret[6 * 0 + 0][index] = ret[index][6 * 0 + 0] = c1[0];
            ret[6 * 0 + 1][index] = ret[index][6 * 0 + 1] = c1[1];
            ret[6 * 0 + 2][index] = ret[index][6 * 0 + 2] = c1[2];
            ret[6 * i + 0][index] = ret[index][6 * i + 0] = -c1[0];
            ret[6 * i + 1][index] = ret[index][6 * i + 1] = -c1[1];
            ret[6 * i + 2][index] = ret[index][6 * i + 2] = -c1[2];

            // Rotation in X
            ret[6 * 0 + 3][index + 1] = ret[index + 1][6 * 0 + 3] = 1;
            ret[6 * i + 3][index + 1] = ret[index + 1][6 * i + 3] = -1;
             

            // Rotation in Y
            ret[6 * 0 + 4][index + 2] = ret[index + 2][6 * 0 + 4] = 1;
            ret[6 * i + 4][index + 2] = ret[index + 2][6 * i + 4] = -1;

        
            // Rotation in Z
            ret[6 * 0 + 5][index + 3] = ret[index + 3][6 * 0 + 5] = 1;
            ret[6 * i + 5][index + 3] = ret[index + 3][6 * i + 5] = -1;

            // Relative Rotation in Y and Z

            double N1 = sqrt( c1[0]*c1[0] + c1[1]*c1[1] );
            double N2 = sqrt( c1[0]*c1[0] + c1[2]*c1[2] );

            if (N1 > N2){
               c2[0] = -c1[1]/N1;
               c2[1] = c1[0]/N1; 
               c2[2] = 0.0;
            }
            else{
               c2[0] = c1[2]/N2;
               c2[1] = 0.0;
               c2[2] = -c1[0]/N2;
            }

            c3[0] = c1[1] * c2[2] - c1[2] * c2[1];
            c3[1] = c1[2] * c2[0] - c1[0] * c2[2];
            c3[2] = c1[0] * c2[1] - c1[1] * c2[0];

            ret[6 * 0 + 0][index + 4] = ret[index + 4][6 * 0 + 0] = c2[0];
            ret[6 * 0 + 1][index + 4] = ret[index + 4][6 * 0 + 1] = c2[1];
            ret[6 * 0 + 2][index + 4] = ret[index + 4][6 * 0 + 2] = c2[2];
            ret[6 * 0 + 3][index + 4] = ret[index + 4][6 * 0 + 3] = -c2[1] * dz + c2[2] * dy; 
            ret[6 * 0 + 4][index + 4] = ret[index + 4][6 * 0 + 4] =  c2[0] * dz - c2[2] * dx; 
            ret[6 * 0 + 5][index + 4] = ret[index + 4][6 * 0 + 5] = -c2[0] * dy + c2[1] * dx; 
            ret[6 * i + 0][index + 4] = ret[index + 4][6 * i + 0] = -c2[0]; 
            ret[6 * i + 1][index + 4] = ret[index + 4][6 * i + 1] = -c2[1]; 
            ret[6 * i + 2][index + 4] = ret[index + 4][6 * i + 2] = -c2[2]; 

            ret[6 * 0 + 0][index + 5] = ret[index + 5][6 * 0 + 0] = c3[0];
            ret[6 * 0 + 1][index + 5] = ret[index + 5][6 * 0 + 1] = c3[1];
            ret[6 * 0 + 2][index + 5] = ret[index + 5][6 * 0 + 2] = c3[2];
            ret[6 * 0 + 3][index + 5] = ret[index + 5][6 * 0 + 3] = -c3[1] * dz + c3[2] * dy; 
            ret[6 * 0 + 4][index + 5] = ret[index + 5][6 * 0 + 4] =  c3[0] * dz - c3[2] * dx; 
            ret[6 * 0 + 5][index + 5] = ret[index + 5][6 * 0 + 5] = -c3[0] * dy + c3[1] * dx; 
            ret[6 * i + 0][index + 5] = ret[index + 5][6 * i + 0] = -c3[0]; 
            ret[6 * i + 1][index + 5] = ret[index + 5][6 * i + 1] = -c3[1]; 
            ret[6 * i + 2][index + 5] = ret[index + 5][6 * i + 2] = -c3[2]; 
            
            index += 6;
        }      
        return ret;
}

int
RigidThreeNodeShell::numNodes()
{
 	return 4;
}

int*
RigidThreeNodeShell::nodes(int *p)
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
RigidThreeNodeShell::numDofs()
{
 	return 30;
}

int*
RigidThreeNodeShell::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[30];

        dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);
        dsa.number(nn[2],DofSet::XYZdisp | DofSet::XYZrot, p+12);
        dsa.number(nn[3],DofSet::Lagrange1, p+18);
        dsa.number(nn[3],DofSet::Lagrange2, p+19);
        dsa.number(nn[3],DofSet::Lagrange3, p+20);
        dsa.number(nn[3],DofSet::Lagrange4, p+21);
        dsa.number(nn[3],DofSet::Lagrange5, p+22);
        dsa.number(nn[3],DofSet::Lagrange6, p+23);
        dsa.number(nn[3],DofSet::Lagrange7, p+24);
        dsa.number(nn[3],DofSet::Lagrange8, p+25);
        dsa.number(nn[3],DofSet::Lagrange9, p+26);
        dsa.number(nn[3],DofSet::Lagrange10, p+27);
        dsa.number(nn[3],DofSet::Lagrange11, p+28);
        dsa.number(nn[3],DofSet::Lagrange12, p+29);
	
        return p;
}

void
RigidThreeNodeShell::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 3,  DofSet::XYZdisp | DofSet::XYZrot);
        for (int i=0; i<18; i++)
            dsa.mark(nn[3], DofSet::AllLagrange[i]);
}

int
RigidThreeNodeShell::getTopNumber()
{
 return 4;
}
