#include	<Element.d/Solid.d/RigidSolid6Dof.h>
#include	<Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <math.h>
#include        <stdio.h>
#include        <stdlib.h>


RigidSolid6Dof::RigidSolid6Dof(int numnode, int *nodenums)
{
        if (numnode < 3){
           fprintf(stderr, " ERROR: minimum number of nodes should be 3\n");
           fprintf(stderr," ... exiting fem program ...\n");
           exit(-1);
        }
        nn = new int[numnode + (numnode-1)];
        nnodes = numnode;
        int i;
        for (i = 0; i < nnodes; ++i)
            nn[i] = nodenums[i];
        for (; i < 2*nnodes-1; ++i)
            nn[i] = -1;
}
        

void
RigidSolid6Dof::renum(int *table)
{
        int index = nnodes + (nnodes-1);
        for (int i = 0; i<index; ++i)
             nn[i] = table[nn[i]];
}

FullSquareMatrix
RigidSolid6Dof::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        FullSquareMatrix ret(6 * (nnodes + (nnodes - 1)), mel);

        ret.zero();

        return ret;

}

FullSquareMatrix
RigidSolid6Dof::stiffness(CoordSet &cs, double *d, int flg)
{
        double x[2], y[2], z[2];
	
        // Node 1
        
        Node &Nd1 = cs.getNode(nn[0]);
        x[0] = Nd1.x; y[0] = Nd1.y; z[0] = Nd1.z;        


        // Initialize the Stiffness Matrix

        FullSquareMatrix ret(6 * (nnodes + (nnodes - 1)), d);
        
        ret.zero();


        // Restraining the rest of the nodes by connecting Node 1 to every 
        // other node by a rigid beam

        int index = 6 * nnodes;
	
	for (int i = 1; i < nnodes; ++i){

           Node &nd2 = cs.getNode(nn[i]);
           x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;        

           double dx = x[1] - x[0];
           double dy = y[1] - y[0];
           double dz = z[1] - z[0];

           double length = sqrt( dx*dx + dy*dy + dz*dz );

           if(length == 0.0) {
             fprintf(stderr,"ERROR: Beam has zero length %d %d\n",nn[0],nn[1]);
             fprintf(stderr," ... exiting fem program ...\n");
             exit(-1);
           }

           double c1[3];
           double c2[3];
           double c3[3];

           c1[0] = dx/length;
           c1[1] = dy/length;
           c1[2] = dz/length;

        
           // Translation in X
           ret[6 * 0 + 0][index] = ret[index][6 * 0 + 0] =  c1[0];
           ret[6 * 0 + 1][index] = ret[index][6 * 0 + 1] =  c1[1];
           ret[6 * 0 + 2][index] = ret[index][6 * 0 + 2] =  c1[2];
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
RigidSolid6Dof::numNodes()
{
 	return (nnodes + (nnodes - 1));
}

int*
RigidSolid6Dof::nodes(int *p)
{
        int ind = nnodes + (nnodes - 1);
 	if(p == 0) p = new int[ind];
        for (int i = 0; i < ind; ++i)
           p[i] = nn[i];

	return p;
}

int
RigidSolid6Dof::numDofs()
{
        int ind = 6 * (nnodes + (nnodes - 1));
 	return ind;
}

int*
RigidSolid6Dof::dofs(DofSetArray &dsa, int *p)
{
        int ind1 = 6 * (nnodes + (nnodes - 1));
 	if(p == 0) p = new int[ind1];

        for (int i = 0; i < nnodes; ++i){
            dsa.number(nn[i], DofSet::XYZdisp | DofSet::XYZrot, p+6*i);
        }

        int ind2 = nnodes + (nnodes - 1);
        for (int j = nnodes; j < ind2; ++j){
            dsa.number(nn[j],DofSet::Lagrange1, p+6*j); 
            dsa.number(nn[j],DofSet::Lagrange2, p+6*j+1); 
            dsa.number(nn[j],DofSet::Lagrange3, p+6*j+2); 
	    dsa.number(nn[j],DofSet::Lagrange4, p+6*j+3);
	    dsa.number(nn[j],DofSet::Lagrange5, p+6*j+4);
	    dsa.number(nn[j],DofSet::Lagrange6, p+6*j+5);
        }

	return p;
}

void
RigidSolid6Dof::markDofs(DofSetArray &dsa)
{
	dsa.mark(nn, nnodes, DofSet::XYZdisp | DofSet::XYZrot);
	dsa.mark(nn+nnodes, nnodes-1, DofSet::Lagrange1 | DofSet::Lagrange2 | DofSet::Lagrange3 |
	                              DofSet::Lagrange4 | DofSet::Lagrange5 | DofSet::Lagrange6 );
}
