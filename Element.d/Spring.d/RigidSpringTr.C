#include	<stdio.h>
#include        <stdlib.h>
#include        <Element.d/Spring.d/RigidSpringTr.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include	<math.h>

RigidSpringTr::RigidSpringTr(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

void
RigidSpringTr::renum(int *table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
}


FullSquareMatrix
RigidSpringTr::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        int index = 0;

        if (prop->A != 0.0)
            index++;
        if (prop->E != 0.0)
            index++;
        if (prop->nu != 0.0)
            index++;

        FullSquareMatrix ret(12 + index, mel);
  
        ret.zero();

        return ret;
}

FullSquareMatrix
RigidSpringTr::stiffness(CoordSet &cs, double *d, int flg)
{
        if (prop->A == 0.0 && prop->E == 0.0 && prop->nu == 0.0){
            fprintf(stderr, "*******************************************\n");
            fprintf(stderr, "*** ERROR: use the appropriate link element\n");
            fprintf(stderr, "***        exiting fem program \n");
            fprintf(stderr, "*******************************************\n");
            exit(-1);
        }

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );

        double c1[3];
       
        if (dx == 0.0)
           c1[0] = 1;
        else{
           c1[0] = dx/length;
        }

        if (dy == 0.0)
           c1[1] = 1;
        else{
           c1[1] = dy/length;
        }

        if (dz == 0.0)
           c1[2] = 1;
        else{
           c1[2] = dz/length;
        }
        

        int index = 0;
        if (prop->A != 0.0)
            index++;
        if (prop->E != 0.0)
            index++;
        if (prop->nu != 0.0)
            index++;

        FullSquareMatrix ret(12 + index , d);

        ret.zero();
       
        int indices = 0;
 
        // Translation in X
        if (prop->A != 0.0){
            indices++;
            ret[0][11 + indices] = ret[11 + indices][0] =  c1[0];
            ret[6][11 + indices] = ret[11 + indices][6] = -c1[0]; 
        }

        // Translation in Y
        if (prop->E != 0.0){
           indices++;
           ret[1][11 + indices] = ret[11 + indices][1] =  c1[1];
           ret[7][11 + indices] = ret[11 + indices][7] = -c1[1];
        }

        // Translation in Z
        if (prop->nu != 0.0){
           indices++;
           ret[2][11 + indices] = ret[11 + indices][2] =  c1[2];
           ret[8][11 + indices] = ret[11 + indices][8] = -c1[2];
        }
        
        return ret;
}

int
RigidSpringTr::numNodes()
{
        return 3;
}

int *
RigidSpringTr::nodes(int *p)
{
        if(p == 0) p = new int[3];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        return p;
}

int
RigidSpringTr::numDofs()
{
        int index = 0;
        if (prop->A != 0.0)
            index++;
        if (prop->E != 0.0)
            index++;
        if (prop->nu != 0.0)
            index++;

        return (12 + index);
}

int *
RigidSpringTr::dofs(DofSetArray &dsa, int *p)
{
        int index = 0;
        if (prop->A != 0.0)
            index++;
        if (prop->E != 0.0)
            index++;
        if (prop->nu != 0.0)
            index++;

        if(p == 0) p = new int[12 + index];

	dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);

        int indices = 0;
        if (prop->A != 0.0){
           dsa.number(nn[2],DofSet::Lagrange1, p+12+indices);
           indices++;
        }
        if (prop->E != 0.0){
           dsa.number(nn[2],DofSet::Lagrange2, p+12+indices);
           indices++;
        }
        if (prop->nu != 0.0){
           dsa.number(nn[2],DofSet::Lagrange3, p+12+indices);
        }

        return p;
}

void
RigidSpringTr::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
        if (prop->A != 0.0)
           dsa.mark(nn[2], DofSet::Lagrange1);
        if (prop->E != 0.0)
           dsa.mark(nn[2], DofSet::Lagrange2);
        if (prop->nu != 0.0)
           dsa.mark(nn[2], DofSet::Lagrange3);
}
