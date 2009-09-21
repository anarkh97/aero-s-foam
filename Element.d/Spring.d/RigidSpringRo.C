#include	<stdio.h>
#include        <stdlib.h>
#include        <Element.d/Spring.d/RigidSpringRo.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include	<math.h>

RigidSpringRo::RigidSpringRo(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

void
RigidSpringRo::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}


FullSquareMatrix
RigidSpringRo::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        int indices = 0;

        if (prop->A != 0.0)
            indices++;
        if (prop->E != 0.0)
            indices++;
        if (prop->nu != 0.0)
            indices++;

        FullSquareMatrix ret(12 + indices, mel);
  
        ret.zero();

        return ret;
}

FullSquareMatrix
RigidSpringRo::stiffness(CoordSet &cs, double *d, int flg)
{
        if (prop->A == 0.0 && prop->E == 0.0 && prop->nu == 0.0){
            fprintf(stderr, "*******************************************\n");
            fprintf(stderr, "*** ERROR: use the appropriate link element\n");
            fprintf(stderr, "***        exiting fem program \n");
            fprintf(stderr, "*******************************************\n");
            exit(-1);
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
        
        // Rotation in X
        if (prop->A != 0.0) {
            indices++; 
            ret[3][11 + indices] = ret[11 + indices][3] =  1.0;
            ret[9][11 + indices] = ret[11 + indices][9] = -1.0; 
        }

        // Rotation in Y
        if (prop->E != 0.0) {
            indices++;
            ret[4][11 + indices] = ret[11 + indices][4] = 1.0;
            ret[10][11 + indices] = ret[11 + indices][10] = -1.0;
        }

        // Rotation in Z
        if (prop->nu != 0.0) {
            indices++;
            ret[5][11 + indices] = ret[11 + indices][5] = 1.0;
            ret[11][11 + indices]= ret[11 + indices][11] = -1.0;
        }

        return ret;
}

int
RigidSpringRo::numNodes()
{
        return 3;
}

int *
RigidSpringRo::nodes(int *p)
{
        if(p == 0) p = new int[3];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        return p;
}

int
RigidSpringRo::numDofs()
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
RigidSpringRo::dofs(DofSetArray &dsa, int *p)
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
RigidSpringRo::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
        if (prop->A != 0.0)
           dsa.mark(nn[2], DofSet::Lagrange1);
        if (prop->E != 0.0)
           dsa.mark(nn[2], DofSet::Lagrange2);
        if (prop->nu != 0.0)
           dsa.mark(nn[2], DofSet::Lagrange3);
}
