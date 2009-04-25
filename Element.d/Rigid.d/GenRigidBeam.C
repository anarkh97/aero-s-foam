#include <stdio.h>
#include <stdlib.h>
#include <Element.d/Rigid.d/GenRigidBeam.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <math.h>

GenRigidBeam::GenRigidBeam(int* nodenums)
{
  // for this constructor, nodenums[0] is the master node, nodenums[1] is the list of constrained dofs,
  // and nodenums[2] is the slave node
  nn[0] = nodenums[0];
  nn[1] = nodenums[2];
  sprintf(cdofs,"%d\0",nodenums[1]+1);
  numcdofs = 0;
  for(int i=0; i<8; ++i) if(cdofs[i] != '\0') numcdofs++;
}

GenRigidBeam::GenRigidBeam(int* _nn, int _numcdofs, char* _cdofs)
{
  nn[0] = _nn[0];
  nn[1] = _nn[1];
  numcdofs = _numcdofs;
  for(int i=0; i<numcdofs; ++i) cdofs[i] = _cdofs[i];
}

void
GenRigidBeam::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

FullSquareMatrix
GenRigidBeam::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
  FullSquareMatrix ret(12+numcdofs, mel);
  ret.zero();
  return ret;
}

FullSquareMatrix
GenRigidBeam::stiffness(CoordSet &cs, double *d, int flg)
{
  Node &node1 = cs.getNode(nn[0]);
  Node &node2 = cs.getNode(nn[1]);
  double l0[3];
  l0[0] = node1.x-node2.x;
  l0[1] = node1.y-node2.y;
  l0[2] = node1.z-node2.z;

  FullSquareMatrix ret(12+numcdofs,d);
  ret.zero();

  for(int i=0; i<numcdofs; ++i) {
    switch(cdofs[i]) {
      case '1' : { // x translation
        // Ux(2)-Ux(1) + Ly*Rz(1)-Lz*Ry(1)
        ret[0][12+i] = ret[12+i][0] = -1.0;
        ret[6][12+i] = ret[12+i][6] = 1.0;
        ret[4][12+i] = ret[12+i][4] = -l0[2];
        ret[5][12+i] = ret[12+i][5] = l0[1];
      }
      break;
      case '2' : { // y translation
        // Ux(2)-Ux(1) + Ly*Rz(1)-Lz*Ry(1)
        ret[1][12+i] = ret[12+i][1] = -1.0;
        ret[7][12+i] = ret[12+i][7] = 1.0;
        ret[3][12+i] = ret[12+i][3] = l0[2];
        ret[5][12+i] = ret[12+i][5] = -l0[0];
      }
      break;
      case '3' : { // z translation
        // Uz(2)-Uz(1) + Lx*Ry(1)-Ly*Rx(1)
        ret[2][12+i] = ret[12+i][2] = -1.0;
        ret[8][12+i] = ret[12+i][8] = 1.0;
        ret[3][12+i] = ret[12+i][3] = -l0[1];
        ret[4][12+i] = ret[12+i][4] = l0[0];
      }
      break;
      case '4' : { // x rotation
        ret[3][12+i] = ret[12+i][3] = 1;
        ret[9][12+i] = ret[12+i][9] = -1; 
      }
      break;
      case '5' : { // y rotation
        ret[4][12+i] = ret[12+i][4] = 1;
        ret[10][12+i] = ret[12+i][10] = -1;
      }
      break;
      case '6' : { // z rotation
        ret[5][12+i] = ret[12+i][5] = 1;
        ret[11][12+i] = ret[12+i][11] = -1;
      }
      break;
      default:
        cerr << " *** WARNING: dof " << cdofs[i] << " is not supported in GenRigidBeam element \n";
      break;
    }
  }
  return ret;
}

int
GenRigidBeam::numNodes()
{
  return 3;
}

int *
GenRigidBeam::nodes(int *p)
{
  if(p == 0) p = new int[3];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  return p;
}

int
GenRigidBeam::numDofs()
{
  return 12+numcdofs;
}

int *
GenRigidBeam::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[12+numcdofs];
  int ldofs[6] = { DofSet::Lagrange1, DofSet::Lagrange2, DofSet::Lagrange3,
                   DofSet::Lagrange4, DofSet::Lagrange5, DofSet::Lagrange6 };

  dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
  dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);
  for(int i=0; i<numcdofs; ++i) 
    dsa.number(nn[2],ldofs[i], p+12+i);

  return p;
}

void
GenRigidBeam::markDofs(DofSetArray &dsa)
{
  int ldofs[6] = { DofSet::Lagrange1, DofSet::Lagrange2, DofSet::Lagrange3,
                   DofSet::Lagrange4, DofSet::Lagrange5, DofSet::Lagrange6 };

  dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
  for(int i=0; i<numcdofs; ++i) dsa.mark(nn[2], ldofs[i]);
}
