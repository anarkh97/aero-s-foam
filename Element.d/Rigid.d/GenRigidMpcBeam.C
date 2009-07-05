#include <Element.d/Rigid.d/GenRigidMpcBeam.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Domain.h>
extern Domain *domain;

GenRigidMpcBeam::GenRigidMpcBeam(int* nodenums)
{
  // for this constructor, nodenums[0] is the master node, nodenums[1] is the list of constrained dofs,
  // and nodenums[2] is the slave node
  nn[0] = nodenums[0];
  nn[1] = nodenums[2];
  sprintf(cdofs,"%d\0",nodenums[1]+1);
  numcdofs = 0;
  for(int i=0; i<8; ++i)
	  if(cdofs[i] != '\0') {
		  numcdofs++;
		  activeDofs |= 1 << (cdofs[i]-'1');
	  }
}

GenRigidMpcBeam::GenRigidMpcBeam(int* _nn, int _numcdofs, char* _cdofs)
{
  nn[0] = _nn[0];
  nn[1] = _nn[1];
  numcdofs = _numcdofs;
  for(int i=0; i<numcdofs; ++i) {
	  cdofs[i] = _cdofs[i];
	  activeDofs |= 1 << (cdofs[i]-'1');
  }
}

void 
GenRigidMpcBeam::computeMPCs(CoordSet &cs, int &lmpcnum)
{
  double rhs = 0.0;
  int i;
  for(i=0; i<numcdofs; ++i) {
    lmpcnum++;
    mpc[i] = new LMPCons(lmpcnum, rhs);
  }

  Node &node1 = cs.getNode(nn[0]);
  Node &node2 = cs.getNode(nn[1]);
  double l0[3];
  l0[0] = node1.x-node2.x;
  l0[1] = node1.y-node2.y;
  l0[2] = node1.z-node2.z;

  for(i=0; i<numcdofs; ++i) {
    switch(cdofs[i]) {
      case '1' : { // x translation
        activeDofs |= DofSet::Xdisp;
        // Ux(2)-Ux(1) + Ly*Rz(1)-Lz*Ry(1)
        LMPCTerm term01(nn[1],0,1.0);
        mpc[i]->addterm(&term01);
        LMPCTerm term02(nn[0],0,-1.0);
        mpc[i]->addterm(&term02);
        LMPCTerm term03(nn[0],5,l0[1]);
        mpc[i]->addterm(&term03);
        LMPCTerm term04(nn[0],4,-l0[2]);
        mpc[i]->addterm(&term04);
      }
      break;
      case '2' : { // y translation
        activeDofs |= DofSet::Ydisp;
        // Uy(2)-Uy(1) + Lz*Rx(1)-Lx*Rz(1)
        LMPCTerm term11(nn[1],1,1.0);
        mpc[i]->addterm(&term11);
        LMPCTerm term12(nn[0],1,-1.0);
        mpc[i]->addterm(&term12);
        LMPCTerm term13(nn[0],3,l0[2]);
        mpc[i]->addterm(&term13);
        LMPCTerm term14(nn[0],5,-l0[0]);
        mpc[i]->addterm(&term14);
      }
      break;
      case '3' : { // z translation
        activeDofs |= DofSet::Zdisp;
        // Uz(2)-Uz(1) + Lx*Ry(1)-Ly*Rx(1)
        LMPCTerm term21(nn[1],2,1.0);
        mpc[i]->addterm(&term21);
        LMPCTerm term22(nn[0],2,-1.0);
        mpc[i]->addterm(&term22);
        LMPCTerm term23(nn[0],4,l0[0]);
        mpc[i]->addterm(&term23);
        LMPCTerm term24(nn[0],3,-l0[1]);
        mpc[i]->addterm(&term24);
      }
      break;
      case '4' : { // x rotation
        activeDofs |= DofSet::Xrot;
        LMPCTerm term1(nn[0],3,1.0);
        mpc[i]->addterm(&term1);
        LMPCTerm term2(nn[1],3,-1.0);
        mpc[i]->addterm(&term2);
      }
      break;
      case '5' : { // y rotation
        activeDofs |= DofSet::Yrot;
        LMPCTerm term1(nn[0],4,1.0);
        mpc[i]->addterm(&term1);
        LMPCTerm term2(nn[1],4,-1.0);
        mpc[i]->addterm(&term2);
      }
      break;
      case '6' : { // z rotation
        activeDofs |= DofSet::Zrot;
        LMPCTerm term1(nn[0],5,1.0);
        mpc[i]->addterm(&term1);
        LMPCTerm term2(nn[1],5,-1.0);
        mpc[i]->addterm(&term2);
      }
      break;
      case '7' : { // temperature
        activeDofs |= DofSet::Temp;
        LMPCTerm term1(nn[0],6,1.0);
        mpc[i]->addterm(&term1);
        LMPCTerm term2(nn[1],6,-1.0);
        mpc[i]->addterm(&term2);
      }
      break;
      case '8' : { // pressure
        activeDofs |= DofSet::Helm;
        LMPCTerm term1(nn[0],7,1.0);
        mpc[i]->addterm(&term1);
        LMPCTerm term2(nn[1],7,-1.0);
        mpc[i]->addterm(&term2);
      }
      break;
      default:
        cerr << " *** WARNING: dof " << cdofs[i] << " is not supported in GenRigidMpcBeam element \n";
      break;
    }
  }

  for(i=0; i<numcdofs; ++i) 
    glMpcNb[i] = domain->addLMPC(mpc[i]);
}

void
GenRigidMpcBeam::updateMPCs(GeomState &gState)
{
  cerr << " *** WARNING: GenRigidMpcBeam::updateMPCs(...) is not implemented \n";
}

void
GenRigidMpcBeam::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

FullSquareMatrix
GenRigidMpcBeam::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix elementMassMatrix(2*numcdofs, mel);
  elementMassMatrix.zero();
  return elementMassMatrix;
}

FullSquareMatrix
GenRigidMpcBeam::stiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(2*numcdofs, k);
  ret.zero();
  return ret;
}

int
GenRigidMpcBeam::numNodes()
{
  return 2;
}

int *
GenRigidMpcBeam::nodes(int *p)
{
  if(p == 0) p = new int[2];
  p[0] = nn[0];
  p[1] = nn[1];
  return p;
}

int
GenRigidMpcBeam::numDofs()
{
  return 2*numcdofs;
}

int *
GenRigidMpcBeam::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[2*numcdofs]; 
  //dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot | DofSet::Temp | DofSet::Helm, p  );
  //dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot | DofSet::Temp | DofSet::Helm, p+numcdofs);
  //dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p);
  //dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+numcdofs);
  cerr << "activeDofs.count() = " << activeDofs.count() << endl;
  dsa.number(nn[0], activeDofs.list(), p);
  dsa.number(nn[1], activeDofs.list(), p+numcdofs);
  return p;
}

void
GenRigidMpcBeam::markDofs(DofSetArray &dsa)
{
  // dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot | DofSet::Temp | DofSet::Helm);
  // dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
  dsa.mark(nn, 2, activeDofs.list());
}

void
GenRigidMpcBeam::getStiffAndForce(GeomState &gState, CoordSet &cs,
                                       FullSquareMatrix &Ktan, double *f)
{
  Ktan.zero();
  for(int i=0; i<2*numcdofs; ++i) f[i] = 0.0;
}

