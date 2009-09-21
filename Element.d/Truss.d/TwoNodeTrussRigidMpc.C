#include <Element.d/Truss.d/TwoNodeTrussRigidMpc.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

TwoNodeTrussRigidMpc::TwoNodeTrussRigidMpc(int* _nn)
 : RigidMpcElement(2, 3, DofSet::XYZdisp, 1, _nn)
{
  /* do nothing */
}

void 
TwoNodeTrussRigidMpc::computeMPCs(CoordSet &cs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double dx = nd2.x - nd1.x;
  double dy = nd2.y - nd1.y;
  double dz = nd2.z - nd1.z;

  double length = sqrt( dx*dx + dy*dy + dz*dz );

  double c1[3], c2[3], c3[3];

  if(length == 0.0) {
    cerr << " *** ERROR: Rigid truss has zero length, nodes: " << nn[0]+1 << " " << nn[1]+1 << ". Exiting...\n";
    exit(-1);
  }
  else {
    mpcs[0] = new LMPCons(0, 0.0);
 
    double c1[3], c2[3], c3[3];
    c1[0] = dx/length;
    c1[1] = dy/length;
    c1[2] = dz/length;

    // translation in x, y, z (1 constraint equation)
    for(int i = 0; i < 3; ++i) {
      mpcs[0]->addterm(new LMPCTerm(nn[0], i, c1[i]));
      mpcs[0]->addterm(new LMPCTerm(nn[1], i, -c1[i]));
    }
  }
}

int
TwoNodeTrussRigidMpc::getTopNumber()
{
  return 101;
}

bool
TwoNodeTrussRigidMpc::isSafe()
{
  return false;
}

