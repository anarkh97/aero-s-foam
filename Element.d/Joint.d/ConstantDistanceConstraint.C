#include <Element.d/Joint.d/ConstantDistanceConstraint.h>

ConstantDistanceConstraint::ConstantDistanceConstraint(int* _nn)
 : MpcElement(2, DofSet::XYZdisp, _nn)
{
}

void 
ConstantDistanceConstraint::buildFrame(CoordSet& cs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double dx = nd2.x - nd1.x;
  double dy = nd2.y - nd1.y;
  double dz = nd2.z - nd1.z;
  double d0[3] = { dx, dy, dz };
  l0 = sqrt( dx*dx + dy*dy + dz*dz );

  if(l0 == 0.0) {
    cerr << " *** ERROR: division by zero in ConstantDistanceConstraint::buildFrame between nodes " << nn[0]+1 << " and " << nn[1]+1 << endl;
    exit(-1);
  }

  for(int i = 0; i < 3; ++i) {
    terms[0+i].coef.r_value = -d0[i]/l0;
    terms[3+i].coef.r_value = d0[i]/l0;
  }

  rhs.r_value = 0;
}

int 
ConstantDistanceConstraint::getTopNumber() 
{ 
  return 106; 
}

void 
ConstantDistanceConstraint::update(GeomState& gState, CoordSet& cs, double)
{
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;
  double d[3] = { dx, dy, dz };
  double l = sqrt(dx*dx + dy*dy + dz*dz);

  if(l == 0.0) {
    cerr << " *** ERROR: division by zero in ConstantDistanceConstraint::update between nodes " << nn[0]+1 << " and " << nn[1]+1 << endl;
    exit(-1);
  }

  // partial derivatives of constraint functions wrt x, y, z 
  for(int i = 0; i < 3; ++i) {
    terms[0+i].coef.r_value = -d[i]/l;
    terms[3+i].coef.r_value = d[i]/l;
  }

  // values of constraint functions
  rhs.r_value = l - l0;
}

void
ConstantDistanceConstraint::getHessian(GeomState& gState, CoordSet& cs, FullSquareMatrix& H)
{
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;
  double d[3] = { dx, dy, dz };
  double l = sqrt(dx*dx + dy*dy + dz*dz);

  double l2 = l*l;
  double l3 = l*l*l;
  H[0][0] = H[3][3] = (l2-dx*dx)/l3;
  H[1][1] = H[4][4] = (l2-dy*dy)/l3;
  H[2][2] = H[5][5] = (l2-dz*dz)/l3;

  H[0][1] = H[1][0] = H[3][4] = H[4][3] = -dx*dy/l3;
  H[0][2] = H[2][0] = H[3][5] = H[5][3] = -dx*dz/l3;
  H[1][2] = H[2][1] = H[4][5] = H[5][4] = -dy*dz/l3;

  H[0][3] = H[3][0] = -H[0][0];
  H[1][4] = H[4][1] = -H[1][1];
  H[2][5] = H[5][2] = -H[2][2];

  H[0][4] = H[4][0] = H[1][3] = H[3][1] = -H[0][1];
  H[0][5] = H[5][0] = H[2][3] = H[3][2] = -H[0][2];
  H[1][5] = H[5][1] = H[2][4] = H[4][2] = -H[1][2];
}

