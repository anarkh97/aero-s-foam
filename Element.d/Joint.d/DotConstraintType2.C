#include <Element.d/Joint.d/DotConstraintType2.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Joint.d/exp-map.h>

DotConstraintType2::DotConstraintType2(int* _nn, int _axis)
 : MpcElement(2, DofSet::XYZdisp | DofSet::XYZrot, _nn)
{
  axis = _axis;
  elemframe = 0;
}

void
DotConstraintType2::setFrame(EFrame *_elemframe)
{
  elemframe = _elemframe;
}

void 
DotConstraintType2::buildFrame(CoordSet& cs)
{
  // build frame if not already defined
  if(false /*elemframe*/) {
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        c0[i][j] = (*elemframe)[i][j];
  }
  else {
    Node &nd1 = cs.getNode(nn[0]);
    Node &nd2 = cs.getNode(nn[1]);

    double dx = nd2.x - nd1.x;
    double dy = nd2.y - nd1.y;
    double dz = nd2.z - nd1.z;

    double l0 = sqrt( dx*dx + dy*dy + dz*dz );

    if(l0 == 0.0) {
      cerr << " *** ERROR: division by zero in DotConstraintType2::buildFrame between nodes " << nn[0]+1 << " and " << nn[1]+1 << endl;
      exit(-1);
    }

    c0[0][0] = dx/l0;
    c0[0][1] = dy/l0;
    c0[0][2] = dz/l0;

    double N1 = sqrt( c0[0][0]*c0[0][0] + c0[0][1]*c0[0][1] );
    double N2 = sqrt( c0[0][0]*c0[0][0] + c0[0][2]*c0[0][2] );

    if (N1 > N2) {
      c0[1][0] = -c0[0][1]/N1;
      c0[1][1] = c0[0][0]/N1;
      c0[1][2] = 0.0;
    }
    else {
      c0[1][0] = c0[0][2]/N2;
      c0[1][1] = 0.0;
      c0[1][2] = -c0[0][0]/N2;
    }

    c0[2][0] = c0[0][1] * c0[1][2] - c0[0][2] * c0[1][1];
    c0[2][1] = c0[0][2] * c0[1][0] - c0[0][0] * c0[1][2];
    c0[2][2] = c0[0][0] * c0[1][1] - c0[0][1] * c0[1][0];
  }

  // fill in the coefficients and rhs (same as update but with initial configuration ie zero rotation at each node)
  double d[3] = { cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z };
  double r[3] = { 0.0, 0.0, 0.0 };

  // partial derivatives of constraint functions wrt x, y, z 
  for(int i = 0; i < 3; ++i) {
    terms[0+i].coef.r_value = -c0[axis][i];
    terms[6+i].coef.r_value = c0[axis][i];
  }

  double dRdvi1[3][3], d1[3];
  for(int i=0; i<3; ++i) {
    // partial derivatives of rotation matrices wrt ith rotation parameters
    Partial_R_Partial_EM3(r, i, dRdvi1);
    // partial derivatives of constraint functions wrt ith rotation parameters
    mat_mult_vec(dRdvi1, c0[axis], d1);
    terms[3+i].coef.r_value = d1[0]*d[0] + d1[1]*d[1] + d1[2]*d[2];
  }

  rhs.r_value = c0[axis][0]*d[0] + c0[axis][1]*d[1] + c0[axis][2]*d[2];
}

int 
DotConstraintType2::getTopNumber() 
{ 
  return 106; 
}

void 
DotConstraintType2::update(GeomState& gState, CoordSet& cs)
{
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;
  double d[3] = { dx, dy, dz };

  // rotated cframe axis
  double c1[3];
  mat_mult_vec(ns1.R, c0[axis], c1);

  // partial derivatives of constraint functions wrt x, y, z 
  for(int i = 0; i < 3; ++i) {
    terms[0+i].coef.r_value = -c1[i];
    terms[6+i].coef.r_value = c1[i];
  }

  // rotation parameters (thetax, thetay, thetaz)
  double r1[3];
  mat_to_vec(ns1.R, r1);

  double dRdvi1[3][3], d1[3];
  for(int i=0; i<3; ++i) {
    // partial derivatives of rotation matrices wrt ith rotation parameters
    Partial_R_Partial_EM3(r1, i, dRdvi1);
    // partial derivatives of constraint functions wrt ith rotation parameters
    mat_mult_vec(dRdvi1, c0[axis], d1);
    terms[3+i].coef.r_value = d1[0]*d[0] + d1[1]*d[1] + d1[2]*d[2];
  }

  rhs.r_value = c1[0]*d[0] + c1[1]*d[1] + c1[2]*d[2];
}

void
DotConstraintType2::getHessian(GeomState& gState, CoordSet& cs, FullSquareMatrix& H)
{
  H.zero();

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;
  double d[3] = { dx, dy, dz };

  // rotated cframes
  double c1[3][3];
  mat_mult_mat(c0, ns1.R, c1, 2);

  // rotation parameters (thetax, thetay, thetaz)
  double r1[3];
  mat_to_vec(ns1.R, r1);

  double d2Rdvidvj1[3][3], dRdvi1[3][3], d1[3][3];
  for(int i=0; i<3; ++i) {
    // second partial derivatives of rotation matrices wrt rotation parameters
    for(int j=i; j<3; ++j) {
      Second_Partial_R_Partial_EM3(r1, i, j, d2Rdvidvj1);
      mat_mult_vec(d2Rdvidvj1, c0[axis], d1[axis]);
      H[3+i][3+j] = H[3+j][3+i] = d1[axis][0]*d[0] + d1[axis][1]*d[1] + d1[axis][2]*d[2];
    }
  }
  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(r1, i, dRdvi1);
    mat_mult_vec(dRdvi1, c0[axis], d1[axis]);
    for(int j=0; j<3; ++j) {
      H[3+i][j] = H[j][3+i] = -d1[axis][j];
      H[3+i][6+j] = H[6+j][3+i] = d1[axis][j];
    }
  }
}

