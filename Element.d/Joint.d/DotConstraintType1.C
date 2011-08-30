#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Joint.d/exp-map.h>

DotConstraintType1::DotConstraintType1(int* _nn, int _axis1, int _axis2)
 : MpcElement(2, DofSet::XYZrot, _nn)
{
  c0 = 0;
  axis1 = _axis1;
  axis2 = _axis2;
}

DotConstraintType1::~DotConstraintType1()
{
  if(c0) delete [] c0;
}

void
DotConstraintType1::setFrame(EFrame *elemframe) 
{ 
  c0 = new double[3][3];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) 
      c0[i][j] = (*elemframe)[i][j]; 
}

void 
DotConstraintType1::buildFrame(CoordSet& cs)
{
  // build frame if not already defined
  if(!c0) {
    c0 = new double[3][3];

    Node &nd1 = cs.getNode(nn[0]);
    Node &nd2 = cs.getNode(nn[1]);

    double dx = nd2.x - nd1.x;
    double dy = nd2.y - nd1.y;
    double dz = nd2.z - nd1.z;

    double l0 = std::sqrt( dx*dx + dy*dy + dz*dz );

    if(l0 == 0.0) {
      cerr << " *** ERROR: division by zero in DotConstraintType1::buildFrame between nodes " << nn[0]+1 << " and " << nn[1]+1 << endl;
      exit(-1);
    }

    c0[0][0] = dx/l0;
    c0[0][1] = dy/l0;
    c0[0][2] = dz/l0;

    double N1 = std::sqrt( c0[0][0]*c0[0][0] + c0[0][1]*c0[0][1] );
    double N2 = std::sqrt( c0[0][0]*c0[0][0] + c0[0][2]*c0[0][2] );

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

  // fill in the coefficients and rhs (same as updateLMPC but with initial configuration ie zero rotation at each node)
  double r[3] = { 0.0, 0.0, 0.0 };
  double dRdvi[3][3], d1[3], d2[3];
  for(int i = 0; i < 3; ++i) {
    // partial derivatives of rotation matrices wrt ith rotation parameters
    Partial_R_Partial_EM3(r, i, dRdvi);
    // partial derivatives of constraint function wrt ith rotation parameters
    mat_mult_vec(dRdvi, c0[axis1], d1);
    mat_mult_vec(dRdvi, c0[axis2], d2);
    terms[0+i].coef.r_value = d1[0]*c0[axis2][0] + d1[1]*c0[axis2][1] + d1[2]*c0[axis2][2];
    terms[3+i].coef.r_value = c0[axis1][0]*d2[0] + c0[axis1][1]*d2[1] + c0[axis1][2]*d2[2];
  }

  // -ve value of constraint function
  rhs.r_value = -(c0[axis1][0]*c0[axis2][0] + c0[axis1][1]*c0[axis2][1] + c0[axis1][2]*c0[axis2][2]);
}

void 
DotConstraintType1::update(GeomState& gState, CoordSet& cs, double)
{
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  // rotated cframes
  double c1[3][3], c2[3][3];
  mat_mult_mat(c0, ns1.R, c1, 2);
  mat_mult_mat(c0, ns2.R, c2, 2);

  // rotation parameters (thetax, thetay, thetaz)
  double r1[3], r2[3];
  mat_to_vec(ns1.R, r1);
  mat_to_vec(ns2.R, r2);

  double dRdvi1[3][3], dRdvi2[3][3], d1[3], d2[3];
  for(int i=0; i<3; ++i) {
    // partial derivatives of rotation matrices wrt ith rotation parameters
    Partial_R_Partial_EM3(r1, i, dRdvi1);
    Partial_R_Partial_EM3(r2, i, dRdvi2);

    // partial derivatives of constraint functions wrt ith rotation parameters
    mat_mult_vec(dRdvi1, c0[axis1], d1);
    mat_mult_vec(dRdvi2, c0[axis2], d2);

    terms[0+i].coef.r_value = d1[0]*c2[axis2][0] + d1[1]*c2[axis2][1] + d1[2]*c2[axis2][2];
    terms[3+i].coef.r_value = c1[axis1][0]*d2[0] + c1[axis1][1]*d2[1] + c1[axis1][2]*d2[2];
  }

  // -ve value of constraint function
  rhs.r_value = -(c1[axis1][0]*c2[axis2][0] + c1[axis1][1]*c2[axis2][1] + c1[axis1][2]*c2[axis2][2]);
}

void
DotConstraintType1::getHessian(GeomState& gState, CoordSet& cs, FullSquareMatrix& H)
{
  H.zero();

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  // rotated cframes
  double c1[3][3], c2[3][3];
  mat_mult_mat(c0, ns1.R, c1, 2);
  mat_mult_mat(c0, ns2.R, c2, 2);

  // rotation parameters (thetax, thetay, thetaz)
  double r1[3], r2[3];
  mat_to_vec(ns1.R, r1);
  mat_to_vec(ns2.R, r2);

  double d2Rdvidvj1[3][3], d2Rdvidvj2[3][3], dRdvi1[3][3], dRdvj2[3][3], d1[3], d2[3];
  for(int i=0; i<3; ++i) {
    // second partial derivatives of rotation matrices wrt rotation parameters
    for(int j=i; j<3; ++j) {
      Second_Partial_R_Partial_EM3(r1, i, j, d2Rdvidvj1);
      Second_Partial_R_Partial_EM3(r2, i, j, d2Rdvidvj2);
      mat_mult_vec(d2Rdvidvj1, c0[axis1], d1);
      mat_mult_vec(d2Rdvidvj2, c0[axis2], d2);

      H[i][j] = H[j][i] = d1[0]*c2[axis2][0] + d1[1]*c2[axis2][1] + d1[2]*c2[axis2][2];
      H[3+i][3+j] = H[3+j][3+i] = c1[axis1][0]*d2[0] + c1[axis1][1]*d2[1] + c1[axis1][2]*d2[2];
    }
  }
  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(r1, i, dRdvi1);
    mat_mult_vec(dRdvi1, c0[axis1], d1);
    for(int j=0; j<3; ++j) {
      Partial_R_Partial_EM3(r2, j, dRdvj2);
      mat_mult_vec(dRdvj2, c0[axis2], d2);

      H[i][3+j] = H[3+j][i] = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2];
    }
  }

}

