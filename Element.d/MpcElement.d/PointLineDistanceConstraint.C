#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/PointLineDistanceConstraint.h>

PointLineDistanceConstraint::PointLineDistanceConstraint(int* _nn)
 : ConstraintFunctionElement<PointLineDistanceConstraintFunction>(1, DofSet::XYZdisp, _nn, 0) 
{
}

void
PointLineDistanceConstraint::setFrame(EFrame *elemframe)
{
  for(int i = 0; i < 3; ++i) {
    x1[i] = (*elemframe)[0][i];
    x2[i] = (*elemframe)[1][i];
  }
}

void
PointLineDistanceConstraint::getConstants(CoordSet& cs, Eigen::Array<double,14,1>& sconst, Eigen::Array<int,1,1>& iconst)
{
  // note: StructProps::relop = -1 --> f(x) <= 0
  //                          = 0  --> f(x) = 0
  //                          = +1 --> f(x) >= 0 in this case ConstraintFunction::operator() returns -f(x)
  //                                             and we enforce -f(x) <= 0
  this->type = (prop) ? int(prop->relop != 0) : 0;

  sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z, x1[0], x1[1], x1[2], x2[0], x2[1], x2[2],
            (prop) ? prop->amplitude : 0, (prop) ? prop->omega : 0, (prop) ? prop->phase : 0,
            (prop) ? prop->B : 1, (prop) ? prop->C : 0;
  iconst << ((prop) ? int(prop->relop == 1) : 0);
}

/*
void 
PointLineDistanceConstraint::buildFrame(CoordSet& cs)
{
  using std::sqrt;
  Node &nd = cs.getNode(nn[0]);
  double x0[3] = { nd.x, nd.y, nd.z };

  double a[3] = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };
  double b[3] = { x1[0]-x0[0], x1[1]-x0[1], x1[2]-x0[2] };
  double c[3] = { a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0] }; // a x b (cross product)
  double numerator = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  double denominator = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  d0 = numerator / denominator;

  // d = sqrt( (a[1]*b[2] - a[2]*b[1])^2 + (a[2]*b[0] - a[0]*b[2])^2 + (a[0]*b[1] - a[1]*b[0])^2 ) / denominator
  //   = sqrt( (a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))^2 +
  //           (a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))^2 +
  //           (a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0]))^2 ) / denominator

  terms[0].coef.r_value = (a[1]*c[2]-a[2]*c[1])/(denominator*numerator);
  terms[1].coef.r_value = (a[2]*c[0]-a[0]*c[2])/(denominator*numerator);
  terms[2].coef.r_value = (a[0]*c[1]-a[1]*c[0])/(denominator*numerator);

  rhs.r_value = 0;
}

void 
PointLineDistanceConstraint::update(GeomState& gState, CoordSet& cs, double)
{
  using std::sqrt;
  // node's current coordinates
  NodeState &nd = gState[nn[0]];
  double x0[3] = { nd.x, nd.y, nd.z };

  double a[3] = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };
  double b[3] = { x1[0]-x0[0], x1[1]-x0[1], x1[2]-x0[2] };
  double c[3] = { a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0] }; // a x b (cross product)
  double numerator = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  double denominator = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  double d = numerator / denominator;

  // d = sqrt( (a[1]*b[2] - a[2]*b[1])^2 + (a[2]*b[0] - a[0]*b[2])^2 + (a[0]*b[1] - a[1]*b[0])^2 ) / denominator
  //   = sqrt( (a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))^2 +
  //           (a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))^2 +
  //           (a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0]))^2 ) / denominator

  // partial derivatives of constraint functions wrt x, y, z 
  terms[0].coef.r_value = (a[1]*c[2]-a[2]*c[1])/(denominator*numerator);
  terms[1].coef.r_value = (a[2]*c[0]-a[0]*c[2])/(denominator*numerator);
  terms[2].coef.r_value = (a[0]*c[1]-a[1]*c[0])/(denominator*numerator);
  
  // -ve value of constraint function
  rhs.r_value = -(d - d0);
}

void
PointLineDistanceConstraint::getHessian(GeomState& gState, CoordSet& cs, FullSquareMatrix& H, double)
{
  using std::sqrt;
  // node's current coordinates
  NodeState &nd = gState[nn[0]];
  double x0[3] = { nd.x, nd.y, nd.z };

  double a[3] = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };
  double b[3] = { x1[0]-x0[0], x1[1]-x0[1], x1[2]-x0[2] };
  double c[3] = { a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0] }; // a x b (cross product)
  double numerator = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  double denominator = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

  // d = sqrt( (a[1]*b[2] - a[2]*b[1])^2 + (a[2]*b[0] - a[0]*b[2])^2 + (a[0]*b[1] - a[1]*b[0])^2 ) / denominator
  //   = sqrt( (a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))^2 +
  //           (a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))^2 +
  //           (a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0]))^2 ) / denominator

  // 1/2 of the first derivatives of numerator^2:
  // a[1]*c[2]-a[2]*c[1] = a[1]*(a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0])) 
  //                      -a[2]*(a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))    
  // a[2]*c[0]-a[0]*c[2] = a[2]*(a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))
  //                      -a[0]*(a[0]*(x1[1]-x0[1]) - a[1]*(x1[0]-x0[0]))
  // a[0]*c[1]-a[1]*c[0] = a[0]*(a[2]*(x1[0]-x0[0]) - a[0]*(x1[2]-x0[2]))
  //                      -a[1]*(a[1]*(x1[2]-x0[2]) - a[2]*(x1[1]-x0[1]))
  double dx = a[1]*c[2]-a[2]*c[1];
  double dy = a[2]*c[0]-a[0]*c[2];
  double dz = a[0]*c[1]-a[1]*c[0];

  // 1/2 of the second derivatives of numerator^2
  double dxx = a[1]*a[1]+a[2]*a[2];
  double dyy = a[0]*a[0]+a[2]*a[2];
  double dzz = a[0]*a[0]+a[1]*a[1];
  double dxy = -a[1]*a[0];
  double dxz = -a[2]*a[0];
  double dyz = -a[2]*a[1];

  // second derivatives
  double n2 = numerator*numerator;
  double n3 = numerator*numerator*numerator;
  H[0][0] = (n2*dxx - dx*dx)/(denominator*n3);
  H[1][1] = (n2*dyy - dy*dy)/(denominator*n3);
  H[2][2] = (n2*dzz - dz*dz)/(denominator*n3);

  H[0][1] = H[1][0] = (n2*dxy - dx*dy)/(denominator*n3);
  H[0][2] = H[2][0] = (n2*dxz - dx*dz)/(denominator*n3);
  H[1][2] = H[2][1] = (n2*dyz - dy*dz)/(denominator*n3);
}
*/
#endif
