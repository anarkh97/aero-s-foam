#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Joint.d/exp-map.h>

NonlinearTorsionalSpring::NonlinearTorsionalSpring(int* _nn, int _axis1, int _axis2)
 : MpcElement(2, DofSet::XYZrot, _nn)
{
  c0 = 0;
  m_axis1 = _axis1;
  m_axis2 = _axis2;
  covariant_derivatives = true;
}

NonlinearTorsionalSpring::~NonlinearTorsionalSpring()
{
  if(c0) delete [] c0;
}

void
NonlinearTorsionalSpring::setProp(StructProp *p, bool _myProp)
{
  prop = (_myProp) ? p : new StructProp(*p); 
  myProp = true;
  if((m_axis1 == 2 && m_axis2 == 1) || (m_axis1 == 1 && m_axis2 == 2)) prop->penalty = prop->kx;
  else if((m_axis1 == 2 && m_axis2 == 0) || (m_axis1 == 0 && m_axis2 == 2)) prop->penalty = prop->ky;
  else prop->penalty = prop->kz;
  prop->lagrangeMult = false;
}

void
NonlinearTorsionalSpring::setFrame(EFrame *elemframe) 
{ 
  if(!c0) c0 = new double[3][3];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) 
      c0[i][j] = (*elemframe)[i][j]; 
}

void 
NonlinearTorsionalSpring::buildFrame(CoordSet& cs)
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
      cerr << " *** ERROR: division by zero in NonlinearTorsionalSpring::buildFrame between nodes " << nn[0]+1 << " and " << nn[1]+1 << endl;
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

  // fill in the coefficients and rhs (same as ::update but with initial configuration ie zero rotation at each node)
  double r[3] = { 0.0, 0.0, 0.0 };
  double dRdvi[3][3], d1[3], d2[3];

  int axis1 = m_axis1, axis2 = m_axis2;
  double cth = (c0[axis1][0]*c0[axis2][0] + c0[axis1][1]*c0[axis2][1] + c0[axis1][2]*c0[axis2][2]);
  double dacosdcth = -1/std::sqrt(1-cth*cth); // derivative of arccos(cth) w.r.t. cth

  for(int i = 0; i < 3; ++i) {
    // partial derivatives of rotation matrices wrt ith rotation parameters
    Partial_R_Partial_EM3(r, i, dRdvi);
    // partial derivatives of constraint function wrt ith rotation parameters
    mat_mult_vec(dRdvi, c0[axis1], d1);
    mat_mult_vec(dRdvi, c0[axis2], d2);
    terms[0+i].coef.r_value = dacosdcth*(d1[0]*c0[axis2][0] + d1[1]*c0[axis2][1] + d1[2]*c0[axis2][2]);
    terms[3+i].coef.r_value = dacosdcth*(c0[axis1][0]*d2[0] + c0[axis1][1]*d2[1] + c0[axis1][2]*d2[2]);
  }

  // -ve value of constraint function
  rhs.r_value = M_PI/2-std::acos(cth);
}

void 
NonlinearTorsionalSpring::update(GeomState& gState, CoordSet& cs, double t)
{
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  // internal states
  updateStates((GeomState *) NULL, gState, cs);

  // rotated cframes
  double c1[3][3], c2[3][3];
  mat_mult_mat(c0, ns1.R, c1, 2);
  mat_mult_mat(c0, ns2.R, c2, 2);

  int axis1 = (quadrant == 0 || quadrant == 2) ? m_axis1 : m_axis2, axis2 = m_axis2;
  double cth = (c1[axis1][0]*c2[axis2][0] + c1[axis1][1]*c2[axis2][1] + c1[axis1][2]*c2[axis2][2]);

  // reparameterization
  if(cth > 0.707106781 || cth < -0.707106781) {
    if(axis1 != axis2) {
      //std::cerr << "reparameterizing #1, cth = " << cth << ", offset = " << offset << std::endl;
      if(quadrant == 0 && cth > 0)      { quadrant = 1; offset2 += M_PI/2; }
      else if(quadrant == 0 && cth < 0) { quadrant = 3; offset2 -= M_PI/2; }
      else if(quadrant == 2 && cth < 0) { quadrant = 3; offset2 += M_PI/2; }
      else if(quadrant == 2 && cth > 0) { quadrant = 1; offset2 -= M_PI/2; }
      else std::cerr << "whoops-a-daisy\n";

      axis1 = m_axis2;

      // NEW: correct for either +ve or -ve force, and arbitrary change of direction
      if(cth > 0) offset = M_PI/2 - offset2;
      else        offset = M_PI/2 + offset2;
    }
    else {
      //std::cerr << "reparameterizing #2, cth = " << cth << ", offset = " << offset << std::endl;
      if(quadrant == 1 && cth < 0)      { quadrant = 2; offset2 += M_PI/2; }
      else if(quadrant == 1 && cth > 0) { quadrant = 0; offset2 -= M_PI/2; }
      else if(quadrant == 3 && cth > 0) { quadrant = 0; offset2 += M_PI/2; }
      else if(quadrant == 3 && cth < 0) { quadrant = 2; offset2 -= M_PI/2; }
      else std::cerr << "whoops-a-daisy\n";

      axis1 = m_axis1;

      // NEW: correct for either +ve or -ve force, and arbitrary change of direction
      if(cth > 0) offset = M_PI/2 + offset2;
      else        offset = M_PI/2 - offset2;
    }
    cth = (c1[axis1][0]*c2[axis2][0] + c1[axis1][1]*c2[axis2][1] + c1[axis1][2]*c2[axis2][2]);
  }

  
  double dacosdcth = -1/std::sqrt(1-cth*cth); // derivative of arccos(cth) w.r.t. cth

  if(covariant_derivatives) {

    // instantaneous rotation parameters (thetax, thetay, thetaz)
    double r[3] = { 0, 0, 0 };

    double dRdvi[3][3], d1[3], d2[3];
    for(int i=0; i<3; ++i) {
      // partial derivatives of instantaneous rotation matrices wrt ith instantaneous rotation parameters
      Partial_R_Partial_EM3(r, i, dRdvi);

      // partial derivatives of constraint functions wrt ith instantaneous rotation parameters
      mat_mult_vec(dRdvi, c1[axis1], d1);
      mat_mult_vec(dRdvi, c2[axis2], d2);

      terms[0+i].coef.r_value = dacosdcth*(d1[0]*c2[axis2][0] + d1[1]*c2[axis2][1] + d1[2]*c2[axis2][2]);
      terms[3+i].coef.r_value = dacosdcth*(c1[axis1][0]*d2[0] + c1[axis1][1]*d2[1] + c1[axis1][2]*d2[2]);
    }
  }
  else {

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

      terms[0+i].coef.r_value = dacosdcth*(d1[0]*c2[axis2][0] + d1[1]*c2[axis2][1] + d1[2]*c2[axis2][2]);
      terms[3+i].coef.r_value = dacosdcth*(c1[axis1][0]*d2[0] + c1[axis1][1]*d2[1] + c1[axis1][2]*d2[2]);
    }
  }

  // -ve value of constraint function
  rhs.r_value = offset-std::acos(cth);
}

void
NonlinearTorsionalSpring::getHessian(GeomState& gState, CoordSet& cs, FullSquareMatrix& H)
{
  H.zero();

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  // rotated cframes
  double c1[3][3], c2[3][3];
  mat_mult_mat(c0, ns1.R, c1, 2);
  mat_mult_mat(c0, ns2.R, c2, 2);

  int axis1 = m_axis1, axis2 = (quadrant == 0 || quadrant == 2) ? m_axis2 : m_axis1;
  double cth = (c1[axis1][0]*c2[axis2][0] + c1[axis1][1]*c2[axis2][1] + c1[axis1][2]*c2[axis2][2]);
  double dacosdcth = -1/std::sqrt(1-cth*cth);
  double d2acosdcth2 = -cth/std::pow(1-cth*cth,1.5);

  if(covariant_derivatives) {

    // instantaneous rotation parameters
    double r[3] = { 0, 0, 0 };

    double d2Rdvidvj[3][3], dRdvi[3][3], dRdvj[3][3], d1[3], d2[3];
    for(int i=0; i<3; ++i) {
      // second partial derivatives of rotation matrices wrt instantaneous rotation parameters
      for(int j=i; j<3; ++j) {
        Second_Partial_R_Partial_EM3(r, i, j, d2Rdvidvj);

        mat_mult_vec(d2Rdvidvj, c1[axis1], d1);
        mat_mult_vec(d2Rdvidvj, c2[axis2], d2);

        H[i][j] = H[j][i] = dacosdcth*(d1[0]*c2[axis2][0] + d1[1]*c2[axis2][1] + d1[2]*c2[axis2][2]) 
                            + d2acosdcth2*(terms[i].coef.r_value/dacosdcth)*(terms[j].coef.r_value/dacosdcth);
        H[3+i][3+j] = H[3+j][3+i] = dacosdcth*(c1[axis1][0]*d2[0] + c1[axis1][1]*d2[1] + c1[axis1][2]*d2[2]) 
                                    + d2acosdcth2*(terms[3+i].coef.r_value/dacosdcth)*(terms[3+j].coef.r_value/dacosdcth);
      }
    }
    for(int i=0; i<3; ++i) {
      Partial_R_Partial_EM3(r, i, dRdvi);
      mat_mult_vec(dRdvi, c1[axis1], d1);
      for(int j=0; j<3; ++j) {
        Partial_R_Partial_EM3(r, j, dRdvj);
        mat_mult_vec(dRdvj, c2[axis2], d2);
        H[i][3+j] = H[3+j][i] = dacosdcth*(d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]) 
                                + d2acosdcth2*(terms[i].coef.r_value/dacosdcth)*(terms[3+j].coef.r_value/dacosdcth);
      }
    }
  }
  else {

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

        H[i][j] = H[j][i] = dacosdcth*(d1[0]*c2[axis2][0] + d1[1]*c2[axis2][1] + d1[2]*c2[axis2][2])
                            + d2acosdcth2*(terms[i].coef.r_value/dacosdcth)*(terms[j].coef.r_value/dacosdcth);
        H[3+i][3+j] = H[3+j][3+i] = dacosdcth*(c1[axis1][0]*d2[0] + c1[axis1][1]*d2[1] + c1[axis1][2]*d2[2])
                                    + d2acosdcth2*(terms[3+i].coef.r_value/dacosdcth)*(terms[3+j].coef.r_value/dacosdcth);
      }
    }
    for(int i=0; i<3; ++i) {
      Partial_R_Partial_EM3(r1, i, dRdvi1);
      mat_mult_vec(dRdvi1, c0[axis1], d1);
      for(int j=0; j<3; ++j) {
        Partial_R_Partial_EM3(r2, j, dRdvj2);
        mat_mult_vec(dRdvj2, c0[axis2], d2);

        H[i][3+j] = H[3+j][i] = dacosdcth*(d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2])
                                + d2acosdcth2*(terms[i].coef.r_value/dacosdcth)*(terms[3+j].coef.r_value/dacosdcth);
      }
    }
  }

}

int
NonlinearTorsionalSpring::numStates()
{
  return 3;
}

void
NonlinearTorsionalSpring::initStates(double *statenp)
{
  statenp[0] = M_PI/2; // offset
  statenp[1] = 0.0;    // offset2
  statenp[2] = 0;      // quadrant
}

void 
NonlinearTorsionalSpring::updateStates(GeomState *, GeomState &gState, CoordSet &)
{
  // TODO: consider if it is better to update the state from the reference state (i.e. the last converged solution)
  //       rather than the current newton iteration, as we do for plasticity
  double *statenp = gState.getElemState(getGlNum()) + stateOffset;
  offset = statenp[0];
  offset2 = statenp[1];
  quadrant = statenp[2];

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  // rotated cframes
  double c1[3][3], c2[3][3];
  mat_mult_mat(c0, ns1.R, c1, 2);
  mat_mult_mat(c0, ns2.R, c2, 2);

  int axis1 = (quadrant == 0 || quadrant == 2) ? m_axis1 : m_axis2, axis2 = m_axis2;
  double cth = (c1[axis1][0]*c2[axis2][0] + c1[axis1][1]*c2[axis2][1] + c1[axis1][2]*c2[axis2][2]);

  // reparameterization
  if(cth > 0.707106781 || cth < -0.707106781) {
    if(axis1 != axis2) {
      //std::cerr << "reparameterizing #1, cth = " << cth << ", offset = " << offset << std::endl;
      if(quadrant == 0 && cth > 0) { 
        //std::cerr << "going from quadrant 0 to quadrant 1\n";
        quadrant = 1; offset2 += M_PI/2; 
      }
      else if(quadrant == 0 && cth < 0) { 
        //std::cerr << "going from quadrant 0 to quadrant 3\n";
        quadrant = 3; offset2 -= M_PI/2;
      }
      else if(quadrant == 2 && cth < 0) {
        //std::cerr << "going from quadrant 2 to quadrant 3\n";
        quadrant = 3; offset2 += M_PI/2;
      }
      else if(quadrant == 2 && cth > 0) {
        //std::cerr << "going from quadrant 2 to quadrant 1\n";
        quadrant = 1; offset2 -= M_PI/2;
      }
      else std::cerr << "whoops-a-daisy\n";

      if(cth > 0) offset = M_PI/2 - offset2;
      else        offset = M_PI/2 + offset2;
    }
    else {
      //std::cerr << "reparameterizing #2, cth = " << cth << ", offset = " << offset << std::endl;
      if(quadrant == 1 && cth < 0) {
        //std::cerr << "going from quadrant 1 to quadrant 2\n";
        quadrant = 2; offset2 += M_PI/2;
      }
      else if(quadrant == 1 && cth > 0) {
        //std::cerr << "going from quadrant 1 to quadrant 0\n";
        quadrant = 0; offset2 -= M_PI/2;
      }
      else if(quadrant == 3 && cth > 0) {
        //std::cerr << "going from quadrant 3 to quadrant 0\n";
        quadrant = 0; offset2 += M_PI/2;
      }
      else if(quadrant == 3 && cth < 0) {
        //std::cerr << "going from quadrant 3 to quadrant 2\n";
        quadrant = 2; offset2 -= M_PI/2;
      }
      else std::cerr << "whoops-a-daisy\n";

      if(cth > 0) offset = M_PI/2 + offset2;
      else        offset = M_PI/2 - offset2;
    }
  }

  statenp[0] = offset;
  statenp[1] = offset2;
  statenp[2] = quadrant;
}
