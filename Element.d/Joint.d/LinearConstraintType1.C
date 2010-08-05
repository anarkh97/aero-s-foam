#include <Element.d/Joint.d/LinearConstraintType1.h>

LinearConstraintType1::LinearConstraintType1(int* _nn, DofSet *nodalDofs)
 : MpcElement(2, nodalDofs, _nn)
{
}

void 
LinearConstraintType1::buildFrame(CoordSet& cs)
{
  terms[0].coef.r_value = 1.0;
  terms[1].coef.r_value = -1.0;

  double v[2];
  for(int i = 0; i < 2; ++i) {
    switch (terms[i].dofnum) {
      case 0 : v[i] = cs.getNode(nn[i]).x; break;
      case 1 : v[i] = cs.getNode(nn[i]).y; break;
      case 2 : v[i] = cs.getNode(nn[i]).z; break;
    }
  }

  rhs.r_value = v[0] - v[1];
}

int 
LinearConstraintType1::getTopNumber() 
{ 
  return 106; 
}

void 
LinearConstraintType1::update(GeomState& gState, CoordSet& cs, double)
{
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  terms[0].coef.r_value = 1.0;
  terms[1].coef.r_value = -1.0;

  double v[2];
  for(int i = 0; i < 2; ++i) {
    switch (terms[i].dofnum) {
      case 0 : v[i] = gState[nn[i]].x; break;
      case 1 : v[i] = gState[nn[i]].y; break;
      case 2 : v[i] = gState[nn[i]].z; break;
    }
  }

  rhs.r_value = v[0] - v[1];
}

