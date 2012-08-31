#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/FlexibleTwoNodeTruss.h>
#include <Element.d/Joint.d/LinearTranslationalSpring.h>

FlexibleTwoNodeTruss::FlexibleTwoNodeTruss(int* _nn)
  : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 1;
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };

  subElems[0] = new LinearTranslationalSpring(indices);
}

void
FlexibleTwoNodeTruss::buildFrame(CoordSet& cs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double x[2], y[2], z[2];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

  double dx = x[1] - x[0];
  double dy = y[1] - y[0];
  double dz = z[1] - z[0];

  l0 = sqrt(dx*dx + dy*dy + dz*dz);

  SuperElement::buildFrame(cs);
}

void
FlexibleTwoNodeTruss::setProp(StructProp *p, bool myProp)
{
  SuperElement::setProp(p, myProp);

  subElems[0]->getProperty()->penalty = p->E*p->A/l0;
}

FullSquareMatrix
FlexibleTwoNodeTruss::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
        FullSquareMatrix elementMassMatrix(6, mel);
        elementMassMatrix.zero();

        const double mass = getMass(cs);

        if (cmflg) {
                // Consistent mass matrix
                const double outDiagMass = mass / 6.0;
                const double diagMass = 2.0 * outDiagMass;

                for (int i = 0; i < 6; ++i) {
                  elementMassMatrix[i][i] = diagMass;
                }
                for (int i = 0; i < 3; ++i) {
                  const int j = i + 3;
                  elementMassMatrix[i][j] = outDiagMass;
                  elementMassMatrix[j][i] = outDiagMass;
                }
        } else {
                // Lumped mass matrix
                const double massPerNode = 0.5 * mass;
                for (int i = 0; i < 6; ++i) {
                  elementMassMatrix[i][i] = massPerNode;
                }
        }

        return elementMassMatrix;
}

double
FlexibleTwoNodeTruss::getMass(CoordSet& cs)
{
        Node &nd1 = cs.getNode( nn[0] );
        Node &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );
        double   mass = (length*prop->A*prop->rho);

        return mass;
}
#endif
