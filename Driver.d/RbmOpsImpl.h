// ****************************************************************************************************
#include "SubDomain.h"
#include <Math.d/SparseSet.h>

template<class Scalar>
void
GenSubDomain<Scalar>::makeZstarAndR(double *centroid)
{
  rigidBodyModesG = std::make_unique<Rbm>(dsa, c_dsa, nodes, sinfo.tolsvd,
                            centroid+3*group, cornerNodes, numCRN, numCRNdof, cornerDofs,
                            numMPC_primal, this->mpc_primal);
}

