#include "VecNodeDof6Conversion.h"

#include <Utils.d/dofset.h>

VecNodeDof6Conversion::VecNodeDof6Conversion(const DofSetArray &dsa) :
  nodeCount_(const_cast<DofSetArray &>(dsa).numNodes()),
  vectorSize_(const_cast<DofSetArray &>(dsa).size()),
  dofLocation_(new int[nodeCount()][6])
{
  DofSetArray &dsa_fix = const_cast<DofSetArray &>(dsa);
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    dofLocation_[iNode][0] = dsa_fix.locate(iNode, DofSet::Xdisp);
    dofLocation_[iNode][1] = dsa_fix.locate(iNode, DofSet::Ydisp);
    dofLocation_[iNode][2] = dsa_fix.locate(iNode, DofSet::Zdisp);
    dofLocation_[iNode][3] = dsa_fix.locate(iNode, DofSet::Xrot);
    dofLocation_[iNode][4] = dsa_fix.locate(iNode, DofSet::Yrot);
    dofLocation_[iNode][5] = dsa_fix.locate(iNode, DofSet::Zrot);
  }
}

VecNodeDof6Conversion::~VecNodeDof6Conversion() {
  delete[] dofLocation_;
}
