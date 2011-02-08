#include "MeshUtils.h"

#include "SimpleBuffer.h"

#include <Driver.d/EFrameData.h>
#include <Driver.d/StructProp.h>

const Elemset &
renumber_nodes(const std::map<int, int> &nodeIndices, Elemset &target) {
  const int indexEnd = nodeIndices.size() > 0 ? nodeIndices.rbegin()->first + 1 : 0;
  SimpleBuffer<int> buffer(indexEnd);
 
  for (std::map<int, int>::const_iterator it = nodeIndices.begin(); it != nodeIndices.end(); ++it) {
    buffer[it->first] = it->second;
  }

  const int iElemEnd = target.last();
  for (int iElem = 0; iElem < iElemEnd; ++iElem) {
    if (Element * e = target[iElem]) {
      e->renum(buffer.array());
    }
  }

  return target;
}

template <>
int BCond::* const
MeshSectionTraits<BCond>::renumbered_slot = &BCond::nnum;

template <>
int Attrib::* const
MeshSectionTraits<Attrib>::renumbered_slot = &Attrib::nele;

template <>
int EFrameData::* const
MeshSectionTraits<EFrameData>::renumbered_slot = &EFrameData::elnum;
