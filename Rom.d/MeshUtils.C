#include "MeshUtils.h"

#include "MeshOutput.h"
#include "RenumberingUtils.h"
#include "SimpleBuffer.h"

#include <Driver.d/EFrameData.h>

#include <Driver.d/Domain.h>

MeshDesc::MeshDesc(Domain *domain,
                   const SPropContainer &prop, const AttribContainer &attrib,
                   const SampledMeshRenumbering &ren) :
  properties_(&prop),
  sampleNodeIds_(ren.reducedSampleNodeIds().begin(), ren.reducedSampleNodeIds().end())
{
  // Nodes
  reduce(ren.reducedNodeIds().begin(), ren.reducedNodeIds().end(), domain->getNodes(), nodes_);

  // Elements & element frames
  reduce(ren.reducedElemIds().begin(), ren.reducedElemIds().end(), ren.nodeRenumbering(), *domain->getEset(), elements_);
  copy_eframes(elements_, std::back_inserter(elemFrames_));

  // Attributes
  reduce(ren.elemRenumbering(), make_value_iterator(attrib.begin()), make_value_iterator(attrib.end()), std::back_inserter(attributes_)); 

  // Boundary conditions
  reduce(ren.nodeRenumbering(), domain->getDBC(), domain->getDBC() + domain->nDirichlet(), std::back_inserter(dirichletBConds_));
  reduce(ren.nodeRenumbering(), domain->getNBC(), domain->getNBC() + domain->nNeumann(), std::back_inserter(neumannBConds_));

  // Initial conditions
  reduce(ren.nodeRenumbering(), domain->getInitDisp(), domain->getInitDisp() + domain->numInitDisp(), std::back_inserter(initDisp_));
  reduce(ren.nodeRenumbering(), domain->getInitDisp6(), domain->getInitDisp6() + domain->numInitDisp6(), std::back_inserter(initDisp_));
  reduce(ren.nodeRenumbering(), domain->getInitVelocity(), domain->getInitVelocity() + domain->numInitVelocity(), std::back_inserter(initVel_));
}

std::ostream &
operator<<(std::ostream &out, const MeshDesc &mesh) {
  out << mesh.nodes();
  out << mesh.elements();
  if (!mesh.elemFrames().empty()) {
    out << make_section(mesh.elemFrames().begin(), mesh.elemFrames().end());
  }

  out << make_section(mesh.attributes().begin(), mesh.attributes().end());
  out << mesh.properties();
  
  out << make_section(mesh.dirichletBConds().begin(), mesh.dirichletBConds().end(), BCond::Displacements);
  out << make_section(mesh.neumannBConds().begin(), mesh.neumannBConds().end(), BCond::Forces);
  
  out << make_section(mesh.initDisp().begin(), mesh.initDisp().end(), BCond::Idisplacements);
  out << make_section(mesh.initVel().begin(), mesh.initVel().end(), BCond::Ivelocities);

  out << make_section(mesh.sampleNodeIds().begin(), mesh.sampleNodeIds().end(), SampleNodeTag());

  return out;
}

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
