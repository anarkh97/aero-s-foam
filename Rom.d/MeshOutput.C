#include "MeshOutput.h"

#include "SimpleBuffer.h"

#include <stdexcept>

namespace Rom {

// Atoms

std::ostream &
operator<<(std::ostream &out, const Attrib &source) {
  out << source.nele + 1 << " "
      << source.attr + 1;
  return out;
}

std::ostream &
operator<<(std::ostream &out, const BCond &source) {
  out << source.nnum   + 1 << " "
      << source.dofnum + 1 << " "
      << source.val;
  return out;
}

std::ostream &
operator<<(std::ostream &out, const EFrameData &source) {
  out << source.elnum + 1 << " ";
  
  const EFrame & eframe = source.frame;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      out << eframe[i][j];
      if (i != 3 || j != 3) {
        out << " ";
      }
    }
  }

  return out;
}

// Sections

std::ostream &
operator<<(std::ostream &out, const CoordSet &source) {
  out << "NODES\n";

  const int size = const_cast<CoordSet &>(source).size();
  for (int i = 0; i < size; ++i) {
    const Node &node = const_cast<CoordSet &>(source).getNode(i);
    out << i + 1
        << " " << node.x
        << " " << node.y
        << " " << node.z;
    out << "\n";
  }

  return out;
}

std::ostream &
operator<<(std::ostream &out, const Elemset &source) {
  out << "TOPOLOGY\n";

  const int size = source.last();
  for (int i = 0; i < size; ++i) {
    Element &ele = *source[i];

    out << i + 1 << " " << ele.getElementType();
    
    SimpleBuffer<int> nodes(ele.numNodes());
    ele.nodes(nodes.array());
    for (int j = 0; j < nodes.size(); ++j) {
      out << " " << nodes[j] + 1;
    }
    out << "\n";
  }

  return out;
}

std::ostream &
operator<<(std::ostream &out, const SPropContainer &source) {
  out << "MATERIALS\n";

  const SPropContainer::const_iterator itEnd = source.end();
  for (SPropContainer::const_iterator it = source.begin(); it != itEnd; ++it) {
    out << it->first + 1 << " ";
    const StructProp &sp = it->second;
    out << sp.A    << " "
        << sp.E    << " "
        << sp.nu   << " "
        << sp.rho  << " "
        << sp.c    << " "
        << sp.k    << " "
        << sp.eh   << " "
        << sp.P    << " "
        << sp.Ta   << " "
        << sp.Q    << " "
        << sp.W    << " "
        << sp.Ixx  << " "
        << sp.Iyy  << " "
        << sp.Izz  << " "
        << sp.ymin << " "
        << sp.ymax << " "
        << sp.zmin << " "
        << sp.zmax << "\n";
  }

  return out;
}

// Headers

template <>
const std::string &
InputFileSectionHelper<int, SampleNodeTag>::header(SampleNodeTag) {
  static const std::string result("SAMPLENODES");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<EFrameData, EmptyTag>::header(EmptyTag) {
  static const std::string result("EFRAMES");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<Attrib, EmptyTag>::header(EmptyTag) {
  static const std::string result("ATTRIBUTES");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<BCond, BCond::BCType>::header(BCond::BCType tag) {
  static const std::string result[] = { "FORCES", "DISPLACEMENTS", "IDISPLACEMENTS", "IVELOCITIES" };
  switch (tag) {
    case BCond::Forces:
      return result[0];
    case BCond::Displacements:
      return result[1];
    case BCond::Idisplacements:
      return result[2];
    case BCond::Ivelocities:
      return result[3];
    default:
      throw std::logic_error("Unknown section tag");
  }
}

} /* end namespace Rom */
