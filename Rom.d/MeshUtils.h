#ifndef ROM_MESHUTILS_H
#define ROM_MESHUTILS_H

#include <Element.d/Element.h>
#include <Driver.d/StructProp.h>

#include <map>
#include <iterator>
#include <cstring>

class SampledMeshRenumbering;

class Domain;

class EFrameData;
class Attrib;
class BCond;

class MeshDesc {
public:
  const CoordSet &nodes() const { return nodes_; }
  const Elemset &elements() const { return elements_; }

  const std::vector<EFrameData> &elemFrames() const { return elemFrames_; }
  const SPropContainer &properties() const { return *properties_; } 
  const std::vector<Attrib> &attributes() const { return attributes_; }

  const std::vector<BCond> &dirichletBConds() const { return dirichletBConds_; }
  const std::vector<BCond> &neumannBConds() const { return neumannBConds_; }

  const std::vector<BCond> &initDisp() const { return initDisp_; }
  const std::vector<BCond> &initVel() const { return initVel_; }

  // Only for a reduced mesh
  const std::vector<int> &sampleNodeIds() const { return sampleNodeIds_; }

  typedef std::map<int, Attrib> AttribContainer;

  // Creates a reduced mesh from the default one
  MeshDesc(Domain *, const SPropContainer &, const AttribContainer &, const SampledMeshRenumbering &);

private:
  CoordSet nodes_; 
  Elemset elements_;

  std::vector<EFrameData> elemFrames_;
  std::vector<Attrib> attributes_;
  const SPropContainer *properties_;

  std::vector<BCond> dirichletBConds_;
  std::vector<BCond> neumannBConds_;
  std::vector<BCond> initDisp_; 
  std::vector<BCond> initVel_;

  const std::vector<int> sampleNodeIds_;

  // Disallow copy & assignment
  MeshDesc(const MeshDesc &);
  MeshDesc &operator=(const MeshDesc &);
};

std::ostream &
operator<<(std::ostream &, const MeshDesc &);

template <typename InputIterator>
const CoordSet &
reduce(InputIterator first, InputIterator last, const CoordSet &origin, CoordSet &target) {
  int index = 0;
  for (InputIterator it = first; it != last; ++it) {
    target.nodeadd(index++, const_cast<CoordSet &>(origin).getNode(*it));
  }
  return target;
}

const Elemset &
renumber_nodes(const std::map<int, int> &nodeIndices, Elemset &target);

template <typename InputIterator>
const Elemset &
reduce(InputIterator first, InputIterator last, const Elemset &origin, Elemset &target) {
  int index = 0;
  for (InputIterator it = first; it != last; ++it) {
    target.elemadd(index++, origin[*it]);
  }
  return target;
}

template <typename InputIterator>
const Elemset &
reduce(InputIterator first, InputIterator last, const std::map<int, int> &nodeIndices, const Elemset &origin, Elemset &target) {
  reduce(first, last, origin, target);
  return renumber_nodes(nodeIndices, target);
}

template <typename ValueType>
struct MeshSectionTraits {
  static int ValueType::* const renumbered_slot;
};

template <typename BCondInputIterator, typename BCondOutputIterator>
BCondOutputIterator
reduce(const std::map<int, int> &indices, BCondInputIterator first, BCondInputIterator last, BCondOutputIterator result) {
  typedef typename std::iterator_traits<BCondInputIterator>::value_type ValueType;
  int ValueType::* const slot = MeshSectionTraits<ValueType>::renumbered_slot;
  
  for (BCondInputIterator source = first; source != last; ++source) {
    const ValueType &value = *source;
    const std::map<int, int>::const_iterator it = indices.find(value.*slot);
    if (it != indices.end()) {
      ValueType newValue(value);
      newValue.*slot = it->second;
      *result++ = newValue;
    }
  }

  return result;
}

// EFrames from Elements

template <typename EFrameDataOutputIterator>
EFrameDataOutputIterator
copy_eframes(const Elemset & elemSet, EFrameDataOutputIterator result) {
  const int iElemEnd = elemSet.last();
  for (int iElem = 0; iElem != iElemEnd; ++iElem) {
    Element *elem = elemSet[iElem];
    if (!elem) {
      continue;
    } 
    const EFrame *eframe = elem->getFrame();
    if (eframe) {
      EFrameData data;
      data.elnum = iElem;
      std::memcpy(&data.frame, eframe, sizeof(data.frame));
      data.next = NULL;
      
      *result++ = data;
    }
  }
}

// Map to sequence conversion

template <typename MapIterator>
class value_iterator {
public:
  typedef typename std::bidirectional_iterator_tag iterator_category;
  typedef typename MapIterator::value_type::second_type value_type;
  typedef const value_type *pointer;
  typedef const value_type &reference;
  typedef typename MapIterator::difference_type difference_type;

  bool operator==(const value_iterator& other) const { return mapIt_ == other.mapIt_; }
  bool operator!=(const value_iterator& other) const { return !(*this == other); }

  const value_type &operator*() const { return mapIt_->second; }
  const value_type *operator->() const { return &mapIt_->second; }

  value_iterator &operator++() { ++mapIt_; return *this; }
  value_iterator operator++(int) { value_iterator temp(*this); ++(*this); return temp; }
  
  value_iterator &operator--() { --mapIt_; return *this; }
  value_iterator operator--(int) { value_iterator temp(*this); --(*this); return temp; }

  explicit value_iterator(MapIterator mapIt) : mapIt_(mapIt) {}

private:
  MapIterator mapIt_;
};

template <typename MapIterator>
value_iterator<MapIterator>
make_value_iterator(MapIterator mapIt) {
  return value_iterator<MapIterator>(mapIt);
}

#endif /* ROM_MESHUTILS_H */
