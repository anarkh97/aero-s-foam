#include "HalfSliceNetworkTopology.h"

namespace Pita {

namespace {

void
buildWriter(const Hs::SeedId & id, std::vector<Hs::SliceId> & writer) {
  switch (id.type()) {
    case Hs::UNDEFINED_SLICE:
      break;
    case Hs::MAIN_SEED:
      writer.push_back(Hs::SliceId(Hs::TAIL_FULL_SLICE, id.rank() - HalfSliceCount(1)));
      break;
    case Hs::LEFT_SEED:
      writer.push_back(Hs::SliceId(Hs::FORWARD_HALF_SLICE, id.rank() - HalfSliceCount(1)));
      break;
    case Hs::RIGHT_SEED:
      writer.push_back(Hs::SliceId(Hs::BACKWARD_HALF_SLICE, id.rank()));
      break;
    default:
      throw Fwk::InternalException(); 
  }
}

void
buildReader(const Hs::SeedId & id, std::vector<Hs::SliceId> & reader) {
  switch (id.type()) {
    case Hs::UNDEFINED_SLICE:
      break;
    case Hs::MAIN_SEED:
      reader.push_back(Hs::SliceId(Hs::BACKWARD_HALF_SLICE, id.rank() - HalfSliceCount(1)));
      reader.push_back(Hs::SliceId(Hs::FORWARD_HALF_SLICE, id.rank()));
      reader.push_back(Hs::SliceId(Hs::HEAD_FULL_SLICE, id.rank()));
      break;
    case Hs::LEFT_SEED:
      reader.push_back(Hs::SliceId(Hs::TAIL_FULL_SLICE, id.rank() - HalfSliceCount(1)));
      break;
    case Hs::RIGHT_SEED:
      reader.push_back(Hs::SliceId(Hs::HEAD_FULL_SLICE, id.rank()));
      break;
    default:
      throw Fwk::InternalException();
  }
}

void
buildSeed(const Hs::SliceId & id, std::vector<Hs::SeedId> & seed) {
  switch (id.type()) {
    case Hs::UNDEFINED_SEED:
      break;
    case Hs::FORWARD_HALF_SLICE:
      seed.push_back(Hs::SeedId(Hs::MAIN_SEED, id.rank()));
      seed.push_back(Hs::SeedId(Hs::LEFT_SEED, id.rank() + HalfSliceCount(1)));
      break;
    case Hs::BACKWARD_HALF_SLICE:
      seed.push_back(Hs::SeedId(Hs::MAIN_SEED, id.rank() + HalfSliceCount(1)));
      seed.push_back(Hs::SeedId(Hs::RIGHT_SEED, id.rank()));
      break;
    case Hs::HEAD_FULL_SLICE:
      seed.push_back(Hs::SeedId(Hs::MAIN_SEED, id.rank()));
      seed.push_back(Hs::SeedId(Hs::RIGHT_SEED, id.rank()));
      break;
    case Hs::TAIL_FULL_SLICE:
      seed.push_back(Hs::SeedId(Hs::MAIN_SEED, id.rank() + HalfSliceCount(1)));
      seed.push_back(Hs::SeedId(Hs::LEFT_SEED, id.rank() + HalfSliceCount(1)));
  }
}

} // end anonymous namespace


template <typename K, typename T>
inline
HalfSliceNetworkTopology::GenIteratorConst<T>
HalfSliceNetworkTopology::createIterator(std::map<K, std::vector<T> > & m, const K & id, 
                                         void (* const buildFunc)(const K &, std::vector<T> &)) const {
  typename std::map<K, std::vector<T> >::const_iterator it = m.find(id);
  if (it != m.end()) {
    return GenIteratorConst<T>(it->second.begin(), it->second.end());
  }
  std::vector<T> tmp;
  buildFunc(id, tmp);
  if (tmp.begin() != tmp.end()) {
    it = m.insert(std::make_pair(id, tmp)).first;
    return GenIteratorConst<T>(it->second.begin(), it->second.end()); 
  } else {
    return GenIteratorConst<T>(tmp.begin(), tmp.end());
  }
}
 
HalfSliceNetworkTopology::SliceIteratorConst
HalfSliceNetworkTopology::writer(const Hs::SeedId & id) const {
  return createIterator(writer_, id, &buildWriter);
}

HalfSliceNetworkTopology::SliceIteratorConst
HalfSliceNetworkTopology::reader(const Hs::SeedId & id) const {
  return createIterator(reader_, id, &buildReader);
}
 
HalfSliceNetworkTopology::SeedIteratorConst
HalfSliceNetworkTopology::seed(const Hs::SliceId & id) const {
  return createIterator(seed_, id, &buildSeed);
}


} // end namespace Pita
