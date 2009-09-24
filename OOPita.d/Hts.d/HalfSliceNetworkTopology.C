#include "HalfSliceNetworkTopology.h"

namespace Pita {

namespace {

void
buildWriter(const Hts::SeedId & id, std::vector<Hts::SliceId> & writer) {
  switch (id.type()) {
    case Hts::UNDEFINED_SLICE:
      break;
    case Hts::MAIN_SEED:
      writer.push_back(Hts::SliceId(Hts::TAIL_FULL_SLICE, id.rank() - HalfSliceCount(1)));
      break;
    case Hts::LEFT_SEED:
      writer.push_back(Hts::SliceId(Hts::FORWARD_HALF_SLICE, id.rank() - HalfSliceCount(1)));
      break;
    case Hts::RIGHT_SEED:
      writer.push_back(Hts::SliceId(Hts::BACKWARD_HALF_SLICE, id.rank()));
      break;
    default:
      throw Fwk::InternalException(); 
  }
}

void
buildReader(const Hts::SeedId & id, std::vector<Hts::SliceId> & reader) {
  switch (id.type()) {
    case Hts::UNDEFINED_SLICE:
      break;
    case Hts::MAIN_SEED:
      reader.push_back(Hts::SliceId(Hts::BACKWARD_HALF_SLICE, id.rank() - HalfSliceCount(1)));
      reader.push_back(Hts::SliceId(Hts::FORWARD_HALF_SLICE, id.rank()));
      reader.push_back(Hts::SliceId(Hts::HEAD_FULL_SLICE, id.rank()));
      break;
    case Hts::LEFT_SEED:
      reader.push_back(Hts::SliceId(Hts::TAIL_FULL_SLICE, id.rank() - HalfSliceCount(1)));
      break;
    case Hts::RIGHT_SEED:
      reader.push_back(Hts::SliceId(Hts::HEAD_FULL_SLICE, id.rank()));
      break;
    default:
      throw Fwk::InternalException();
  }
}

void
buildSeed(const Hts::SliceId & id, std::vector<Hts::SeedId> & seed) {
  switch (id.type()) {
    case Hts::UNDEFINED_SEED:
      break;
    case Hts::FORWARD_HALF_SLICE:
      seed.push_back(Hts::SeedId(Hts::MAIN_SEED, id.rank()));
      seed.push_back(Hts::SeedId(Hts::LEFT_SEED, id.rank() + HalfSliceCount(1)));
      break;
    case Hts::BACKWARD_HALF_SLICE:
      seed.push_back(Hts::SeedId(Hts::MAIN_SEED, id.rank() + HalfSliceCount(1)));
      seed.push_back(Hts::SeedId(Hts::RIGHT_SEED, id.rank()));
      break;
    case Hts::HEAD_FULL_SLICE:
      seed.push_back(Hts::SeedId(Hts::MAIN_SEED, id.rank()));
      seed.push_back(Hts::SeedId(Hts::RIGHT_SEED, id.rank()));
      break;
    case Hts::TAIL_FULL_SLICE:
      seed.push_back(Hts::SeedId(Hts::MAIN_SEED, id.rank() + HalfSliceCount(1)));
      seed.push_back(Hts::SeedId(Hts::LEFT_SEED, id.rank() + HalfSliceCount(1)));
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
HalfSliceNetworkTopology::writer(const Hts::SeedId & id) const {
  return createIterator(writer_, id, &buildWriter);
}

HalfSliceNetworkTopology::SliceIteratorConst
HalfSliceNetworkTopology::reader(const Hts::SeedId & id) const {
  return createIterator(reader_, id, &buildReader);
}
 
HalfSliceNetworkTopology::SeedIteratorConst
HalfSliceNetworkTopology::seed(const Hts::SliceId & id) const {
  return createIterator(seed_, id, &buildSeed);
}


} // end namespace Pita
