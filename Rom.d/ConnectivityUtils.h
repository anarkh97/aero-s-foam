#ifndef ROM_CONNECTIVITYUTILS_H
#define ROM_CONNECTIVITYUTILS_H

#include <Utils.d/Connectivity.h>

#include <set>
#include <algorithm>
#include <cassert>

template <typename InputIterator, typename OutputIterator>
OutputIterator
neighborhood(const Connectivity &connectivity,
             InputIterator first, InputIterator last,
             OutputIterator result) {
  Connectivity & conn_fix = const_cast<Connectivity &>(connectivity);
  std::set<int> found;

  for (InputIterator it = first; it != last; ++it) {
    const int key = *it;
    assert(key < conn_fix.csize()); 
    found.insert(key);

    const int *begin = conn_fix[key];
    found.insert(begin, begin + conn_fix.num(key));
  }

  return std::copy(found.begin(), found.end(), result);
}

template <typename InputIterator, typename OutputIterator>
OutputIterator
connections(const Connectivity &connectivity,
             InputIterator first, InputIterator last,
             OutputIterator result) {
  Connectivity & conn_fix = const_cast<Connectivity &>(connectivity);
  std::set<int> found;

  for (InputIterator it = first; it != last; ++it) {
    const int key = *it;
    assert(key < conn_fix.csize()); 

    const int *begin = conn_fix[key];
    found.insert(begin, begin + conn_fix.num(key));
  }

  return std::copy(found.begin(), found.end(), result);
}

#endif /* ROM_CONNECTIVITYUTILS_H */
