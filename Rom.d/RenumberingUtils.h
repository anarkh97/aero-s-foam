#ifndef ROM_RENUMBERINGUTILS_H
#define ROM_RENUMBERINGUTILS_H

#include <utility>
#include <map>

template <typename InputIterator, typename OutputIterator>
OutputIterator
inverse_numbering(InputIterator first, InputIterator last,
                  OutputIterator result, int firstIndex = 0) {
  int index = firstIndex;

  for (InputIterator it = first; it != last; ++it) {
    *result++ = std::pair<const int, int>(*it, index++);
  }

  return result;
}

template <typename InputIterator, typename OutputIterator>
OutputIterator
renumber(const std::map<int, int> &renum, InputIterator first, InputIterator last, OutputIterator result) {
  for (InputIterator it = first; it != last; ++it) {
    std::map<int, int>::const_iterator r = renum.find(*it);
    if (r != renum.end()) {
      *result++ = r->second;
    }
  }

  return result;
}

#endif /* ROM_RENUMBERINGUTILS_H */
