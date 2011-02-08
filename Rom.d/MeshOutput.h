#ifndef ROM_MESHOUTPUT_H
#define ROM_MESHOUTPUT_H

#include <Driver.d/StructProp.h>
#include <Driver.d/EFrameData.h>
#include <Element.d/Element.h>

#include <ostream>
#include <iterator>
#include <string>

// Specific sections

std::ostream &
operator<<(std::ostream &, const CoordSet &);

std::ostream &
operator<<(std::ostream &, const Elemset &);

std::ostream &
operator<<(std::ostream &, const SPropContainer &);

// Atoms

std::ostream &
operator<<(std::ostream &, const Attrib &);

std::ostream &
operator<<(std::ostream &, const BCond &);

std::ostream &
operator<<(std::ostream &, const EFrameData &);

// Sections from atoms

template <typename ValueType, typename TagType> struct InputFileSectionHelper {
  static const std::string &header(TagType);
  static ValueType transformation(const ValueType &v) { return v; }
};

template <typename InputIterator, typename TagType>
class InputFileSection {
public:
  typedef typename std::iterator_traits<InputIterator>::value_type ValueType;

  const std::string &header() const {
    return InputFileSectionHelper<ValueType, TagType>::header(tag_);
  }

  InputIterator begin() const { return first_; }
  InputIterator end()   const { return last_;  }

  InputFileSection(InputIterator first, InputIterator last, TagType tag) :
    first_(first), last_(last), tag_(tag)
  {}

private:
  InputIterator first_, last_;
  TagType tag_;
};

template <typename InputIterator, typename TagType>
std::ostream &
operator<<(std::ostream &out, const InputFileSection<InputIterator, TagType> &source) {
  out << source.header() << "\n";
  InputIterator itEnd = source.end();
  for (InputIterator it = source.begin(); it != itEnd; ++it) {
    typedef typename InputFileSection<InputIterator, TagType>::ValueType ValueType;
    out << InputFileSectionHelper<ValueType, TagType>::transformation(*it) << "\n";
  }
}

struct EmptyTag {};

struct SampleNodeTag {};

template <>
inline
int
InputFileSectionHelper<int, SampleNodeTag>::transformation(const int &v) {
  return v + 1;
}

// Convenience functions

template <typename InputIterator>
InputFileSection<InputIterator, EmptyTag>
make_section(InputIterator first, InputIterator last) {
  return InputFileSection<InputIterator, EmptyTag>(first, last, EmptyTag());
}

template <typename InputIterator, typename TagType>
InputFileSection<InputIterator, TagType>
make_section(InputIterator first, InputIterator last, TagType tag) {
  return InputFileSection<InputIterator, TagType>(first, last, tag);
}


#endif /* ROM_MESHOUTPUT_H */
