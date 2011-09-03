#ifndef ROM_MESHOUTPUT_H
#define ROM_MESHOUTPUT_H

#include <Driver.d/StructProp.h>
#include <Driver.d/EFrameData.h>
#include <Element.d/Element.h>

class NLMaterial;
class ExpMat;

#include <ostream>
#include <iterator>
#include <string>
#include <sstream>
#include <utility>
#include <stdexcept>

namespace Rom {

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

std::ostream &
operator<<(std::ostream &, const ExpMat &);

// Sections from atoms

template <typename ValueType, typename TagType>
struct InputFileSectionHelper {
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

  return out;
}


struct EmptyTag {};

struct SampleNodeTag {};

template <>
inline
int
InputFileSectionHelper<int, SampleNodeTag>::transformation(const int &v) {
  return v + 1;
}

template <typename TagType>
struct InputFileSectionHelper<std::pair<const int, typename TagType::SecondType>, TagType> {
  typedef std::pair<const int, typename TagType::SecondType> ValueType;
  static const std::string &header(TagType);
  static std::string transformation(const ValueType &);
};

template <typename TagType>
std::string
InputFileSectionHelper<std::pair<const int, typename TagType::SecondType>, TagType>::transformation(const ValueType &p) {
  std::ostringstream result;
  result << p.first + 1 << " " << TagType::valueTransformation(p.second);
  return result.str();
}

struct ElementPressureTag {
  typedef double SecondType;
  static double valueTransformation(double x) { return x; }
};

struct MatUsageTag {
  typedef int SecondType;
  static int valueTransformation(int i) { return i + 1; }
};

struct MatLawTag {
  typedef NLMaterial* SecondType;
  static const ExpMat &valueTransformation(const NLMaterial *);
};


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

} /* end namespace Rom */

#endif /* ROM_MESHOUTPUT_H */
