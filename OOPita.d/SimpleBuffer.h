#ifndef PITA_SIMPLEBUFFER_H
#define PITA_SIMPLEBUFFER_H

namespace Pita {

template <typename Scalar>
class SimpleBuffer {
public:
  typedef Scalar DataType;
  
  explicit SimpleBuffer(size_t numValues = 0);
  ~SimpleBuffer();
  
  size_t size() const { return numValues_; }
  void sizeIs(size_t numValues);

  Scalar * array() { return values_; }
  const Scalar * array() const { return values_; }
  
  const Scalar & operator[](size_t n) const { return values_[n]; }
  Scalar & operator[](size_t n) { return values_[n]; }
  
private:
  size_t numValues_;
  Scalar * values_;
};

template <typename Scalar>
inline
SimpleBuffer<Scalar>::SimpleBuffer(size_t numValues) :
  numValues_(numValues),
  values_(numValues > 0 ? new Scalar[numValues] : NULL)
{}

template <typename Scalar>
inline
SimpleBuffer<Scalar>::~SimpleBuffer() {
  delete[] values_;
}

template <typename Scalar>
inline
void SimpleBuffer<Scalar>::sizeIs(size_t numValues) {
  if (numValues != numValues_) {
    Scalar * v = numValues > 0 ? new Scalar[numValues] : NULL;
    delete[] values_;
    values_ = v;
    numValues_ = numValues;
  }
}

} // end namespace Pita

#endif /* PITA_SIMPLEBUFFER_H */
