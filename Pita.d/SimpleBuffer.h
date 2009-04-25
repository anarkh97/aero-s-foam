#include <cstdlib>

template <typename Scalar>
class SimpleBuffer
{
  public:
  typedef Scalar DataType;
  
  SimpleBuffer();
  ~SimpleBuffer();

  Scalar * array() { return values_; }
  const Scalar * array() const { return values_; }
  const Scalar & operator[](size_t n) const { return values_[n]; }
  Scalar & operator[](size_t n) { return values_[n]; }
  size_t size() const { return numValues_; }
  void size(size_t numValues);
  void empty();

  private:
  static const size_t sizeOfScalar_ = sizeof(Scalar);

  Scalar * values_;
  size_t numValues_;
};

template <typename Scalar>
SimpleBuffer<Scalar>::SimpleBuffer() : values_(NULL), numValues_(0) {}

template <typename Scalar>
SimpleBuffer<Scalar>::~SimpleBuffer()
{
  std::free(values_);
}

template <typename Scalar>
void SimpleBuffer<Scalar>::size(size_t numValues)
{
  if (numValues == numValues_)
    return;
  if (numValues > 0)
  {
    numValues_ = numValues;
    values_ = (Scalar*) std::realloc(values_, numValues * sizeOfScalar_);
  }
  else
  {
    empty();
  }
}

template <typename Scalar>
void SimpleBuffer<Scalar>::empty()
{
  std::free(values_);
  values_ = NULL;
  numValues_ = 0;
}
