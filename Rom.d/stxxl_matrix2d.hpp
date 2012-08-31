#ifndef _STXXL_MATRIX2D_
#define _STXXL_MATRIX2D_

#include <stxxl/vector>

template <typename T>
class stxxl_matrix2d
{
    stxxl::vector<T> *v;
    stxxl::unsigned_type width, height;
    bool delete_v;

public:
    stxxl_matrix2d(stxxl::unsigned_type height, stxxl::unsigned_type width) : height(height), width(width)
    {
        v = new stxxl::vector<T>();
        delete_v = true;
        v->resize(width * height);
    }

    stxxl_matrix2d(stxxl::vector<T> *_v, stxxl::unsigned_type height, stxxl::unsigned_type width) : height(height), width(width), v(_v)
    {
        delete_v = false;
    }

    ~stxxl_matrix2d()
    {
      if(delete_v) delete v;
    }

    stxxl::unsigned_type cols() const
    {
        return width;
    }

    stxxl::unsigned_type rows() const
    {
        return height;
    }


    T & coeffRef(stxxl::unsigned_type rowId, stxxl::unsigned_type colId)
    {
        //col-major
        return (*v)[colId * height + rowId];
    }

    const T & coeff(stxxl::unsigned_type rowId, stxxl::unsigned_type colId) const
    {
        //col-major
        return (*v)[colId * height + rowId];
    }
};
#endif
