#ifndef _STXXL_MATRIX2D_
#define _STXXL_MATRIX2D_

// simple wrapper to view stxxl::vector (or std::vector) as a 2d matrix

template <typename vector_type>
class stxxl_matrix2d
{
    vector_type *v;
    typename vector_type::size_type width, height;
    bool delete_v;

public:
    stxxl_matrix2d(typename vector_type::size_type height, typename vector_type::size_type width) : height(height), width(width)
    {
        v = new vector_type();
        delete_v = true;
        v->resize(width * height);
    }

    stxxl_matrix2d(vector_type *_v, typename vector_type::size_type height, typename vector_type::size_type width) : height(height), width(width), v(_v)
    {
        delete_v = false;
    }

    ~stxxl_matrix2d()
    {
      if(delete_v) delete v;
    }

    typename vector_type::size_type cols() const
    {
        return width;
    }

    typename vector_type::size_type rows() const
    {
        return height;
    }

    typename vector_type::reference coeffRef(typename vector_type::size_type rowId, typename vector_type::size_type colId)
    {
        //col-major
        return (*v)[colId * height + rowId];
    }

    typename vector_type::const_reference coeff(typename vector_type::size_type rowId, typename vector_type::size_type colId) const
    {
        //col-major
        return (*v)[colId * height + rowId];
    }
};
#endif
