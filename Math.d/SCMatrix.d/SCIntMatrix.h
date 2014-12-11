#ifndef SCINTMATRIX_H_
#define SCINTMATRIX_H_

#ifdef NNLS_DEV
#include "SCBaseMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCBaseMatrix.h"
#endif

#include <string>

class SCIntMatrix : public SCBaseMatrix {
    public:
        SCIntMatrix(int context, int m, int n, int mb, int nb, bool pvec=false);
        SCIntMatrix(const SCIntMatrix& matrix);
        ~SCIntMatrix();

        int pivot(int *ip, int *desc_ip);
        int setMatrixColumn(int j, int *col);
        int zero();
        void write(std::string, int m=0, int n=0);
        void writeLocal(std::string filename);
        int identityPermutation();
        int * getMatrix() {return _matrix;};
        void setElement(int i, int j, int value);
        int getElementLocal(int i) {return _matrix[i];}
        int getElement(int i, int j);
        int getLocalElements(int *elems);
        int permute(char direc, char rowcol, SCIntMatrix &ip, int m=0, int n=0);
        int reorder(SCIntMatrix& order);
        void swap(int i, int j);
        int distributeVector();

    private:
        int * _matrix;     // Local _mlocal X _nlocal matrix
        bool _pvec;        // Flag to determine if this is a permutation vector. 
                           // Needs more memory allocated it if is.
        void init();
};

#endif // SCINTMATRIX_H_
