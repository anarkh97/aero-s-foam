#ifndef	_FULLSQUAREMATRIX_H_
#define	_FULLSQUAREMATRIX_H_
#include <Utils.d/NodeSpaceArray.h>
#include <stdio.h>
#include <Utils.d/MyComplex.h>
#include <Math.d/Vector.h>

template<class Scalar>
class GenFullSquareMatrix {
private:
	int	size;
	int     myval;
	Scalar*	value;
        int     length;
public:
        GenFullSquareMatrix(int, Scalar *l=0);
	GenFullSquareMatrix(); 
        GenFullSquareMatrix(GenFullSquareMatrix<Scalar> &m, Scalar s);
	~GenFullSquareMatrix();

        void setSize(int); 
        void changeSize(int i, int numMax); 
        void setMyval(int _myval) { myval = _myval; }

	template<class Scalar1> void operator = (const GenFullSquareMatrix<Scalar1> &M2);

        GenFullSquareMatrix<Scalar> operator * (Scalar v);
        GenFullSquareMatrix<Scalar> operator / (Scalar v);
        GenFullSquareMatrix<Scalar> operator + (const GenFullSquareMatrix<Scalar> &M2);
        GenFullSquareMatrix<Scalar> operator - (const GenFullSquareMatrix<Scalar> &M2);
        GenFullSquareMatrix<Scalar> &operator *= (Scalar v);
        GenFullSquareMatrix<Scalar> &operator /= (Scalar v);
        GenFullSquareMatrix<Scalar> &operator += (const GenFullSquareMatrix<Scalar> &M2); 
        GenFullSquareMatrix<Scalar> &operator -= (const GenFullSquareMatrix<Scalar> &M2);
        GenFullSquareMatrix<Scalar> &operator += (const Tensor_d2s0 &t);
	Scalar *operator[] (int row);
        int dim()    { return size; }
        int numRow() { return size; }
        int numCol() { return size; }
        void multiply(GenFullSquareMatrix<Scalar> &res, double d);
	void zero();
	void unitary();
	void unitaryDiag();
        void copy(GenFullSquareMatrix<Scalar> &m);

	void symmetrize();
        void print(const char *msg = "", const char *msg2="");
        void printDiagonals();

        Scalar* data() { return value; }
        const Scalar* data() const { return value; }
        Scalar getValue(int i) const { return value[i]; }
        void setValue(int i, Scalar s) { value[i] = s; }
        void copy(Scalar *d);

        void add(GenFullSquareMatrix<Scalar> &m, int *rc);
	template<class Scalar1, class Scalar2, class Scalar3> void multiply(GenVector<Scalar1>& a, GenVector<Scalar2>& b, Scalar3 c=1.0);
        void multiply(GenFullSquareMatrix<Scalar> &M2, GenFullSquareMatrix<Scalar> &result);
        void eigenVals(Scalar*);
        void eigenV(Scalar*);
        Scalar trace() {
        	Scalar res = 0;
        	for(int i = 0; i < dim(); ++i)
        		res += (*this)[i][i];
        	return res;
        }
        //void invert(GenFullSquareMatrix<Scalar>);
};

template<class Scalar>
inline
Scalar *
GenFullSquareMatrix<Scalar>::operator[](int row)
 { return value+size*row; }

typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;

#ifdef _TEMPLATE_FIX_
#include <Math.d/FullSquareMatrix.C>
#endif

#endif

