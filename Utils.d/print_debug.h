#ifndef _PRINT_DEBUG_H_
#define _PRINT_DEBUG_H_

#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>

template<class Scalar>
void print_debug(GenDistrVector<Scalar> *vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec->size(); ++i) cerr << i << "," << vec->data()[i] << " "; cerr << endl;
 }
 else {
   for(int i=0; i<vec->size(); ++i) cerr << vec->data()[i] << " "; cerr << endl;
 }
}

template<class Scalar>
void print_debug(GenDistrVector<Scalar> &vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec.size(); ++i) cerr << i << "," << vec.data()[i] << " "; cerr << endl;
 }
 else {
   for(int i=0; i<vec.size(); ++i) cerr << vec.data()[i] << " "; cerr << endl;
 }
}

template<class Scalar>
void print_debug(GenVector<Scalar> *vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec->size(); ++i) cerr << i << "," << vec->data()[i] << " "; cerr << endl; 
 }
 else {
   for(int i=0; i<vec->size(); ++i) cerr << vec->data()[i] << " "; cerr << endl;
 }
}

template<class Scalar>
void print_debug(GenVector<Scalar> &vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec.size(); ++i) cerr << i << "," << vec.data()[i] << " "; cerr << endl;
 }
 else {
   for(int i=0; i<vec.size(); ++i) cerr << vec.data()[i] << " "; cerr << endl;
 }
}

template<class Scalar>
void print_debug(Scalar *vec, int len)
{
 for(int i=0; i<len; ++i) cerr << vec[i] << " "; cerr << endl;
}

#endif
