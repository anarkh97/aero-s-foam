// $Id$

#include "Contact_Defines.h"

extern "C" {

#ifdef WIN32

  void __stdcall FORTRAN(contact_make_rank)
                        ( int&, int&, int&, Real*, int*, int*, int* );

  void __stdcall FORTRAN(contact_get_bound)
                        ( int&, Real*, int*, int&, Real*, Real*, 
                          int&, int*, int*, int&, Real* );

  void __stdcall FORTRAN(contact_make_list)
                        ( int&, int*, int*, int*, int*, 
                          int&, int*, int&, int&, int& );

  void __stdcall FORTRAN(contact_make_rank_t)
                        ( int&, int&, int&, Real*, Real*, int*, int*, int* ,
                          int* ,int& );

  void __stdcall FORTRAN(contact_get_bound_t)
                        ( int&, Real*, int*, int&, Real*, Real*, 
                          int&, int*, int*, int&, Real*, int& );

 void __stdcall FORTRAN(contact_indexx)
                        ( int&, Real*, int*, int& );

 void __stdcall FORTRAN(contact_indexx_float)
                        ( int&, float*, int*, int& );

  void __stdcall FORTRAN(contact_rank)
                        ( int&, int*, int*, int& );

#else

  void FORTRAN(contact_make_rank)( int&, int&, int&, Real*, int*, int*, int*);

  void FORTRAN(contact_get_bound)( int&, Real*, int*, int&, Real*, Real*, int&, 
				   int*, int*, int&, Real*);

  void FORTRAN(contact_make_list)( int&, int*, int*, int*, int*, int&, int*, 
				   int&, int&, int& );
				    
  void FORTRAN(contact_make_rank_t)( int&, int&, int&, Real*, Real*, int*, int*, int*,
                                     int*,int&  );

  void FORTRAN(contact_get_bound_t)( int&, Real*, int*, int&, Real*, Real*, int&, 
				     int*, int*, int&, Real*, int& );

  void FORTRAN(contact_indexx)
                        ( const int&, Real*, int*, const int& );

  void FORTRAN(contact_indexx_float)
                        ( const int&, float*, int*, const int& );

  void FORTRAN(contact_rank)
                        ( const int&, int*, int*, const int& );

#endif

}
