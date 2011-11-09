// $Id$

#include "Contact_Defines.h"

extern "C" {

#ifdef _WIN32

  void __stdcall FORTRAN(cnodetriangle_cpproj)
                        ( int&, Real*, Real*, Real*, int*,
                          Real*, Real*, Real&, Real&, Real& );

  void __stdcall FORTRAN(cnodetriangle_movsrch)
                        ( int&,  Real*, Real*, Real*, 
                          Real*, Real*, Real*, Real*, Real& );

  void __stdcall FORTRAN(cnodetriangle_cpproj_aug)
                        ( int&, 
                          Real*, Real*, Real*,
                          Real*, Real*, Real*,
                          int*, Real*, Real*, Real&, Real&, Real& );

  void __stdcall FORTRAN(cnodetriangle_movsrch_aug)
                        ( int&,  
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*,
                          Real*, Real&, Real&, Real&);

  void __stdcall FORTRAN(cnodetriangle_cpproj_aug_noauto)
                        ( int&, 
                          Real*, Real*, Real*,
                          Real*, Real*, Real*,
                          int*, Real*, Real*, Real&, Real&, Real& );

  void __stdcall FORTRAN(cnodetri_movsrch_aug_noauto)
                        ( int&,  
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*, 
                          Real*, Real*, Real*,
                          Real*, Real&, Real& );

  void __stdcall FORTRANNO_(oq4h8)
                           ( Real*, Real*, int*, int*, 
		             int*, int*, Real*, Real* );

#else

  void FORTRAN(cnodetriangle_cpproj)( int&, Real*, Real*, Real*, int*,
				      Real*, Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetriangle_movsrch)( int&, Real*, Real*, Real*, 
				       Real*, Real*, Real*, Real*, Real& );

  void FORTRAN(cnodetriangle_cpproj_aug)( int&, 
                                          Real*, Real*, Real*,
                                          Real*, Real*, Real*,
                                          int*, Real*, Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetriangle_movsrch_aug)( int&,  
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*,
                                           Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetriangle_cpproj_aug_noauto)( int&, 
                                          Real*, Real*, Real*,
                                          Real*, Real*, Real*,
                                          int*, Real*, Real*, Real&, Real&, Real& );

  void FORTRAN(cnodetri_movsrch_aug_noauto)( int&,  
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*, 
                                           Real*, Real*, Real*,
                                           Real*, Real&, Real& );

  void FORTRAN(cnodeline_cpproj)( int&, Real*, Real*, Real*, int*,
				  Real*, Real*, Real&, Real&, Real& );

  void FORTRANNO_(oq4h8)( Real*, Real*, int*, int*, int*, int*, Real*, Real* );

#endif

}

