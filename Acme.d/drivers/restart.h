// $Id: restart.h,v 2002.1 2002/01/04 16:24:32 khbrown Exp $

#ifndef _restart_h_
#define _restart_h_

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif
#include "Contact_Defines.h"

extern "C" {
  
  void FORTRAN(write_restart)( int&, Real*, int&, int&, int& );
  void FORTRAN(read_restart_buffer_size)( int&, int&, int&, int& );
  void FORTRAN(read_restart_data)( Real*, int&, int&, int& );
  void FORTRAN(write_enf_restart)( int&, Real*, int&, int&, int& );
  void FORTRAN(read_enf_restart_buffer_size)( int&, int&, int&, int& );
  void FORTRAN(read_enf_restart_data)( Real*, int&, int&, int& );

}

#endif
