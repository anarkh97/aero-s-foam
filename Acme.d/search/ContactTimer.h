// $Id$

#ifndef contact_timer_h_
#define contact_timer_h_

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif

class CString;

class ContactTimer {

 public:

  ContactTimer(
#ifndef CONTACT_NO_MPI
	       MPI_Comm& communicator
#else
	       int communicator
#endif
	       );

  ~ContactTimer();

  int Register_Timer( const CString );

  void Start_Timer( int timer_handle );
  void Stop_Timer( int timer_handle );

  void Output_Processor_Timings();
  void Output_Timings();

 private:

  // Don't allow copying 
  ContactTimer( const ContactTimer& );
  ContactTimer operator=(const ContactTimer& );

#ifndef CONTACT_NO_MPI
  MPI_Comm communicator;
#else
  int communicator;
#endif

  int allocated_size;
  int number_timers;
  double* times;
  double* start_times;
  CString** names;
  int* timing_calls;

  void Increase_Size();
  
};

#endif
