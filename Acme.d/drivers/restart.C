// $Id: restart.C,v 2002.2 2002/12/23 23:31:20 khbrown Exp $

#include "restart.h"
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include "Contact_Defines.h"

#include <stdio.h>

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "Contact_Communication.h"
#endif

extern "C"
void FORTRAN(write_restart)( int& isize, Real* data, int& icomm,
                             int& my_proc_id, int& num_procs ) {

  char current_dir[1024];
  char trans[81];

#ifndef CONTACT_NO_MPI
  int  len_current_dir;
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#endif

  // get current working directory
  if (my_proc_id == 0) {
    char * dir_pointer = getenv("CONTACT_CWD");
    if ( NULL == dir_pointer ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "FATAL ERROR: the environment variable CONTACT_CWD" << endl
  	   << "   must be set to the test directory." << endl;
#endif
      exit(1);
    }
    strcpy(current_dir,dir_pointer);
#ifndef CONTACT_NO_MPI
    len_current_dir = strlen(current_dir)+1;
#endif
  }

#ifndef CONTACT_NO_MPI
  contact_broadcast(&len_current_dir,1,0,comm);
  contact_broadcast(current_dir,len_current_dir,0,comm);
#endif

  // open up restart file to write to
  strcat (current_dir,"/restart.data.out");
  if (num_procs>1) {
    sprintf(trans,".%d.%d",num_procs,my_proc_id);
    strcat(current_dir,trans);
  }
  ofstream data_file(current_dir);
  
  data_file << isize << endl;;
  for (int i = 0; i < isize; i++) {
    data_file << data[i] << endl;
  }
  data_file.close();
}

void FORTRAN(write_enf_restart)( int& isize, Real* data, int& icomm,
				 int& my_proc_id, int& num_procs ) {

  char current_dir[1024];
  char trans[81];

#ifndef CONTACT_NO_MPI
  int  len_current_dir;
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#endif

  // get current working directory
  if (my_proc_id == 0) {
    char * dir_pointer = getenv("CONTACT_CWD");
    if ( NULL == dir_pointer ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "FATAL ERROR: the environment variable CONTACT_CWD" << endl
  	   << "   must be set to the test directory." << endl;
#endif
      exit(1);
    }
    strcpy(current_dir,dir_pointer);
#ifndef CONTACT_NO_MPI
    len_current_dir = strlen(current_dir)+1;
#endif
  }

#ifndef CONTACT_NO_MPI
  contact_broadcast(&len_current_dir,1,0,comm);
  contact_broadcast(current_dir,len_current_dir,0,comm);
#endif

  // open up restart file to write to
  strcat (current_dir,"/enf.restart.data.out");
  if (num_procs>1) {
    sprintf(trans,".%d.%d",num_procs,my_proc_id);
    strcat(current_dir,trans);
  }
  ofstream data_file(current_dir);
  
  data_file << isize << endl;;
  for (int i = 0; i < isize; i++) {
    data_file << data[i] << endl;
  }
  data_file.close();
}

extern "C"
void FORTRAN(read_restart_buffer_size)( int& isize, int& icomm,
                                        int& my_proc_id, 
                                        int& num_procs ) {

  char current_dir[1024];
  char trans[81];
  
#ifndef CONTACT_NO_MPI
  int  len_current_dir;
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#endif
  // get current working directory
  if (my_proc_id==0) {
    char * dir_pointer = getenv("CONTACT_CWD");
    if ( NULL == dir_pointer ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "FATAL ERROR: the environment variable CONTACT_CWD" << endl
	   << "   must be set to the test directory." << endl;
#endif
      exit(1);
    }
    strcpy(current_dir,dir_pointer);
#ifndef CONTACT_NO_MPI
    len_current_dir = strlen(current_dir)+1;
#endif
  }

#ifndef CONTACT_NO_MPI
  contact_broadcast(&len_current_dir,1,0,comm);
  contact_broadcast(current_dir,len_current_dir,0,comm);
#endif

  // open restart file to read from
  strcat (current_dir,"/restart.data.in");
  if (num_procs>1) {
    sprintf(trans,".%d.%d",num_procs,my_proc_id);
    strcat(current_dir,trans);
  }
  ifstream data_file(current_dir);
  data_file >> isize;
  data_file.close();

}

extern "C"
void FORTRAN(read_enf_restart_buffer_size)( int& isize, int& icomm,
					    int& my_proc_id, 
					    int& num_procs ) {

  char current_dir[1024];
  char trans[81];
  
#ifndef CONTACT_NO_MPI
  int  len_current_dir;
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#endif

  // get current working directory
  if (my_proc_id==0) {
    char * dir_pointer = getenv("CONTACT_CWD");
    if ( NULL == dir_pointer ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "FATAL ERROR: the environment variable CONTACT_CWD" << endl
	   << "   must be set to the test directory." << endl;
#endif
      exit(1);
    }
    strcpy(current_dir,dir_pointer);
#ifndef CONTACT_NO_MPI
    len_current_dir = strlen(current_dir)+1;
#endif
  }

#ifndef CONTACT_NO_MPI
  contact_broadcast(&len_current_dir,1,0,comm);
  contact_broadcast(current_dir,len_current_dir,0,comm);
#endif

  // open restart file to read from
  strcat (current_dir,"/enf.restart.data.in");
  if (num_procs>1) {
    sprintf(trans,".%d.%d",num_procs,my_proc_id);
    strcat(current_dir,trans);
  }
  ifstream data_file(current_dir);
  data_file >> isize;
  data_file.close();

}

extern "C"
void FORTRAN(read_restart_data)( Real* data, int& icomm,
                                 int& my_proc_id, 
                                 int& num_procs ) {

  char current_dir[1024];
  char trans[81];
  
#ifndef CONTACT_NO_MPI
  int  len_current_dir;
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#endif
  // get current working directory
  if (my_proc_id==0) {
    char * dir_pointer = getenv("CONTACT_CWD");
    if ( NULL == dir_pointer ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "FATAL ERROR: the environment variable CONTACT_CWD" << endl
	   << "   must be set to the test directory." << endl;
#endif
      exit(1);
    }
    strcpy(current_dir,dir_pointer);
#ifndef CONTACT_NO_MPI
    len_current_dir = strlen(current_dir)+1;
#endif
  }

#ifndef CONTACT_NO_MPI
  contact_broadcast(&len_current_dir,1,0,comm);
  contact_broadcast(current_dir,len_current_dir,0,comm);
#endif

  // open restart file to read from
  strcat (current_dir,"/restart.data.in");
  if (num_procs>1) {
    sprintf(trans,".%d.%d",num_procs,my_proc_id);
    strcat(current_dir,trans);
  }
  ifstream data_file(current_dir);
  int isize;
  data_file >> isize;
  for (int i = 0; i < isize; i++) {
    data_file >> data[i];
  }
  data_file.close();

}

extern "C"
void FORTRAN(read_enf_restart_data)( Real* data, int& icomm,
				     int& my_proc_id, 
				     int& num_procs ) {

  char current_dir[1024];
  char trans[81];
  
#ifndef CONTACT_NO_MPI
  int  len_current_dir;
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#endif

  // get current working directory
  if (my_proc_id==0) {
    char * dir_pointer = getenv("CONTACT_CWD");
    if ( NULL == dir_pointer ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "FATAL ERROR: the environment variable CONTACT_CWD" << endl
	   << "   must be set to the test directory." << endl;
#endif
      exit(1);
    }
    strcpy(current_dir,dir_pointer);
#ifndef CONTACT_NO_MPI
    len_current_dir = strlen(current_dir)+1;
#endif
  }

#ifndef CONTACT_NO_MPI
  contact_broadcast(&len_current_dir,1,0,comm);
  contact_broadcast(current_dir,len_current_dir,0,comm);
#endif

  // open restart file to read from
  strcat (current_dir,"/enf.restart.data.in");
  if (num_procs>1) {
    sprintf(trans,".%d.%d",num_procs,my_proc_id);
    strcat(current_dir,trans);
  }
  ifstream data_file(current_dir);
  int isize;
  data_file >> isize;
  for (int i = 0; i < isize; i++) {
    data_file >> data[i];
  }
  data_file.close();

}


