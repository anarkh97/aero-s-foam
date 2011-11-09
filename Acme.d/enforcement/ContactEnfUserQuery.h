// $Id$

#include "Contact_Defines.h"


extern "C" {
  
  void FORTRAN(userquery_table_last_abscissa)(void *obj, int *id, Real* value);
  void FORTRAN(userquery_table_interpolate)(void *obj, int *id, 
                                            Real* abscissa, Real* ordinate);
  void FORTRAN(userquery_number_of_nodes)(void *obj, int *nnodes);
  void FORTRAN(userquery_node_state_data)(void *obj, int* id, int* node, 
                                          int* offset, Real* value);
  void FORTRAN(userset_node_state_data)(void *obj, int* id, int* node,
                                        int* offset, Real* value);
  void FORTRAN(userset_nfi_failure_tied)(void *enf_obj, int* model_id, 
                                         void *nfi_obj);

}
