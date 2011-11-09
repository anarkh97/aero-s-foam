// $Id: cxxdriver.C,v 2002.101 2004/06/18 17:47:04 mwglass Exp $

// a test comment to force a commit
// another test comment to force a commit
// yet another test comment to force a commit
// another yet another test comment to force a commit

//
// NOTE: THESE DEFINES MUST AGREE WITH THE DEFINES IN read_mesh.F
//
#include "params.h"

#include "ContactSearch.h"
#include "Contact_Defines.h"
#include "output_interactions.h"
#include <stdlib.h>
#include <unistd.h>
#include <iostream.h>
#include <fstream.h>
#include "contact_assert.h"
#ifndef CONTACT_NO_EXODUS_OUTPUT
#include "exodusII.h"
#endif
#include "restart.h"
#include <string.h>
#include <assert.h>
#include "read_data.h"
#include "test_init.h"

#ifndef CONTACT_NO_MPI
#include <stdio.h>
#include "mpi.h"
#include "Contact_Communication.h"
#endif

#include "ContactTDEnforcement.h"
#include "ContactGapRemoval.h"
#include "ContactTiedKinematics.h"
#include "ContactVolumeTransfer.h"
#include "ContactTDFaceFaceEnf.h"
#include "ContactMPCs.h"
#include "ContactTDUserSubTypes.h"
#include "usersubs.h"

extern "C" int MPI_Init( int*, char*** );
extern "C" void FORTRAN(register_user_subs)();

int main( int argc, char* argv[] )
{
  // set list of nodes for debugging
  int num_debug_nodes = 2;
  int debug_nodes[2];
  debug_nodes[0] = 11;
  debug_nodes[1] = 10;
  //debug_nodes[2] = 1;
  // if this is compiled with MPI, then start up MPI
#ifndef CONTACT_NO_MPI
  // initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_proc_id = contact_processor_number(comm);
  int num_procs = contact_number_of_processors(comm);
  MPI_INT_TYPE icomm = MPI_COMM_C2F(comm);

  // vars for parallel
  int  num_comm_partners(0);
  int* comm_proc_ids = new int[MAX_PROCS];
  int* number_nodes_to_partner = new int[MAX_PROCS];
  int* comm_node = new int[MAX_PROCS*MAX_NODES];

#else
  
  int comm       = 0;
  int icomm      = 0;
  int num_procs  = 1;
  int my_proc_id = 0;

  // parallel vars, but in serial
  int* comm_proc_ids = new int[1];
  int* number_nodes_to_partner = new int[1];
  int* comm_node = new int[1];
  int  num_comm_partners     = 0;
  comm_proc_ids[0]           = 0;
  number_nodes_to_partner[0] = 0;
  comm_node[0]               = 0;

#endif

  if( my_proc_id == 0 ){
    cout << "ACME Version: " << ACME_Version() << endl;
    cout << "ACME Date:    " << ACME_VersionDate() << endl;
    cout << flush;
  }

  int * elem_ignore_list = new int[MAX_ELEMENTS*(MAX_STEPS+1)];
  memset(elem_ignore_list,0,sizeof(int)*MAX_ELEMENTS*(MAX_STEPS+1));
  int * num_elem_ignore = new int[MAX_STEPS];
  memset(num_elem_ignore,0,sizeof(int)*MAX_STEPS);
  int * face_ignore_list = new int[MAX_ELEMENTS*(MAX_STEPS+1)];
  memset(face_ignore_list,0,sizeof(int)*MAX_ELEMENTS*(MAX_STEPS+1));
  int * num_face_ignore = new int[MAX_STEPS];
  memset(num_face_ignore,0,sizeof(int)*MAX_STEPS);
  int * node_ignore_list = new int[MAX_NODES*(MAX_STEPS+1)];
  memset(node_ignore_list,0,sizeof(int)*MAX_NODES*(MAX_STEPS+1));
  int * num_node_ignore = new int[MAX_STEPS];
  memset(num_node_ignore,0,sizeof(int)*MAX_STEPS);
  
  int* face_conn = new int[8*MAX_ELEMENTS];
  int* element_conn = new int[8*MAX_ELEMENTS];
  int* node_ids = new int [2*MAX_NODES];
  int* node_eids = new int[MAX_NODES];
  int* face_ids = new int [2*MAX_ELEMENTS];
  int* element_ids = new int[2*MAX_ELEMENTS];
  int node_block_ids[MAX_BLOCKS];
  int node_block_types[MAX_BLOCKS];
  int face_block_types[MAX_BLOCKS];
  int element_block_types[MAX_BLOCKS];
  int faces_per_block[MAX_BLOCKS];
  int nodes_per_block[MAX_BLOCKS];
  int elements_per_block[MAX_BLOCKS];
  
  int num_node_deaths_per_block[MAX_BLOCKS];
  int node_deaths_global_ids[2*MAX_NODES];
  int num_face_deaths_per_block[MAX_BLOCKS];
  int face_deaths_global_ids[2*MAX_ELEMENTS];
  int num_element_deaths_per_block[MAX_BLOCKS];
  int element_deaths_global_ids[2*MAX_ELEMENTS];
  int num_node_births_per_block[MAX_BLOCKS];
  int node_births_exodus_ids[MAX_NODES];
  int node_births_global_ids[2*MAX_NODES];
  int number_face_births_per_block[MAX_BLOCKS];
  int face_births_global_ids[2*MAX_ELEMENTS];
  int face_births_connectivity[2*8*MAX_ELEMENTS];
  int number_element_births_per_block[MAX_BLOCKS];
  int element_births_global_ids[2*MAX_ELEMENTS];
  int element_births_connectivity[2*8*MAX_ELEMENTS];
  int num_node_exports=0;
  int node_export_gids[2*MAX_NODES];
  int node_export_pids[MAX_NODES];
  int num_face_exports=0;
  int face_export_gids[2*MAX_ELEMENTS];
  int face_export_pids[MAX_ELEMENTS];
  int num_elem_exports=0;
  int elem_export_gids[2*MAX_ELEMENTS];
  int elem_export_pids[MAX_ELEMENTS];
  
  double* initial_position = new double[3*MAX_NODES];
  double* position1 = new double[3*MAX_NODES];
  double* position2 = new double[3*MAX_NODES];
  double* nodal_vars = new double[MAX_NODES*MAX_NODAL_VARS];
  double* element_vars = new double[MAX_ELEMENTS*MAX_ELEMENT_VARS];
  
  int as_type[MAX_ANALYTIC_SURFACES];
  double as_data[8*MAX_ANALYTIC_SURFACES];
  double search_data[ContactSearch::NSIZSD*MAX_BLOCKS*MAX_BLOCKS];
  double enforcement_data[5*MAX_BLOCKS*MAX_BLOCKS];
  double enforcement_data_vars[6];
  double* enforcement_vars = new double[3*MAX_NODES];
  double* attributes = new double[MAX_ELEMENTS*MAX_ATTRIBUTES];
  int* num_attributes = new int[MAX_BLOCKS];
  double* node_radius = new double[MAX_NODES];
  double* shell_thickness = new double[MAX_ELEMENTS];
  double* shell_offset = new double[MAX_ELEMENTS];
  
  int number_of_tables = 0;
  int table_ids[MAX_TABLES];
  int table_num_points[MAX_TABLES];
  Real table_abscissas[MAX_TABLES*MAX_TABLE_POINTS];
  Real table_ordinates[MAX_TABLES*MAX_TABLE_POINTS];
  
  int multiple_interaction_status = 0;
  int compute_node_areas = 0;
  int restart_flag = 0;
  int normal_smoothing_status = 0;
  double smoothing_data[3];
  double interaction_data[1];
  
  ContactSearch::ContactErrorCode error;

  int num_enf_models = 0;
  int enf_model_ids[MAX_ENFORCEMENT_MODELS];
  int enf_model_types[MAX_ENFORCEMENT_MODELS];
  int enf_model_sub_num[MAX_ENFORCEMENT_MODELS];
  int enf_model_sub_ids[MAX_ENFORCEMENT_MODELS*MAX_ENF_MODEL_SUBS];
  char enf_model_sub_names[MAX_ENFORCEMENT_MODELS*MAX_ENF_MODEL_SUBS*64];
  Real* enf_model_data = NULL;

  int idim = 0;
  int search_type = -1;
  int number_node_blocks = 0;
  int number_face_blocks = 0;
  int num_as = 0;
  int number_element_blocks = 0;
  int exodus_out = 0;

  int i,j;
  for( i=0 ; i<MAX_BLOCKS ; i++ ) node_block_ids[i] = i+1;

  int ierror = ACME_MPI_Compatibility(MPI_COMPILE);
  if( ierror ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    cerr << "Incompatible MPI compilation" << endl;
#endif
    exit(1);
  }
  FORTRAN(register_user_subs)();

  int num_enf_type_1_vars = 3;
  int num_nodal_vars   = MAX_NODAL_VARS;
  int num_element_vars = MAX_ELEMENT_VARS;
  int enforcement_type = 0;
  int num_load_steps = 0;
  int use_config = 0;
  int start_step = 1;
  int topology_update = 0;
  enf_model_data = new Real[MAX_ENFORCEMENT_MODELS*MAX_ENF_MODEL_DATA_VARS];
  memset( enf_model_data, 0, 
	  MAX_ENFORCEMENT_MODELS*MAX_ENF_MODEL_DATA_VARS*sizeof(Real) );
  // read in the mesh and contact data
  FORTRAN(test_init) (idim, 
	              num_load_steps,
	              use_config,
	              search_type, 
	              restart_flag, 
	              start_step,
	              multiple_interaction_status, 
	              interaction_data,
	              normal_smoothing_status, 
	              smoothing_data,
	              compute_node_areas,
	              num_as, 
	              as_type, 
	              as_data, 
	              search_data, 
	              number_of_tables,
	              table_ids,
	              table_num_points,
	              table_abscissas,
	              table_ordinates,
	              enforcement_type, 
	              &(enforcement_data[0]),
	              &(enforcement_data_vars[0]),
	              num_enf_models, 
	              enf_model_types, 
	              enf_model_ids,
  	              enf_model_sub_num,
  	              enf_model_sub_ids,
 	              enf_model_sub_names,
  	              enf_model_data, 
	              num_node_ignore, 
	              node_ignore_list, 
	              num_face_ignore, 
	              face_ignore_list, 
	              num_elem_ignore, 
	              elem_ignore_list, 
	              topology_update,
	              my_proc_id, 
	              num_procs, 
                      icomm);


  // modify the original topology for any ignored
  // faces, this will be the starting topology
  FORTRAN(start_topology)(*num_node_ignore,
			  node_ignore_list,
			  *num_face_ignore,
			  face_ignore_list,
			  *num_elem_ignore,
			  elem_ignore_list,
                          
			  number_node_blocks,
			  node_block_types,
			  nodes_per_block,
			  node_eids,
			  node_ids,
                          
			  number_face_blocks,
			  face_block_types,
			  faces_per_block,
			  face_ids,
			  face_conn,
                          
			  number_element_blocks,
			  element_block_types,
			  elements_per_block,
			  element_ids,
			  element_conn,
                          
                          initial_position,
                          nodal_vars,
                          num_nodal_vars,
			  element_vars,
			  num_element_vars,
			  attributes,
			  num_attributes,
                          node_radius,
                          shell_thickness,
                          shell_offset,
                          
			  num_comm_partners,
			  comm_proc_ids,
			  number_nodes_to_partner,
			  comm_node,
                          
			  my_proc_id,
			  num_procs);
   
  // get total number of nodes, faces, and elements on this processor
  int tot_num_nodes = 0;
  for (i = 0; i < number_node_blocks; i++)
    tot_num_nodes += nodes_per_block[i];
  int tot_num_faces = 0;
  for (i = 0; i < number_face_blocks; i++) 
    tot_num_faces += faces_per_block[i];
  int tot_num_elems = 0;
  for (i = 0; i < number_element_blocks; i++ )
    tot_num_elems += elements_per_block[i];

  ContactSearch* search;
  ContactTDEnforcement* tdenf;
  ContactGapRemoval*    gr;
  ContactTiedKinematics* tiedkin;
  ContactVolumeTransfer* voltrans;
  ContactMPCs* mpcenf;
#ifdef CONTACT_TD_FACE_FACE_ENF
  ContactTDFaceFaceEnf* tdffenf;
#endif
  ContactEnforcement* enf_object = NULL;

  if( restart_flag != 1 && restart_flag != 3) {
    
    // I'm going to do an explicit cast of int to enum.  This works only if
    // the compiler doesn't use a smaller size for the enum (current none
    // do).  The following preconditions will fail if we ever get a 
    // compiler that does this.
    PRECONDITION( sizeof(int) == sizeof(ContactSearch::ContactNode_Type) );
    PRECONDITION( sizeof(int) == sizeof(ContactSearch::ContactFace_Type) );
    PRECONDITION( sizeof(int) == sizeof(ContactSearch::ContactEdge_Type) );
    PRECONDITION( sizeof(int) == sizeof(ContactSearch::AnalyticSurface_Type) );

    // construct search    
    int num_states = 1;
    search = new ContactSearch(idim, num_states, num_as,
			       number_node_blocks,
			       (ContactSearch::ContactNode_Type*)
			       node_block_types,
			       nodes_per_block,&(node_eids[0]),
			       &(node_ids[0]),
			       initial_position,
			       number_face_blocks,
			       (ContactSearch::ContactFace_Type*)
			       face_block_types,  
			       &(faces_per_block[0]),&(face_ids[0]),
                               &(face_conn[0]),
			       shell_offset,
			       number_element_blocks,
			       (ContactSearch::ContactElement_Type*)
			       element_block_types,
			       elements_per_block,
                               &(element_ids[0]),
                               &(element_conn[0]),
			       num_comm_partners,comm_proc_ids, 
			       number_nodes_to_partner, 
			       comm_node, comm, error);
    
    if( error ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "Error creating Search " << i << " Error Code = "
	   << error << endl;
      for( i=0 ; i<search->Number_of_Errors() ; i++ )
        cerr << search->Error_Message(i+1) << endl;
#endif
      delete search;
      exit(1);
    }

    // setup analytic surfaces
    for( i=0 ; i<num_as ; i++ ) {
      ContactSearch::AnalyticSurface_Type astype =
	(ContactSearch::AnalyticSurface_Type) as_type[i];
      error = search->Add_Analytic_Surface( astype, &as_data[8*i] );
      if( error ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	cerr << "Error adding ANALYTICAL SURFACE " << i << " Error Code = "
	     << error << endl;
	for( i=0 ; i<search->Number_of_Errors() ; i++ )
	  cerr << search->Error_Message(i+1) << endl;
#endif
	delete search;
	exit(1);
      }
    }
    
    // Set multiple interaction option
    if( search->Set_Search_Option( ContactSearch::MULTIPLE_INTERACTIONS,
				   (ContactSearch::Search_Option_Status) 
				   multiple_interaction_status,
				   interaction_data ) ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      for( i=0 ; i<search->Number_of_Errors() ; i++ )
	cerr << search->Error_Message(i+1) << endl;
#endif
      delete search;
      exit(1);
    }
    
    // Set normal smoothing option
    if( search->Set_Search_Option( ContactSearch::NORMAL_SMOOTHING,
				   (ContactSearch::Search_Option_Status) 
				   normal_smoothing_status,
				   smoothing_data ) ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      for( i=0 ; i<search->Number_of_Errors() ; i++ )
	cerr << search->Error_Message(i+1) << endl;
#endif
      delete search;
      exit(1);
    }
    
    // Set compute_node_areas option
    if( search->Set_Search_Option( ContactSearch::COMPUTE_NODE_AREAS,
				   (ContactSearch::Search_Option_Status) 
				   compute_node_areas,
				   NULL ) ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      for( i=0 ; i<search->Number_of_Errors() ; i++ )
	cerr << search->Error_Message(i+1) << endl;
#endif
      delete search;
      exit(1);
    }
    
    // set search data for all contact entities
    int num_entity_keys = num_as +
                          number_node_blocks + 
                          number_face_blocks + 
                          number_element_blocks;
    if (number_face_blocks+number_element_blocks) num_entity_keys--;
    if( search->Check_Search_Data_Size(ContactSearch::NSIZSD,num_entity_keys) ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      for( i=0 ; i<search->Number_of_Errors() ; i++ )
	cerr << search->Error_Message(i+1) << endl;
#endif
      delete search;
      exit(1);
    }
    search->Set_Search_Data(search_data);

    // Add the tables
    int table_data_offset = 0;
    for( i=0 ; i<number_of_tables ; i++ ){
      search->Add_Table( table_ids[i], table_num_points[i],
			 &(table_abscissas[table_data_offset]),
			 &(table_ordinates[table_data_offset]) );
      table_data_offset += table_num_points[i];
    }

    // create enforcement
    if( enforcement_type ){
      switch( enforcement_type ){
      case(1):{
	// Transient Dynamic Enforcement
	tdenf = new ContactTDEnforcement( enforcement_data, search, error );
	enf_object = tdenf;
	break;
      }
      case(2):{
	gr = new ContactGapRemoval( enforcement_data, search, error );
	enf_object = gr;
	break;
      }
      case(3):{
        tiedkin = new ContactTiedKinematics( enforcement_data, search, error );
        enf_object = tiedkin;
        break;
      }
      case(4):{
        voltrans = new ContactVolumeTransfer( enforcement_data, search, error );
        enf_object = voltrans;
        break;
      }
      case(5):{
	mpcenf = new ContactMPCs( enforcement_data, search, error );
	enf_object = mpcenf;
	break;
      }
#ifdef CONTACT_TD_FACE_FACE_ENF
      case(6):{
	tdffenf = new ContactTDFaceFaceEnf( enforcement_data, search, error );
	enf_object = tdffenf;
	break;
      }
#endif
      }
    }
    Real real_data[MAX_ENF_MODEL_DATA_VARS];
    int int_data[MAX_ENF_MODEL_DATA_VARS];
    for( i=0 ; i<num_enf_models ; i++ ){
      switch( enf_model_types[i] ){
      case ContactEnforcement::TD_SPOT_WELD:{  // 4
	real_data[0] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	real_data[1] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+1];
	int_data[0]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+2];
	int_data[1]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+3];
	enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
	break;
      }
      case ContactTDEnforcement::TD_SPRING_WELD:{ // 7
	int_data[0]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	int_data[1]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+1];
	real_data[0] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+2];
	int_data[2]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+3];
	int_data[3]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+4];
	enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
	break;
      }
      case ContactTDEnforcement::TD_ADHESION:{ // 8
	int_data[0]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
	break;
      }
      case ContactTDEnforcement::TD_COHESIVE_ZONE:{ // 9
	real_data[0] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	real_data[1] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+1];
	int_data[0]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+2];
	enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
	break;
      }
      case ContactTDEnforcement::TD_JUNCTION:{ // 10
	int_data[0]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	int_data[1]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+1];
	real_data[0] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+2];
	enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
	break;
      }
      case ContactTDEnforcement::TD_THREADED:{ // 11
	int_data[0]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	int_data[1]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+1];
	int_data[2]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+2];
	int_data[3]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+3];
	int_data[4]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+4];
	real_data[0]  = enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+5];
	real_data[1]  = enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+6];
	real_data[2]  = enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+7];
	enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
	break;
      }
      case ContactEnforcement::TD_AREA_WELD:{  // 13
	real_data[0] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	real_data[1] =       enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+1];
	int_data[0]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+2];
	int_data[1]  = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+3];
	enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
	break;
      }
      case ContactTDEnforcement::TD_USER:{ // 14
        int ii;
        int n=4;
	int_data[0] = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0];
	int_data[1] = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+1];
	int_data[2] = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+2];
	int_data[3] = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+3];
        for (ii=0; ii<int_data[0]; ii++, n++) {
	  int_data[ii+4] = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+n];
	}
        for (ii=0; ii<int_data[1]; ii++, n++) {
	  real_data[ii] = (int) enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+n];
	}
        enf_object->Add_Enforcement_Model( 
	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], &(int_data[0]), &(real_data[0]) );
        for (ii=0; ii<enf_model_sub_num[i]; ii++) {
          int j = i*MAX_ENF_MODEL_SUBS+ii;
          int jj = j*64;
          int size = strlen(&enf_model_sub_names[jj]);
          switch (enf_model_sub_ids[j]) {
          case 0: {
            CONTACT_INIT_MODEL_FN * fn;
            FORTRAN(get_usersub)( (void**)&fn ,
                                  &enf_model_sub_names[jj],
                                  size );
            enf_object->User_Initialize_Model_Fn(enf_model_ids[i],fn);
            } break;
          case 1: {
            CONTACT_INIT_TIME_STEP_FN * fn;
            FORTRAN(get_usersub)( (void**)&fn ,
                                  &enf_model_sub_names[jj],
                                  size );
            enf_object->User_Initialize_Time_Step_Fn(enf_model_ids[i],fn);
            } break;
          case 2: {
            CONTACT_INIT_NODE_STATE_DATA_FN * fn;
            FORTRAN(get_usersub)( (void**)&fn ,
                                  &enf_model_sub_names[jj],
                                  size );
            enf_object->User_Initialize_Node_State_Data_Fn(enf_model_ids[i],fn);
            } break;
          case 3: {
            CONTACT_INTERACTION_TYPE_FN * fn;
            FORTRAN(get_usersub)( (void**)&fn ,
                                  &enf_model_sub_names[jj],
                                  size );
            enf_object->User_Interaction_Type_Fn(enf_model_ids[i],fn);
            } break;
          case 4: {
            CONTACT_INTERACTION_ACTIVE_FN * fn;
            FORTRAN(get_usersub)( (void**)&fn ,
                                  &enf_model_sub_names[jj],
                                  size );
            enf_object->User_Active_Fn(enf_model_ids[i],fn);
            } break;
          case 5: {
            CONTACT_LIMIT_FORCE_FN * fn;
            FORTRAN(get_usersub)( (void**)&fn ,
                                  &enf_model_sub_names[jj],
                                  size );
            enf_object->User_Limit_Force_Fn(enf_model_ids[i],fn);
            } break;
          }
        }
	break;
      }
      default:
	enf_object->Add_Enforcement_Model( 
 	  (ContactEnforcement::Enforcement_Model_Types) enf_model_types[i],
	  enf_model_ids[i], NULL, 
	  &(enf_model_data[MAX_ENF_MODEL_DATA_VARS*i+0]) );
	break;
      }
    }
  } else {
    // original binary restart. new (SIERRA-friendly) restart is below
    int rsize;
    FORTRAN(read_restart_buffer_size)( rsize,icomm,my_proc_id,num_procs );
    Real* rbuf = new Real[rsize];
    FORTRAN(read_restart_data)( rbuf,icomm,my_proc_id,num_procs );
    search = new ContactSearch( rbuf, &(node_ids[0]), &(face_ids[0]),
                                &(element_ids[0]), comm, error );
    delete [] rbuf;
    if( enforcement_type ){
      FORTRAN(read_enf_restart_buffer_size)( rsize,icomm,my_proc_id,num_procs );
      rbuf = new Real[rsize];
      FORTRAN(read_enf_restart_data)( rbuf,icomm,my_proc_id,num_procs );
      switch( enforcement_type ){
      case(1):{
        // Transient Dynamic Enforcement
        tdenf = new ContactTDEnforcement( search, rbuf, error );
        break;
      }
      case(2):{
        gr = new ContactGapRemoval( search, rbuf, error );
        break;
      }
      case(3):{
        tiedkin = new ContactTiedKinematics( search, rbuf, error );
        break;
      }
      case(4):{
        voltrans = new ContactVolumeTransfer( search, rbuf, error );
        break;
      }
      case(5):{
	mpcenf = new ContactMPCs( search, rbuf, error );
	break;
      }
#ifdef CONTACT_TD_FACE_FACE_ENF
      case(6):{
	tdffenf = new ContactTDFaceFaceEnf( search, rbuf, error );
	break;
      }
#endif
      }
      delete [] rbuf;
    }
  }
  
  // setup debug nodes
  for ( i = 0; i < num_debug_nodes; i++ ) 
    search->Add_Debug_Node(debug_nodes[i]);


  // now run the search and enforcement. If the number of load steps from
  // the input file is greater than 1, make sure the input file has enough
  // data for that number of load steps.
  int num_steps = 0;
  int num_pos_vars = num_nodal_vars;
  if (enforcement_type == 1 ) num_pos_vars -= num_enf_type_1_vars;
  if (num_pos_vars%idim != 0) {
    cerr << "FATAL ERROR: number of nodal vars in mesh must be a" << endl
	 << "  multiple of the spatial dimension, plus enforcement" << endl
	 << "  variables if enforcement has been requested." << endl;
    exit(1);
  }
  if (search_type == 0) num_steps = num_pos_vars/idim+1;
  else num_steps = num_pos_vars/idim;
  if (num_steps < num_load_steps) {
    cerr << "FATAL ERROR:  number of load steps specified in input" << endl
	 << "  file is greater than number of load steps defined in" << endl
	 << "  mesh file." << endl;
    exit(1);
  }

  // if a load step to start on was provided, make sure it is less than
  // the maximum number of steps
  if ( start_step > 1 ) {
    if ( start_step > num_load_steps ) {
      cerr << "FATAL ERROR: start step of " << start_step << " is larger "
	   << endl << "  than the number of steps specified, " 
	   << num_steps << endl;
      exit(1);
    }
  }

  // copy enforcement data
  if (enforcement_type == 1) {
    int offset = tot_num_nodes*num_pos_vars;
    for (i = 0; i < tot_num_nodes*num_enf_type_1_vars; i++)
      enforcement_vars[i] = nodal_vars[offset + i];
  }

  // loop over the number of load steps we plan to run
  for (int current_step = start_step;
       current_step <= num_load_steps;
       current_step ++){
    
    
    if( my_proc_id == 0 )
      cout << ">> Executing step " << current_step << endl;
    
    if (current_step>1) {
      // if a new decomposition has been computed, then read the
      // new mesh data
      if (topology_update == 2 || topology_update == 3 ) {
#ifdef DOTHIS
        int step = current_step-1;
        FORTRAN(new_dlb) (icomm,
                          idim,
                          number_node_blocks,
			  node_block_types,
                          nodes_per_block,
                          node_eids,
                          node_ids,
                          nodal_vars,
                          num_nodal_vars,
                          number_face_blocks,
                          face_block_types,
			  shell_block_offsets,
                          faces_per_block,
                          face_ids,
                          face_conn,
			  face_ignore_list,
			  num_face_ignore,
			  num_load_steps,
                          number_element_blocks,
                          element_block_types,
                          elements_per_block,
                          element_ids,
                          element_conn,
			  element_vars,
			  num_element_vars,
                          initial_position,
			  attributes,
			  num_attributes,
			  shell_layer_offsets,
                          num_comm_partners,
                          comm_proc_ids,
                          number_nodes_to_partner,
                          comm_node,
                          search_type,
                          enforcement_type,
                          restart_flag,
                          step,
                          my_proc_id,
                          num_procs);
#endif
      }

      int node_ignore_offset = 0;
      int face_ignore_offset = 0;
      int elem_ignore_offset = 0;
      for (i = 1; i < current_step; i++ ) {
	node_ignore_offset += num_node_ignore[i-1];
	face_ignore_offset += num_face_ignore[i-1];
	elem_ignore_offset += num_elem_ignore[i-1];
      }
      
      // modify the topology for any ignored faces
      if (topology_update) {
        FORTRAN(get_topology_birthdeath)
          ( num_node_ignore[current_step-1],
            &node_ignore_list[node_ignore_offset],
            num_face_ignore[current_step-1],
            &face_ignore_list[face_ignore_offset],
            num_elem_ignore[current_step-1],
            &elem_ignore_list[elem_ignore_offset],
            num_node_deaths_per_block, 
	    node_deaths_global_ids,
            num_face_deaths_per_block, 
	    face_deaths_global_ids,
            num_element_deaths_per_block, 
	    element_deaths_global_ids,
            num_node_births_per_block, 
	    node_births_exodus_ids,
	    node_births_global_ids,
	    number_face_births_per_block, 
	    face_births_global_ids,
	    face_births_connectivity,
	    number_element_births_per_block, 
	    element_births_global_ids,
	    element_births_connectivity,
            my_proc_id,
            num_procs );
        FORTRAN(get_topology_dlb)
          ( num_node_exports, 
	    &node_export_gids[0],
	    &node_export_pids[0],
	    num_face_exports, 
	    &face_export_gids[0],
	    &face_export_pids[0],
	    num_elem_exports, 
	    &elem_export_gids[0],
            &elem_export_pids[0]);
        FORTRAN(get_topology_comm_plan)
          ( num_comm_partners,
            comm_proc_ids,
            number_nodes_to_partner,
            comm_node,
            my_proc_id,
            num_procs );
        FORTRAN(get_topology_host_ids)
          ( tot_num_nodes,
            node_eids,
            node_ids,
            tot_num_faces,
            face_ids,
            tot_num_elems,
            element_ids );
        FORTRAN(get_topology_variables)
          ( initial_position,
            nodal_vars,
            num_nodal_vars,
	    element_vars,
	    num_element_vars,
	    attributes,
	    num_attributes,
	    node_radius,
	    shell_thickness,
            shell_offset );
      }
    }

    // set up the current position
    if ( current_step > 1 && use_config ) {
      // if this is not the first step, compute the current config
      if( use_config == 1 ){
	// Compute using the given displacements for the previous step
	int var_index = (current_step-2)*idim*tot_num_nodes;
	for (i = 0; i < tot_num_nodes; i++){
	  position1[idim*i] = 
	    initial_position[idim*i] + nodal_vars[var_index + i];
	  position1[idim*i+1] = 
	    initial_position[idim*i+1] + 
	    nodal_vars[var_index+tot_num_nodes+i];
	  if (idim == 3 ) {
	    position1[idim*i+2] = 
	      initial_position[idim*i+2] + 
	      nodal_vars[var_index+2*tot_num_nodes+i];
          }
        }
      } else if( use_config == 2 ){
	// Simply read it from the database
	int var_index = (2*(current_step-1)-1)*idim*tot_num_nodes;
	for( i=0 ; i<tot_num_nodes ; i++ ){
	  position1[idim*i  ] = nodal_vars[var_index + i];
	  position1[idim*i+1] = nodal_vars[var_index+tot_num_nodes+i];
	  if (idim == 3 ) {
	    position1[idim*i+2] = nodal_vars[var_index+2*tot_num_nodes+i];
	  }
	}
      }
    } else {
      // this is first step. set first position as initial position
      memcpy (position1, initial_position, tot_num_nodes*idim*sizeof(Real));
    }
      
    if (current_step>1) {
      // now call the routine to update the 
      // toplogy if the topology is different
      if (topology_update) {
	ContactSearch::ContactErrorCode error;
	search->UpdateSearch(num_node_deaths_per_block, 
		             node_deaths_global_ids,
                             num_face_deaths_per_block, 
	                     face_deaths_global_ids,
                             num_element_deaths_per_block, 
	                     element_deaths_global_ids,
                             
                             num_node_births_per_block, 
	                     node_births_exodus_ids,
	                     node_births_global_ids,
	                     number_face_births_per_block, 
	                     face_births_global_ids,
	                     face_births_connectivity,
	                     number_element_births_per_block, 
	                     element_births_global_ids,
	                     element_births_connectivity,
                             
	                     num_node_exports,
                             node_export_gids,
                             node_export_pids,
	                     num_face_exports,
                             face_export_gids,
                             face_export_pids,
	                     num_elem_exports,
                             elem_export_gids,
                             elem_export_pids,
                             
                             tot_num_nodes,
                             node_ids,
                             tot_num_faces,
                             face_ids,
                             tot_num_elems,
                             element_ids,
                             
	                     num_comm_partners,
	                     comm_proc_ids,
	                     number_nodes_to_partner,
	                     comm_node,
                             
	                     error);
      }
    }
    
    
    // setup nodal configuration data
    for (j=0, i=0; i<number_node_blocks; j+=nodes_per_block[i]*idim, i++) {
      search->Set_Node_Block_Configuration(ContactSearch::CURRENT_CONFIG,
					   node_block_ids[i], &position1[j]);
    }

    // if we have shell blocks, set the thickness and lofting
    int index = 0;
    for( i=0 ; i<number_face_blocks ; i++ ){
      if( face_block_types[i] == 5 || face_block_types[i] == 6 ){
	search->Set_Face_Block_Attributes( ContactSearch::SHELL_THICKNESS,
					   i+1, &shell_thickness[index] );
	search->Set_Face_Block_Attributes( ContactSearch::LOFTING_FACTOR,
			                   i+1, &shell_offset[index] );
        index += faces_per_block[i];
      }
    }

    // if we have sph blocks, set the radius
    index = 0;
    for( i=0 ; i<number_node_blocks ; i++ ){
      if( node_block_types[i] == 2 ){
	search->Set_Node_Block_Attributes( ContactSearch::RADIUS,
					   i+1, &node_radius[index] );
        index += nodes_per_block[i];
      }
    }

    // if this is restart from new (SIERRA-friendly) restart, 
    // then read and install the new restart data
    if (restart_flag == 4 || restart_flag == 6){
      // Note: We currently do NOT store edge restart variables, which
      // is currently o.k. since there aren't any. When we do, we will
      // need a method of querying the search object to how many edges
      // there are, since the host code doen't know.
      int rsize;
      FORTRAN(read_restart_buffer_size)( rsize,icomm,my_proc_id,num_procs );
      Real* restart_buffer = new Real[rsize];
      memset (restart_buffer, 0, rsize*sizeof(Real));
      for (i = 0; i < rsize; i++) restart_buffer[i]=0;
      FORTRAN(read_restart_data)( restart_buffer,icomm,
				  my_proc_id,num_procs );
      int num_gen_search_vars      = (int)restart_buffer[0];
      int num_nodal_search_vars    = (int)restart_buffer[1];
      int num_edge_search_vars     = (int)restart_buffer[2];
      int num_face_search_vars     = (int)restart_buffer[3];
      int num_element_search_vars  = (int)restart_buffer[4];
      int num_gen_enforce_vars     = (int)restart_buffer[5];
      int num_nodal_enforce_vars   = (int)restart_buffer[6];
      int num_edge_enforce_vars    = (int)restart_buffer[7];
      int num_face_enforce_vars    = (int)restart_buffer[8];
      int num_element_enforce_vars = (int)restart_buffer[9];
      Real* rbuffer                = &restart_buffer[10];
      
      
      search->Implant_General_Restart_Variable( rbuffer);
      rbuffer += num_gen_search_vars;
      for (i = 1; i <= num_nodal_search_vars; i++)
	search->Implant_Nodal_Restart_Variable( i, &rbuffer[(i-1)*tot_num_nodes]);
      rbuffer += num_nodal_search_vars*tot_num_nodes;
      // here we process edge variables -- this will fail if 
      // num_edge_search_vars != 0; we need a way to get what the
      // edges are (the host code doesn't know)
      for (i = 1; i <= num_edge_search_vars; i++)
	search->Implant_Edge_Restart_Variable( i, NULL);
      for (i = 1; i <= num_face_search_vars; i++)
	search->Implant_Face_Restart_Variable( i, &rbuffer[(i-1)*tot_num_faces]);
      rbuffer += num_face_search_vars*tot_num_faces;
      for (i = 1; i <= num_element_search_vars; i++)
	search->Implant_Element_Restart_Variable( i, &rbuffer[(i-1)*tot_num_elems]);
      rbuffer += num_element_search_vars*tot_num_elems;
      search->Complete_Restart();
      if( enforcement_type ){
	switch( enforcement_type ){
	case(1):{
	  tdenf->Implant_General_Restart_Variable( rbuffer);
	  rbuffer += num_gen_enforce_vars;
	  for (i = 1; i <= num_nodal_enforce_vars; i++)
	    tdenf->Implant_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes]);
          rbuffer += num_nodal_enforce_vars*tot_num_nodes;
	  // here we process edge variables -- this will fail if 
	  // num_edge_enforce_vars != 0; we need a way to get what the
	  // edges are (the host code doesn't know)
	  for (i = 1; i <= num_edge_enforce_vars; i++)
	    tdenf->Implant_Edge_Restart_Variable( i, NULL);
	  for (i = 1; i <= num_face_enforce_vars; i++)
	    tdenf->Implant_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces]);
          rbuffer += num_face_enforce_vars*tot_num_faces;
	  for (i = 1; i <= num_element_enforce_vars; i++)
	    tdenf->Implant_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems]);
          rbuffer += num_element_enforce_vars*tot_num_elems;
	  break;
	}
	case(2):{
	  gr->Implant_General_Restart_Variable( rbuffer);
	  rbuffer += num_gen_enforce_vars;
	  for (i = 1; i <= num_nodal_enforce_vars; i++)
	    gr->Implant_Nodal_Restart_Variable(i, &restart_buffer[(i-1)*tot_num_nodes]);
          rbuffer += num_nodal_enforce_vars*tot_num_nodes;
	  // here we process edge variables -- this will fail if 
	  // num_edge_enforce_vars != 0; we need a way to get what the
	  // edges are (the host code doesn't know)
	  for (i = 1; i <= num_edge_enforce_vars; i++)
	    gr->Implant_Edge_Restart_Variable( i, NULL);
	  for (i = 1; i <= num_face_enforce_vars; i++)
	    gr->Implant_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces]);
          rbuffer += num_face_enforce_vars*tot_num_faces;
	  for (i = 1; i <= num_element_enforce_vars; i++)
	    gr->Implant_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems]);
          rbuffer += num_element_enforce_vars*tot_num_elems;
	  break;
	}
        case(3):{
	  tiedkin->Implant_General_Restart_Variable(rbuffer);
	  rbuffer += num_gen_enforce_vars;
          for (i = 1; i <= num_nodal_enforce_vars; i++)
            tiedkin->Implant_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes]);
          rbuffer += num_nodal_enforce_vars*tot_num_nodes;
          // here we process edge variables -- this will fail if
          // num_edge_enforce_vars != 0; we need a way to get what the
          // edges are (the host code doesn't know)
          for (i = 1; i <= num_edge_enforce_vars; i++)
            tiedkin->Implant_Edge_Restart_Variable( i, NULL);
          for (i = 1; i <= num_face_enforce_vars; i++)
            tiedkin->Implant_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces]);
          rbuffer += num_face_enforce_vars*tot_num_faces;
	  for (i = 1; i <= num_element_enforce_vars; i++)
	    tiedkin->Implant_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems]);
          rbuffer += num_element_enforce_vars*tot_num_elems;
          break;
        }
        case(5):{
	  mpcenf->Implant_General_Restart_Variable(rbuffer);
	  rbuffer += num_gen_enforce_vars;
          for (i = 1; i <= num_nodal_enforce_vars; i++)
            mpcenf->Implant_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes]);
          rbuffer += num_nodal_enforce_vars*tot_num_nodes;
          // here we process edge variables -- this will fail if
          // num_edge_enforce_vars != 0; we need a way to get what the
          // edges are (the host code doesn't know)
          for (i = 1; i <= num_edge_enforce_vars; i++)
            mpcenf->Implant_Edge_Restart_Variable( i, NULL);
          for (i = 1; i <= num_face_enforce_vars; i++)
            mpcenf->Implant_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces]);
          rbuffer += num_face_enforce_vars*tot_num_faces;
	  for (i = 1; i <= num_element_enforce_vars; i++)
	    mpcenf->Implant_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems]);
          rbuffer += num_element_enforce_vars*tot_num_elems;
          break;
        }
#ifdef CONTACT_TD_FACE_FACE_ENF
        case(6):{
	  tdffenf->Implant_General_Restart_Variable(rbuffer);
	  rbuffer += num_gen_enforce_vars;
          for (i = 1; i <= num_nodal_enforce_vars; i++)
            tdffenf->Implant_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes]);
          rbuffer += num_nodal_enforce_vars*tot_num_nodes;
          // here we process edge variables -- this will fail if
          // num_edge_enforce_vars != 0; we need a way to get what the
          // edges are (the host code doesn't know)
          for (i = 1; i <= num_edge_enforce_vars; i++)
            tdffenf->Implant_Edge_Restart_Variable( i, NULL);
          for (i = 1; i <= num_face_enforce_vars; i++)
            tdffenf->Implant_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces]);
          rbuffer += num_face_enforce_vars*tot_num_faces;
	  for (i = 1; i <= num_element_enforce_vars; i++)
	    tdffenf->Implant_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems]);
          rbuffer += num_element_enforce_vars*tot_num_elems;
          break;
        }
#endif
	}
      }
      delete [] restart_buffer;
    }

    // Perform the search
    if( search_type==0){
      error = search->Static_Search_1_Configuration();
      if( error ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	cerr << "Error in Static_Search: Error Code = "
	     << error << endl;
	for( i=0 ; i<search->Number_of_Errors() ; i++ )
	  cerr << search->Error_Message(i+1) << endl;
#endif
	delete search;
	exit(1);
      }
    } else {
      if( use_config != 2 ){
	for (i = 0; i < tot_num_nodes; i++){
	  int var_index = (current_step-1)*idim*tot_num_nodes;
	  position2[idim*i] = 
	    initial_position[idim*i] + nodal_vars[var_index + i];
	  position2[idim*i+1] = 
	    initial_position[idim*i+1] + 
	    nodal_vars[var_index+tot_num_nodes+i];
	  if (idim == 3 ) 
	    position2[idim*i+2] = 
	      initial_position[idim*i+2] + 
	      nodal_vars[var_index+2*tot_num_nodes+i];
	}
      } else {
	for (i = 0; i < tot_num_nodes; i++){
	  int var_index = (current_step-1)*2*idim*tot_num_nodes;
	  position2[idim*i  ] = position1[idim*i] + nodal_vars[var_index + i];
	  position2[idim*i+1] = position1[idim*i+1] + 
	    nodal_vars[var_index+tot_num_nodes+i];
	  if (idim == 3 ) {
	    position2[idim*i+2] = position1[idim*i+2] + 
	      nodal_vars[var_index+2*tot_num_nodes+i];
	  }
        }	
      }

      for (j=0, i=0; i<number_node_blocks; j+=nodes_per_block[i]*idim, i++) {
	error = 
	  search->Set_Node_Block_Configuration(ContactSearch::PREDICTED_CONFIG,
					       node_block_ids[i], 
					       &position2[j]);
	if( error ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	  cerr << "Error in Set_NodeBlk_Configuration: Error code = "
	       << error << endl;
	  for( i=0 ; i<search->Number_of_Errors() ; i++ )
	    cerr << search->Error_Message(i+1) << endl;
#endif
	  delete search;
	  exit(1);
	}
      }

      // Now execute the search
      if( search_type==1 ){
	error = search->Static_Search_2_Configuration();
	if( error ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	  cerr << "Error in Static_Search_2: Error Code = "
	       << error << endl;
	  for( i=0 ; i<search->Number_of_Errors() ; i++ )
	    cerr << search->Error_Message(i+1) << endl;
#endif
	  delete search;
	  exit(1);
	}
      } else if( search_type==2 ){
	error = search->Dynamic_Search_2_Configuration();
	if( error ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	  cerr << "Error in Dynamic_Search : Error Code = "
	       << error << endl;
	  for( i=0 ; i<search->Number_of_Errors() ; i++ )
	    cerr << search->Error_Message(i+1) << endl;
#endif
	  delete search;
	  exit(1);
	}
      } else if( search_type==3 ){
	Real* mass = new Real[tot_num_nodes];
	for( i=0; i<tot_num_nodes; i++ ) {
	  mass[i] = enforcement_vars[i];
	}
	error = search->Dynamic_Search_Augmented_2_Configuration(mass,
								 enforcement_data_vars[0],
								 enforcement_data_vars[1] );
	delete [] mass;
	if( error ){
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	  cerr << "Error in Dynamic_Search_Augmented : Error Code = "
	       << error << endl;
	  for( i=0 ; i<search->Number_of_Errors() ; i++ )
	    cerr << search->Error_Message(i+1) << endl;
#endif
	  delete search;
	  exit(1);
	}
      }
    }
    
    search->Display();
    
    // Get the node-node interactions
    int num_interactions,data_size;
    search->Size_NodeNode_Interactions(num_interactions,data_size);
    int* slave_node_block_id = NULL;
    int* slave_node_indexes_in_block = NULL;
    int* master_node_block_id = NULL;
    int* master_node_indexes_in_block = NULL;
    int* master_node_proc = NULL;
    Real* interactions_data;
    if( num_interactions ){
      slave_node_block_id         = new int[num_interactions];
      slave_node_indexes_in_block = new int[num_interactions];
      master_node_block_id	 = new int[num_interactions];
      master_node_indexes_in_block = new int[num_interactions];
      master_node_proc = new int[num_interactions];
      interactions_data = new Real[num_interactions*data_size];
      search->Get_NodeNode_Interactions(slave_node_block_id,
                                        slave_node_indexes_in_block,
					master_node_block_id,
					master_node_indexes_in_block,
					master_node_proc,
					interactions_data);
    }
    // This must be called in sync
    FORTRAN(output_nodenode_interactions)( num_interactions,
                                           slave_node_block_id,
					   slave_node_indexes_in_block,
					   master_node_block_id,
					   master_node_indexes_in_block,
					   master_node_proc,
					   interactions_data, 
                                           data_size );
    if( num_interactions ){
      delete [] slave_node_block_id;
      delete [] slave_node_indexes_in_block;
      delete [] master_node_block_id;
      delete [] master_node_indexes_in_block;
      delete [] master_node_proc;
      delete [] interactions_data;
    }
    
    // Get the node-face interactions
    search->Size_NodeFace_Interactions(num_interactions,data_size);
    int* node_block_id = NULL;
    int* node_indexes_in_block = NULL;
    int* node_entity_keys = NULL;
    int* face_block_id = NULL;
    int* face_indexes_in_block = NULL;
    int* face_proc = NULL;
    if( num_interactions ){
      node_block_id         = new int[num_interactions];
      node_indexes_in_block = new int[num_interactions];
      node_entity_keys      = new int[num_interactions];
      face_block_id         = new int[num_interactions];
      face_indexes_in_block = new int[num_interactions];
      face_proc = new int[num_interactions];
      interactions_data = new Real[num_interactions*data_size];
      search->Get_NodeFace_Interactions(node_block_id,node_indexes_in_block,
					node_entity_keys,
					face_block_id,face_indexes_in_block,
					face_proc,interactions_data);
    }
    // This must be called in sync
    FORTRAN(output_nodeface_interactions)( num_interactions,node_block_id,
					   node_indexes_in_block,
					   face_block_id,face_indexes_in_block,
					   face_proc,interactions_data, 
                                           data_size );
    if( num_interactions ){
      delete [] node_block_id;
      delete [] node_indexes_in_block;
      delete [] node_entity_keys;
      delete [] face_block_id;
      delete [] face_indexes_in_block;
      delete [] face_proc;
      delete [] interactions_data;
    }
    
    // Get the node-surface interactions
    search->Size_NodeSurface_Interactions( num_interactions,data_size );
    int* n_blk_id = NULL;
    int* n_ind = NULL;
    int* surf_id = NULL;
    Real* data = NULL;
    if( num_interactions ){
      n_blk_id = new int[num_interactions];
      n_ind    = new int[num_interactions];
      surf_id  = new int[num_interactions];
      data     = new Real[num_interactions*data_size];
      search->Get_NodeSurface_Interactions( n_blk_id, n_ind, surf_id, data );
    }
    // This must be called in sync
    FORTRAN(output_nodesurface_interactions)(num_interactions,n_blk_id, n_ind,
					     surf_id, data, &data_size );
    if( num_interactions ){
      delete [] n_blk_id;
      delete [] n_ind;
      delete [] surf_id;
      delete [] data;
    }
    
    // Enforcement
    if( enforcement_type ){
      switch( enforcement_type ){
      case(1):{
	Real* force     = new Real[3*tot_num_nodes];
	Real* mass      = new Real[tot_num_nodes];
	Real* density   = new Real[tot_num_nodes];
	Real* wavespeed = new Real[tot_num_nodes];
	for( i=0; i<tot_num_nodes; i++ ) {
	  mass[i]      = enforcement_vars[i];
	  density[i]   = enforcement_vars[tot_num_nodes + i];
	  wavespeed[i] = enforcement_vars[tot_num_nodes*2 + i];
	}
	tdenf->Set_Number_of_Iterations( (int) enforcement_data_vars[2] );
        if( (int) enforcement_data_vars[3] == 1 )
          tdenf->Enforce_Symmetry_on_Nodes( (int) enforcement_data_vars[4],
                                            (int) enforcement_data_vars[5] );
	if( tdenf->Compute_Contact_Force( enforcement_data_vars[0],
					  enforcement_data_vars[1],
					  mass, density,
					  wavespeed, force ) ){
	  int ne = tdenf->Number_of_Errors();
	  for( i=0 ; i<ne ; i++ ){
	    cerr << tdenf->Error_Message(i+1) << endl;
	  }
	}
	delete [] force;
	delete [] mass;
	delete [] density;
	delete [] wavespeed;
	break;
      }
      case(2):{
	Real* dispc = new Real[3*tot_num_nodes];
	int max_iterations = (int) enforcement_data_vars[0];
	double trivial_gap = enforcement_data_vars[1];
	if( gr->Compute_Gap_Removal( max_iterations, trivial_gap, dispc ) ){
	  int ne = gr->Number_of_Errors();
	  for( i=0 ; i<ne ; i++ ){
	    cerr << tdenf->Error_Message(i+1) << endl;
	  }
	}
	delete [] dispc;
	break;
      }
      case(3):{
        Real* pos = new Real[3*tot_num_nodes];
        if( tiedkin->Compute_Position( pos ) ){
          int ne = tiedkin->Number_of_Errors();
          for( i=0 ; i<ne ; i++ ){
            cerr << tiedkin->Error_Message(i+1) << endl;
	  }
	}
	delete [] pos;
	break;
      }
      case(4):{
	// transfer all element and nodal variables from genesis input 
	int num_node_vars   = num_nodal_vars;
	int num_elem_vars   = num_element_vars;  
	Real* donor_node_vars = new Real[num_node_vars*tot_num_nodes];
	Real* donor_elem_vars = new Real[num_elem_vars*tot_num_elems];
	Real* receiver_node_vars = new Real[num_node_vars*tot_num_nodes];
	Real* receiver_elem_vars = new Real[num_elem_vars*tot_num_elems];
	Real* volume_fraction = new Real[tot_num_elems];

	// set node variables for transfer
	int node_count = 0;
	for( i=0; i<tot_num_nodes; i++ ) {
	  for( j=0; j<num_node_vars; j++ ) {
	    donor_node_vars[node_count++] = nodal_vars[j*tot_num_nodes+i];
	  }
	}
	
	// set element variables for transfer
	int elem_count = 0;
	for( i=0; i<tot_num_elems; i++ ) {
	  for( j=0; j<num_elem_vars; j++ ) {
	    donor_elem_vars[elem_count++] = element_vars[j*tot_num_elems+i];  
	  }
	}
     
	ContactSearch::ContactErrorCode ec = 
	  voltrans->Compute_Volume_Transfer( num_node_vars,
					     num_elem_vars,
					     donor_node_vars,
					     donor_elem_vars,
					     receiver_node_vars,
					     receiver_elem_vars,
					     volume_fraction );

	delete [] donor_node_vars;
	delete [] donor_elem_vars;
	delete [] receiver_node_vars;
	delete [] receiver_elem_vars;
	delete [] volume_fraction;
	break;
      }
      case(5):{
	int num_mpcs = 0;
	ContactMPCs::Id_Numbering_Scheme mpc_id_num_scheme = ContactMPCs::ACME_LOCAL_ID;
	mpcenf->Compute_MPCs( (ContactMPCs::Id_Numbering_Scheme) mpc_id_num_scheme );
	ContactSearch::ContactErrorCode ec = 
	  mpcenf->Number_of_MPC_Equations(num_mpcs);
	int* snode_pid   = new int[num_mpcs];
	int* snode_lid   = new int[num_mpcs];
	int* mface_pid   = new int[num_mpcs];
	int* mface_lid   = new int[num_mpcs];
	int* nface_nodes = new int[num_mpcs];
	int* fnode_pid   = new int[8*num_mpcs];
	int* fnode_lid   = new int[8*num_mpcs];
	Real* fnode_coefs = new Real[8*num_mpcs];
	mpcenf->Get_MPC_Equations(num_mpcs,
				  snode_pid,
				  snode_lid,
				  mface_pid,
				  mface_lid,
				  nface_nodes,
				  fnode_pid,
				  fnode_lid,
				  fnode_coefs );

	delete [] snode_pid;
	delete [] snode_lid;
	delete [] mface_pid;
	delete [] mface_lid;
	delete [] nface_nodes;
	delete [] fnode_pid;
	delete [] fnode_lid;
	delete [] fnode_coefs;

	break;
      }
#ifdef CONTACT_TD_FACE_FACE_ENF
      case(6):{
	Real* mass = new Real[tot_num_nodes];
	for( i=0; i<tot_num_nodes; i++ ) {
	  mass[i] = enforcement_vars[i];
	}
	Real* force = new Real[3*tot_num_nodes];
	if( tdffenf->Compute_Forces( mass, force ) ){
	  int ne = tdffenf->Number_of_Errors();
	  for( i=0 ; i<ne ; i++ ){
	    cerr << tdffenf->Error_Message(i+1) << endl;
	  }
	}
	delete [] mass;
	delete [] force;
	break;
      }
#endif
      }
    }

#ifndef CONTACT_NO_EXODUS_OUTPUT
    Real time = current_step;
    int plot_step = current_step;
    FORTRAN(create_plot_file)( icomm, my_proc_id, plot_step, exodus_out );
    search->Exodus_Output( exodus_out, time );
    ex_close( exodus_out );
#endif
    
  } // end of loop on nubmer of steps

  // create restart tape if requested
  if( restart_flag == 2 || restart_flag == 3 ){
    int restart_size = search->Restart_Size();
    Real* restart_buffer = new Real[restart_size];
    search->Extract_Restart_Data( restart_buffer );
    FORTRAN(write_restart)(restart_size,restart_buffer,
			   icomm,my_proc_id,num_procs);
    delete [] restart_buffer;
    ContactEnforcement* enf;
    if( enforcement_type ){
      switch( enforcement_type ){
      case(1):
	enf = tdenf;
	break;
      case(2):
	enf = gr;
	break;
      case(3):
	enf = tiedkin;
	break; 
      case(4):
	enf = voltrans;
	break; 
      case(5):
	enf = mpcenf;
#ifdef CONTACT_TD_FACE_FACE_ENF
      case(6):
	enf = tdffenf;
	break;
#endif
      }
      restart_size = enf->Restart_Size();
      restart_buffer = new Real[restart_size];
      enf->Extract_Restart_Data( restart_buffer );
      FORTRAN(write_enf_restart)(restart_size, restart_buffer,
				 icomm, my_proc_id, num_procs );
      delete [] restart_buffer;
    }
  } else if( restart_flag == 5 || restart_flag == 6 ) {
    // Note: We currently do NOT store edge restart variables, which
    // is currently o.k. since there aren't any. When we do, we will
    // need a method of querying the search object to how many edges
    // there are, since the host code doen't know.
    int num_gen_search_vars = search->Number_General_Restart_Variables();
    int num_nodal_search_vars = search->Number_Nodal_Restart_Variables();
    int num_edge_search_vars = search->Number_Edge_Restart_Variables();
    assert ( num_edge_search_vars == 0);
    int num_face_search_vars = search->Number_Face_Restart_Variables();
    int num_element_search_vars = search->Number_Element_Restart_Variables();
    int num_gen_enforce_vars(0);
    int num_nodal_enforce_vars(0);
    int num_edge_enforce_vars(0);
    int num_face_enforce_vars(0);
    int num_element_enforce_vars(0);
    if (enforcement_type) {
      switch( enforcement_type ){
      case(1):{
	num_gen_enforce_vars = tdenf->Number_General_Restart_Variables();
	num_nodal_enforce_vars = tdenf->Number_Nodal_Restart_Variables();
	num_edge_enforce_vars = tdenf->Number_Edge_Restart_Variables();
	num_face_enforce_vars = tdenf->Number_Face_Restart_Variables();
	num_element_enforce_vars = tdenf->Number_Element_Restart_Variables();
	break;
      }
      case(2):{
	num_gen_enforce_vars = gr->Number_General_Restart_Variables();
	num_nodal_enforce_vars = gr->Number_Nodal_Restart_Variables();
	num_edge_enforce_vars = gr->Number_Edge_Restart_Variables();
	num_face_enforce_vars = gr->Number_Face_Restart_Variables();
	num_element_enforce_vars = gr->Number_Element_Restart_Variables();
	break;
      }
      case(3):{
	num_gen_enforce_vars = tiedkin->Number_General_Restart_Variables();
	num_nodal_enforce_vars = tiedkin->Number_Nodal_Restart_Variables();
	num_edge_enforce_vars = tiedkin->Number_Edge_Restart_Variables();
	num_face_enforce_vars = tiedkin->Number_Face_Restart_Variables();
	num_element_enforce_vars = tiedkin->Number_Element_Restart_Variables();
	break;
      }
      case(4):{
	num_gen_enforce_vars = voltrans->Number_General_Restart_Variables();
	num_nodal_enforce_vars = voltrans->Number_Nodal_Restart_Variables();
	num_edge_enforce_vars = voltrans->Number_Edge_Restart_Variables();
	num_face_enforce_vars = voltrans->Number_Face_Restart_Variables();
	num_element_enforce_vars = voltrans->Number_Element_Restart_Variables();
	break;
      }
      case(5):{
	num_gen_enforce_vars = mpcenf->Number_General_Restart_Variables();
	num_nodal_enforce_vars = mpcenf->Number_Nodal_Restart_Variables();
	num_edge_enforce_vars = mpcenf->Number_Edge_Restart_Variables();
	num_face_enforce_vars = mpcenf->Number_Face_Restart_Variables();
	num_element_enforce_vars = mpcenf->Number_Element_Restart_Variables();
	break;
      }
#ifdef CONTACT_TD_FACE_FACE_ENF
      case(6):{
	num_gen_enforce_vars = tdffenf->Number_General_Restart_Variables();
	num_nodal_enforce_vars = tdffenf->Number_Nodal_Restart_Variables();
	num_edge_enforce_vars = tdffenf->Number_Edge_Restart_Variables();
	num_face_enforce_vars = tdffenf->Number_Face_Restart_Variables();
	num_element_enforce_vars = tdffenf->Number_Element_Restart_Variables();
	break;
      }
#endif
      }
      assert ( num_edge_enforce_vars == 0);
    }
    int restart_size = num_gen_search_vars
                     + num_nodal_search_vars * tot_num_nodes
                     + num_face_search_vars * tot_num_faces
                     + num_element_search_vars * tot_num_elems
                     + num_gen_enforce_vars
                     + num_nodal_enforce_vars * tot_num_nodes
                     + num_face_enforce_vars * tot_num_faces
                     + num_element_enforce_vars * tot_num_elems
                     + 10;
    Real* restart_buffer = new Real[restart_size];
    memset (restart_buffer, 0, restart_size*sizeof(Real));
    
    restart_buffer[0] = num_gen_search_vars;
    restart_buffer[1] = num_nodal_search_vars;
    restart_buffer[2] = num_edge_search_vars;
    restart_buffer[3] = num_face_search_vars;
    restart_buffer[4] = num_element_search_vars;
    restart_buffer[5] = num_gen_enforce_vars;
    restart_buffer[6] = num_nodal_enforce_vars;
    restart_buffer[7] = num_edge_enforce_vars;
    restart_buffer[8] = num_face_enforce_vars;
    restart_buffer[9] = num_element_enforce_vars;
    Real* rbuffer     = &restart_buffer[10];

    // general variables    
    search->Extract_General_Restart_Variable(rbuffer);
    rbuffer += num_gen_search_vars;

    // nodal variables
    for (i = 1; i <= num_nodal_search_vars; i++)
      search->Extract_Nodal_Restart_Variable( i, &rbuffer[(i-1)*tot_num_nodes]);
    rbuffer += num_nodal_search_vars * tot_num_nodes;
    
    // here we process edge variables -- this will fail if 
    // num_edge_search_vars != 0; we need a way to get what the
    // edges are (the host code doesn't know)
    for (i = 1; i <= num_edge_search_vars; i++)
      search->Extract_Edge_Restart_Variable( i, NULL);
      
    for (i = 1; i <= num_face_search_vars; i++)
      search->Extract_Face_Restart_Variable( i, &rbuffer[(i-1)*tot_num_faces]);
    rbuffer += num_face_search_vars * tot_num_faces;
                     
    for (i = 1; i <= num_element_search_vars; i++)
      search->Extract_Element_Restart_Variable( i, &rbuffer[(i-1)*tot_num_elems]);
    rbuffer += num_element_search_vars * tot_num_elems;
      
    if( enforcement_type ){
      switch( enforcement_type ){
      case(1):{
	tdenf->Extract_General_Restart_Variable(rbuffer );
	rbuffer += num_gen_enforce_vars;
	for (i = 1; i <= num_nodal_enforce_vars; i++)
	  tdenf->Extract_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes] );
	rbuffer += num_nodal_enforce_vars * tot_num_nodes;
	// here we process edge variables -- this will fail if 
	// num_edge_enforce_vars != 0; we need a way to get what the
	// edges are (the host code doesn't know)
	for (i = 1; i <= num_edge_enforce_vars; i++)
	  tdenf->Extract_Edge_Restart_Variable( i, NULL);
	for (i = 1; i <= num_face_enforce_vars; i++)
	  tdenf->Extract_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces] );
	rbuffer += num_face_enforce_vars * tot_num_faces;
	for (i = 1; i <= num_element_enforce_vars; i++)
	  tdenf->Extract_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems] );
	rbuffer += num_element_enforce_vars * tot_num_elems;
	break;
      }
      case(2):{
	gr->Extract_General_Restart_Variable(rbuffer);
	rbuffer += num_gen_enforce_vars;
	for (i = 1; i <= num_nodal_enforce_vars; i++)
	  gr->Extract_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes] );
	rbuffer += num_nodal_enforce_vars * tot_num_nodes;
	// here we process edge variables -- this will fail if 
	// num_edge_enforce_vars != 0; we need a way to get what the
	// edges are (the host code doesn't know)
	for (i = 1; i <= num_edge_enforce_vars; i++)
	  gr->Extract_Edge_Restart_Variable( i, NULL);
	for (i = 1; i <= num_face_enforce_vars; i++)
	  gr->Extract_Face_Restart_Variable(i, &restart_buffer[(i-1)*tot_num_faces] );
	rbuffer += num_face_enforce_vars * tot_num_faces;
	for (i = 1; i <= num_element_enforce_vars; i++)
	  gr->Extract_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems] );
	rbuffer += num_element_enforce_vars * tot_num_elems;
	break;
      }
      case(3):{
	tiedkin->Extract_General_Restart_Variable(rbuffer);
	rbuffer += num_gen_enforce_vars;
	for (i = 1; i <= num_nodal_enforce_vars; i++)
	  tiedkin->Extract_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes] );
	rbuffer += num_nodal_enforce_vars * tot_num_nodes;
	// here we process edge variables -- this will fail if
	// num_edge_enforce_vars != 0; we need a way to get what the
	// edges are (the host code doesn't know)
	for (i = 1; i <= num_edge_enforce_vars; i++)
	  tiedkin->Extract_Edge_Restart_Variable( i, NULL);
	for (i = 1; i <= num_face_enforce_vars; i++)
	  tiedkin->Extract_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces] );
	rbuffer += num_face_enforce_vars * tot_num_faces;
	for (i = 1; i <= num_element_enforce_vars; i++)
	  tiedkin->Extract_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems] );
	rbuffer += num_element_enforce_vars * tot_num_elems;
	break;
      }
      case(4):{
	voltrans->Extract_General_Restart_Variable(rbuffer);
	rbuffer += num_gen_enforce_vars;
	for (i = 1; i <= num_nodal_enforce_vars; i++)
	  voltrans->Extract_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes] );
	rbuffer += num_nodal_enforce_vars * tot_num_nodes;
	// here we process edge variables -- this will fail if
	// num_edge_enforce_vars != 0; we need a way to get what the
	// edges are (the host code doesn't know)
	for (i = 1; i <= num_edge_enforce_vars; i++)
	  voltrans->Extract_Edge_Restart_Variable( i, NULL);
	for (i = 1; i <= num_face_enforce_vars; i++)
	  voltrans->Extract_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces] );
	rbuffer += num_face_enforce_vars * tot_num_faces;
	for (i = 1; i <= num_element_enforce_vars; i++)
	  voltrans->Extract_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems] );
	rbuffer += num_element_enforce_vars * tot_num_elems;
	break;
      }
      case(5):{
	mpcenf->Extract_General_Restart_Variable(rbuffer);
	rbuffer += num_gen_enforce_vars;
	for (i = 1; i <= num_nodal_enforce_vars; i++)
	  mpcenf->Extract_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes] );
	rbuffer += num_nodal_enforce_vars * tot_num_nodes;
	// here we process edge variables -- this will fail if
	// num_edge_enforce_vars != 0; we need a way to get what the
	// edges are (the host code doesn't know)
	for (i = 1; i <= num_edge_enforce_vars; i++)
	  mpcenf->Extract_Edge_Restart_Variable( i, NULL);
	for (i = 1; i <= num_face_enforce_vars; i++)
	  mpcenf->Extract_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces] );
	rbuffer += num_face_enforce_vars * tot_num_faces;
	for (i = 1; i <= num_element_enforce_vars; i++)
	  mpcenf->Extract_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems] );
	rbuffer += num_element_enforce_vars * tot_num_elems;
	break;
      }
#ifdef CONTACT_TD_FACE_FACE_ENF
      case(6):{
	tdffenf->Extract_General_Restart_Variable(rbuffer);
	rbuffer += num_gen_enforce_vars;
	for (i = 1; i <= num_nodal_enforce_vars; i++)
	  tdffenf->Extract_Nodal_Restart_Variable(i, &rbuffer[(i-1)*tot_num_nodes] );
	rbuffer += num_nodal_enforce_vars * tot_num_nodes;
	// here we process edge variables -- this will fail if
	// num_edge_enforce_vars != 0; we need a way to get what the
	// edges are (the host code doesn't know)
	for (i = 1; i <= num_edge_enforce_vars; i++)
	  tdffenf->Extract_Edge_Restart_Variable( i, NULL);
	for (i = 1; i <= num_face_enforce_vars; i++)
	  tdffenf->Extract_Face_Restart_Variable(i, &rbuffer[(i-1)*tot_num_faces] );
	rbuffer += num_face_enforce_vars * tot_num_faces;
	for (i = 1; i <= num_element_enforce_vars; i++)
	  tdffenf->Extract_Element_Restart_Variable(i, &rbuffer[(i-1)*tot_num_elems] );
	rbuffer += num_element_enforce_vars * tot_num_elems;
	break;
      }
#endif
      }
    }
    FORTRAN(write_restart)(restart_size,restart_buffer,
			   icomm,my_proc_id,num_procs);
    delete [] restart_buffer;
  }
  
  // now delete objects
  if( enforcement_type ){
    switch( enforcement_type ){
    case(1):{
      delete tdenf;
      break;
    }
    case(2):{
      delete gr;
      break;
    }
    case(3):{
      delete tiedkin;
      break;
    }
    case(4):{
      delete voltrans;
      break;
    }
    case(5):{
      delete mpcenf;
      break;
    }
#ifdef CONTACT_TD_FACE_FACE_ENF
    case(6):{
      delete tdffenf;
      break;
    }
#endif
    }

  }

  delete search;

  delete [] enf_model_data;
  delete [] face_conn;
  delete [] element_conn;
  delete [] node_eids;
  delete [] node_ids;
  delete [] face_ids;
  delete [] element_ids;
  delete [] initial_position;
  delete [] position1;
  delete [] position2;
  delete [] nodal_vars;
  delete [] enforcement_vars;
  delete [] num_node_ignore;
  delete [] node_ignore_list;
  delete [] num_face_ignore;
  delete [] face_ignore_list;
  delete [] num_elem_ignore;
  delete [] elem_ignore_list;
  delete [] comm_proc_ids;
  delete [] number_nodes_to_partner;
  delete [] comm_node;
  delete [] element_vars;
  delete [] num_attributes;
  delete [] attributes;
  delete [] node_radius;
  delete [] shell_thickness;
  delete [] shell_offset;
#ifndef CONTACT_NO_MPI
  //contact_global_sync(comm);
  MPI_Finalize();
#endif
}

  
