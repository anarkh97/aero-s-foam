// $Id: read_data.C,v 2002.44 2004/06/22 17:00:00 mwglass Exp $

#include "stdlib.h"
#include "iostream.h"
#include "fstream.h"
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactTDEnforcement.h"
#include "ContactGapRemoval.h"
#include "ContactTiedKinematics.h"
#include "ContactTDFaceFaceEnf.h"
#include "ContactVolumeTransfer.h"
#include "read_data.h"

void FORTRAN(read_data) ( int & num_load_steps,
			  int & use_config,
			  int & search_type,
			  int & restart_flag,
			  int & start_step,
			  int & multiple_interaction_status,
			  Real* interaction_data,
			  int & normal_smoothing_status,
			  Real* smoothing_data,
			  int & compute_node_areas,
			  int & number_of_analytical_surfaces, 
			  int * as_type, 
			  Real* as_data,
			  int & number_face_blocks,
			  int * face_block_types,
			  Real* shell_block_offset,
			  Real* search_data,
			  int&  number_of_tables,
			  int*  table_ids,
			  int*  table_num_points,
			  Real* table_abscissas,
			  Real* table_ordinates,
			  char * filename,
			  int& enforcement_type,
			  Real* enforcement_data,
			  Real* enforcement_data_vars,
			  int& num_enf_models,
			  int* enf_model_types,
			  int* enf_model_ids, 
  	                  int* enf_model_sub_num,
  	                  int* enf_model_sub_ids,
 	                  char* enf_model_sub_names, 
			  Real* enf_model_data,
			  int& update_topology, 
			  int* face_ignore_list, 
			  int* num_face_ignore)
{
  // open file
  ifstream data_file;
  data_file.open(filename);
  if ( ! data_file ) {
    cerr << "ERROR: Cannot open " << filename << endl;
    exit(1);
  }
  char line[1024];
  int i,j,ikeys;

  use_config = 0;
  search_type = 0;
  restart_flag = 0;
  normal_smoothing_status = 0;
  compute_node_areas = 0;
  number_of_tables = 0;
  enforcement_type = 0;
  update_topology = 0;

  // read number of load steps
  data_file.getline(line,sizeof(line));
  data_file >> num_load_steps;
  data_file.get();  // parse the newline special character
  if( num_load_steps > 1 ){
    // Get the flag to control using the initial configuration
    data_file.getline(line,sizeof(line));
    data_file >> use_config;
    data_file.get();
  }

  // read search type
  data_file.getline(line,sizeof(line));
  data_file >> search_type;
  data_file.get();  // parse the newline special character

  // read the restart flag
  data_file.getline(line,sizeof(line));
  data_file >> restart_flag;
  data_file.get(); // parse the newline special character
  if(restart_flag == 1 || restart_flag == 3 ||
     restart_flag == 4 || restart_flag == 6 ) {
    // if restart, get the flag to control using the initial configuration
    data_file.getline(line,sizeof(line));
    data_file >> use_config;
    data_file.get();
  }

  if( restart_flag != 1 && restart_flag != 3 ){
    
    // get the load step at which to restart the analysis if 
    // restart_flag = 4 or 6 (sierra-friendly restart)
    if( restart_flag == 4 || restart_flag == 6 ){
      data_file.getline(line,sizeof(line));
      data_file >> start_step;
      data_file.get();
      if (start_step > 1 ) num_load_steps = num_load_steps + start_step - 1;
    }
    
    //  Read the number of analytic surfaces
    data_file.getline(line,sizeof(line));
    data_file >> number_of_analytical_surfaces;
    data_file.get();
    
    // read in all the analytic surface descriptions
    if( number_of_analytical_surfaces >0 ){
      for (i = 0; i< number_of_analytical_surfaces; i++) {
	data_file.getline(line,sizeof(line));
	data_file >> as_type[i];
	data_file.get();
	
	if ( as_type[i] == 1) {
	  data_file.getline(line,sizeof(line));
	  data_file >> as_data[8*i];
	  data_file >> as_data[8*i+1];
	  data_file >> as_data[8*i+2];
	  data_file >> as_data[8*i+3];
	  data_file >> as_data[8*i+4];
	  data_file >> as_data[8*i+5];
	  data_file.get();
	}
	else if ( as_type[i] == 2) {
	  data_file.getline(line,sizeof(line));
	  data_file >> as_data[8*i];
	  data_file >> as_data[8*i+1];
	  data_file >> as_data[8*i+2];
	  data_file >> as_data[8*i+3];
	  data_file.get();
	}
	else if ( as_type[i] == 3 || as_type[i] == 4 ) {
	  data_file.getline(line,sizeof(line));
	  data_file >> as_data[8*i];
	  data_file >> as_data[8*i+1];
	  data_file >> as_data[8*i+2];
	  data_file >> as_data[8*i+3];
	  data_file >> as_data[8*i+4];
	  data_file >> as_data[8*i+5];
	  data_file >> as_data[8*i+6];
	  data_file >> as_data[8*i+7];
	  data_file.get();
	}
      }
    }
    
    // Read the number of face blocks
    data_file.getline(line,sizeof(line));
    data_file >> number_face_blocks;
    data_file.get();

    // Read the face block types
    if( number_face_blocks ){
      data_file.getline(line,sizeof(line));
      for( i=0 ; i<number_face_blocks ; i++ ){
	data_file >> face_block_types[i];
	if( face_block_types[i] == 5 || face_block_types[i] == 6 ){
	  data_file.get();
	  data_file.getline(line,sizeof(line));
	  data_file >> shell_block_offset[i];
	} else
	  shell_block_offset[i] = -1.0;
      }
      data_file.get();
    }

    // read the number of entity keys
    data_file.getline(line,sizeof(line));
    data_file >> ikeys;
    data_file.get();

    // read search data
    data_file.getline(line,sizeof(line));
    for ( i = 0; i < ikeys*ikeys*ContactSearch::NSIZSD; i++){
      data_file >> search_data[i];
    }
    
    // read the interaction data
    data_file.get();
    data_file.getline(line,sizeof(line));
    data_file >> multiple_interaction_status;
    data_file.get();
    if( multiple_interaction_status == 1 ){
      data_file.getline(line,sizeof(line));
      data_file >> interaction_data[0];
      data_file.get();
    }
    
    // read the normal smoothing data
    data_file.getline(line,sizeof(line));
    data_file >> normal_smoothing_status;
    data_file.get();
    if( normal_smoothing_status == 1){
      data_file.getline(line,sizeof(line));
      data_file >> smoothing_data[0];
      data_file.get();
      data_file.getline(line,sizeof(line));
      data_file >> smoothing_data[1];
      data_file.get();
      data_file.getline(line,sizeof(line));
      data_file >> smoothing_data[2];
      data_file.get();
    }

    // read compute_node_area
    data_file.getline(line,sizeof(line));
    data_file >> compute_node_areas;
    data_file.get();

    // read the table data
    int table_data_offset = 0;
    data_file.getline(line,sizeof(line));
    data_file >> number_of_tables;
    data_file.get();
    for( i=0 ; i<number_of_tables ; i++ ){
      data_file.getline(line,sizeof(line));
      data_file >> table_ids[i];
      data_file.get();
      data_file.getline(line,sizeof(line));
      data_file >> table_num_points[i];
      data_file.get();
      for( j=0 ; j<table_num_points[i] ; j++ ){
	data_file >> table_abscissas[table_data_offset];
	data_file >> table_ordinates[table_data_offset];
	data_file.get();
	table_data_offset++;
      }
    }
  } else {

    // get the load step at which to restart the analysis
    data_file.getline(line,sizeof(line));
    data_file >> start_step;
    data_file.get();
    if (start_step > 1 ) num_load_steps = num_load_steps + start_step - 1;

    // Read the number of face blocks
    data_file.getline(line,sizeof(line));
    data_file >> number_face_blocks;
    data_file.get();

    // Read the face block types
    if( number_face_blocks > 0 ){
      data_file.getline(line,sizeof(line));
      for( i=0 ; i<number_face_blocks ; i++ ){
	data_file >> face_block_types[i];
	if( face_block_types[i] == 5 || face_block_types[i] == 6 ){
	  data_file.get();
	  data_file.getline(line,sizeof(line));
	  data_file >> shell_block_offset[i];
	} else
	  shell_block_offset[i] = -1.0;
      }
      data_file.get();
    }

    // read the number of entity keys
    data_file.getline(line,sizeof(line));
    data_file >> ikeys;
    data_file.get();
  }
  data_file.getline(line,sizeof(line));
  data_file >> enforcement_type;
  data_file.get();
  if( enforcement_type ){
    if( restart_flag != 1 && restart_flag != 3 ){
      switch( enforcement_type ){
      case 1:{
	data_file.getline(line,sizeof(line));
	for( i=0 ; i<ikeys*ikeys*ContactTDEnforcement::NSIZED ; i++ )
	  data_file >> enforcement_data[i];
	data_file.get();
	data_file.getline(line,sizeof(line));
	data_file >> enforcement_data_vars[0];
	data_file.get();
	data_file >> enforcement_data_vars[1];
	data_file.get();
	data_file >> enforcement_data_vars[2];
	data_file.get();
        data_file >> enforcement_data_vars[3];
	data_file.get();
	if( (int) enforcement_data_vars[3] == 1 ){
          data_file >> enforcement_data_vars[4];
          data_file.get();
          data_file >> enforcement_data_vars[5];
          data_file.get();
        }
	break;
      }
      case 2:{
	data_file.getline(line,sizeof(line));
	for( i=0 ; i<ikeys*ikeys*ContactGapRemoval::NSIZED ; i++ )
	  data_file >> enforcement_data[i];
	data_file.get();
	data_file.getline(line,sizeof(line));
	data_file >> enforcement_data_vars[0];
	data_file.get();
	data_file >> enforcement_data_vars[1];
	data_file.get();
	break;
      }
      case 3:{
	data_file.getline(line,sizeof(line));
	for( i=0 ; i<ikeys*ikeys*ContactTiedKinematics::NSIZED ; i++ )
	  data_file >> enforcement_data[i];
	data_file.get();
        break;
      }
      case 4:{
	data_file.getline(line,sizeof(line));
	for( i=0 ; i<ikeys*ikeys*ContactVolumeTransfer::NSIZED ; i++ )
	  data_file >> enforcement_data[i];
	data_file.get();
        break;
      }
      case 5:{
	// currently no enforcement data for MPCs
	break;
      }
#ifdef CONTACT_TD_FACE_FACE_ENF
      case 6:{
	data_file.getline(line,sizeof(line));
	for( i=0 ; i<ikeys*ikeys*ContactTDFaceFaceEnf::NSIZED ; i++ )
	  data_file >> enforcement_data[i];
	data_file.get();
        break;
      }
#endif
      }
      // read number of enforcement models
      data_file.getline(line,sizeof(line));
      data_file >> num_enf_models;
      data_file.get();  // parse the newline special character
      
      for( i=0 ; i<num_enf_models ; i++ ){
	// Get the type and id
	data_file.getline(line,sizeof(line));
	data_file >> enf_model_types[i];
	data_file.get();
	data_file.getline(line,sizeof(line));
	data_file >> enf_model_ids[i];
	data_file.get();  // parse the newline special character
	switch( enf_model_types[i] ){
	case ContactEnforcement::TD_FRICTIONLESS: {
	  data_file.getline(line,sizeof(line));
	  break;
	}
	case ContactEnforcement::TD_CONSTANT_FRICTION: { 
	  data_file.getline(line,sizeof(line));
	  data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+0];
	  data_file.get();  // parse the newline special character
	  break;
	}
	case ContactEnforcement::TD_TIED: { 
	  data_file.getline(line,sizeof(line));
	  break;
	}
	case ContactEnforcement::TD_SPOT_WELD: {
	  for( j=0 ; j<4 ; j++ ){
	    data_file.getline(line,sizeof(line));
	    data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
	    data_file.get();   // parse the newline special character
	  }
	  break;
	}
	case ContactEnforcement::TD_PRESSURE_DEPENDENT: {
	  for( j=0 ; j<4 ; j++ ){
	    data_file.getline(line,sizeof(line));
	    data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
	    data_file.get();   // parse the newline special character
	  }	  
	  break;
	}
	case ContactEnforcement::TD_VELOCITY_DEPENDENT: {
	  for( j=0 ; j<3 ; j++ ){
	    data_file.getline(line,sizeof(line));
	    data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
	    data_file.get();   // parse the newline special character
	  }	  
	  break;
	}
	case ContactEnforcement::TD_PV_DEPENDENT: {
	  for( j=0 ; j<6 ; j++ ){
	    data_file.getline(line,sizeof(line));
	    data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
	    data_file.get();   // parse the newline special character
	  }	  
	  break;
	}
        case ContactEnforcement::TD_SPRING_WELD: {
          for( j=0 ; j<5 ; j++ ){
            data_file.getline(line,sizeof(line));
            data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
            data_file.get();   // parse the newline special character
          }      
          break;
        }
        case ContactEnforcement::TD_THREADED: {
          for( j=0 ; j<8 ; j++ ){
            data_file.getline(line,sizeof(line));
            data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
            data_file.get();   // parse the newline special character
          }      
          break;
        }
        case ContactEnforcement::TD_ADHESION: {
            data_file.getline(line,sizeof(line));
            data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS];
            data_file.get();   // parse the newline special character
          break;
        }
        case ContactEnforcement::TD_JUNCTION: {
          for( j=0 ; j<3 ; j++ ){
            data_file.getline(line,sizeof(line));
            data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
            data_file.get();   // parse the newline special character
          }      
          break;
        }
        case ContactEnforcement::TD_COHESIVE_ZONE: {
          for( j=0 ; j<3 ; j++ ){
            data_file.getline(line,sizeof(line));
            data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
            data_file.get();   // parse the newline special character
          }      
          break;
        }
	case ContactEnforcement::TD_AREA_WELD: {
	  for( j=0 ; j<4 ; j++ ){
	    data_file.getline(line,sizeof(line));
	    data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+j];
	    data_file.get();   // parse the newline special character
	  }
	  break;
	}
        case ContactEnforcement::TD_USER: {
          data_file.getline(line,sizeof(line));
          data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS];
          data_file.get();   // parse the newline special character
          data_file.getline(line,sizeof(line));
          data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+1];
          data_file.get();   // parse the newline special character
          int n = (int) (enf_model_data[i*MAX_ENF_MODEL_DATA_VARS]+
                         enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+1]+2);
	  for( j=0 ; j<n ; j++ ){
	    data_file.getline(line,sizeof(line));
	    data_file >> enf_model_data[i*MAX_ENF_MODEL_DATA_VARS+2+j];
	    data_file.get();   // parse the newline special character
	  }
          data_file.getline(line,sizeof(line));
          data_file >> enf_model_sub_num[i];
          data_file.get();
          data_file.getline(line,sizeof(line));
	  for( j=0 ; j<enf_model_sub_num[i] ; j++ ){
            char id_name[64];
	    data_file >> id_name;
            if (strcmp(id_name,"INIT_MODEL_FN")==0) {
	      enf_model_sub_ids[i*MAX_ENF_MODEL_SUBS+j] = 0;
            } else if (strcmp(id_name,"INIT_TIME_STEP_FN")==0) {
	      enf_model_sub_ids[i*MAX_ENF_MODEL_SUBS+j] = 1;
            } else if (strcmp(id_name,"INIT_NODE_STATE_DATA_FN")==0) {
	      enf_model_sub_ids[i*MAX_ENF_MODEL_SUBS+j] = 2;
            } else if (strcmp(id_name,"INTERACTION_TYPE_FN")==0) {
	      enf_model_sub_ids[i*MAX_ENF_MODEL_SUBS+j] = 3;
            } else if (strcmp(id_name,"INTERACTION_ACTIVE_FN")==0) {
	      enf_model_sub_ids[i*MAX_ENF_MODEL_SUBS+j] = 4;
            } else if (strcmp(id_name,"LIMIT_FORCE_FN")==0) {
	      enf_model_sub_ids[i*MAX_ENF_MODEL_SUBS+j] = 5;
            } else {
              cerr << "Warning: model name : " << id_name << ", not found\n";
	      enf_model_sub_ids[i*MAX_ENF_MODEL_SUBS+j] = -1;
            }
	    data_file >> id_name;
            strcpy(&enf_model_sub_names[(i*MAX_ENF_MODEL_SUBS+j)*64],id_name);
            data_file.get();
          }
          break;
        }

	default:
	  cerr << "Unrecognized enforcement model" << endl;
	}
      }
    } else {
      switch( enforcement_type ){
      case 1:{
	data_file.getline(line,sizeof(line));
	data_file >> enforcement_data_vars[0];
	data_file.get();
	data_file >> enforcement_data_vars[1];
	data_file.get();
 	data_file >> enforcement_data_vars[2];
	data_file.get();
        data_file >> enforcement_data_vars[3];
	data_file.get();
 	break;
      }
      case 2:
      case 3:
      case 4:
      case 5:
	break;
      }
    }
  } else
    num_enf_models = 0;

  // read flag for topology changes
  data_file.getline(line,sizeof(line));
  data_file >> update_topology;
  data_file.get();  // parse the newline special character

  if ( update_topology == 1) {
    // if update_topology = 1, then we have a list of faces to ignore for
    // each load step, plus another at the beginning (for testing element
    // birth)

    int count = 0;
    for (i = 0; i < num_load_steps; i++) {
      data_file.getline(line,sizeof(line));
      data_file >> num_face_ignore[i];
      data_file.get();  // parse the newline special character

      data_file.getline(line,sizeof(line));
      for (j = 0; j < num_face_ignore[i]; j++) {
	data_file >> face_ignore_list[count];
	data_file.get();  // parse the newline special character
	count++;
      }
    }
  } else {
    for (i = 0; i < num_load_steps; i++) {
      num_face_ignore[i] = 0;
    }
  }

  // read flag for dlb changes
  /*
  int* do_dlb_at_this_step = new int[num_load_steps];
  int do_dlb = 0;
  data_file.getline(line,sizeof(line));
  data_file >> do_dlb;
  data_file.get();  // parse the newline special character

  for (i = 0; i < num_load_steps; i++) {
    do_dlb_at_this_step[i] = 0;
  }
  if ( do_dlb == 1) {
    // if do_dlb = 1, then we have a list of load steps at which dlb occurs
    int count = 0;
    data_file.getline(line,sizeof(line));
    data_file >> count;
    data_file.get();  // parse the newline special character
    data_file.getline(line,sizeof(line));
    for (i = 0; i < count; i++) {
      int step = 0;
      data_file >> step;
      data_file.get();  // parse the newline special character
      do_dlb_at_this_step[step] = 1;
    }
  }
  delete [] do_dlb_at_this_step;
  */

  data_file.close();
}

