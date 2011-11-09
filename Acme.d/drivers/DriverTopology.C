// $Id: DriverTopology.C,v 2002.5 2004/06/24 14:51:35 mwglass Exp $

#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "params.h"

#include "exodusII.h"
#include "DriverTopology.h"
#include "AcmeNode.h"
#include "AcmeFace.h"
#include "AcmeElem.h"
#include "AcmeEntityBlock.h"
#include "AcmeBlockEntityList.h"
#include "AcmeZoltanID.h"
#include "AcmeZoltanComm.h"
#include "AcmeZoltanCallbacks.h"

#include "contact_assert.h"
#include <strings.h>
#include <string.h>
#include "Contact_Communication.h"

#ifndef CONTACT_NO_MPI
#ifndef CONTACT_NO_EXODUS_INPUT
#include "ne_nemesisI.h"
#endif
#endif

DriverTopology::DriverTopology(char* file, 
                           int& dim,
                           int* face_type, 
                           double* shell_block_offset,
                           int my_proc_id, 
                           int num_procs,
#ifndef CONTACT_NO_MPI
			   MPI_Comm& comm) : postream( comm )
#else
			   int& comm) : postream( comm )
#endif
{
  int i,j,k;
  float version;
  int comp_ws(8), io_ws(0);
  double* nodal_vars   = NULL;
  int* exo_truth_table = NULL;
  
  dimension            = 0;
  number_of_nodes      = 0;
  number_of_faces      = 0;
  number_of_elems      = 0;
  number_node_blocks   = 0;
  number_face_blocks   = 0;
  number_elem_blocks   = 0;
  node_list            = NULL;
  face_list            = NULL;
  elem_list            = NULL;
  node_block           = NULL;
  face_block           = NULL;
  elem_block           = NULL;
  num_node_vars        = 0;
  num_elem_vars        = 0;
  NumCommProcs         = 0;
  CommProcIds          = NULL;
  NumNodesToProc       = NULL;
  CommNodes            = NULL;
  ActiveNumCommProcs   = 0;
  ActiveCommProcIds    = NULL;
  ActiveNumNodesToProc = NULL;
  ActiveCommNodes      = NULL;
#ifndef CONTACT_NO_MPI
  zoltan               = NULL;
#endif
  AcmeComm             = comm;
  
  int idexo = ex_open (file, EX_READ, &comp_ws, &io_ws, &version);
  if ( idexo<0 ) {
    cerr << "ERROR: Cannot open " << file << endl;
    exit(1);
  }
  
  char exo_title[81];
  int exo_nodes    = 0;
  int exo_elements = 0;
  int exo_blocks   = 0;
  int exo_nodesets = 0;
  int exo_sidesets = 0;
  ex_get_init( idexo, exo_title, &dimension, &exo_nodes, &exo_elements,
	       &exo_blocks, &exo_nodesets, &exo_sidesets );
  if ( exo_nodes > MAX_NODES ) {
    cout << "FATAL ERROR: number of nods in mesh file" << endl
	 << "   exceeds the max number for which the driver was " << endl
	 << "   compiled. Recompile driver with MAX_NODES greater" << endl
	 << "   than or equal to " << exo_nodes << "." << endl;
    exit(1);
  }
  if ( exo_elements > MAX_ELEMENTS ) {
    cout << "FATAL ERROR: number of elements in mesh file" << endl
	 << "   exceeds the max number for which the driver was " << endl
	 << "   compiled. Recompile driver with MAX_ELEMENTS greater" << endl
	 << "   than or equal to " << exo_elements << "." << endl;
    exit(1);
  }
  if (exo_blocks  > MAX_BLOCKS ) {
    cout << "FATAL ERROR: number of blocks in mesh file" << endl
	 << "   exceeds the max number for which the driver was " << endl
	 << "   compiled. Recompile driver with MAX_BLOCKS greater" << endl
	 << "   than or equal to " << exo_blocks << "." << endl;
    exit(1);
  }
  number_of_nodes = exo_nodes;
  dim = dimension;
  
  int   nsteps = 0;
  float fdum;
  char  cdum;
  ex_inquire (idexo, EX_INQ_TIME, &nsteps, &fdum, &cdum);
  if (nsteps>0) {
    ex_get_var_param (idexo, (char*) "n", &num_node_vars);
    if ( num_node_vars > MAX_NODAL_VARS ) {
      cout << "FATAL ERROR: number of nodal variables in mesh file" << endl
	   << "   exceeds the max number for which the driver was " << endl
	   << "   compiled. Recompile driver with max vars greater" << endl
	   << "   than or equal to " << num_node_vars << "." << endl;
      exit(1);
    }
    ex_get_var_param (idexo, (char*) "e", &num_elem_vars);
    if ( num_elem_vars > MAX_ELEMENT_VARS ) {
      cout << "FATAL ERROR: number of element variables in mesh file" << endl
	   << "   exceeds the max number for which the driver was " << endl
	   << "   compiled. Recompile driver with MAX_ELEMENT_VARS greater" << endl
	   << "   than or equal to " << num_elem_vars << "." << endl;
      exit(1);
    }
    if (num_elem_vars>0) {
      exo_truth_table = new int[num_elem_vars*exo_blocks];
      ex_get_elem_var_tab(idexo,exo_blocks,num_elem_vars,exo_truth_table);
    }
    // read in all the nodal vars -- it will be the responsibility of
    // the search/enforcement to interpret what these variables mean
    if (num_node_vars>0) {
      nodal_vars = new double[num_node_vars*number_of_nodes];
      for ( i = 0; i < num_node_vars; i++) {
        ex_get_nodal_var( idexo, 1, i+1, number_of_nodes, 
                          &nodal_vars[i*number_of_nodes] );
      }
    }
  }
  
  int* node_ids = new int [number_of_nodes];
  ex_get_node_num_map(idexo, node_ids);
  Real* x_coord = new Real[number_of_nodes];
  Real* y_coord = new Real[number_of_nodes];
  Real* z_coord = new Real[number_of_nodes];
  ex_get_coord( idexo, x_coord, y_coord, z_coord );

  number_node_blocks = 0;
  number_face_blocks = 0;
  number_elem_blocks = 0;
  number_of_faces = 0;
  number_of_elems = 0;
  char elem_name[81];
  AcmeEntity::DerivedEntityType* elem_type = new AcmeEntity::DerivedEntityType[exo_blocks];
  for( i=0 ; i<exo_blocks ; i++ ){
    int num_nodes_per_elem, num_attr, entities_per_block;
    ex_get_elem_block( idexo, i+1, elem_name, &entities_per_block, 
		       &num_nodes_per_elem, &num_attr );
    // 
    // FACE BLOCKS
    //
    if( strncasecmp(elem_name,"SHELL",5) == 0 ||
	strncasecmp(elem_name,"QUAD",4)  == 0){
      if( num_nodes_per_elem == 4 ){
        elem_type[i] = AcmeEntity::ACME_QUADFACEL4;
      }
      else if( num_nodes_per_elem == 8 ){
        elem_type[i] = AcmeEntity::ACME_QUADFACEQ8;
      }
    }
    else if( strncasecmp(elem_name,"TRI",3)==0 && num_nodes_per_elem == 3 )
      elem_type[i] = AcmeEntity::ACME_TRIFACEL3;
    else if( strncasecmp(elem_name,"TRI",3)==0 && num_nodes_per_elem == 6 )
      elem_type[i] = AcmeEntity::ACME_TRIFACEQ6;
    // 
    // ELEMENT BLOCKS
    //
    else if( strncasecmp(elem_name,"HEX",3)==0 )
      elem_type[i] = AcmeEntity::ACME_HEXELEML8;
    // 
    // NODE BLOCKS
    //
    else if( strncasecmp(elem_name,"SPH",3)==0 ){
      elem_type[i] = AcmeEntity::ACME_POINT;
    }
    // 
    // UNKNOWN (i.e. an empty block in parallel)
    //
    else if( strncasecmp(elem_name,"NULL",4)==0 )
      elem_type[i] = AcmeEntity::ACME_UNKNOWN;
    elem_type[i] = (AcmeEntity::DerivedEntityType)contact_global_maximum( (int)elem_type[i], comm );
    switch (elem_type[i]) {
    case AcmeEntity::ACME_POINT:
      number_node_blocks++;
      break;
    case AcmeEntity::ACME_QUADFACEL4:
    case AcmeEntity::ACME_QUADFACEQ8:
    case AcmeEntity::ACME_TRIFACEL3:
    case AcmeEntity::ACME_TRIFACEQ6:
    case AcmeEntity::ACME_LINEFACEL2:
    case AcmeEntity::ACME_LINEFACEQ3:
      number_of_faces += entities_per_block;
      if (face_type[number_face_blocks]==5) number_of_faces += entities_per_block;
      if (face_type[number_face_blocks]==6) number_of_faces += entities_per_block;
      number_face_blocks++;
      break;
    case AcmeEntity::ACME_HEXELEML8:
      number_of_elems += entities_per_block;
      number_elem_blocks++;
      break;
    }
  }
  
  // Put nodes that are in a node set into a new node block
  number_node_blocks += exo_nodesets;
  if (number_face_blocks+number_elem_blocks) number_node_blocks++;
  
  int* elem_ids = new int [exo_elements];
  ex_get_elem_num_map(idexo, elem_ids);

  node_block = new AcmeEntityBlock*[number_node_blocks];
  face_block = new AcmeEntityBlock*[number_face_blocks];
  elem_block = new AcmeEntityBlock*[number_elem_blocks];
  
  int n_block = 0;
  int f_block = 0;
  int e_block = 0;
  if (number_face_blocks+number_elem_blocks) n_block++;
  n_block   += exo_nodesets;
  int  index = 0;
  int* node_connectivity = new int[number_of_nodes];
  int  econn[50];
  for( i=0 ; i<exo_blocks ; i++ ){
    int num_nodes_per_elem, num_attr, entities_per_block;
    ex_get_elem_block( idexo, i+1, elem_name, &entities_per_block, 
		       &num_nodes_per_elem, &num_attr );
    switch (elem_type[i]) {
    case AcmeEntity::ACME_POINT: {
      int* exo_ids=NULL;
      if( entities_per_block>0 ) {
        exo_ids = new int[entities_per_block];
        ex_get_elem_conn( idexo, i+1, node_connectivity );
        for (j=0; j<entities_per_block; j++) {
          exo_ids[j] = node_ids[node_connectivity[j]-1];
          node_ids[node_connectivity[j]-1] = -1;
        }
      }
      node_block[n_block] = new AcmeEntityBlock( AcmeEntity::CT_NODE,
                                                 elem_type[i], 
					         n_block,
					         entities_per_block,
                                                 exo_ids );
      node_block[n_block]->NodeType(ContactSearch::POINT);
      if (exo_ids!=NULL) delete [] exo_ids;
      if( entities_per_block>0 ) {
        int* connectivity = new int[num_nodes_per_elem*entities_per_block];
        ex_get_elem_conn( idexo, i+1, connectivity );
        AcmeBlockEntityList* entity_list = node_block[n_block]->EntityList();
        if (num_attr>0) {
          node_block[n_block]->NumAttributes(num_attr);
          Real* attr = new double[num_attr*entities_per_block];
          ex_get_elem_attr (idexo, i+1, attr);
          j = 0;
          entity_list->IteratorStart();
          while (AcmeEntity* entity = entity_list->IteratorForward()) {
            entity->NumAttributes(num_attr);
            entity->Attributes(&attr[j]);
            j += num_attr;
          }
          delete [] attr;
        }
        if (num_node_vars>0) {
          Real* data = new Real[num_node_vars];
          j = 0;
          entity_list->IteratorStart();
          while (AcmeEntity* entity = entity_list->IteratorForward()) {
            for (k=0; k<num_node_vars; k++) data[k] = nodal_vars[k*number_of_nodes+connectivity[j]-1];
            entity->NumData(num_node_vars);
            entity->Data(data);
            j++;
          }
          delete [] data;
        }
        j = 0;
        entity_list->IteratorStart();
        while (AcmeEntity* entity = entity_list->IteratorForward()) {
          AcmeNode* node = static_cast<AcmeNode*>(entity);
          k = connectivity[j]-1;
          node->Coordinates(x_coord[k],y_coord[k],z_coord[k]);
          j++;
        }
        delete [] connectivity;
      }
      n_block++;
      }break;
    case AcmeEntity::ACME_QUADFACEL4:
    case AcmeEntity::ACME_QUADFACEQ8:
    case AcmeEntity::ACME_TRIFACEL3:
    case AcmeEntity::ACME_TRIFACEQ6:
    case AcmeEntity::ACME_LINEFACEL2:
    case AcmeEntity::ACME_LINEFACEQ3: {
      if (face_type[f_block]==5 || face_type[f_block]==6) {
        // special shell case for shells
        int* exo_ids=NULL;
        if( entities_per_block>0 ) {
          exo_ids = new int[2*entities_per_block];
          for (j=0; j<entities_per_block; j++) {
            exo_ids[                   j] = 2*elem_ids[index+j];
            exo_ids[entities_per_block+j] = 2*elem_ids[index+j]+1;
          }
        }
        face_block[f_block] = new AcmeEntityBlock( AcmeEntity::CT_FACE,
                                                   elem_type[i], 
					           f_block,
					           2*entities_per_block,
                                                   exo_ids );
        face_block[f_block]->FaceType((ContactSearch::ContactFace_Type)face_type[f_block]);
        face_block[f_block]->IsShell(true);
        face_block[f_block]->ShellOffset(shell_block_offset[f_block]);
        if (exo_ids!=NULL) delete [] exo_ids;
        if( entities_per_block>0 ) {
          AcmeBlockEntityList* entity_list = face_block[f_block]->EntityList();
          int* connectivity = new int[num_nodes_per_elem*2*entities_per_block];
          ex_get_elem_conn( idexo, i+1, connectivity );
          int* c0 = &connectivity[0];
          int* c1 = &connectivity[num_nodes_per_elem*entities_per_block];
          for (j=0; j<entities_per_block; j++) {
            c1[0] = c0[0];
            for (k=1; k<num_nodes_per_elem; k++) {
	      c1[num_nodes_per_elem-k] = c0[k];
            }
            c0 += num_nodes_per_elem;
            c1 += num_nodes_per_elem;
          }
          j = 0;
          entity_list->IteratorStart();
          while (AcmeEntity* entity = entity_list->IteratorForward()) {
            for (k=0; k<num_nodes_per_elem; k++) {
              econn[k] = node_ids[connectivity[j+k]-1];
            }
            AcmeFace* face = static_cast<AcmeFace*>(entity);
            face->NodeIds(&econn[0]);
            j += num_nodes_per_elem;
          }
          delete [] connectivity;
          if (num_attr>0) {
            face_block[f_block]->NumAttributes(num_attr);
            Real* attr = new double[num_attr*2*entities_per_block];
            ex_get_elem_attr (idexo, i+1, &attr[0]);
            memcpy( &attr[num_attr*entities_per_block], attr,
                    num_attr*entities_per_block*sizeof(Real));
            j = 0;
            entity_list->IteratorStart();
            while (AcmeEntity* entity = entity_list->IteratorForward()) {
              entity->NumAttributes(num_attr);
              entity->Attributes(&attr[j]);
              j += num_attr;
            }
            delete [] attr;
          }
          if (num_elem_vars>0) {
            Real* data0 = new Real[num_elem_vars];
            Real* data1 = new double[num_elem_vars*entities_per_block];
            for (j=0; j<num_elem_vars; j++) {
              if (exo_truth_table[num_elem_vars*i+j]) {
                ex_get_elem_var(idexo,1,j+1,i+1,entities_per_block,
                                &data1[j*entities_per_block]);
                memcpy( &data1[j*2*entities_per_block+entities_per_block], 
                        &data1[j*2*entities_per_block],
                        entities_per_block*sizeof(double));
              }
            }
            j = 0;
            entity_list->IteratorStart();
            while (AcmeEntity* entity = entity_list->IteratorForward()) {
              for (k=0; k<num_elem_vars; k++) data0[k] = data1[k*entities_per_block+j];
              entity->NumData(num_elem_vars);
              entity->Data(data0);
              j++;
            }
            delete [] data0;
            delete [] data1;
          }
        }
      } else {
        int* exo_ids=NULL;
        if( entities_per_block>0 ) {
          exo_ids = new int[entities_per_block];
          for (j=0; j<entities_per_block; j++) {
            exo_ids[j] = 2*elem_ids[index+j];
          }
        }
        face_block[f_block] = new AcmeEntityBlock( AcmeEntity::CT_FACE,
                                                   elem_type[i], 
					           f_block,
					           entities_per_block,
                                                   exo_ids );
        face_block[f_block]->FaceType((ContactSearch::ContactFace_Type)face_type[f_block]);
        if (exo_ids!=NULL) delete [] exo_ids;
        if( entities_per_block>0 ) {
          AcmeBlockEntityList* entity_list = face_block[f_block]->EntityList();
          int* connectivity = new int[num_nodes_per_elem*entities_per_block];
          ex_get_elem_conn( idexo, i+1, connectivity );
          j = 0;
          entity_list->IteratorStart();
          while (AcmeEntity* entity = entity_list->IteratorForward()) {
            for (k=0; k<num_nodes_per_elem; k++) {
              econn[k] = node_ids[connectivity[j+k]-1];
            }
            AcmeFace* face = static_cast<AcmeFace*>(entity);
            face->NodeIds(&econn[0]);
            j += num_nodes_per_elem;
          }
          delete [] connectivity;
          if (num_attr>0) {
            face_block[f_block]->NumAttributes(num_attr);
            Real* attr = new double[num_attr*entities_per_block];
            ex_get_elem_attr (idexo, i+1, attr);
            j = 0;
            entity_list->IteratorStart();
            while (AcmeEntity* entity = entity_list->IteratorForward()) {
              entity->NumAttributes(num_attr);
              entity->Attributes(&attr[j]);
              j += num_attr;
            }
            delete [] attr;
          }
          if (num_elem_vars>0) {
            Real* data0 = new Real[num_elem_vars];
            Real* data1 = new double[num_elem_vars*entities_per_block];
            for (j=0; j<num_elem_vars; j++) {
              if (exo_truth_table[num_elem_vars*i+j]) {
                ex_get_elem_var(idexo,1,j+1,i+1,entities_per_block,
                                &data1[j*entities_per_block]);
              }
            }
            j = 0;
            entity_list->IteratorStart();
            while (AcmeEntity* entity = entity_list->IteratorForward()) {
              for (k=0; k<num_elem_vars; k++) data0[k] = data1[k*entities_per_block+j];
              entity->NumData(num_elem_vars);
              entity->Data(data0);
              j++;
            }
            delete [] data0;
            delete [] data1;
          }
        }
      }
      f_block++;
      }break;
    case AcmeEntity::ACME_HEXELEML8: {
      int* exo_ids=NULL;
      if( entities_per_block>0 ) {
        exo_ids = new int[entities_per_block];
        for (j=0; j<entities_per_block; j++) {
          exo_ids[j] = 2*elem_ids[index+j];
        }
      }
      elem_block[e_block] = new AcmeEntityBlock( AcmeEntity::CT_ELEM,
                                                 elem_type[i], 
					         e_block,
					         entities_per_block,
                                                 exo_ids );
      if( strncasecmp(elem_name,"HEX8CAR",8)==0 ) {
        elem_block[e_block]->ElemType(ContactSearch::CARTESIANHEXELEMENTL8);
      } else {
        elem_block[e_block]->ElemType(ContactSearch::HEXELEMENTL8);
      }
      if (exo_ids!=NULL) delete [] exo_ids;
      if( entities_per_block>0 ) {
        AcmeBlockEntityList* entity_list = elem_block[e_block]->EntityList();
        int* connectivity = new int[num_nodes_per_elem*entities_per_block];
        ex_get_elem_conn( idexo, i+1, connectivity );
        j = 0;
        entity_list->IteratorStart();
        while (AcmeEntity* entity = entity_list->IteratorForward()) {
            for (k=0; k<num_nodes_per_elem; k++) {
              econn[k] = node_ids[connectivity[j+k]-1];
            }
            AcmeElem* elem = static_cast<AcmeElem*>(entity);
            elem->NodeIds(&econn[0]);
            j += num_nodes_per_elem;
        }
        delete [] connectivity;
        if (num_attr>0) {
          elem_block[e_block]->NumAttributes(num_attr);
          Real* attr = new double[num_attr*entities_per_block];
          ex_get_elem_attr (idexo, i+1, attr);
          j = 0;
          entity_list->IteratorStart();
          while (AcmeEntity* entity = entity_list->IteratorForward()) {
            entity->NumAttributes(num_attr);
	    entity->Attributes(&attr[j]);
            j += num_attr;
          }
          delete [] attr;
        }
        if (num_elem_vars>0) {
          Real* data0 = new Real[num_elem_vars];
          Real* data1 = new double[num_elem_vars*entities_per_block];
          for (j=0; j<num_elem_vars; j++) {
            if (exo_truth_table[num_elem_vars*i+j]) {
              ex_get_elem_var(idexo,1,j+1,i+1,entities_per_block,
                              &data1[j*entities_per_block]);
            }
          }
          j = 0;
          entity_list->IteratorStart();
          while (AcmeEntity* entity = entity_list->IteratorForward()) {
            for (k=0; k<num_elem_vars; k++) data0[k] = data1[k*entities_per_block+j];
            entity->NumData(num_elem_vars);
	    entity->Data(data0);
            j++;
          }
          delete [] data0;
          delete [] data1;
        }
      }
      e_block++;
      break;
      }
    }
    index += entities_per_block;
  }
  delete [] elem_type;
  delete [] elem_ids;
  delete [] node_connectivity;
  if (exo_truth_table) delete [] exo_truth_table;

  if( exo_nodesets > 0 ){
    int* nset_ids = new int[exo_nodesets];
    ex_get_node_set_ids( idexo, nset_ids );
    int num_nset_nodes = 0;
    for( i=0 ; i<exo_nodesets ; i++ ){
      int num_nodes_in_set;
      int num_dist_fact;
      ex_get_node_set_param( idexo, nset_ids[i], &num_nodes_in_set,
			     &num_dist_fact );
      num_nset_nodes += num_nodes_in_set;
    }
    int* list_of_nodes = new int[number_of_nodes];
    for( i=0 ; i<exo_nodesets ; i++ ){
      int num_nodes_in_set;
      int num_dist_fact;
      ex_get_node_set( idexo, nset_ids[i], list_of_nodes );
      ex_get_node_set_param( idexo, nset_ids[i], &num_nodes_in_set,
			     &num_dist_fact );
      int* exo_ids = new int[num_nodes_in_set];
      for (j=0; j<num_nodes_in_set; j++) {
        exo_ids[j] = node_ids[list_of_nodes[j]-1];
        node_ids[list_of_nodes[j]-1] = -1;
      }                  
      node_block[i+1] = new AcmeEntityBlock( AcmeEntity::CT_NODE, 
                                             AcmeEntity::ACME_NODE, 
					     i+1,
					     num_nodes_in_set,
                                             exo_ids );
      node_block[i+1]->NodeType(ContactSearch::NODE);
      delete [] exo_ids;
      AcmeBlockEntityList* entity_list = node_block[i+1]->EntityList();
      if (num_node_vars>0) {
        j = 0;
        Real* data = new Real[num_node_vars];
        entity_list->IteratorStart();
        while (AcmeEntity* entity = entity_list->IteratorForward()) {
          for (k=0; k<num_node_vars; k++) data[k] = nodal_vars[k*number_of_nodes+list_of_nodes[j]-1];
          entity->NumData(num_node_vars);
          entity->Data(data);
          j++;
        }
        delete [] data;
      }
      j = 0;
      entity_list->IteratorStart();
      while (AcmeEntity* entity = entity_list->IteratorForward()) {
        AcmeNode* node = static_cast<AcmeNode*>(entity);
        k = list_of_nodes[j]-1;
        node->Coordinates(x_coord[k],y_coord[k],z_coord[k]);
        j++;
      }
    }
    delete [] list_of_nodes;
    delete [] nset_ids;
  }  
  
  if (number_face_blocks+number_elem_blocks) {
    int nnodes = number_of_nodes;
    for( i=1 ; i<number_node_blocks ; i++ ) {
      nnodes -= node_block[i]->Number_of_Entities();
    }
    int* node_map = new int [number_of_nodes];
    int* exo_ids = new int[nnodes];
    for (j=0, i=0; i<number_of_nodes; i++) {
      if (node_ids[i]>0) {
        exo_ids[j]    = node_ids[i];
        node_map[j++] = i;
      }
    }
    POSTCONDITION(j==nnodes);
    node_block[0] = new AcmeEntityBlock( AcmeEntity::CT_NODE, 
                                         AcmeEntity::ACME_NODE, 
                                         0,
                                         nnodes,
                                         exo_ids );
    node_block[0]->NodeType(ContactSearch::NODE);
    delete [] exo_ids;
    AcmeBlockEntityList* entity_list = node_block[0]->EntityList();
    if (num_node_vars>0) {
      j = 0;
      Real* data = new Real[num_node_vars];
      entity_list->IteratorStart();
      while (AcmeEntity* entity = entity_list->IteratorForward()) {
        for (k=0; k<num_node_vars; k++) {
          data[k] = nodal_vars[k*number_of_nodes+node_map[j]];
        }
        entity->NumData(num_node_vars);
        entity->Data(data);
        j++;
      }
      delete [] data;
    }
    j = 0;
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      k = node_map[j];
      node->Coordinates(x_coord[k],y_coord[k],z_coord[k]);
      j++;
    }
    delete [] node_map;
  }
  for( i=0 ; i<number_node_blocks ; i++ ) {
    node_block[i]->EntityList()->Sort();
  }
  
  node_list = new AcmeTopologyEntityList();
  node_list->BuildList(node_block, number_node_blocks);
  POSTCONDITION(number_of_nodes==node_list->NumEntities());
  POSTCONDITION(number_node_blocks==node_list->Num_Blocks());
  face_list = new AcmeTopologyEntityList();
  face_list->BuildList(face_block, number_face_blocks);
  POSTCONDITION(number_of_faces==face_list->NumEntities());
  POSTCONDITION(number_face_blocks==face_list->Num_Blocks());
  elem_list = new AcmeTopologyEntityList();
  elem_list->BuildList(elem_block, number_elem_blocks);
  POSTCONDITION(number_of_elems==elem_list->NumEntities());
  POSTCONDITION(number_elem_blocks==elem_list->Num_Blocks());
  
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      face->ConnectNodes(node_list);
    }
  }
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      elem->ConnectNodes(node_list);
    }
  }
  
#ifndef CONTACT_NO_MPI
  if( num_procs > 1 ){
    create_zoltan_object();
    int  num_internal_nodes(0), num_border_nodes(0), num_external_nodes(0);
    int  num_internal_elems(0), num_border_elems(0), num_elem_comm_procs(0);
    ne_get_loadbal_param( idexo, &num_internal_nodes, &num_border_nodes,
			  &num_external_nodes, &num_internal_elems,
			  &num_border_elems, &NumCommProcs,
			  &num_elem_comm_procs, my_proc_id );
    PRECONDITION( MAX_PROCS >= num_procs );
    PRECONDITION( MAX_NODES >= number_of_nodes );
    int* comm_elem_proc_ids = new int[num_procs];
    int* num_elems_to_proc  = new int[num_procs];
    int* comm_ids           = new int[num_procs * number_of_nodes];
    CommProcIds             = new int[MAX_PROCS];
    NumNodesToProc          = new int[MAX_PROCS];
    CommNodes               = new int[MAX_PROCS * MAX_NODES];
    ActiveCommProcIds       = new int[MAX_PROCS];
    ActiveNumNodesToProc    = new int[MAX_PROCS];
    ActiveCommNodes         = new int[MAX_PROCS * MAX_NODES];
    ne_get_cmap_params( idexo, CommProcIds, NumNodesToProc,
			comm_elem_proc_ids, num_elems_to_proc, my_proc_id );

    int icount = 0;
    for( i=0 ; i<NumCommProcs ; i++ ){
      ne_get_node_cmap( idexo, CommProcIds[i], &CommNodes[icount],
		        &comm_ids[icount], my_proc_id );
      icount += NumNodesToProc[i];
    }
    // renumber comm_nodes to reflect the renumbering of 
    // nodes that occurs when we have multiple node blocks
    ex_get_node_num_map(idexo, node_ids);
    for( i=0 ; i<icount ; i++ ) {
      CommNodes[i] = node_list->Find(node_ids[CommNodes[i]-1])->ProcArrayIndex()+1;
    }
    delete [] comm_elem_proc_ids;
    delete [] num_elems_to_proc;
    delete [] comm_ids;
  }
#endif
  Set_Ownership();
  number_of_active_nodes = number_of_nodes;
  number_of_active_faces = number_of_faces;
  number_of_active_elems = number_of_elems;
  if (nodal_vars) delete [] nodal_vars;
  if (x_coord) delete [] x_coord;
  if (y_coord) delete [] y_coord;
  if (z_coord) delete [] z_coord;
  if (node_ids) delete [] node_ids;
}

void
DriverTopology::GetInitialTopology(int& num_node_ignore,
                                 int* node_ignore_list,
				 int& num_face_ignore,
                                 int* face_ignore_list,
				 int& num_elem_ignore,
                                 int* elem_ignore_list,
                                 
                                 int& num_node_blocks, 
				 int* node_types,
				 int* nodes_per_block,
				 int* node_eids,
                                 int* node_ids,
                                 
				 int& num_face_blocks, 
                                 int* face_types,
				 int* faces_per_block, 
                                 int* face_ids,
                                 int* face_connectivity,
                                 
				 int& num_element_blocks,
				 int* element_types,
				 int* elements_per_block, 
                                 int* element_ids,
                                 int* element_connectivity,
                                 
				 double* position,
				 double* nodal_vars, 
                                 int&    num_nodal_vars,
				 double* element_vars,
				 int&    num_element_vars,
				 double* attributes,
				 int*    num_attributes,
				 double* node_radius,
				 double* shell_thickness,
				 double* shell_offset,
                                 
				 int& num_comm_procs, 
                                 int* comm_proc_ids,
				 int* num_nodes_to_proc, 
                                 int* comm_nodes,
                                 
                                 int  my_proc_id,
                                 int  num_procs)
{
  int  i, k, n, nn;
  int* conn;
  
  num_nodal_vars   = num_node_vars;
  num_element_vars = num_elem_vars;
  double* nvars    = nodal_vars;
  double* evars    = element_vars;
  double* attr     = attributes;
  
  // Check for nodes to make inactive
  for( i=1 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      int base_id = node->Exodus_ID();
      for (k=0; k<num_node_ignore; k++) {
        if (base_id==node_ignore_list[k]) {
          node->State(AcmeEntity::INACTIVE);
          break;
        }
      }
    }
  }
  
  // Check for faces to make inactive and flag
  // their connected nodes as being inactive
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      int base_id = face->Exodus_ID()/2;
      for (k=0; k<num_face_ignore; k++) {
        if (base_id==face_ignore_list[k]) {
          face->State(AcmeEntity::INACTIVE);
          for (n=0; n<face->Nodes_Per_Face(); n++) {
            AcmeNode* node = face->Node(n);
            face->Node(n)->State(AcmeEntity::INACTIVE);
          }
          break;
        }
      }
    }
  }
  
  // Check for elements to make inactive and flag
  // their connected nodes as being inactive
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      int base_id = elem->Exodus_ID();
      for (k=0; k<num_elem_ignore; k++) {
        if (base_id==elem_ignore_list[k]) {
          elem->State(AcmeEntity::INACTIVE);
          for (n=0; n<elem->Nodes_Per_Element(); n++) {
            elem->Node(n)->State(AcmeEntity::INACTIVE);
          }
          break;
        }
      }
    }
  }
  
  // Make sure shared nodes on faces are marked as active
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        AcmeFace* face = static_cast<AcmeFace*>(entity);
        for (n=0; n<face->Nodes_Per_Face(); n++) {
          face->Node(n)->State(AcmeEntity::ACTIVE);
        }
      }
    }
  }
  
  // Make sure shared nodes on elements are marked as active
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        AcmeElem* elem = static_cast<AcmeElem*>(entity);
        for (n=0; n<elem->Nodes_Per_Element(); n++) {
          elem->Node(n)->State(AcmeEntity::ACTIVE);
        }
      }
    }
  }
  
  // Determine the total number of active nodes
  number_of_active_nodes = 0;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        entity->HostArrayIndex(number_of_active_nodes++);
      } else {
        entity->HostArrayIndex(-1);
      }
    }
  }
  
  // Determine the total number of active faces
  number_of_active_faces = 0;
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        entity->HostArrayIndex(number_of_active_faces++);
      } else {
        entity->HostArrayIndex(-1);
      }
    }
  }
  
  // Determine the total number of active elements
  number_of_active_elems = 0;
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        entity->HostArrayIndex(number_of_active_elems++);
      } else {
        entity->HostArrayIndex(-1);
      }
    }
  }
  
  for( i=0 ; i<number_of_active_faces ; i++ ) {
    shell_offset[i]    = -1.0;
    shell_thickness[i] =  0.0;
  }
  for( i=0 ; i<number_of_active_nodes ; i++ ) {
    node_radius[i] = 0.0;
  }
  
  double* n_radius = node_radius;
  double* coords   = position;
  num_node_blocks  = number_node_blocks;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    node_types[i] = node_block[i]->NodeType();
    nodes_per_block[i] = 0;
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      if (node->State() == AcmeEntity::ACTIVE) {
        n = node->HostArrayIndex();
	node_eids[n]    = node->Exodus_ID();
        node_ids[2*n  ] = 0;
        node_ids[2*n+1] = node->Exodus_ID();
        for (k=0; k<dimension; k++) {
          coords[k] = node->Coordinates()[k];
        }
        for (k=0; k<num_node_vars; k++) {
          nvars[number_of_active_nodes*k+n] = node->Data(k);
        }
        if( node_types[i] == 2 && node->NumAttributes()>0){
          *n_radius++ = node->Attributes(0);
        }
        coords += dimension;
        nodes_per_block[i]++;
      }
    }
  }
  
  double* s_thickness = shell_thickness;
  double* s_offset    = shell_offset;
  conn = face_connectivity;
  num_face_blocks = number_face_blocks;
  for( i=0 ; i<number_face_blocks ; i++ ) {
    num_attributes[i] = face_block[i]->NumAttributes();
    face_types[i] = face_block[i]->FaceType();
    faces_per_block[i] = 0;
    nn = 0;
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      if (face->State() == AcmeEntity::ACTIVE) {
        n = face->HostArrayIndex();
        face_ids[2*n  ] = 0;
        face_ids[2*n+1] = face->Exodus_ID();
        faces_per_block[i]++;
        for (k=0; k<face->Nodes_Per_Face(); k++) {
          conn[k] = face->Node(k)->HostArrayIndex()+1;
        }
        conn += face->Nodes_Per_Face();
        for (k=0; k<num_elem_vars; k++) {
          evars[(number_of_active_faces+number_of_active_elems)*k+n] = face->Data(k);
        }
        for (k=0; k<face->NumAttributes(); k++) {
          attr[face_block[i]->Number_of_Entities()*k+nn] = face->Attributes(k);
        }
        if( face_block[i]->IsShell() ){
          if (face->Exodus_ID()%2) {
            *s_offset++ = 1.0-face_block[i]->ShellOffset();
          } else {
            *s_offset++ = face_block[i]->ShellOffset();
          }
          *s_thickness++ = face->Attributes(0);
        }
        nn++;
      }
    }
    attr += face_block[i]->NumAttributes()*faces_per_block[i];
  }
  
  conn = element_connectivity;
  num_element_blocks = number_elem_blocks;
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    num_attributes[i] = elem_block[i]->NumAttributes();
    element_types[i]  = elem_block[i]->ElemType();
    elements_per_block[i] = 0;
    nn = 0;
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      if (elem->State() == AcmeEntity::ACTIVE) {
        n = elem->HostArrayIndex();
        element_ids[2*n  ] = 0;
        element_ids[2*n+1] = elem->Exodus_ID();
        elements_per_block[i]++;
        for (k=0; k<elem->Nodes_Per_Element(); k++) {
          conn[k] = elem->Node(k)->HostArrayIndex()+1;
        }
        conn += elem->Nodes_Per_Element();
        for (k=0; k<num_elem_vars; k++) {
          evars[(number_of_active_faces+number_of_active_elems)*k+n] = 
            elem->Data(k);
        }
        for (k=0; k<elem->NumAttributes(); k++) {
          attr[elem_block[i]->Number_of_Entities()*k+nn] = 
            elem->Attributes(k);
        }
        nn++;
      }
    }
    attr += num_attributes[i]*elements_per_block[i];
  }
  
#ifndef CONTACT_NO_MPI
  if (num_node_ignore+num_face_ignore+num_elem_ignore) {
    ComputeNodalCommPlan();
  }
  GetCommPlan( num_comm_procs,
               comm_proc_ids,
               num_nodes_to_proc,
               comm_nodes,
               my_proc_id,
               num_procs);
#endif

#ifdef DEBUG
  postream<<"GetInitialTopology()\n";
  if (num_face_ignore>0) {
    postream << "face ignore list: ";
    for (k=0; k<num_face_ignore; k++) {
      postream << " "<<face_ignore_list[k];
    }
    postream<<"\n";
  }
  DisplayByList();
  DisplayCommPlan(num_comm_procs, comm_proc_ids,
		  num_nodes_to_proc, comm_nodes);
#endif
}

void
DriverTopology::GetModsForBirthDeath(int& num_node_ignore,
                                   int* node_ignore_list,
			           int& num_face_ignore,
                                   int* face_ignore_list,
			           int& num_elem_ignore,
                                   int* elem_ignore_list,
                              
                                   int* num_node_deaths_per_block, 
	                           int* node_deaths_global_ids,
                                   int* num_face_deaths_per_block, 
	                           int* face_deaths_global_ids,
                                   int* num_elem_deaths_per_block, 
	                           int* elem_deaths_global_ids,
                                   
                                   int* num_node_births_per_block, 
	                           int* node_births_exodus_ids,
	                           int* node_births_global_ids,
	                           int* num_face_births_per_block, 
	                           int* face_births_global_ids,
	                           int* face_births_connectivity,
	                           int* num_elem_births_per_block, 
	                           int* elem_births_global_ids,
	                           int* elem_births_connectivity,
                                   
                                   int  my_proc_id,
                                   int  num_procs)
{
  int  i, k, n;
  int  n_birth, n_death;
  int* conn;
  
  // MARK UNCONNECTED NODES FOR DEATH
  for( i=1 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      int base_id = node->Exodus_ID();
      node->Action(AcmeEntity::NONE);
      k = num_node_ignore;
      for (k=0; k<num_node_ignore; k++) {
        if (base_id==node_ignore_list[k]) {
          node->Action(AcmeEntity::MARK_FOR_DEATH);
          break;
        }
      }
    }
  }
  // MARK UNCONNECTED NODES FOR BIRTH
  for( i=1 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      if (node->Action() == AcmeEntity::NONE) {
        node->Action(AcmeEntity::MARK_FOR_BIRTH);
      }
    }
  }
  
  // MARK FACES AND THEIR CONNECTED NODES FOR DEATH
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      int base_id = face->Exodus_ID()/2;
      face->Action(AcmeEntity::NONE);
      k = num_face_ignore;
      for (k=0; k<num_face_ignore; k++) {
        if (base_id==face_ignore_list[k]) {
          face->Action(AcmeEntity::MARK_FOR_DEATH);
          for (n=0; n<face->Nodes_Per_Face(); n++) {
            face->Node(n)->Action(AcmeEntity::MARK_FOR_DEATH);
          }
          break;
        }
      }
    }
  }
  // MARK FACES AND THEIR CONNECTED NODES FOR BIRTH
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      if (face->Action() == AcmeEntity::NONE) {
        face->Action(AcmeEntity::MARK_FOR_BIRTH);
        for (n=0; n<face->Nodes_Per_Face(); n++) {
          face->Node(n)->Action(AcmeEntity::MARK_FOR_BIRTH);
        }
      }
    }
  }
  
  // MARK ELEMENTS AND THEIR CONNECTED NODES FOR BIRTH
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      int base_id = elem->Exodus_ID();
      elem->Action(AcmeEntity::NONE);
      for (k=0; k<num_elem_ignore; k++) {
        if (base_id==elem_ignore_list[k]) {
          elem->Action(AcmeEntity::MARK_FOR_DEATH);
          for (n=0; n<elem->Nodes_Per_Element(); n++) {
            elem->Node(n)->Action(AcmeEntity::MARK_FOR_DEATH);
          }
          break;
        }
      }
    }
  }
  // MARK ELEMENTS AND THEIR CONNECTED NODES FOR DEATH
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      if (elem->Action() == AcmeEntity::NONE) {
        elem->Action(AcmeEntity::MARK_FOR_BIRTH);
        for (n=0; n<elem->Nodes_Per_Element(); n++) {
          elem->Node(n)->Action(AcmeEntity::MARK_FOR_BIRTH);
        }
      }
    }
  }
  
  //----------------------------------------------------------------
  
  n       = 0;
  n_birth = 0;
  n_death = 0;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    num_node_deaths_per_block[i] = 0;
    num_node_births_per_block[i] = 0;
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      if (node->State()==AcmeEntity::ACTIVE && 
          (node->Action()==AcmeEntity::MARK_FOR_DEATH ||
           node->Action()==AcmeEntity::MARK_FOR_DEATH1)) {
        node->State(AcmeEntity::INACTIVE);
        if (node->Action()==AcmeEntity::MARK_FOR_DEATH1) {
          node_deaths_global_ids[2*n_death  ] = 0;
          node_deaths_global_ids[2*n_death+1] = node->Exodus_ID();
          num_node_deaths_per_block[i]++;
          n_death++;
        }
      }
      if (node->State()==AcmeEntity::INACTIVE && node->Action()==AcmeEntity::MARK_FOR_BIRTH) {
        node->State(AcmeEntity::ACTIVE);
        node_births_exodus_ids[  n_birth  ] = node->Exodus_ID();
        node_births_global_ids[2*n_birth  ] = 0;
        node_births_global_ids[2*n_birth+1] = node->Exodus_ID();
        num_node_births_per_block[i]++;
        n_birth++;
      }
      if (node->State()==AcmeEntity::ACTIVE) {
        node->HostArrayIndex(n++);
      } else {
        node->HostArrayIndex(-1);
      }
      node->Action(AcmeEntity::NONE);
    }
  }
  
  n_birth = 0;
  n_death = 0;
  conn    = face_births_connectivity;
  for( i=0 ; i<number_face_blocks ; i++ ) {
    num_face_deaths_per_block[i] = 0;
    num_face_births_per_block[i] = 0;
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      if (face->State()==AcmeEntity::ACTIVE && face->Action()==AcmeEntity::MARK_FOR_DEATH) {
        face->State(AcmeEntity::INACTIVE);
        face_deaths_global_ids[2*n_death  ] = 0;
        face_deaths_global_ids[2*n_death+1] = face->Exodus_ID();
        num_face_deaths_per_block[i]++;
        n_death++;
      }
      if (face->State()==AcmeEntity::INACTIVE && face->Action()==AcmeEntity::MARK_FOR_BIRTH) {
        face->State(AcmeEntity::ACTIVE);
        face_births_global_ids[2*n_birth  ] = 0;
        face_births_global_ids[2*n_birth+1] = face->Exodus_ID();
        for (k=0; k<face->Nodes_Per_Face(); k++) {
          *conn++ = 0;
          *conn++ = face->Node(k)->Exodus_ID();
        }
        num_face_births_per_block[i]++;
        n_birth++;
      }
      face->Action(AcmeEntity::NONE);
    }
  }
  
  n_birth = 0;
  n_death = 0;
  conn    = elem_births_connectivity;
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    num_elem_deaths_per_block[i] = 0;
    num_elem_births_per_block[i] = 0;
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      if (elem->State()==AcmeEntity::ACTIVE && elem->Action()==AcmeEntity::MARK_FOR_DEATH) {
        elem->State(AcmeEntity::INACTIVE);
        elem_deaths_global_ids[2*n_death  ] = 0;
        elem_deaths_global_ids[2*n_death+1] = elem->Exodus_ID();
        num_elem_deaths_per_block[i]++;
        n_death++;
      }
      if (elem->State()==AcmeEntity::INACTIVE && elem->Action()==AcmeEntity::MARK_FOR_BIRTH) {
        elem->State(AcmeEntity::ACTIVE);
        elem_births_global_ids[2*n_birth  ] = 0;
        elem_births_global_ids[2*n_birth+1] = elem->Exodus_ID();
        for (k=0; k<elem->Nodes_Per_Element(); k++) {
          *conn++ = 0;
          *conn++ = elem->Node(k)->Exodus_ID();
        }
        num_elem_births_per_block[i]++;
        n_birth++;
      }
      elem->Action(AcmeEntity::NONE);
    }
  }
  
  // Update state of shared nodes
  //UpdateSharedNodeState();
  
  //---------------------------------------------------------------
    
  // Determine the total number of active nodes
  number_of_active_nodes = 0;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        entity->HostArrayIndex(number_of_active_nodes++);
      } else {
        entity->HostArrayIndex(-1);
      }
    }
  }
  
  // Determine the total number of active faces
  number_of_active_faces = 0;
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        number_of_active_faces++;
      }
    }
  }
  
  // Determine the total number of active elements
  number_of_active_elems = 0;
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        number_of_active_elems++;
      }
    }
  }

#ifdef DEBUG
  postream<<"GetModsForBirthDeath()\n";
  postream<<"  ignoring "<<num_node_ignore<<" nodes\n";
  postream<<"  ignoring "<<num_face_ignore<<" faces\n";
  postream<<"  ignoring "<<num_elem_ignore<<" elems\n";
  postream<<"Before Birth/Death...\n";
  postream<<"  num active nodes = "<<tmp_active_nodes<<"\n";
  postream<<"  num active faces = "<<tmp_active_faces<<"\n";
  postream<<"  num active elems = "<<tmp_active_elems<<"\n";
  postream<<"After Birth/Death...\n";
  postream<<"  num active nodes = "<<number_of_active_nodes<<"\n";
  postream<<"  num active faces = "<<number_of_active_faces<<"\n";
  postream<<"  num active elems = "<<number_of_active_elems<<"\n";
  for( n=0, i=0 ; i<number_node_blocks ; i++ ) {
    postream<<"killing "<<num_node_deaths_per_block[i]<<" nodes in block "<<i<<"\n";
    for (j=0; j<num_node_deaths_per_block[i]; j++) {
      postream<<"  node ("<<node_deaths_global_ids[2*n]
              <<", "<<node_deaths_global_ids[2*n+1]<<")\n";
      n++;
    }
  }
  for( n=0, i=0 ; i<number_face_blocks ; i++ ) {
    postream<<"killing "<<num_face_deaths_per_block[i]<<" faces in block "<<i<<"\n";
    for (j=0; j<num_face_deaths_per_block[i]; j++) {
      postream<<"  face ("<<face_deaths_global_ids[2*n]
              <<", "<<face_deaths_global_ids[2*n+1]<<")\n";
      n++;
    }
  }
  for( n=0, i=0 ; i<number_elem_blocks ; i++ ) {
    postream<<"killing "<<num_elem_deaths_per_block[i]<<" elems in block "<<i<<"\n";
    for (j=0; j<num_elem_deaths_per_block[i]; j++) {
      postream<<"  elem ("<<elem_deaths_global_ids[2*n]
              <<", "<<elem_deaths_global_ids[2*n+1]<<")\n";
      n++;
    }
  }
  postream.flush();
  //DisplayByList();
#endif
  ComputeNodalCommPlan();
}

void
DriverTopology::GetModsForDLB( int& num_node_exports,
                       	     int* node_export_gids,
                             int* node_export_pids,
		       	     int& num_face_exports,
                       	     int* face_export_gids,
                             int* face_export_pids,
		       	     int& num_elem_exports,
                       	     int* elem_export_gids,
                             int* elem_export_pids )
{
#ifndef CONTACT_NO_MPI
  int i;
  
  int my_proc_id  = contact_processor_number(AcmeComm);
  int total_procs = contact_number_of_processors(AcmeComm);
  
  if (total_procs==1) return;
  
  AcmeNode** nodes = 
    REINTERPRET_CAST(AcmeNode**)(node_list->EntityList());
  AcmeFace** faces = 
    REINTERPRET_CAST(AcmeFace**)(face_list->EntityList());
  AcmeElem** elems = 
    REINTERPRET_CAST(AcmeElem**)(elem_list->EntityList());
  
  zoltan->Set_GeomCallBacks(this);
  zoltan->Balance();

  LB_ID_TYPE zoltan_lid[ACME_ZOLTAN_LID_SIZE];
  LB_ID_TYPE zoltan_gid[ACME_ZOLTAN_GID_SIZE];
  int        zoltan_pid;

  for (i=0; i<number_node_blocks; i++) {
    int n     = node_block[i]->Number_of_Entities();
    int sum   = contact_global_sum(n, AcmeComm );
    int guess = (int)(1.5*(double)sum/(double)total_procs)+1;
    if (n == 0) node_block[i]->EntityList()->SetupHash(guess);
  }
  for (i=0; i<number_face_blocks; i++) {
    int n     = face_block[i]->Number_of_Entities();
    int sum   = contact_global_sum(n, AcmeComm );
    int guess = (int)(1.5*(double)sum/(double)total_procs)+1;
    if (n == 0) face_block[i]->EntityList()->SetupHash(guess);
  }
  for (i=0; i<number_elem_blocks; i++) {
    int n     = elem_block[i]->Number_of_Entities();
    int sum   = contact_global_sum(n, AcmeComm );
    int guess = (int)(1.5*(double)sum/(double)total_procs)+1;
    if (n == 0) elem_block[i]->EntityList()->SetupHash(guess);
  }
  int zoltancomm_size = 12*number_of_nodes +
                        12*number_of_faces + 
                        12*number_of_elems;
  AcmeZoltanComm Toplevel_ZoltanComm( AcmeZoltanComm::ZOLTAN_EXPORT, 
                                      zoltancomm_size );
                                         
  Assign_New_Ownership();
  
  //=========================================================================
  //  F A C E S
  //=========================================================================
  for (i=0; i<number_of_faces; i++) {
    AcmeFace* face = faces[i];
    if (face->temp_tag != my_proc_id) {
      // off processor - add to communication object
      zoltan_pid = face->temp_tag;
      face->ZoltanLID(zoltan_lid);
      face->ZoltanGID(zoltan_gid);
      Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
      // add the nodes too, duplicates will
      // be culled out as they're added
      for(int j=0 ; j<face->Nodes_Per_Face() ; j++ ){
        face->Node(j)->ZoltanLID(zoltan_lid);
        face->Node(j)->ZoltanGID(zoltan_gid);
        Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
      }
    }
  }
  
  //=========================================================================
  //  E L E M E N T S
  //=========================================================================
  for (i=0; i<number_of_elems; i++) {
    AcmeElem* elem = elems[i];
    if (elem->temp_tag != my_proc_id) {
      // off processor - add to communication object
      zoltan_pid = elem->temp_tag;
      elem->ZoltanLID(zoltan_lid);
      elem->ZoltanGID(zoltan_gid);
      Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
      // add the nodes too, duplicates will
      // be culled out as they're added
      for(int j=0 ; j<elem->Nodes_Per_Element() ; j++ ){
        elem->Node(j)->ZoltanLID(zoltan_lid);
        elem->Node(j)->ZoltanGID(zoltan_gid);
        Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
      }
    }
  }
  // In this block, we are only moving "POINTS" not "NODES"
  int base = 0;
  if (number_face_blocks+number_elem_blocks > 0) base=1;
  for(i=base; i<number_node_blocks; i++) {
    if( node_block[i]->DerivedType() == AcmeEntity::ACME_POINT ){
      int nnodes = node_list->BlockNumEntities(i);
      AcmeNode** nodes  = 
        REINTERPRET_CAST(AcmeNode**)(node_list->BlockEntityList(i));
      for (int j=0; j<nnodes; j++) {
        AcmeNode* node = nodes[j];
        if (node->temp_tag != my_proc_id) {
          // off processor - add to communication object
          zoltan_pid = node->temp_tag;
          node->ZoltanLID(zoltan_lid);
          node->ZoltanGID(zoltan_gid);
          Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
        }
      }
    }
  }
  
  int       num_toplevel_import   = -1;
  LB_ID_PTR import_toplevel_gids  = NULL;
  LB_ID_PTR import_toplevel_lids  = NULL;
  int*      import_toplevel_procs = NULL;
  int       num_toplevel_export   = Toplevel_ZoltanComm.Num_Export();
  LB_ID_PTR export_toplevel_gids  = Toplevel_ZoltanComm.Export_GIDS();
  LB_ID_PTR export_toplevel_lids  = Toplevel_ZoltanComm.Export_LIDS();
  int*      export_toplevel_procs = Toplevel_ZoltanComm.Export_Procs();

  //================================================
  // Migrate all the off processor toplevel objects
  //================================================
  zoltan->Set_DLBmultiCallBacks(this);
  zoltan->Help_Migrate(num_toplevel_import,  import_toplevel_gids, 
                       import_toplevel_lids, import_toplevel_procs, 
                       num_toplevel_export,  export_toplevel_gids,
		       export_toplevel_lids, export_toplevel_procs);

#ifdef DEBUG
  postream<<"After DLB Migration...\n";
  for( i=0 ; i<number_face_blocks ; i++ ) {
    postream<<"Block "<<i<<" has "<<face_block[i]->EntityList()->NumEntities()<<" faces\n";
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    int m = 0;
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      postream << "  face "<<m++<<"\n";
      postream << "    exoid  = "<<entity->Exodus_ID()<<"\n";
      postream << "    owner  = "<<entity->Owner()<<"\n";
      postream << "    active = "<<entity->State()<<"\n";
      postream << "    hindex = "<<entity->HostArrayIndex()<<"\n";
      postream << "    pindex = "<<entity->ProcArrayIndex()<<"\n";
      postream << "    tag    = "<<entity->temp_tag<<"\n";
    }
  }
  if (export_toplevel_gids!=NULL && 
      export_toplevel_lids!=NULL && 
      export_toplevel_procs!=NULL) {
    postream << "DLB EXPORTS...\n";
    for (i=0; i< num_toplevel_export; i++){
      postream << "  " << i << ":  GID = ("
               << AcmeZoltanGID::Type(&export_toplevel_gids[i*ACME_ZOLTAN_GID_SIZE]) << ", "
               << AcmeZoltanGID::ID  (&export_toplevel_gids[i*ACME_ZOLTAN_GID_SIZE]) << ")"
               << "    LID = ("
               << AcmeZoltanLID::Type (&export_toplevel_lids[i*ACME_ZOLTAN_LID_SIZE]) << ", "
               << AcmeZoltanLID::Index(&export_toplevel_lids[i*ACME_ZOLTAN_LID_SIZE]) << ")"
               << "    PROC = "
               << export_toplevel_procs[i]<< "\n";
    }
  }
  if (import_toplevel_gids!=NULL && 
      import_toplevel_lids!=NULL && 
      import_toplevel_procs!=NULL) {
    postream << "DLB IMPORTS...\n";
    for (i=0; i< num_toplevel_import; i++){
      postream << "  " << i << ":  GID = ("
               << AcmeZoltanGID::Type(&import_toplevel_gids[i*ACME_ZOLTAN_GID_SIZE]) << ", "
               << AcmeZoltanGID::ID  (&import_toplevel_gids[i*ACME_ZOLTAN_GID_SIZE]) << ")"
               << "    LID = ("
               << AcmeZoltanLID::Type (&import_toplevel_lids[i*ACME_ZOLTAN_LID_SIZE]) << ", "
               << AcmeZoltanLID::Index(&import_toplevel_lids[i*ACME_ZOLTAN_LID_SIZE]) << ")"
               << "    PROC = "
               << import_toplevel_procs[i]<< "\n";
    }
  }
  postream.flush();
#endif
  
  postream<<"Deleting exported faces/elements...\n";
  num_node_exports = 0;
  num_face_exports = 0;
  num_elem_exports = 0;
  for (i=0; i<num_toplevel_export; i++) {
    ZOLTAN_ID_PTR lid = &export_toplevel_lids[i*ACME_ZOLTAN_LID_SIZE];
    ZOLTAN_ID_PTR gid = &export_toplevel_gids[i*ACME_ZOLTAN_GID_SIZE];
    int entity_type  = AcmeZoltanGID::Type(gid);
    int entity_id    = AcmeZoltanGID::ID(gid); 
    int entity_index = AcmeZoltanLID::Index(lid); 
    switch (entity_type) {
    case AcmeEntity::CT_NODE: {
      AcmeNode* node = nodes[entity_index];
      if (node->State()==AcmeEntity::ACTIVE) {
        node_export_gids[2*num_node_exports  ] = 0;
        node_export_gids[2*num_node_exports+1] = entity_id;
        node_export_pids[  num_node_exports  ] = export_toplevel_procs[i];
        num_node_exports++;
      }
      } break;
    case AcmeEntity::CT_FACE: {
      AcmeFace* face = faces[entity_index];
      if (face->State()==AcmeEntity::ACTIVE) {
        face_export_gids[2*num_face_exports  ] = 0;
        face_export_gids[2*num_face_exports+1] = entity_id;
        face_export_pids[  num_face_exports  ] = export_toplevel_procs[i];
        num_face_exports++;
      }
      postream<<"Deleting face id "<<face->Exodus_ID()<<" from block "<<face->BlockID()<<"\n";
      face_block[face->BlockID()]->EntityList()->Delete(face);
      } break;
    case AcmeEntity::CT_ELEM: {
      AcmeElem* elem = elems[entity_index];
      if (elem->State()==AcmeEntity::ACTIVE) {
        elem_export_gids[2*num_elem_exports  ] = 0;
        elem_export_gids[2*num_elem_exports+1] = entity_id;
        elem_export_pids[  num_elem_exports  ] = export_toplevel_procs[i];
        num_elem_exports++;
      }
      elem_block[elem->BlockID()]->EntityList()->Delete(elem);
      } break;
    }
  } 
  for( i=0 ; i<number_face_blocks ; i++ ) {
    postream<<"Block "<<i<<" has "<<face_block[i]->EntityList()->NumEntities()<<" faces\n";
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    int m = 0;
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      postream << "  face "<<m++<<"\n";
      postream << "    exoid  = "<<entity->Exodus_ID()<<"\n";
      postream << "    owner  = "<<entity->Owner()<<"\n";
      postream << "    active = "<<entity->State()<<"\n";
      postream << "    hindex = "<<entity->HostArrayIndex()<<"\n";
      postream << "    pindex = "<<entity->ProcArrayIndex()<<"\n";
      postream << "    tag    = "<<entity->temp_tag<<"\n";
    }
  }
  postream.flush();
  
  node_list->BuildList(node_block, number_node_blocks);
  
  // need to delete any hanging nodes
  AcmeBlockEntityList* entity_list = node_block[0]->EntityList();
  entity_list->IteratorStart();
  while (AcmeEntity* entity = entity_list->IteratorForward()) {
    entity->temp_tag = 0;
  }
  for( i=1 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      entity->temp_tag = 1;
    }
  }
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      face->ConnectNodes(node_list);
      for (int n=0; n<face->Nodes_Per_Face(); n++) {
        face->Node(n)->temp_tag = 1;
      }
    }
  }
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      elem->ConnectNodes(node_list);
      for (int n=0; n<elem->Nodes_Per_Element(); n++) {
        elem->Node(n)->temp_tag = 1;
      }
    }
  }
  int Knt=0;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->temp_tag==0) entity_list->Delete(entity);
    }
    Knt += node_block[i]->EntityList()->NumEntities();
  }
  
  entity_list = node_block[0]->EntityList();
  entity_list->IteratorStart();
  while (AcmeEntity* entity = entity_list->IteratorForward()) {
    entity->State(AcmeEntity::INACTIVE);
  }
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      if (face->State()==AcmeEntity::ACTIVE) {
        for (int n=0; n<face->Nodes_Per_Face(); n++) {
          face->Node(n)->State(AcmeEntity::ACTIVE);
        }
      }
    }
  }
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      if (elem->State()==AcmeEntity::ACTIVE) {
        for (int n=0; n<elem->Nodes_Per_Element(); n++) {
          elem->Node(n)->State(AcmeEntity::ACTIVE);
        }
      }
    }
  }
  
  node_list->BuildList(node_block, number_node_blocks);
  POSTCONDITION(number_node_blocks==node_list->Num_Blocks());
  POSTCONDITION(Knt==node_list->NumEntities());
  face_list->BuildList(face_block, number_face_blocks);
  POSTCONDITION(number_face_blocks==face_list->Num_Blocks());
  elem_list->BuildList(elem_block, number_elem_blocks);
  POSTCONDITION(number_elem_blocks==elem_list->Num_Blocks());
  
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      face->ConnectNodes(node_list);
    }
  }
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      elem->ConnectNodes(node_list);
    }
  }

  // Determine the total number of nodes
  number_of_nodes = 0;
  number_of_active_nodes = 0;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    number_of_nodes += node_block[i]->EntityList()->NumEntities();
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        entity->HostArrayIndex(number_of_active_nodes++);
      } else {
        entity->HostArrayIndex(-1);
      }
    }
  }
  POSTCONDITION(number_of_nodes==node_list->NumEntities());
  // Determine the total number of faces 
  number_of_faces = 0; 
  number_of_active_faces = 0;
  for( i=0 ; i<number_face_blocks ; i++ ) {
    number_of_faces += face_block[i]->EntityList()->NumEntities();
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        entity->HostArrayIndex(number_of_active_faces++);
      } else {
        entity->HostArrayIndex(-1);
      }
    }
  }
  POSTCONDITION(number_of_faces==face_list->NumEntities());
  // Determine the total number of elements
  number_of_elems = 0; 
  number_of_active_elems = 0;
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    number_of_elems += elem_block[i]->EntityList()->NumEntities();
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        entity->HostArrayIndex(number_of_active_elems++);
      } else {
        entity->HostArrayIndex(-1);
      }
    }
  }
  POSTCONDITION(number_of_elems==elem_list->NumEntities());
  
  ComputeNodalCommPlan();
  //ComputeActiveNodalCommPlan();
  Set_Ownership();
   
#ifdef DEBUG
  postream<<"After DLB...\n";
  postream<<"  number_of_nodes        = "<<number_of_nodes<<"\n";
  postream<<"  number_of_faces        = "<<number_of_faces<<"\n";
  postream<<"  number_of_elems        = "<<number_of_elems<<"\n";
  postream<<"  number_of_active_nodes = "<<number_of_active_nodes<<"\n";
  postream<<"  number_of_active_faces = "<<number_of_active_faces<<"\n";
  postream<<"  number_of_active_elems = "<<number_of_active_elems<<"\n";
  postream<<"  num_node_exports       = "<<num_node_exports<<"\n";
  postream<<"  num_face_exports       = "<<num_face_exports<<"\n";
  postream<<"  num_elem_exports       = "<<num_elem_exports<<"\n";
  //DisplayByList();
  DisplayCommPlan(NumCommProcs,CommProcIds,
                  NumNodesToProc,CommNodes);
#endif
#endif
}

#ifndef CONTACT_NO_MPI
void DriverTopology::Assign_New_Ownership()
{
  int i;
  
  AcmeNode** Nodes = 
    REINTERPRET_CAST(AcmeNode**)(node_list->EntityList());
  for (i=0; i<number_of_nodes; i++) {
    int proc_num;
    Real *position = Nodes[i]->Coordinates();
    zoltan->Point_Assign(position,&proc_num);
    Nodes[i]->temp_tag = proc_num;
    postream<<"node "<<i<<" ("<<Nodes[i]->State()<<") assigned to proc "<<proc_num<<"\n";
  }

  AcmeFace** Faces = 
    REINTERPRET_CAST(AcmeFace**)(face_list->EntityList());
  for (i=0; i<number_of_faces; i++) {
    AcmeFace* face = Faces[i];
    postream<<"face "<<i<<" ("<<face->State()<<"):\n";
    face->temp_tag = -1;
    for(int j=0 ; j<face->Nodes_Per_Face() ; j++ ){
      postream<<"  node "<<j<<" ("<<face->Node(j)->Ownership()<<") assigned to proc "<<face->Node(j)->temp_tag<<"\n";
      if( face->Node(j)->Ownership() == AcmeEntity::OWNED ){
        if (face->temp_tag<0) {
          face->temp_tag = face->Node(j)->temp_tag;
        } else {
          face->temp_tag = MIN(face->temp_tag,face->Node(j)->temp_tag);
        }
      }
    }
    postream<<"  face assigned to proc "<<face->temp_tag<<"\n";
  }

  AcmeElem** Elems = 
    REINTERPRET_CAST(AcmeElem**)(elem_list->EntityList());
  for (i=0; i<number_of_elems; i++) {
    AcmeElem* elem = Elems[i];
    postream<<"elem "<<i<<" ("<<elem->State()<<"):\n";
    elem->temp_tag = -1;
    for(int j=0 ; j<elem->Nodes_Per_Element() ; j++ ){
      postream<<"  node "<<j<<" ("<<elem->Node(j)->Ownership()<<") assigned to proc "<<elem->Node(j)->temp_tag<<"\n";
      if( elem->Node(j)->Ownership() == AcmeEntity::OWNED ){
        if (elem->temp_tag<0) {
          elem->temp_tag = elem->Node(j)->temp_tag;
        } else {
          elem->temp_tag = MIN(elem->temp_tag,elem->Node(j)->temp_tag);
        }
      }
    }
    postream<<"  elem assigned to proc "<<elem->temp_tag<<"\n";
  }
  postream.flush();			       
}
#endif

void
DriverTopology::GetCommPlan(int& num_comm_procs, 
                          int* comm_proc_ids,
		          int* num_nodes_to_proc, 
                          int* comm_nodes,
                          int  my_proc_id,
		          int  num_procs)
{
#ifndef CONTACT_NO_MPI

  if (num_procs==1) return;
  num_comm_procs = NumCommProcs;
  for ( int ii = 0; ii < MAX_PROCS; ii++) {
    comm_proc_ids[ii] = CommProcIds[ii] ;
    num_nodes_to_proc[ii] = NumNodesToProc[ii];
    for ( int jj = 0; jj < MAX_NODES; jj++) {
      comm_nodes[ii*MAX_NODES+jj] = CommNodes[ii*MAX_NODES+jj];
    }
  }
#ifdef DEBUG
  postream << ">>> Here is the nodal comm plan\n";
  int cnt = 0;
  for ( int i = 0;  i < num_comm_procs; i++ ) {
    postream << my_proc_id << ":  -> partner " << i+1 << " is P" 
             << comm_proc_ids[i] << " with " 
             << num_nodes_to_proc[i] << " nodes.\n";
    for (int j = 0; j < num_nodes_to_proc[i]; j++ ) {
      postream << my_proc_id << ":   -> node " << j << ", id " 
               << node_list->Entity(comm_nodes[cnt]-1)->Exodus_ID() << "\n";
      cnt++;
    }
  }
  postream.flush();
#endif
#endif
}

void
DriverTopology::ComputeNodalCommPlan()
{
#ifndef CONTACT_NO_MPI
  int i, j;
  int total_procs = contact_number_of_processors(AcmeComm);
  int my_proc_id  = contact_processor_number(AcmeComm);
  if (total_procs==1) return;
#ifdef DEBUG_IT
  postream << ">>> Here is the original nodal comm plan\n";
  int cnt = 0;
  for ( i = 0;  i < NumCommProcs; i++ ) {
    postream << my_proc_id << ":  -> partner " << i+1 << " is P" 
             << CommProcIds[i] << " with " << NumNodesToProc[i]
             << " nodes.\n";
    for (j = 0; j < NumNodesToProc[i]; j++ ) {
      postream << my_proc_id << ":   -> node " << j << ", id " 
               << node_list->Entity(CommNodes[cnt]-1)->Exodus_ID() << "\n";
      cnt++;
    }
  }
  postream.flush();
#endif
  NumCommProcs = 0;
  for (i=0; i<MAX_PROCS; i++) {
    CommProcIds[i]    = -1;
    NumNodesToProc[i] = 0;
  }
  //===================================================
  // Need to compute the new node communication plan.
  // Just brute force it for now.  Not elegant but for
  // these small regression tests, should be sufficiant.  
  //===================================================
  int max_nodes=0;
  contact_global_maximum(&number_of_active_nodes, &max_nodes, 1, AcmeComm );
  if (max_nodes>0) {
    int mesg=1009;
    int* tmp_nodes = new int[max_nodes];
    RequestHandle rh;
    AcmeEntity** entities = node_list->EntityList();
    for (int master_proc=0; master_proc<total_procs; master_proc++) {
      int nnodes = number_of_active_nodes;
      contact_broadcast( &nnodes, 1, master_proc, AcmeComm );
      if (nnodes>0) {
        if (master_proc==my_proc_id) {
          for( j=0, i=0 ; i<number_of_nodes ; i++ ) {
            if (entities[i]->State()==AcmeEntity::ACTIVE) {
              tmp_nodes[j++] = entities[i]->Exodus_ID();
            }
          }
        }
        contact_broadcast( tmp_nodes, nnodes, master_proc, AcmeComm );
        if (my_proc_id!=master_proc) {
          for (i=0; i<nnodes; i++) {
            AcmeEntity* entity = node_list->Find(tmp_nodes[i]);
            if (entity==NULL) {
              tmp_nodes[i] = -1;
            } else {
              if (entity->State()==AcmeEntity::INACTIVE) tmp_nodes[i] = -1;
            }
          }
          contact_blocking_send( mesg, tmp_nodes, nnodes, 
                                 master_proc, AcmeComm );
        } else {
          //contact_global_sync(AcmeComm);
          int* nptr = CommNodes;
          for (int slave_proc=0; slave_proc<total_procs; slave_proc++) {
            if (slave_proc==master_proc) continue;
            rh = contact_nonblocking_receive( mesg, tmp_nodes, nnodes, 
                                              slave_proc, AcmeComm );
            contact_wait_msg_done( rh );
            for (i=0; i<nnodes; i++) {
              if (tmp_nodes[i]>=0) {
                AcmeEntity* entity = node_list->Find(tmp_nodes[i]);
                if (CommProcIds[NumCommProcs]<0) {
                  CommProcIds[NumCommProcs] = slave_proc;
                }
                NumNodesToProc[NumCommProcs]++;
                *nptr++ = entity->HostArrayIndex()+1;
              }
            }
            if (CommProcIds[NumCommProcs]>=0) NumCommProcs++;
          }
        }
      }
    }
    int knt = 0;
    for (int i = 0;  i < NumCommProcs; i++ ) {
      for (int j = 0; j < NumNodesToProc[i]; j++ ) {
        int k;
        for (k=0; k<number_of_nodes; k++) {
          if (entities[k]->HostArrayIndex()==CommNodes[knt+j]-1) break;
        }
        tmp_nodes[j] = node_list->Entity(k)->Exodus_ID();
      }
      SortCommNodes(NumNodesToProc[i],&CommNodes[knt],tmp_nodes);
      knt += NumNodesToProc[i];
    }
    delete [] tmp_nodes;
  }
#ifdef DEBUG
  postream << ">>> Here is the updated nodal comm plan\n";
  int cnt = 0;
  for ( i = 0;  i < NumCommProcs; i++ ) {
    postream << my_proc_id << ":  -> partner " << i+1 << " is P" 
             << CommProcIds[i] << " with " << NumNodesToProc[i]
             << " nodes.\n";
    for (j = 0; j < NumNodesToProc[i]; j++ ) {
      int k;
      for (k=0; k<number_of_nodes; k++) {
        if (entities[k]->HostArrayIndex()==CommNodes[cnt]-1) break;
      }
      postream << my_proc_id << ":   -> node " << j << ", id " 
               << node_list->Entity(k)->Exodus_ID() << "\n";
      cnt++;
    }
  }
  postream.flush();
#endif
#endif
}

void
DriverTopology::ComputeActiveNodalCommPlan()
{
#ifndef CONTACT_NO_MPI
  int total_procs = contact_number_of_processors(AcmeComm);
  int my_proc_id  = contact_processor_number(AcmeComm);
  if (total_procs==1) return;
  int i, j, len, mesg=1008;

  // count total number of nodes to communicate
  int total_nodes_to_comm = 0; 
  for (i = 0; i < NumCommProcs; i++) {
    total_nodes_to_comm += NumNodesToProc[i];
  }

  // allocate buffers for communication
  int * send_buf = new int[MAX(total_nodes_to_comm,1)];
  int * recv_buf = new int[MAX(total_nodes_to_comm,1)];
  RequestHandle * recv_handles = new RequestHandle[MAX(NumCommProcs,1)];
  RequestHandle * rh = recv_handles;

  // Gather my data into send_buf
  int count = 0;
  for( i=0 ; i<NumCommProcs ; i++ ){
    for( j=0 ; j<NumNodesToProc[i] ; j++ ){
      int node_num = CommNodes[count]-1;
      POSTCONDITION ( node_num >= 0);
      send_buf[count] = node_list->Entity(node_num)->State();
      count++;
      POSTCONDITION(count <= total_nodes_to_comm);
    }
  }
  
  // Post Receives
  count = 0;
  for( i=0 ; i<NumCommProcs ; i++ ){
    len = NumNodesToProc[i];
    rh[i] = contact_nonblocking_receive( mesg, &recv_buf[count], len, 
                                         CommProcIds[i], AcmeComm );
    count += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( AcmeComm );

  if( NumCommProcs > 0 ) {
    // Send my data
    count = 0;
    for( i=0 ; i<NumCommProcs ; i++ ){
      len = NumNodesToProc[i];
      contact_blocking_send( mesg, &send_buf[count], len, 
                             CommProcIds[i], AcmeComm );
      count += len;
    }
    // Wait till all of my messages have arrived
    for( i=0 ; i<NumCommProcs ; i++ ){
      contact_wait_msg_done( rh[i] );
    }
  }
  
  ActiveNumCommProcs = 0;
  int old_comm_count = 0;
  int new_comm_count = 0;
  for (i = 0; i < NumCommProcs; i++) {
    int num_for_partner = 0;
    for (j = 0; j < NumNodesToProc[i]; j++) {
      int node = CommNodes[old_comm_count] - 1;//-1: fortran->c
      if (node_list->Entity(node)->State() == AcmeEntity::ACTIVE) {
        if (recv_buf[old_comm_count] > 0 ){
          ActiveCommNodes[new_comm_count] = node_list->Entity(node)->HostArrayIndex()+1;//+1: c->fortran
          new_comm_count++;
          num_for_partner++;
        }
      }
      old_comm_count ++;
    }
    if (num_for_partner > 0) { 
      ActiveCommProcIds[ActiveNumCommProcs]    = CommProcIds[i] ;
      ActiveNumNodesToProc[ActiveNumCommProcs] = num_for_partner;
      ActiveNumCommProcs++;
    }
  }

  // cleanup
  delete [] send_buf;
  delete [] recv_buf;
  delete [] recv_handles;
  
#ifdef DEBUG
  // ---- Code used for debugging the driver ---
  postream << ">>> Here is original nodal comm plan\n";
  int cnt = 0;
  for ( i = 0;  i < NumCommProcs; i++ ) {
    postream << my_proc_id << ":  -> partner " << i+1 << " is P" 
             << CommProcIds[i] << " with " << NumNodesToProc[i]
             << " nodes.\n";
    for (j = 0; j < NumNodesToProc[i]; j++ ) {
      postream << my_proc_id << ":   -> node " << j << ", id " 
               << node_list->Entity(CommNodes[cnt]-1)->Exodus_ID() << "\n";
      cnt++;
    }
  }
  postream << ">>> Here is active nodal comm plan\n";
  cnt = 0;
  for ( i = 0;  i < ActiveNumCommProcs; i++ ) {
    postream << my_proc_id << ":  -> partner " << i+1 << " is P" 
             << ActiveCommProcIds[i] << " with " << ActiveNumNodesToProc[i]
             << " nodes.\n";
    for (j = 0; j < ActiveNumNodesToProc[i]; j++ ) {
      postream << my_proc_id << ":   -> node " << j << ", id " 
               << node_list->Entity(ActiveCommNodes[cnt]-1)->Exodus_ID() << "\n";
      cnt++;
    }
  }
  postream.flush();
#endif

#endif
}

void
DriverTopology::GetHostIDs(int& num_nodes,
                         int* node_exo_ids,
                         int* node_host_ids,
                         int& num_faces,
                         int* face_host_ids,
                         int& num_elems,
                         int* elem_host_ids)
{
  int i;
  
  num_nodes = 0;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State() == AcmeEntity::ACTIVE) {
	node_exo_ids[num_nodes]      = entity->Exodus_ID();
        node_host_ids[2*num_nodes  ] = 0;
        node_host_ids[2*num_nodes+1] = entity->Exodus_ID();
        num_nodes++;
      }
    }
  }
  
  num_faces = 0;
  for( i=0 ; i<number_face_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        face_host_ids[2*num_faces  ] = 0;
        face_host_ids[2*num_faces+1] = entity->Exodus_ID();
        num_faces++;
      }
    }
  }
  
  num_elems = 0;
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        elem_host_ids[2*num_elems  ] = 0;
        elem_host_ids[2*num_elems+1] = entity->Exodus_ID();
        num_elems++;
      }
    }
  }
}

void
DriverTopology::GetVariables(double* position,
		           double* nodal_vars, 
                           int&    num_nodal_vars,
		           double* element_vars,
		           int&    num_element_vars,
		           double* attributes,
		           int*    num_attributes,
		           double* node_radius,
		           double* shell_thickness,
		           double* shell_offset)
{
  int i, k, n, nn;
  
  for( i=0 ; i<number_of_active_faces ; i++ ) {
    shell_offset[i]    = -1.0;
    shell_thickness[i] =  0.0;
  }
  for( i=0 ; i<number_of_active_nodes ; i++ ) {
    node_radius[i] = 0.0;
  }
  
  num_nodal_vars = num_node_vars;
  double* nvars  = nodal_vars;
  double* coords = position;
  for( n=0, i=0 ; i<number_node_blocks ; i++ ) {
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        AcmeNode* node = static_cast<AcmeNode*>(entity);
        for (k=0; k<dimension; k++) {
          coords[k] = node->Coordinates()[k];
        }
        for (k=0; k<num_node_vars; k++) {
          nvars[number_of_active_nodes*k+n] = node->Data(k);
        }
        coords += dimension;
        n++;
      }
    }
  }
  
  num_element_vars = num_elem_vars;
  double* evars    = element_vars;
  double* attr     = attributes;
  
  for( n=0, i=0 ; i<number_face_blocks ; i++ ) {
    nn = 0;
    num_attributes[i] = face_block[i]->NumAttributes();
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        AcmeFace* face = static_cast<AcmeFace*>(entity);
        for (k=0; k<num_elem_vars; k++) {
          evars[(number_of_active_faces+number_of_active_elems)*k+n] = 
            face->Data(k);
        }
        for (k=0; k<face->NumAttributes(); k++) {
          attr[face_block[i]->Number_of_Entities()*k+nn] = 
            face->Attributes(k);
        }
        if( face_block[i]->FaceType() == 5 || face_block[i]->FaceType() == 6 ){
          if (face->Exodus_ID()%2) {
            shell_offset[n] = 1.0-face_block[i]->ShellOffset();
          } else {
            shell_offset[n] = face_block[i]->ShellOffset();
          }
        }
	else {
	  shell_offset[n] = -1.0;
	}
        n++; nn++;
      }
    }
    attr += face_block[i]->NumAttributes()*face_block[i]->Number_of_Entities();
  }
  
  for( n=0, i=0 ; i<number_elem_blocks ; i++ ) {
    num_attributes[i] = elem_block[i]->NumAttributes();
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) {
        AcmeElem* elem = static_cast<AcmeElem*>(entity);
        for (k=0; k<num_elem_vars; k++) {
          evars[(number_of_active_faces+number_of_active_elems)*k+n] = 
            elem->Data(k);
        }
        for (k=0; k<elem->NumAttributes(); k++) {
          attr[elem_block[i]->Number_of_Entities()*k+nn] = 
            elem->Attributes(k);
        }
        n++;
      }
    }
    attr += num_attributes[i]*elem_block[i]->Number_of_Entities();
  }
}

#ifndef CONTACT_NO_MPI
void DriverTopology::create_zoltan_object()
{
  //========================
  // Create a zoltan object 
  //========================
  int zoltan_error;
  zoltan = new AcmeZoltan(AcmeComm, zoltan_error);
  if( zoltan_error == LB_FATAL ){
    cerr<<"Error Creating Zoltan Object";
    exit(1);
  }
  zoltan->Set_Method( (char*)"RCB" );
  zoltan->Set_Param( (char*)"RCB_REUSE",       (char*)"1" );
  // Check_Geom causes Zoltan to do unnecessary checking, so turn it off
  // For us this is
  //   1) check dot weights > 0 (we don't specify them)
  //   2) Make sure #dots_in = #dots_out
  //   3) Check that the dots are loadbalanced within a tolerance
  //   4) Check all points are within the RCB boxes
  zoltan->Set_Param( (char*)"CHECK_GEOM",      (char*)"0" );
  zoltan->Set_Param( (char*)"KEEP_CUTS",       (char*)"1" );
  zoltan->Set_Param( (char*)"DEBUG_LEVEL",     (char*)"0" );
#if CONTACT_DEBUG_PRINT_LEVEL>=10
  zoltan->Set_Param( (char*)"RETURN_LISTS",    (char*)"ALL");
#else
  zoltan->Set_Param( (char*)"RETURN_LISTS",    (char*)"NONE");
#endif
  char num_lid_entries[2];
  char num_gid_entries[2];
  sprintf(num_lid_entries,"%d",ACME_ZOLTAN_LID_SIZE);
  sprintf(num_gid_entries,"%d",ACME_ZOLTAN_GID_SIZE);
  zoltan->Set_Param( (char*)"NUM_LID_ENTRIES", num_lid_entries );
  zoltan->Set_Param( (char*)"NUM_GID_ENTRIES", num_gid_entries );
#ifdef __PUMAGON__
  zoltan->Set_Param( (char*)"TFLOPS_SPECIAL",  (char*)"1" );
#endif
}
#endif

void DriverTopology::Set_Ownership()
{
  int i;
  int my_proc_id = contact_processor_number(AcmeComm);
  int num_procs = contact_number_of_processors(AcmeComm);
  
  AcmeNode** Nodes = 
    REINTERPRET_CAST(AcmeNode**)(node_list->EntityList());
  for (i=0; i<number_of_nodes; i++) {
    AcmeNode* node = Nodes[i];
    node->Ownership(AcmeEntity::OWNED);
    node->Owner(my_proc_id);
  }

  AcmeFace** Faces = 
    REINTERPRET_CAST(AcmeFace**)(face_list->EntityList());
  for (i=0; i<number_of_faces; i++) {
    AcmeFace* face = Faces[i];
    face->Ownership(AcmeEntity::OWNED);
    face->Owner(my_proc_id);
  }

  AcmeElem** Elems = 
    REINTERPRET_CAST(AcmeElem**)(elem_list->EntityList());
  for (i=0; i<number_of_elems; i++) {
    AcmeElem* elem = Elems[i];
    elem->Ownership(AcmeEntity::OWNED);
    elem->Owner(my_proc_id);
  }
  
#ifndef CONTACT_NO_MPI
  if (num_procs>1) {
    int j;
    
    //-------------------------------------------
    // compute the owner of all entities in the node list
    
    // first loop over all the nodes and set temp tag to
    // the local processor id
    for (i=0; i<number_of_nodes; i++) {
      Nodes[i]->temp_tag = my_proc_id;
    }
    
    // loop over the communication list and set the 
    // owner to the lowest numbered processor
    int cnt = 0;
    for (i=0; i<NumCommProcs; i++) {
      int proc_id = CommProcIds[i];
      if (proc_id<my_proc_id) {
        for ( j = 0; j < NumNodesToProc[i]; j++ ) {
          if (Nodes[CommNodes[cnt+j]-1]->temp_tag>proc_id) {
            Nodes[CommNodes[cnt+j]-1]->temp_tag = proc_id;
          }
        }
      }
      cnt += NumNodesToProc[i];
    }
    
    // loop over all entities in the entity list and set the owner flag
    for (i=0; i<number_of_nodes; i++) {
      AcmeNode* node = Nodes[i];
      if ( my_proc_id != node->temp_tag ){
        node->Ownership(AcmeEntity::NOT_OWNED);
        node->Owner(node->temp_tag);
      }
    }
  }
#endif
}

void DriverTopology::DisplayByList()
{
  int i, k;
  
  postream<<"num nodes (total/active) = "<<number_of_nodes<<" / "<<number_of_active_nodes<<"\n";
  postream<<"num faces (total/active) = "<<number_of_faces<<" / "<<number_of_active_faces<<"\n";
  postream<<"num elems (total/active) = "<<number_of_elems<<" / "<<number_of_active_elems<<"\n";
  postream<<"num node blocks  = "<<number_node_blocks<<"\n";
  for( i=0 ; i<number_node_blocks ; i++ ) {
    int ntotal = node_block[i]->EntityList()->NumEntities();
    int nactive = 0;
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) nactive++;
    }
    postream<<"  total/active = "<<ntotal<<" / "<<nactive<<"\n";
  }
  postream<<"num face blocks  = "<<number_face_blocks<<"\n";
  for( i=0 ; i<number_face_blocks ; i++ ) {
    int ntotal = face_block[i]->EntityList()->NumEntities();
    int nactive = 0;
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) nactive++;
    }
    postream<<"  total/active = "<<ntotal<<" / "<<nactive<<"\n";
  }
  postream<<"num elem blocks  = "<<number_elem_blocks<<"\n";
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    int ntotal = elem_block[i]->EntityList()->NumEntities();
    int nactive = 0;
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) nactive++;
    }
    postream<<"  total/active = "<<ntotal<<" / "<<nactive<<"\n";
  }
  postream<<"\n";
  
  AcmeNode** Nodes = 
    REINTERPRET_CAST(AcmeNode**)(node_list->EntityList());
  for (i=0; i<number_of_nodes; i++) {
    AcmeNode* node = Nodes[i];
    postream << "  Node "<<i<<" ----------\n";
    postream << "    owner  = "<<node->Owner()<<"\n";
    postream << "    active = "<<node->State()<<"\n";
    postream << "    exoid  = "<<node->Exodus_ID()<<"\n";
    postream << "    hindex = "<<node->HostArrayIndex()<<"\n";
    postream << "    pindex = "<<node->ProcArrayIndex()<<"\n";
    postream << "    block  = "<<node->BlockID()<<"\n";
    postream << "    coords = ";
    for (k=0; k<dimension; k++) {
      postream<<node->Coordinates()[k]<<" ";
    }
    postream<<"\n";
    postream << "    data   = ";
    for (k=0; k<node->NumData(); k++) {
      postream<<node->Data(k)<<" ";
    }
    postream<<"\n";
    postream << "    attr   = ";
    for (k=0; k<node->NumAttributes(); k++) {
      postream<<node->Attributes(k)<<" ";
    }
    postream<<"\n";
  }
  postream<<"\n";
  
  AcmeFace** Faces = 
    REINTERPRET_CAST(AcmeFace**)(face_list->EntityList());
  for( i=0 ; i<number_of_faces ; i++ ) {
    AcmeFace* face = Faces[i];
    postream << "  Face "<<i<<" ----------\n";
    postream << "    owner  = "<<face->Owner()<<"\n";
    postream << "    active = "<<face->State()<<"\n";
    postream << "    exoid  = "<<face->Exodus_ID()<<"\n";
    postream << "    hindex = "<<face->HostArrayIndex()<<"\n";
    postream << "    pindex = "<<face->ProcArrayIndex()<<"\n";
    postream << "    block  = "<<face->BlockID()<<"\n";
    postream << "    conn0  = ";
    for (k=0; k<face->Nodes_Per_Face(); k++) {
      postream<<face->NodeId(k)<<" ";
    }
    postream<<"\n";
    postream << "    conn1  = ";
    for (k=0; k<face->Nodes_Per_Face(); k++) {
      postream<<face->Node(k)->HostArrayIndex()+1<<" ";
    }
    postream<<"\n";
    postream << "    data   = ";
    for (k=0; k<face->NumData(); k++) {
      postream<<face->Data(k)<<" ";
    }
    postream<<"\n";
    postream << "    attr   = ";
    for (k=0; k<face->NumAttributes(); k++) {
      postream<<face->Attributes(k)<<" ";
    }
    postream<<"\n";
    if( face_block[face->BlockID()]->IsShell() ){
      postream << "    offset = "<<face_block[face->BlockID()]->ShellOffset()<<"\n";
    }
  }
  postream<<"\n";
  
  AcmeElem** Elems = 
    REINTERPRET_CAST(AcmeElem**)(elem_list->EntityList());
  for( i=0 ; i<number_of_elems ; i++ ) {
    AcmeElem* elem = Elems[i];
    postream << "  Elem "<<i<<" ----------\n";
    postream << "    owner  = "<<elem->Owner()<<"\n";
    postream << "    active = "<<elem->State()<<"\n";
    postream << "    exoid  = "<<elem->Exodus_ID()<<"\n";
    postream << "    hindex = "<<elem->HostArrayIndex()<<"\n";
    postream << "    pindex = "<<elem->ProcArrayIndex()<<"\n";
    postream << "    block  = "<<elem->BlockID()<<"\n";
    postream << "    conn0  = ";
    for (k=0; k<elem->Nodes_Per_Element(); k++) {
      postream<<elem->NodeId(k)<<" ";
    }
    postream<<"\n";
    postream << "    conn1  = ";
    for (k=0; k<elem->Nodes_Per_Element(); k++) {
      postream<<elem->Node(k)->HostArrayIndex()+1<<" ";
    }
    postream<<"\n";
    postream << "    data   = ";
    for (k=0; k<elem->NumData(); k++) {
      postream<<elem->Data(k)<<" ";
    }
    postream<<"\n";
    postream << "    attr   = ";
    for (k=0; k<elem->NumAttributes(); k++) {
      postream<<elem->Attributes(k)<<" ";
    }
    postream<<"\n";
  }
  postream<<"\n";
  postream.flush();
}

void DriverTopology::DisplayByBlock()
{
  int i, k;
  
  postream<<"num nodes (total/active) = "<<number_of_nodes<<" / "<<number_of_active_nodes<<"\n";
  postream<<"num faces (total/active) = "<<number_of_faces<<" / "<<number_of_active_faces<<"\n";
  postream<<"num elems (total/active) = "<<number_of_elems<<" / "<<number_of_active_elems<<"\n";
  postream<<"num node blocks  = "<<number_node_blocks<<"\n";
  for( i=0 ; i<number_node_blocks ; i++ ) {
    int ntotal = node_block[i]->EntityList()->NumEntities();
    int nactive = 0;
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) nactive++;
    }
    postream<<"  total/active = "<<ntotal<<" / "<<nactive<<"\n";
  }
  postream<<"num face blocks  = "<<number_face_blocks<<"\n";
  for( i=0 ; i<number_face_blocks ; i++ ) {
    int ntotal = face_block[i]->EntityList()->NumEntities();
    int nactive = 0;
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) nactive++;
    }
    postream<<"  total/active = "<<ntotal<<" / "<<nactive<<"\n";
  }
  postream<<"num elem blocks  = "<<number_elem_blocks<<"\n";
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    int ntotal = elem_block[i]->EntityList()->NumEntities();
    int nactive = 0;
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      if (entity->State()==AcmeEntity::ACTIVE) nactive++;
    }
    postream<<"  total/active = "<<ntotal<<" / "<<nactive<<"\n";
  }
  postream<<"\n";
  
  int n=0;
  for( i=0 ; i<number_node_blocks ; i++ ) {
    postream << "Node Block "<<i<<" ----------\n";
    postream << "  Type is "<<node_block[i]->NodeType()<<"\n";
    AcmeBlockEntityList* entity_list = node_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      postream << "  Node "<<n++<<" ----------\n";
      postream << "    owner  = "<<node->Owner()<<"\n";
      postream << "    active = "<<node->State()<<"\n";
      postream << "    exoid  = "<<node->Exodus_ID()<<"\n";
      postream << "    hindex = "<<node->HostArrayIndex()<<"\n";
      postream << "    pindex = "<<node->ProcArrayIndex()<<"\n";
      postream << "    coords = ";
      for (k=0; k<dimension; k++) {
        postream<<node->Coordinates()[k]<<" ";
      }
      postream<<"\n";
      postream << "    data   = ";
      for (k=0; k<node->NumData(); k++) {
        postream<<node->Data(k)<<" ";
      }
      postream<<"\n";
      postream << "    attr   = ";
      for (k=0; k<node->NumAttributes(); k++) {
        postream<<node->Attributes(k)<<" ";
      }
      postream<<"\n";
    }
  }
  postream<<"\n";
  
  n=0;
  for( i=0 ; i<number_face_blocks ; i++ ) {
    postream << "Face Block "<<i<<" ----------\n";
    postream << "  Type is "<<face_block[i]->FaceType()<<"\n";
    AcmeBlockEntityList* entity_list = face_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      postream << "  Face "<<n++<<" ----------\n";
      postream << "    exoid  = "<<face->Exodus_ID()<<"\n";
      postream << "    owner  = "<<face->Owner()<<"\n";
      postream << "    active = "<<face->State()<<"\n";
      postream << "    hindex = "<<face->HostArrayIndex()<<"\n";
      postream << "    pindex = "<<face->ProcArrayIndex()<<"\n";
      postream << "    conn0  = ";
      for (k=0; k<face->Nodes_Per_Face(); k++) {
        postream<<face->NodeId(k)<<" ";
      }
      postream<<"\n";
      postream << "    conn1  = ";
      for (k=0; k<face->Nodes_Per_Face(); k++) {
        postream<<face->Node(k)->HostArrayIndex()+1<<" ";
      }
      postream<<"\n";
      postream << "    data   = ";
      for (k=0; k<face->NumData(); k++) {
        postream<<face->Data(k)<<" ";
      }
      postream<<"\n";
      postream << "    attr   = ";
      for (k=0; k<face->NumAttributes(); k++) {
        postream<<face->Attributes(k)<<" ";
      }
      postream<<"\n";
      if( face_block[i]->IsShell() ){
        postream << "    offset = "<<face_block[i]->ShellOffset()<<"\n";
      }
    }
  }
  postream<<"\n";
  
  n=0;
  for( i=0 ; i<number_elem_blocks ; i++ ) {
    postream << "ELem Block "<<i<<" ----------\n";
    postream << "  Type is "<<elem_block[i]->ElemType()<<"\n";
    AcmeBlockEntityList* entity_list = elem_block[i]->EntityList();
    entity_list->IteratorStart();
    while (AcmeEntity* entity = entity_list->IteratorForward()) {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      postream << "  Elem "<<n++<<" ----------\n";
      postream << "    exoid  = "<<elem->Exodus_ID()<<"\n";
      postream << "    owner  = "<<elem->Owner()<<"\n";
      postream << "    active = "<<elem->State()<<"\n";
      postream << "    hindex = "<<elem->HostArrayIndex()<<"\n";
      postream << "    pindex = "<<elem->ProcArrayIndex()<<"\n";
      postream << "    conn0  = ";
      for (k=0; k<elem->Nodes_Per_Element(); k++) {
        postream<<elem->NodeId(k)<<" ";
      }
      postream<<"\n";
      postream << "    conn1  = ";
      for (k=0; k<elem->Nodes_Per_Element(); k++) {
        postream<<elem->Node(k)->HostArrayIndex()+1<<" ";
      }
      postream<<"\n";
      postream << "    data   = ";
      for (k=0; k<elem->NumData(); k++) {
        postream<<elem->Data(k)<<" ";
      }
      postream<<"\n";
      postream << "    attr   = ";
      for (k=0; k<elem->NumAttributes(); k++) {
        postream<<elem->Attributes(k)<<" ";
      }
      postream<<"\n";
    }
  }
  postream<<"\n";
  postream.flush();
}

void DriverTopology::DisplayCommPlan(int  num_comm_procs, 
                                   int* comm_proc_ids,
		                   int* num_nodes_to_proc, 
                                   int* comm_nodes)
{
#ifndef CONTACT_NO_MPI
  if (contact_number_of_processors(AcmeComm)>1) {
    postream<<"\nNodal comm plan:\n";
    int cnt = 0;
    int my_proc_id = contact_processor_number(AcmeComm);
    for (int i = 0;  i < num_comm_procs; i++ ) {
      postream << "    P" << my_proc_id << " <-> P" << comm_proc_ids[i] 
           << " with " << num_nodes_to_proc[i] << " nodes.\n";
      for (int j = 0; j < num_nodes_to_proc[i]; j++ ) {
        postream << "      " << j << ":  index = " << comm_nodes[cnt]-1 
             << "  exo_id = " << node_list->Entity(comm_nodes[cnt]-1)->Exodus_ID() << "\n";
        cnt++;
      }
    }
    postream<<"\n";
    postream.flush();
  }
#endif
}

void DriverTopology::SortCommNodes(int cnt, int* index, int* ids)
{
  int i, j, k, n, id, indx;

  if (cnt>1) {
    k = (cnt>>1)+1;
    n = cnt;
    for (;;) {
      if (k>1) {
        --k;
        indx = index[k-1];
        id   = ids[k-1];
      } else {
        indx       = index[n-1];
        index[n-1] = index[0];
        id         = ids[n-1];
        ids[n-1]   = ids[0];
        if (--n == 1) {
          index[0] = indx;
          ids[0]   = id;
          break;
        }
      }
      i = k;
      j = k<<1;
      while (j<=n) {
        if ((j<n) && (ids[j-1]<ids[j])) ++j;
        if (id<ids[j-1]) {
          index[i-1] = index[j-1];
          ids[i-1]   = ids[j-1];
          j += (i=j);
        } else {
          j = n+1;
        }
      }
      index[i-1] = indx;
      ids[i-1]   = id;
    }
  }
}

void DriverTopology::UpdateSharedNodeState()
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( AcmeComm ) == 1 ) return;
  int mesg = 1001;
  int i,j,len;

  int size = 0;
  for( i=0 ; i<NumCommProcs ; i++ ){
    size += NumNodesToProc[i];
  }
  char* buffer   = new char[2*sizeof(int)*size];
  char* csnd_buf = buffer;
  char* crcv_buf = buffer+sizeof(int)*size;
  int *sb, *send_buf, *rb, *recv_buf;
  sb = send_buf = REINTERPRET_CAST(int*) (csnd_buf);
  rb = recv_buf = REINTERPRET_CAST(int*) (crcv_buf);
  RequestHandle* rh = new RequestHandle[NumCommProcs];

  // Gather my data into send_buf
  int knt = 0;
  int num_procs = NumCommProcs;
  for( i=0 ; i<num_procs ; i++ ){
    for( j=0 ; j<NumNodesToProc[i] ; j++ ){
      *sb++ = node_list->Entity(CommNodes[knt+j]-1)->State();
    }
    knt += NumNodesToProc[i];
  }
  
  // Post Receives
  for( i=0 ; i<num_procs ; i++ ){
    len   = NumNodesToProc[i];
    rh[i] = contact_nonblocking_receive( mesg, rb, len, CommProcIds[i], AcmeComm );
    rb   += len;
  }

  // Sync to ensure all receives have been posted
  contact_global_sync( AcmeComm );

  if( num_procs == 0 )return;

  // Send my data
  sb = send_buf;
  for( i=0 ; i<num_procs ; i++ ){
    len = NumNodesToProc[i];
    contact_blocking_send( mesg, sb, len, CommProcIds[i], AcmeComm );
    sb += len;
  }

  // Gather my data into send_buf again but zero the data to avoid double
  // counting later.
  knt = 0;
  sb  = send_buf;
  for( i=0 ; i<num_procs ; i++ ){
    for( j=0 ; j<NumNodesToProc[i] ; j++ ){
      *sb++ = node_list->Entity(CommNodes[knt+j]-1)->State();
    }
    knt += NumNodesToProc[i];
  }  

  // Wait till all of my messages have arrived
  for( i=0 ; i<num_procs ; i++ ){
    contact_wait_msg_done( rh[i] );
  }

  knt = 0;
  rb  = recv_buf;
  for( i=0 ; i<num_procs ; i++ ){
    for( j=0 ; j<NumNodesToProc[i] ; j++ ){
      if ((AcmeEntity::EntityState)(*rb++)==AcmeEntity::ACTIVE) {
        node_list->Entity(CommNodes[knt+j]-1)->State(AcmeEntity::ACTIVE);
      }
    }
    knt += NumNodesToProc[i];
  }
  
  delete [] buffer;
  delete [] rh;
#endif
}
