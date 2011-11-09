// $Id: DriverTopology.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _DriverTopology_h_
#define _DriverTopology_h_

#include "AcmeEntityBlock.h"
#include "AcmeTopologyEntityList.h"
#include "AcmeZoltan.h"
#include "Contact_Defines.h"
#include "ContactParOStream.h"
#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif

class DriverTopology {

  public:
    
    DriverTopology( char*, int&, int*, double*, int, int,
#ifndef CONTACT_NO_MPI
                  MPI_Comm& communicator);
#else
                  int& communicator);
#endif
    ~DriverTopology();
    
    int Dimensionality(){return dimension;};
    int NumNodes()      {return number_of_nodes;};
    int NumFaces()      {return number_of_faces;};
    int NumElems()      {return number_of_elems;};
    int NumNodeBlocks() {return number_node_blocks;};
    int NumFaceBlocks() {return number_face_blocks;};
    int NumElemBlocks() {return number_elem_blocks;};
    int NumNodeVars()   {return num_node_vars;};
    int NumElemVars()   {return num_elem_vars;};
    AcmeEntityBlock* NodeBlock(int i) { return node_block[i]; };
    AcmeEntityBlock* FaceBlock(int i) { return face_block[i]; };
    AcmeEntityBlock* ElemBlock(int i) { return elem_block[i]; };
    AcmeTopologyEntityList* NodeList() { return node_list; };
    AcmeTopologyEntityList* FaceList() { return face_list; };
    AcmeTopologyEntityList* ElemList() { return elem_list; };
    
#ifndef CONTACT_NO_MPI
  inline MPI_Comm Get_Comm() {return AcmeComm;};
  inline AcmeZoltan* Get_Zoltan() {return zoltan;};
#else
  inline int Get_Comm() {return AcmeComm;};
#endif
    
    void GetInitialTopology(int& num_node_ignore,
                            int* node_ignore_list,
			    int& num_face_ignore,
                            int* face_ignore_list,
			    int& num_elem_ignore,
                            int* elem_ignore_list,
                            
                            int& number_node_blocks, 
			    int* node_block_types,
			    int* nodes_per_block,
			    int* node_eids,
                            int* node_ids,
                            
			    int& number_face_blocks, 
                            int* face_types,
			    int* faces_per_block, 
                            int* face_ids,
                            int* face_connectivity,
                            
			    int& number_element_blocks,
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
			    int  num_procs);

    void GetModsForBirthDeath(int& num_node_ignore,
                              int* node_ignore_list,
			      int& num_face_ignore,
                              int* face_ignore_list,
			      int& num_elem_ignore,
                              int* elem_ignore_list,
                         
                              int* num_node_deaths_per_block, 
	                      int* node_deaths_global_ids,
                              int* num_face_deaths_per_block, 
	                      int* face_deaths_global_ids,
                              int* num_element_deaths_per_block, 
	                      int* element_deaths_global_ids,
                              
                              int* num_node_births_per_block, 
	                      int* node_births_exodus_ids,
	                      int* node_births_global_ids,
	                      int* number_face_births_per_block, 
	                      int* face_births_global_ids,
	                      int* face_births_connectivity,
	                      int* number_element_births_per_block, 
	                      int* element_births_global_ids,
	                      int* element_births_connectivity,
                              
                              int  my_proc_id,
                              int  num_procs);

    void GetModsForDLB(int& num_node_exports,
                       int* node_export_id_list,
                       int* node_export_pids,
		       int& num_face_exports,
                       int* face_export_id_list,
                       int* face_export_pids,
		       int& num_elem_exports,
                       int* elem_export_id_list,
                       int* elem_export_pids);
    
    void GetCommPlan(int& num_comm_procs, 
                     int* comm_proc_ids,
		     int* num_nodes_to_proc, 
                     int* comm_nodes,
                     int  my_proc_id,
		     int  num_procs);

    void GetHostIDs(int& num_nodes,
                    int* node_exo_ids,
                    int* node_host_ids,
                    int& num_faces,
                    int* face_host_ids,
                    int& num_elements,
                    int* element_host_ids);
    
    void GetVariables(double* position,
		      double* nodal_vars, 
                      int&    num_nodal_vars,
		      double* element_vars,
		      int&    num_element_vars,
		      double* attributes,
		      int*    num_attributes,
		      double* node_radius,
		      double* shell_thickness,
		      double* shell_offset);
                     
    void Set_Ownership(); 
#ifndef CONTACT_NO_MPI
    void Assign_New_Ownership();
#endif

    void SortCommNodes(int, int*, int*);

    void DisplayByList();
    void DisplayByBlock();
    void DisplayCommPlan(int  num_comm_procs, 
                         int* comm_proc_ids,
		         int* num_nodes_to_proc, 
                         int* comm_nodes);
    
  private:
    
    int        dimension;
    int        number_of_nodes;
    int        number_of_faces;
    int        number_of_elems;
    int        number_of_active_nodes;
    int        number_of_active_faces;
    int        number_of_active_elems;
    int        number_node_blocks;
    int        number_face_blocks;
    int        number_elem_blocks;
    AcmeEntityBlock** node_block;
    AcmeEntityBlock** face_block;
    AcmeEntityBlock** elem_block;
    AcmeTopologyEntityList* node_list;
    AcmeTopologyEntityList* face_list;
    AcmeTopologyEntityList* elem_list;
    int        num_node_vars;
    int        num_elem_vars;
    int        NumCommProcs;
    int*       CommProcIds;
    int*       NumNodesToProc;
    int*       CommNodes;
    int        ActiveNumCommProcs;
    int*       ActiveCommProcIds;
    int*       ActiveNumNodesToProc;
    int*       ActiveCommNodes;
#ifndef CONTACT_NO_MPI
  MPI_Comm AcmeComm;
  AcmeZoltan* zoltan;
  void create_zoltan_object();
#else
  int AcmeComm;
#endif
  void UpdateSharedNodeState();
  void ComputeNodalCommPlan();
  void ComputeActiveNodalCommPlan();
  ContactParOStream postream;
};

#endif
