// STL
#include <algorithm>
#include <vector>
#include <map>

// FEM headers
#include <Element.d/Element.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/DistHelper.h>

#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>
#include <Mortar.d/NodalMortarShapeFct.d/NodalMortarShapeFct.h>

// Locally define flags
#include <Mortar.d/MortarDriver.d/MortarHandlerDefs.h>

#if defined(HB_THREAD_FFI_M_N) || defined(HB_THREAD_NODALMORTAR)
#include <Threads.d/PHelper.h>
#endif

extern MortarElement* CreateMortarElement(FaceElement*, CoordSet&, bool);

template<typename Scalar>
void
MortarHandler::ComputeOneFFIMandN(int iFFI, CoordSet &cs, std::vector<FFIPolygon<Scalar> > &CtcPolygons)
{
  Connectivity* SlaveACMEBlocksMap = PtrSlaveEntity->GetPtrACMEBlocksMap();
  int SlaveBlkId = Slave_face_block_id[iFFI]-1;
  int slave_face = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
  int iMortar = ActiveSlaveFacesToMortarEl[slave_face];
  MortarElement* MortarEl = MortarEls[iMortar];
  if(InteractionType==MortarHandler::FSI) {
    CtcPolygons[iFFI].ComputeNormalN(MortarEl, cs, MortarIntegrationRule);
  } 
  else if(InteractionType==MortarHandler::CTC) {
    CoordSet& csm = PtrMasterEntity->GetNodeSet();
    double offset = PtrMasterEntity->GetShellThickness()/2 + PtrSlaveEntity->GetShellThickness()/2;
    CtcPolygons[iFFI].ComputeNormalGeoGap(MortarEl, cs, csm, MortarIntegrationRule, offset);
    CtcPolygons[iFFI].ComputeGradNormalM(MortarEl, cs, csm, MortarIntegrationRule);
    CtcPolygons[iFFI].ComputeGradNormalN(MortarEl, cs, csm, MortarIntegrationRule);
  }
  else { // TIED
    CtcPolygons[iFFI].ComputeM(MortarEl, cs, MortarIntegrationRule);
    CtcPolygons[iFFI].ComputeN(MortarEl, cs, MortarIntegrationRule);
  }
}

template<typename Scalar>
void
MortarHandler::MakeOneNodalMortarLMPC(int i, std::vector<FFIPolygon<Scalar> > &CtcPolygons, bool Dual)
{
  NodalMortars[i].SetRefData(ActiveSlaveNodes[i], MortarScaling);
  NodalMortars[i].MakeSlaveLink(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, Dual);
  NodalMortars[i].MakeMasterLink(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
  if(InteractionType==MortarHandler::FSI)
    NodalMortars[i].BuildWetFSICoupling(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
  else if(InteractionType==MortarHandler::CTC)
    NodalMortars[i].BuildMortarCtcLMPC(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
  else // TIED
    NodalMortars[i].BuildMortarLMPC(ActiveSlaveNodeToElem, &ActiveSlaveElemSet, SlaveFaceToFFIConnect, &CtcPolygons[0]);
}

template<class Scalar>
void
MortarHandler::CreateFFIPolygon()
{
  // clear data from previous call to this function
  MortarEls.clear();
  NodalMortars.clear();
  ActiveSlaveNodes.clear();
  ActiveSlaveElemSet.deleteElems();
  ActiveMasterElemSet.deleteElems();
  if(ActiveSlaveNodeToElem) { delete ActiveSlaveNodeToElem; ActiveSlaveNodeToElem = 0; }
  if(ActiveMasterNodeToElem) { delete ActiveMasterNodeToElem; ActiveMasterNodeToElem = 0; }
  if(SlaveFaceToFFIConnect) { delete SlaveFaceToFFIConnect; SlaveFaceToFFIConnect = 0; }
  ActiveSlaveFacesToMortarEl.clear();

  // BRUTE FORCE TEST
  CoordSet& cs = PtrSlaveEntity->GetNodeSet();

#ifdef HB_MORTAR_TIMER
  double t0,t00,dt0,dt1,dt2,dt3,dt4;
  t00 = getTime(); t0 = t00;
#endif  
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * Create the FFI polygon \n");
#endif  
  // Allocate CtcPolygons array
  std::vector<FFIPolygon<Scalar> > CtcPolygons(nFFI);

  FaceElemSet* MasterElemSet = GetPtrMasterFaceElemSet();
  FaceElemSet* SlaveElemSet  = GetPtrSlaveFaceElemSet();
  Connectivity* SlaveACMEBlocksMap = PtrSlaveEntity->GetPtrACMEBlocksMap();
  Connectivity* MasterACMEBlocksMap= PtrMasterEntity->GetPtrACMEBlocksMap();

  // Create each FFIPolygon
  // !!! ACME indices in Fortran indexing !!!
  //filePrint(stderr," ### Write ACMEFFI_data array offset\n");
  //FILE* FFIDataFile = fopen("FFIDataOffset.txt","w");
  std::vector<int> IndexActiveSlaveFaces ; IndexActiveSlaveFaces.reserve(512);
  std::vector<int> IndexActiveMasterFaces; IndexActiveMasterFaces.reserve(512);
  int SlaveBlkId, MasterBlkId;
  for(int iFFI=0; iFFI < nFFI; iFFI++) {
    SlaveBlkId  = Slave_face_block_id[iFFI]-1;
    if(Slave_face_index_in_block[iFFI] <= 0) cerr << "error here in MortarHandler::CreateFFIPolygon, iFFI = " << iFFI
                                                  << ", Slave_face_index_in_block[iFFI] = " << Slave_face_index_in_block[iFFI] << endl;
    int slave_face  = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
    FaceElement *SlaveFaceEl  = (*SlaveElemSet)[slave_face];

    if(SelfContact) { // TODO
      MasterBlkId = Master_face_block_id[iFFI]-1;
      MasterACMEBlocksMap = SlaveACMEBlocksMap;
      MasterElemSet = SlaveElemSet;
      PtrMasterEntity = PtrSlaveEntity;
    }
    else {
      MasterBlkId = Master_face_block_id[iFFI]-1 - 4; // in all the ACME blocks, the master blocks are last
    }
    if(Master_face_index_in_block[iFFI] <= 0) cerr << "error here in MortarHandler::CreateFFIPolygon, iFFI = " << iFFI 
                                                << ", Master_face_index_in_block[iFFI] = " << Master_face_index_in_block[iFFI] << endl;
    int master_face = (*MasterACMEBlocksMap)[MasterBlkId][Master_face_index_in_block[iFFI]-1];
    FaceElement *MasterFaceEl = (*MasterElemSet)[master_face];

    int FFIDataOffset = ACMEFFI_index[iFFI];
    int nVertices = (int) ACMEFFI_data[FFIDataOffset];
    double* ACME_FFI_Data = &ACMEFFI_data[FFIDataOffset];

    if(GeomType == MortarHandler::NON_MATCHING) {
      CtcPolygons[iFFI].SetFFIPolygon(MasterFaceEl, SlaveFaceEl, nVertices, ACME_FFI_Data,
                                      MAX_FFI_DERIVATIVES, MAX_FFI_SECOND_DERIVATIVES);
    }
    else {
      CtcPolygons[iFFI].SetFFIPolygon(MasterFaceEl, SlaveFaceEl, nVertices, ACME_FFI_Data, 0, 0);
    }
    
    // Add ptr to FFI in slave face element 
    //SlaveFaceEl->AddPtrFFI(&CtcPolygons[iFFI]); 
#ifdef FFI_MORTAR_DEBUG
    filePrint(stderr," * -------------------------------- \n");
    filePrint(stderr,"   -> create FFIPolygon %d\n", iFFI+1);
    std::cerr << "   -> MasterBlkId = " << MasterBlkId <<" master_face = "<<master_face<< std::endl;
    std::cerr << "   -> SlaveBlkId  = " << SlaveBlkId<<" slave_face = "<<slave_face << std::endl;
    filePrint(stderr,"   -> FFIDataOffset = %d\n",FFIDataOffset);
    std::cerr<<"    -> slave  face element:\n"; SlaveFaceEl->print();
    std::cerr<<"    -> master face element:\n"; MasterFaceEl->print();
    CtcPolygons[iFFI].Print();
#endif

    // Get indices of active slave/master face elements 
    IndexActiveSlaveFaces.push_back(slave_face);
    IndexActiveMasterFaces.push_back(master_face);

    //filePrint(FFIDataFile," %6d, nVertices = %3d\n",FFIDataOffset,nVertices);
  }
  //fclose(FFIDataFile);
#ifdef HB_MORTAR_TIMER
  dt0 = (getTime()-t0)/1000.;
#endif
  // Eliminates duplicate face indices
  std::sort(IndexActiveSlaveFaces.begin() , IndexActiveSlaveFaces.end());
  std::sort(IndexActiveMasterFaces.begin(), IndexActiveMasterFaces.end());

  IndexActiveSlaveFaces.erase(unique(IndexActiveSlaveFaces.begin(),IndexActiveSlaveFaces.end()), 
                               IndexActiveSlaveFaces.end());

  IndexActiveMasterFaces.erase(unique(IndexActiveMasterFaces.begin(),IndexActiveMasterFaces.end()),
                               IndexActiveMasterFaces.end());

  int nActiveSlaveFaces = IndexActiveSlaveFaces.size();
  
  // Create active master & slave face element set
  // & associated connectivity object
#ifdef MORTAR_DEBUG
  int nActiveMasterFaces= IndexActiveMasterFaces.size();
  filePrint(stderr,"   * number of active slave face elem. = %6d\n", nActiveSlaveFaces);
  filePrint(stderr,"   * number of active master face elem.= %6d\n", nActiveMasterFaces);
  filePrint(stderr,"   * active slave & master face element sets\n"); 
#endif
#ifdef HB_MORTAR_TIMER
  t0 = getTime();
#endif
  std::map<int,int> IndexActiveSlaveFacesToActiveSlaveElemSetMap; 
  for(int i=0; i<int(IndexActiveSlaveFaces.size()); ++i) {
    ActiveSlaveElemSet.elemadd(i, (*SlaveElemSet)[IndexActiveSlaveFaces[i]]); 
    IndexActiveSlaveFacesToActiveSlaveElemSetMap[IndexActiveSlaveFaces[i]] = i;
  }  
  Connectivity ActiveSlaveElemToNode(&ActiveSlaveElemSet);
  ActiveSlaveNodeToElem = ActiveSlaveElemToNode.reverse();

  for(int i=0; i<int(IndexActiveMasterFaces.size()); ++i) {
    ActiveMasterElemSet.elemadd(i, (*MasterElemSet)[IndexActiveMasterFaces[i]]); 
  }
  Connectivity ActiveMasterElemToNode(&ActiveMasterElemSet);
  ActiveMasterNodeToElem = ActiveMasterElemToNode.reverse();

  ActiveSlaveNodes.clear();
  ActiveSlaveNodes.reserve(128);
  for(int iel=0; iel<int(IndexActiveSlaveFaces.size()); iel++) {
    size_t old_size(ActiveSlaveNodes.size());
    ActiveSlaveNodes.insert(ActiveSlaveNodes.end(), ActiveSlaveElemSet[iel]->nNodes(), 0);
    ActiveSlaveElemSet[iel]->GetNodes(&ActiveSlaveNodes[old_size]);
  }
  //filePrint(stderr,"   * ActiveSlaveNodes.size() = %d\n", (int) ActiveSlaveNodes.size());
  std::sort(ActiveSlaveNodes.begin(),ActiveSlaveNodes.end());
  ActiveSlaveNodes.erase(std::unique(ActiveSlaveNodes.begin(),ActiveSlaveNodes.end()),
                         ActiveSlaveNodes.end());

  int nActiveSlaveNodes = ActiveSlaveNodes.size();
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * nActiveSlaveNodes = %d\n", nActiveSlaveNodes);
#endif

  // Create the map ActiveSlaveFaces to FFIs
  // 1) count the number of FFIs per ActiveSlaveFaces
  int* SlaveFaceToFFIpointer = new int[nActiveSlaveFaces+1];
  int* SlaveFaceToFFItarget = new int[nFFI];
  { 
    std::vector<int> nFFIperSlaveFace(nActiveSlaveFaces, 0);
    for(int iFFI=0; iFFI<nFFI; iFFI++) {
      int SlaveBlkId = Slave_face_block_id[iFFI]-1;
      int slave_face = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
      nFFIperSlaveFace[IndexActiveSlaveFacesToActiveSlaveElemSetMap[slave_face]]++;  
    }
    // 2) set the map pointer array
    SlaveFaceToFFIpointer[0] = 0;
    for(int i=0; i<nActiveSlaveFaces; i++)
      SlaveFaceToFFIpointer[i+1] = SlaveFaceToFFIpointer[i]+nFFIperSlaveFace[i];

    // 3) fill the target (map) array
    std::vector<int> Offset(nActiveSlaveFaces, 0);
    for(int iFFI=0; iFFI<nFFI; iFFI++){
      int SlaveBlkId = Slave_face_block_id[iFFI]-1;
      int slave_face = (*SlaveACMEBlocksMap)[SlaveBlkId][Slave_face_index_in_block[iFFI]-1];
      int i = IndexActiveSlaveFacesToActiveSlaveElemSetMap[slave_face];
      SlaveFaceToFFItarget[SlaveFaceToFFIpointer[i]+Offset[i]] = iFFI; 
      Offset[i]++;
    }
  }
  SlaveFaceToFFIConnect = new Connectivity(nActiveSlaveFaces, SlaveFaceToFFIpointer, SlaveFaceToFFItarget);

  //SlaveFaceToFFIConnect.print();
#ifdef HB_MORTAR_TIMER
  dt1 = (getTime()-t0)/1000.;
  t0 = getTime();
#endif
  // Create Mortar element
  bool Dual = false;
  if(MortarType==MortarHandler::DUAL) { Dual = true; }
#ifdef MORTAR_DEBUG
  if(Dual) filePrint(stderr,"   * Create Dual Mortar elements \n");
  else     filePrint(stderr,"   * Create Std Mortar elements \n"); 
#endif
  MortarEls.assign(nActiveSlaveFaces, (MortarElement*) 0);

  for(int i=0; i<nActiveSlaveFaces; ++i) {
    FaceElement* SlaveFaceEl  = (*SlaveElemSet)[IndexActiveSlaveFaces[i]];
    MortarEls[i] = CreateMortarElement(SlaveFaceEl, cs, Dual); 
    ActiveSlaveFacesToMortarEl[IndexActiveSlaveFaces[i]] = i;
  }
#ifdef HB_MORTAR_TIMER
  dt2 = (getTime()-t0)/1000.;
  t0 = getTime();
#endif
  // Integrate shape fcts product
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * Compute FFI contributions to the shape fct products\n"); 
#endif
#ifdef HB_THREAD_FFI_M_N
  execParal2R(nFFI, this, &MortarHandler::ComputeOneFFIMandN<Scalar>, cs, CtcPolygons); // TODO
#else
  for(int i=0; i<nFFI; ++i) ComputeOneFFIMandN(i, cs, CtcPolygons);
#endif
#ifdef HB_MORTAR_TIMER
  dt3 = (getTime()-t0)/1000.;
  t0 = getTime();
#endif
  // Create nodal Mortar shape fcts (NO BOUNDARY MODIFICATION) 
#ifdef MORTAR_DEBUG
  filePrint(stderr,"   * Create nodal Mortar shape fcts\n"); 
  filePrint(stderr,"     -> NO BOUNDARY MODIFICATION\n"); 
  filePrint(stderr,"     -> NodalMortars.size() = %d\n",NodalMortars.size()); 
#endif
  NodalMortars.assign(nActiveSlaveNodes, NodalMortarShapeFct());
#ifdef HB_THREAD_NODALMORTAR
  execParal2R(nActiveSlaveNodes, this, &MortarHandler::MakeOneNodalMortarLMPC, CtcPolygons, Dual); // TODO
#else
  for(int i=0; i<nActiveSlaveNodes; ++i) MakeOneNodalMortarLMPC(i, CtcPolygons, Dual);
#endif
#ifdef HB_MORTAR_TIMER
  dt4 = (getTime()-t0)/1000.;
  double dttot = dt0+dt1+dt2+dt3+dt4;
  filePrint(stderr,"   * CPU time statistic of MortarHandler::CreateFFIPolygon(...):\n");
  filePrint(stderr,"    # making the FFI contact polygons (seq): %e s (%2.4f %%)\n",dt0,100.*dt0/dttot);
  filePrint(stderr,"    # making mapping & connectivities (seq): %e s (%2.4f %%)\n",dt1,100.*dt1/dttot);
  filePrint(stderr,"    # creating the mortar elements    (seq): %e s (%2.4f %%)\n",dt2,100.*dt2/dttot);
#ifdef HB_THREAD_FFI_M_N
  filePrint(stderr,"    # computing FFI M & N contribution(// ): %e s (%2.4f %%)\n",dt3,100.*dt3/dttot);
#else
  filePrint(stderr,"    # computing FFI M & N contribution(seq): %e s (%2.4f %%)\n",dt3,100.*dt3/dttot);
#endif
#ifdef HB_THREAD_NODALMORTAR
  filePrint(stderr,"    # assembling nodal mortar LMPCs   (// ): %e s (%2.4f %%)\n",dt4,100.*dt4/dttot);
#else
  filePrint(stderr,"    # assembling nodal mortar LMPCs   (seq): %e s (%2.4f %%)\n",dt4,100.*dt4/dttot);
#endif
  filePrint(stderr,"    # total CPU time                       : %e s\n",dttot);
#endif
}

