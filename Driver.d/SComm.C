#include <Driver.d/SComm.h>
#include <Driver.d/SubDomain.h>

SComm::~SComm()
{
  if(exchangeData) { delete [] exchangeData; exchangeData=0; }
  //if(sharedDOFs) { delete sharedDOFs; sharedDOFs=0; }
  if(neighb) { delete [] neighb; neighb = 0; }
  if(sharedNodes) { delete sharedNodes; sharedNodes = 0; }
  if(subNums/* && !salinasFlag*/) { delete [] subNums; subNums = 0; }
  if(remoteId) { delete [] remoteId; remoteId = 0; }
  if(isEdgeNeighb) { delete [] isEdgeNeighb; isEdgeNeighb = 0; }
  // don't delete glSubToLocal (DistDom)
  if(NumNeighb) { delete [] NumNeighb; NumNeighb = 0; }
  if(SubNums) { for(int i=0; i<numDofType; ++i) { if(SubNums[i]) delete [] SubNums[i]; } delete [] SubNums; }
  if(TypeMap) { for(int i=0; i<numDofType; ++i) { if(TypeMap[i]) delete [] TypeMap[i]; } delete [] TypeMap; }
  if(SharedDOFs) { for(int i=0; i<numDofType; ++i) { if(SharedDOFs[i]) delete SharedDOFs[i]; } delete [] SharedDOFs; }
  if(sharedDOFsPlus) delete sharedDOFsPlus;
}

void 
SComm::deleteTypeSpecificList(DofType type)
{
  if(SubNums && SubNums[type]) { delete [] SubNums[type]; SubNums[type] = 0; }
  if(TypeMap && TypeMap[type]) { delete [] TypeMap[type]; TypeMap[type] = 0; }
  if(SharedDOFs && SharedDOFs[type]) { delete SharedDOFs[type]; SharedDOFs[type] = new Connectivity(); }
  NumNeighb[type] = 0;
}

void
SComm::setEdgeNeighb(int _numEdgeNeighb, bool *_isEdgeNeighb)
{
  numEdgeNeighb = _numEdgeNeighb;
  if(isEdgeNeighb) delete [] isEdgeNeighb;
  isEdgeNeighb = _isEdgeNeighb;
}

void *
SComm::getExchangeData(int iSub)
{
#ifdef DISTRIBUTED
  cerr << " *** ERROR: SComm::getExchangeData is not implemented for distributed memory \n";
  return 0;
#else
  return neighb[iSub]->getExchangePointer(remoteId[iSub]);
#endif
}

void
SComm::setExchangeData(int iSub, void *data)
{
#ifdef DISTRIBUTED
  cerr << " *** ERROR: SComm::setExchangeData is not implemented for distributed memory \n";
#else
  exchangeData[iSub] = data;
#endif
}

bool
SComm::isLocal(int iSub)
{
  // PJSA 6-3-04: can use this function to determine how & when to delete exchange data memory
  // see sendExpDOFList and gatherDOFList for example
#ifdef DISTRIBUTED
  if(glSubToLocal != 0 && (glSubToLocal[subNums[iSub]] < 0))
    return false;
  else
#endif
    return true;
}

void 
SComm::setTypeSpecificList(DofType type, int *_subNums, Connectivity *_sharedDOFs)
{
  if(!SharedDOFs) {
    NumNeighb = new int[numDofType];
    SubNums = new int * [numDofType];
    SharedDOFs = new Connectivity * [numDofType];
    for(int i=0; i<numDofType; ++i) { 
      NumNeighb[i] = 0;
      SubNums[i] = 0;
      // SharedDOFs[i] = 0;
      SharedDOFs[i] = new Connectivity();
    }
  }
  NumNeighb[type] = _sharedDOFs->csize();
  if(SubNums[type]) delete [] SubNums[type]; 
  SubNums[type] = _subNums;
  if(SharedDOFs[type]) delete SharedDOFs[type]; 
  SharedDOFs[type] = _sharedDOFs;
}

int *
SComm::mergeTypeSpecificLists()
{
  // build combined list of all types 0, 1 and 2 shared dofs: 
  // also make **TypeMap and *boundDofFlag (returned)
  // update **neighb and *remoteId 
  // resize **exchangeData but i don't think it is necessary to add "virtual nodes"
  // to sharedNodes list. however, check in code where sharedNodes is used and convert to "std"
  // NumNeighb and SubNums etc.
  if(TypeMap) { for(int i=0; i<numDofType; ++i) { if(TypeMap[i]) delete [] TypeMap[i]; } delete [] TypeMap; }

  int i,j,k,l;
  // step 1. initialize combined list to type 0
  Connectivity *allSharedDOFs = SharedDOFs[0]->copy();
  int *allSubNums = new int[NumNeighb[0]];
  for(i=0; i<NumNeighb[0]; ++i) allSubNums[i] = SubNums[0][i];
  TypeMap = new int * [numDofType];
  TypeMap[0] = new int[SharedDOFs[0]->numConnect()]; 

  // step 2. loop over all other types and add to combined list
  for(i=1; i<int(all); ++i) {
    if(NumNeighb[i] > 0) {
      allSharedDOFs->combine(SharedDOFs[i], allSubNums, SubNums[i]); 
      TypeMap[i] = new int[SharedDOFs[i]->numConnect()];
    }
    else {
      TypeMap[i] = 0;
    }
  }
  //int allNumNeighb = allSharedDOFs->csize();

  // step 3. make boundDofFlag and TypeMaps
  int *boundDofFlag = new int[allSharedDOFs->numConnect()];
  int count = 0;
  for(i = 0; i < allSharedDOFs->csize(); ++i) {
    int neighb = allSubNums[i];
    for(j = 0; j < int(all); ++j) {
      for(k = 0; k < SharedDOFs[j]->csize(); ++k) {
        if(SubNums[j][k] == neighb) {
          for(l = 0; l < SharedDOFs[j]->num(k); ++l) {
            TypeMap[j][SharedDOFs[j]->offset(k)+l] = count;
            boundDofFlag[count++] = j;
          }
          break;
        }
      }
    }
  } 
  //numNeighb = allNumNeighb;
  //subNums = allSubNums;
  //sharedDOFs = allSharedDOFs;
  setTypeSpecificList(all, allSubNums, allSharedDOFs);
  sharedDOFs = allSharedDOFs; // tmp fix to compile MDAxi

  for(i = int(all); i<numDofType; ++i) TypeMap[i] = 0;

  return boundDofFlag;
}

void 
SComm::setTypeMap(DofType t, int *map)
{
  if(TypeMap[t]) delete [] TypeMap[t];
  TypeMap[t] = map;
}
