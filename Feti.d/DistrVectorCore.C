#include <Feti.d/DistrVector.h>

void
DistrInfo::initialize()
{
  domLen = 0;
  subLen = 0; 
  subOffset = 0;
  threadOffset = 0; 
  threadLen = 0;
  masterFlag = 0;
}

DistrInfo::DistrInfo(int ns)
{
 initialize();
 numDom = ns;
 domLen = new int[numDom];
}

DistrInfo::~DistrInfo()
{
  if(domLen) { delete [] domLen; domLen = 0; }  // union with subLen
  if(masterFlag) { delete [] masterFlag; masterFlag = 0; }
  if(subOffset) { delete [] subOffset; subOffset = 0; }
  if(threadOffset) { delete [] threadOffset; threadOffset = 0; }
  if(threadLen) { delete [] threadLen; threadLen = 0; }
}

void DistrInfo::setMasterFlag()
{
  // this version for internal DistrInfo's or if interface values are weighted
  masterFlag = new bool[len];
  int i;
  for(i = 0; i < len; ++i)
    masterFlag[i] = true;
}

void
DistrInfo::computeOffsets()
{
  if(subOffset)
    return;
  int iSub, iThread;

  subOffset = new int[numLocSub];
  int cOffset = 0;
  for(iSub = 0; iSub < numLocSub; ++iSub) {
    subOffset[iSub] = cOffset;
    cOffset += subLen[iSub];
  }
  len = cOffset;
 
  // XML NEED TO UPDATE FOR OPENMP
  numLocThreads = threadManager->numThr();
  threadOffset = new int[numLocThreads+1];
  threadLen = new int[numLocThreads+1]; // PJSA: for DistVec class in Math.d
  int nRemain = len%numLocThreads;
  int npt = len/numLocThreads;
  for(iThread = 0; iThread < numLocThreads+1; ++iThread) {
    threadOffset[iThread] = (iThread < nRemain ? iThread : nRemain)
       +iThread*npt;
    threadLen[iThread] = (iThread < nRemain ? 1 : 0) + npt; // PJSA: for DistVec class in Math.d
  }
  threadOffset[0] = 0;
  
}

void
DistrInfo::recomputeOffsets()
{
  if(subOffset) { delete [] subOffset; subOffset = 0; }
  if(threadOffset) { delete [] threadOffset; threadOffset = 0; }
  if(threadLen) { delete [] threadLen; threadLen = 0; }

  computeOffsets();
}

