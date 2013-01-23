#ifndef _FS_COMMUNICATOR_H_
#define _FS_COMMUNICATOR_H_
#include <fstream>
#include <Utils.d/Connectivity.h>
#include <complex>

#include <algorithm>
using namespace std;
//#ifdef SUNOS
//#include <algorithm.cc>
//#else
//#include <algo.h>

#ifdef USE_MPI
#ifdef RS6000_SYS
#include "/usr/lpp/ppe.poe/include/mpi.h"
#else
#include <iostream>
using namespace std;
#ifndef MPI_NO_CPPBIND
#define MPI_NO_CPPBIND
#endif
#include <mpi.h>
#endif // USE_MPI
#endif // RS6000_SYS

#include <map>
#include <vector>
#include <iterator>

struct
FSRecInfo {
  int cpu, len;
};

class Communicator;

class FSCommunicator {
     int thisCPU;
     int numCPU;
#ifdef USE_MPI
     MPI_Comm comm;
#endif
   public:
     FSCommunicator();
#ifdef USE_MPI
     FSCommunicator(Communicator *);
     MPI_Comm& getComm() { return comm; }
#endif
     int cpuNum() { return thisCPU; }
     int size() { return numCPU; }

     void exchange(int tag, int numNeighb, int *cpus,
      int *sndPtr, int *sndData, int *rcvPtr, int *rcvData);
     void exchange(int tag, int numNeighb, int *cpus,
       int *sndLen, int **sndData, int *crvLen, int **rcvData); 
     void exchange(int tag, int numNeighb, int *cpus,
       int *sndLen, double **sndData, int *crvLen, double **rcvData); 
     void exchange(int tag, int numNeighb, int *cpus,
       int *sndLen, complex<double> **sndData, int *crvLen, complex<double> **rcvData);

     void abort(int errorCode);
     void sync() {
#ifdef USE_MPI
       MPI_Barrier(comm);
#endif
     }
    
    // PJSA: templated functions
    template <class Type>
       Type globalSum(Type);
    template <class Type>
       Type globalMax(Type);
    template <class Type>
       Type globalMin(Type);
    template <class Type>
       void globalMax(int, Type*);
    template <class Type>
       void globalMin(int, Type*);
    template <class Type>
       void globalSum(int, Type*);
    template <class Type>
       void sendTo(int cpu, int tag, Type *data, int len);
    template <class Type>
       FSRecInfo recFrom(int tag, Type *data, int len);
     
    template <class Type>
       void allGather(Type *send_data, int send_count, Type *recv_data,
                      int recv_count);
    template <class Type>
       void allGatherv(Type *send_data, int send_count, Type *recv_data,
                      int recv_counts[], int displacements[]);

#ifdef USE_MPI
    template <class Type>
       void globalMpiOp(int num, Type *data, MPI_Op mpi_op);
#endif

     template <class Type>
       void reduce(int num, Type*data, int root = 0);

     template <class Type>
       void broadcast(int num, Type* data, int root = 0);

};

void initComm(int, char *[]);
void closeComm();
// extern FSCommunicator *communicator;

class Connectivity;

class Triplet { 
  public:
    int from, to, cpuID;
    Triplet(int _from, int _to, int _cpuID) {
      from = _from; to = _to; cpuID = _cpuID;
    }
    bool operator < (const Triplet &x) const { // required by the STL sort algorithm
      // we want to order first by cpuID (with local coms first), then by origin and then by destination
      return cpuID < x.cpuID || (cpuID  == x.cpuID &&
        (from < x.from || (from == x.from && to < x.to))); 
    }
    bool operator == (const Triplet &x) const {
      return cpuID == x.cpuID && from == x.from && to == x.to; 
    }
};

template <class T>
class FSSubRecInfo {
  public:
    T *data;
    int len;
    int leadDim;
    int nvec;
    bool isSend;
    FSSubRecInfo() { 
      data = 0; len = 0; leadDim = 0; nvec = 0; isSend = false; 
    }
    FSSubRecInfo(int _len, int _leadDim, int _nvec) {
      data = 0; len = _len; leadDim = _leadDim; nvec = _nvec; isSend = true;
    }
    FSSubRecInfo(int _len, int _leadDim, int _nvec, bool _isSend) {
      data = 0; len = _len; leadDim = _leadDim; nvec = _nvec; isSend = _isSend;
    }
    FSSubRecInfo(const FSSubRecInfo &c) {
      data = c.data; len = c.len; leadDim = c.leadDim; nvec = c.nvec; isSend = c.isSend;
    }

};

typedef map<Triplet, int> MapType;
typedef MapType::value_type ValuePair;
typedef MapType::iterator MapIter;

/*****************************************************************************
   FSCommPattern represent a communication pattern. 
     Communication is based on the model that a message from one
     subdomain to another subdomain is made of a number of vectors.
     Such vectors are stored in a matrix form. That matrix may have
     larger columns than the actual message length. This allows a subdomain
     to send different subparts of a given matrix to different subdomains.

   A communication patterns contains the following information:
     - Length of a vector from one subdomain to its neighbors
     - Number of vectors in each message
     - Leading dimension of the matrix containing the vectors
 *****************************************************************************/
template <class T>
class FSCommPattern {
  public:
     enum Mode { Share, CopyOnSend };
     enum Symmetry { Sym, NonSym };
  protected:
     FSCommunicator *communicator;  // PJSA
     Mode mode;
     Symmetry sym;
     int myCPU;
     int numChannels;
     MapType channelMap;
     vector< FSSubRecInfo<T> > sRecInfo;
     Connectivity *subToCPU;
     int numCPUs;
     int *reverseChannel; // corresponding reverse channel
     int *channelToCPU; // this CPU is marked as -1
     int numNeighbCPUs;
     int *neighbCPUs;
     Connectivity *sendConnect;
     Connectivity *recvConnect;
     // Buffers for local copy-communication
     T *localDBuffer;
     // Buffers for cross-memory communication on a cpu per cpu basis
     int *crossSendLen, *crossRecvLen;
     T **crossSendBuffer;
     T **crossRecvBuffer;
     void makeSendAndRecvConnectivities();

  public:
     FSCommPattern(FSCommunicator *communicator, Connectivity *cpuToSub, int myCPU = 0, Mode = Share, Symmetry = Sym);
     // FSCommPattern() {};
     ~FSCommPattern();
     void setLen(int glFrom, int glTo, int size, int lDim = -1, int nv = 1);
     void finalize(); // complete the internal setup
     FSSubRecInfo<T> recData(int glFrom, int glTo);
     void sendData(int glFrom, int glTo, T *);
     FSSubRecInfo<T> getSendBuffer(int glFrom, int glTo);
     void exchange();

     Mode getMode() { return mode; }
     Symmetry getSym() {return sym; }
     int getMyCPU() { return myCPU; }
     int getNumChannels() { return numChannels; }
     int *getChannelToCPU() { return channelToCPU; }
     MapType getChannelMap() { return channelMap; }
     Connectivity *getSubToCPU() { return subToCPU; }
};


template <class T>
FSCommPattern<T>::FSCommPattern(FSCommunicator *_communicator, Connectivity *cpuToSub, int _myCPU, Mode _mode, Symmetry _sym) 
{
  communicator = _communicator;
  myCPU = _myCPU;
  mode = _mode;
  sym = _sym;
  numChannels = 0;
  numCPUs = communicator->size();
  subToCPU = cpuToSub->reverse();
  reverseChannel = 0;
  neighbCPUs = 0;
  channelToCPU = 0;
  crossSendBuffer = 0;
  crossRecvBuffer = 0;
  crossSendLen = 0;
  crossRecvLen = 0;
  localDBuffer = 0;
  sendConnect = 0;
  recvConnect = 0;
}

template <class T>
FSCommPattern<T>::~FSCommPattern()
{
  if(reverseChannel) { delete [] reverseChannel; reverseChannel = 0; }
  if(neighbCPUs) { delete [] neighbCPUs; neighbCPUs = 0; }
  if(channelToCPU) { delete [] channelToCPU; channelToCPU = 0; }
  if(crossSendBuffer) { 
     delete [] crossSendBuffer[0]; crossSendBuffer[0] = 0;
     delete [] crossSendBuffer; crossSendBuffer = 0; 
  }
  if(crossRecvBuffer) { delete [] crossRecvBuffer; crossRecvBuffer = 0; }
  if(crossSendLen) { delete [] crossSendLen; crossSendLen = 0; }
  if(crossRecvLen) { delete [] crossRecvLen; crossRecvLen = 0; }
  if(localDBuffer) { delete [] localDBuffer; localDBuffer = 0; }
  if(subToCPU) { delete subToCPU; subToCPU = 0; }
  if(sendConnect) { delete sendConnect; sendConnect = 0; }
  if(recvConnect) { delete recvConnect; recvConnect = 0; }
}

/*  // PJSA: this copy constructor doesn't seem to be functioning correctly
template <class T>
template <class TB>
FSCommPattern<T>::FSCommPattern(FSCommPattern<TB> &cp, Mode _mode)
{
  myCPU = cp.getMyCPU();
  mode = _mode;
  sym = (FSCommPattern<T>::Symmetry) cp.getSym();
  numChannels = cp.getNumChannels();
  numCPUs = communicator->size();
  subToCPU = cp.getSubToCPU();

  MapType tmpMap = cp.getChannelMap();
  MapIter iter = tmpMap.begin();
  while(iter != tmpMap.end()) {
    Triplet key = (*iter).first;
    int channelID = (*iter).second;
    channelMap.insert(ValuePair(Triplet(key.from, key.to, key.cpuID), channelID));
    FSSubRecInfo<TB> rInfo = cp.getSendBuffer(key.from, key.to);
    sRecInfo.insert(sRecInfo.end(), FSSubRecInfo<T>(rInfo.len, rInfo.leadDim, rInfo.nvec, rInfo.isSend));
    ++iter;
  }
  finalize();
}
*/

template <class T>
void
FSCommPattern<T>::setLen(int glFrom, int glTo, int len, int leadDim, int nvec)
{
  // serial function

  int cpuID = (*subToCPU)[glTo][0];
  if(leadDim < 0) leadDim = len;

  if(cpuID == myCPU) {
    channelMap.insert(ValuePair(Triplet(glFrom, glTo, -1), numChannels));
    sRecInfo.insert(sRecInfo.end(), FSSubRecInfo<T>(len, leadDim, nvec));
    numChannels++;
  }
  else {
    // make send channel
    channelMap.insert(ValuePair(Triplet(glFrom, glTo, cpuID), numChannels));
    sRecInfo.insert(sRecInfo.end(), FSSubRecInfo<T>(len, leadDim, nvec));
    numChannels++;
    // make receive channel
    channelMap.insert(ValuePair(Triplet(glTo, glFrom, cpuID), numChannels));
    sRecInfo.insert(sRecInfo.end(), FSSubRecInfo<T>());  // empty for now, see ::finalize()
    numChannels++;
  }
}

template <class T>
void
FSCommPattern<T>::finalize()
{
  int i, j;
  // Step 1. build isSend, reverseChannel and channelToCPU
  reverseChannel = new int[numChannels];
  channelToCPU = new int[numChannels];
  MapIter iter = channelMap.begin();
  while(iter != channelMap.end()) {
    Triplet key = (*iter).first;
    i = (*iter).second; // channelID
    if(sRecInfo[i].isSend) reverseChannel[i] = i+1;
    else reverseChannel[i] = i-1;
    channelToCPU[i] = key.cpuID;
    ++iter;
  }    

  // Step 2. make send and recv connectivities
  makeSendAndRecvConnectivities();

  // Step 3. make reverse channel FSSubRecInfo's
  if(numCPUs > 1) {
    if(sym == Sym) {
      for(i = 0; i < numChannels; ++i) {
        if(channelToCPU[i] >= 0 && sRecInfo[i].isSend) {
          sRecInfo[reverseChannel[i]].len  = sRecInfo[i].len;
          sRecInfo[reverseChannel[i]].leadDim = sRecInfo[i].leadDim;
          sRecInfo[reverseChannel[i]].nvec = sRecInfo[i].nvec;
        }
      }
    }
    else {
      // send and receive the message length and number of vectors
      int *sendMsg = new int[2 * sendConnect->numConnect()];
      int *recvMsg = new int[2 * recvConnect->numConnect()];
      int **sendPtr = new int*[numNeighbCPUs];
      int **recvPtr = new int*[numNeighbCPUs];
      int *sendLen = new int[numNeighbCPUs];
      int *recvLen = new int[numNeighbCPUs];
      int offset = 0;
      int totLen = 0;
      for(i = 0; i < numNeighbCPUs; ++i) {
        sendPtr[i] = sendMsg + offset;
        recvPtr[i] = recvMsg + offset;
        sendLen[i] = 2 * sendConnect->num(i);
        recvLen[i] = 2 * sendConnect->num(i);
        for(j = 0; j < sendConnect->num(i); ++j) {
          sendPtr[i][2*j] = sRecInfo[(*sendConnect)[i][j]].len;
          sendPtr[i][2*j+1] = sRecInfo[(*sendConnect)[i][j]].nvec;
          totLen += sendMsg[offset] * sendMsg[offset+1];
          offset += 2;
        }
      }
      communicator->exchange(101, numNeighbCPUs, neighbCPUs, sendLen, sendPtr,
                             recvLen, recvPtr);
      offset = 0;
      totLen = 0;
      for(i = 0; i < numNeighbCPUs; ++i)
        for(j = 0; j < sendConnect->num(i); ++j) {
          sRecInfo[(*recvConnect)[i][j]].len = recvMsg[offset];
          sRecInfo[(*recvConnect)[i][j]].leadDim = recvMsg[offset];
          sRecInfo[(*recvConnect)[i][j]].nvec = recvMsg[offset+1];
          totLen += recvMsg[offset]*recvMsg[offset+1];
          offset += 2;
        }
      delete [] sendMsg;
      delete [] recvMsg;
      delete [] sendPtr;
      delete [] recvPtr;
      delete [] sendLen;
      delete [] recvLen;
    }
  }

  // Step 4. allocate the cross-memory communication buffers (if necessary)
  if((mode == Share) && (numCPUs == 1)) return;
  if((numCPUs != 1) && (numNeighbCPUs > 0)) {
    int totLen = 0;
    crossSendBuffer = new T*[numNeighbCPUs];
    crossRecvBuffer  = new T*[numNeighbCPUs];
    crossSendLen = new int[numNeighbCPUs];
    crossRecvLen  = new int[numNeighbCPUs];
    for(i = 0; i < numNeighbCPUs; ++i)
      crossSendLen[i] = crossRecvLen[i] = 0;
    for(i = 0; i < numNeighbCPUs; ++i) {
      for(j = 0; j < sendConnect->num(i); ++j) {
        int channel = (*sendConnect)[i][j];
        int msgLen = sRecInfo[channel].len * sRecInfo[channel].nvec;
        crossSendLen[i] += msgLen;
        totLen += msgLen;
      }
      for(j = 0; j < recvConnect->num(i); ++j) {
        int channel = (*recvConnect)[i][j];
        int msgLen = sRecInfo[channel].len * sRecInfo[channel].nvec;
        crossRecvLen[i] += msgLen;
        totLen += msgLen;
      }
    }
    T *crBuff = new T[totLen];
    for(i = 0; i < numNeighbCPUs; ++i) {
      crossSendBuffer[i] = crBuff;
      for(j = 0; j < sendConnect->num(i); ++j) {
        int channel = (*sendConnect)[i][j];
        int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
        sRecInfo[channel].data = crBuff;
        crBuff += msgLen;
      }
      crossRecvBuffer[i] = crBuff;
      for(j = 0; j < recvConnect->num(i); ++j) {
        int channel = (*recvConnect)[i][j];
        int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
        sRecInfo[channel].data = crBuff;
        crBuff += msgLen;
      }
    }
  }
  if(mode == Share) return;
        
  // allocate buffer for Copy On Send mode
  int len = 0;
  for(i = 0; i < numChannels; ++i)
    if((numCPUs == 1) || (channelToCPU[i] < 0))
      len += sRecInfo[i].len * sRecInfo[i].nvec;
  localDBuffer = new T[len];
  T *cBuf = localDBuffer;
  for(i = 0; i < numChannels; ++i)
    if(numCPUs == 1 || channelToCPU[i] < 0) {
      sRecInfo[i].data = cBuf;
      cBuf += sRecInfo[i].len * sRecInfo[i].nvec;
    }
}

template <class T>
void
FSCommPattern<T>::exchange()
{
  //if(numCPUs == 1) return; // PJSA 11-7-05 
  communicator->exchange(101, numNeighbCPUs, neighbCPUs, crossSendLen,
                         crossSendBuffer, crossRecvLen, crossRecvBuffer);
}


template <class T>
FSSubRecInfo<T>
FSCommPattern<T>::recData(int glFrom, int glTo)
{
  int cpuID = (*subToCPU)[glFrom][0];
  if(cpuID == myCPU) cpuID = -1;
  int channel = channelMap.find(Triplet(glFrom, glTo, cpuID))->second;
  FSSubRecInfo<T> ret = sRecInfo[channel];
  if(mode != Share) ret.leadDim = ret.len;
  return ret;
}

template <class T>
void
FSCommPattern<T>::sendData(int glFrom, int glTo, T *data)
{
  int cpuID = (*subToCPU)[glTo][0];
  if(cpuID == myCPU) cpuID = -1;
  int channel = channelMap.find(Triplet(glFrom, glTo, cpuID))->second;

  if((mode == Share) && ((numCPUs == 1) || (channelToCPU[channel] < 0))) {
    sRecInfo[channel].data = data;
    return;
  }
  // For CopyOnSend, or non-local communication copy the data into the right buffer
  int i, iVec;
  FSSubRecInfo<T> chObj = sRecInfo[channel];
  for(iVec = 0; iVec < chObj.nvec; ++iVec)
    for(i = 0; i < chObj.len; ++ i)
      chObj.data[iVec*chObj.len + i] = data[iVec*chObj.leadDim + i];
}

template <class T>
FSSubRecInfo<T>
FSCommPattern<T>::getSendBuffer(int glFrom, int glTo)
{
  int cpuID = (*subToCPU)[glTo][0];
  if(cpuID == myCPU) cpuID = -1;
  int channel = channelMap.find(Triplet(glFrom, glTo, cpuID))->second;
  FSSubRecInfo<T> chObj = sRecInfo[channel];
  if(mode == Share) // Check if we have a send buffer. If not, return NULL
    chObj.data = NULL;
  return chObj;
}

template <class T>
void
FSCommPattern<T>::makeSendAndRecvConnectivities()
{
 int i;
 // First step, find all the neighboring CPUs
 int *glToLocCPU = new int[numCPUs];
 for(i = 0; i < numCPUs; ++i)
   glToLocCPU[i] = -1;

 numNeighbCPUs = 0;
 for(i = 0; i < numChannels; ++i) {
    int cpuID;
    if((cpuID = channelToCPU[i]) >= 0) {
      if(glToLocCPU[cpuID] < 0) glToLocCPU[cpuID] = numNeighbCPUs++;
    }
 }
 neighbCPUs = new int[numNeighbCPUs];
 for(i = 0; i < numCPUs; ++i)
   if(glToLocCPU[i] >= 0)
     neighbCPUs[glToLocCPU[i]] = i;

 // Now count the send and receives (actually they should be the same for
 // a symmetric communication
 int *sendPtr = new int[numNeighbCPUs+1];
 int *recvPtr = new int[numNeighbCPUs+1];
 for(i = 0; i < numNeighbCPUs+1; ++i)
   sendPtr[i] = recvPtr[i] = 0;

 MapIter iter = channelMap.begin();
 while(iter != channelMap.end()) {
   Triplet key = (*iter).first;
   int cpuID = key.cpuID;
   if(cpuID >= 0)
     if((*subToCPU)[key.from][0] == myCPU) // send pair
       sendPtr[glToLocCPU[cpuID]]++;
     else // receive
       recvPtr[glToLocCPU[cpuID]]++;
   ++iter;
 }

 // Make the actual pointers (shifted for easy fill in)
 for(i = 0; i < numNeighbCPUs; ++i) {
   recvPtr[i+1] += recvPtr[i];
   sendPtr[i+1] += recvPtr[i];
 }

 int *sendTrgt = new int[sendPtr[numNeighbCPUs]];
 int *recvTrgt = new int[recvPtr[numNeighbCPUs]];
 // Final fill in
 iter = channelMap.begin();
 while(iter != channelMap.end()) {
   Triplet key = (*iter).first;
   int channelID = (*iter).second;
   int cpuID = key.cpuID;
   if(cpuID >= 0) {
     if((*subToCPU)[key.from][0] == myCPU) {
       sendPtr[glToLocCPU[cpuID]]--;
       sendTrgt[sendPtr[glToLocCPU[cpuID]]] = channelID;
     }
     else { 
       recvPtr[glToLocCPU[cpuID]]--;
       recvTrgt[recvPtr[glToLocCPU[cpuID]]] = channelID;
     }
   }
   ++iter;
 }

 sendConnect = new Connectivity(numNeighbCPUs, sendPtr, sendTrgt);
 recvConnect = new Connectivity(numNeighbCPUs, recvPtr, recvTrgt);

 delete [] glToLocCPU;
// delete [] sendPtr; 
// delete [] recvPtr;
// delete [] sendTrgt;
// delete [] sendPtr;

}

#ifdef _TEMPLATE_FIX_
#include <Driver.d/Communicator.C>
#endif

#endif
