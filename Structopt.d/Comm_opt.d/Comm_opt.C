#include <cassert>
#include <limits.h>
#include <Structopt.d/Comm_opt.d/Comm_opt.h>
#include <Utils.d/Connectivity.h>

//------------------------------------------------------------------------------
template<class T, class Scalar>
void DistVecReducer<T, Scalar>::reduce(DistVec<T>& v)
{
  const int maxSub = 512;
  T* subs[maxSub];
  for(int iNode=0; iNode < nodeToSub->csize(); ++iNode)
    {
      const int numSub = nodeToSub->num(iNode);
      int mySubs = 0;
      if(numSub <= 1) { continue; }
      assert(numSub <= maxSub);
      for(int iSub=0; iSub<numSub; ++iSub)
	{
	  const int glbSubNo = nodeToSub[iNode][iSub];
	  const int locSubNo = glSubToLocal? glSubToLocal[glbSubNo] : glbSubNo;
	  if(locSubNo<0) { continue; }
	  const int locNode = subDomains[locSubNo]->globalToLocal(iNode);
	  if(locNode < 0) { continue; }
	  subs[mySubs] = v.subData(locSubNo) + locNode;
	  ++mySubs;
	}
      reduceOp(mySubs, subs);
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void MaxDistVecReducer<Scalar>::reduceOp(int n, double** ptr)
{
  double buffer = MINDOUBLE;
  assert(ptr != 0);
  for(int i=0; i<n; ++i)
    { buffer = std::max(*(ptr[i]), buffer); }
  this->fsComm->globalMax(1, &buffer);
  for(int i=0; i<n; ++i)
    { *(ptr[i]) = buffer; }
  return;
}

//------------------------------------------------------------------------------
template<class T, class Scalar>
void SumDistVecReducer<T, Scalar>::reduceOp(int n, T** ptr)
{
  T buffer = 0;
  assert(ptr != 0);
  for(int i=0; i<n; ++i)
    { buffer += *(ptr[i]); }
  this->fsComm->globalSum(1, &buffer);
  for(int i=0; i<n; ++i)
    { *(ptr[i]) = buffer; }
}

//------------------------------------------------------------------------------
template<class Scalar>
void PowSumDistVecReducer<Scalar>::reduceOp(int n, double** ptr)
{
  double buffer = 0;
  assert(ptr != 0);
  for(int i=0; i<n; ++i)
    { 
      double dummy = *(ptr[i])>0;
      if(dummy > 0)
	{ buffer += pow(dummy, pw); }
    }
  this->fsComm->globalSum(1, &buffer);
  buffer = (buffer>0.0)? pow(buffer, 1.0/pw) : 0.0;
  for(int i=0; i<n; ++i)
    { *(ptr[i]) = buffer;  }
}

