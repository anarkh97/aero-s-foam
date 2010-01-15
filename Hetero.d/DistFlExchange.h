#ifndef _DISTFLEXCHANGE_H_
#define _DISTFLEXCHANGE_H_

#include <Element.d/Element.h>
#include <Utils.d/OutputInfo.h>
#include <Hetero.d/FlExchange.h>

#include <map>
#include <algorithm>
using namespace std;

class CoordSet;
class State;
class Connectivity;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class VecType> class SysState;

#define FL_NEGOT 10000

// This is an inerpolation point. The element elemNum uses x and y to
// compute the interpolated displacements/velocities

/*
struct InterpPoint {
    int subNumber;
    int elemNum;
    double xy[2];
    int *dofs;
};
*/

typedef map<int, InterpPoint> MatchMap;

class DofSetArray;

class DistFlExchanger {
  double *buffer;
  //double *buff;
  int bufferLen;
  //int buffLen;

  //double *pArray;
  //int     pArrayLen;

  //int *nbData;
  //int *senderId;

  int numFluidNeighbors;   // number of fluid mpi's w/matches in this mpi
  int *idSendTo;  	   // list of fluid mpi's to sendTo
  int *nbSendTo;	   // num of match data per fluid neighbor
  //int *consOrigin; // reverse table of idSendTo
  InterpPoint **sndTable;  // match data by local subdomain in mpi
  //InterpPoint *globMatches;// all match data from fluid mpi's

  //int numWetElements;
  Elemset **eset;
  DofSetArray **cdsa;
  DofSetArray **dsa;

  double *localF;
  double aforce[4];
  double aflux;
  int rcvParity, sndParity;
  OutputInfo *oinfo;
  int isCollocated;
  double alpha[2];
  double alph[2];
  double dt;
  double dtemp;
  DistrVector *tmpDisp;
  Vector *dsp;
  Vector *vel;
  Vector *acc; 
  Vector *pVel; 
  State *newState;

  CoordSet **cs;      	     // nodes in this mpi process

public:

  DistFlExchanger(CoordSet **, Elemset **, DofSetArray **, 
		DofSetArray **, OutputInfo *oinfo = 0);
  MatchMap getMatchData();
  void negotiate();
  void sendDisplacements(SysState<DistrVector> &, double **, double **, int = -1);
  void sendParam(int, double, double, int, int, double a[2]);
  double getFluidLoad(DistrVector &, int, double,
                      double alphaf, int& iscollocated);

  void setBufferLength(int size)  { bufferLen = size; }
 

  int cmdCom(int);

  void initSndParity(int pinit) { sndParity = pinit; }
  void initRcvParity(int pinit) { rcvParity = pinit; }
  void flipRcvParity() { if(rcvParity >= 0) rcvParity = 1-rcvParity; }
  void flipSndParity() { if(sndParity >= 0) sndParity = 1-sndParity; }
} ;

#define FLTOSTMT 1000
#define STTOFLMT 2000
#define STTOSTMT 4000
#define FLTOSTHEAT 5000
#define STTOFLHEAT 6000
#define STCMDMSG 8000
#define FLCMDMSG 9000
#define OPTPARMSG 8100
#define OPTRESMSG 9100
#define NBPRESSDATAMAX 7

#endif

