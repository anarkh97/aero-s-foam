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
class DistrGeomState;

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
  double *buffer, *buff;
  int bufferLen, buffLen;

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

  CoordSet **cs;      	     // nodes in this mpi process

public:

  DistFlExchanger(CoordSet **, Elemset **, DofSetArray **, 
		DofSetArray **, OutputInfo *oinfo = 0);
  MatchMap* getMatchData();
  void negotiate();
  void thermoread(int);

  void sendDisplacements(SysState<DistrVector>&, double**, double**, int = -1, DistrGeomState* = 0);
  void sendTemperature(SysState<DistrVector>&);
  void sendStrucTemp(DistrVector&);
  
  double getFluidLoad(DistrVector&, int, double, double, int&, DistrGeomState* = 0);
  double getFluidFlux(DistrVector& flux, double time, double &bflux);
  void getStrucTemp(double*);

  void sendParam(int, double, double, int, int, double a[2]);
  void sendTempParam(int algnum, double step, double totaltime,
                     int rstinc, double alphat[2]);

  void setBufferLength(int size)  { bufferLen = size; }
 
  void sendModeFreq(double *modFrq, int numFrq);
  void sendModeShapes(int numFrq, int nNodes, double (**)[6],
                      SysState<DistrVector>&, double factor = 1.0);

  void initSndParity(int pinit) { sndParity = pinit; }
  void initRcvParity(int pinit) { rcvParity = pinit; }
  void flipRcvParity() { if(rcvParity >= 0) rcvParity = 1-rcvParity; }
  void flipSndParity() { if(sndParity >= 0) sndParity = 1-sndParity; }

 int cmdCom(int);

};

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
#define FL_NEGOT 10000

#endif

